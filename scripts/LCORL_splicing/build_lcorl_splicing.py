#!/usr/bin/env python3
"""Build a coordinate-validated LCORL splicing animation on Bunya.

Python >=3.8; standard library only. All input coordinates are 1-based inclusive.
No 3D model, GUI, GPU, browser, network, or rendering library is needed to build.
The result is a self-contained offline HTML/SVG animation, not a slide sequence.
"""
import argparse
import csv
import gzip
import hashlib
import io
import json
import platform
import re
import shutil
import sys
import zipfile
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import unquote

PROJECT = Path('/QRISdata/Q9774/lcorl_project')
SOURCE = PROJECT / 'results/LCORL_11SNP_isoform_comparison_20261025_231527'
GFF = PROJECT / 'LCORL_full_bam_analysis/annotation/Bos_taurus.ARS-UCD1.2.105.chr.gff3'
COMP = str.maketrans('ACGT', 'TGCA')
CODE = dict(zip(
    (a+b+c for a in 'TCAG' for b in 'TCAG' for c in 'TCAG'),
    'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'))


def require(ok, message):
    if not ok:
        raise ValueError(message)


def stable(value):
    return re.sub(r'\.\d+$', '', value.split(':')[-1])


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024*1024), b''):
            h.update(block)
    return h.hexdigest()


def fasta(path):
    result, header = {}, None
    for line in path.read_text().splitlines():
        if line.startswith('>'):
            header = line[1:].split('|')[0]
            require(header not in result, 'Duplicate FASTA header: '+header)
            result[header] = ''
        elif line.strip():
            require(header is not None, 'Sequence before FASTA header')
            result[header] += line.strip().upper()
    return result


def read_gff(path, transcript_ids):
    """Retain exact rows used; no guessed exon count or inferred intron boundaries."""
    tx = {t: {'exons': [], 'cds': []} for t in transcript_ids}
    retained = []
    opener = gzip.open if path.suffix == '.gz' else open
    with opener(path, 'rt') as stream:
        for line in stream:
            if line.startswith('##FASTA'):
                break
            if line.startswith('#') or not line.strip():
                continue
            f = line.rstrip('\n').split('\t')
            if len(f) != 9 or f[2] not in ('mRNA', 'transcript', 'exon', 'CDS'):
                continue
            a = {unquote(k): unquote(v) for k, v in
                 (x.split('=', 1) for x in f[8].split(';') if '=' in x)}
            targets = [stable(a.get('ID', ''))] if f[2] in ('mRNA', 'transcript') else [
                stable(p) for p in a.get('Parent', '').split(',')]
            targets = [t for t in targets if t in tx]
            if not targets:
                continue
            retained.append(line.rstrip('\n'))
            for t in targets:
                record = tx[t]
                require(f[6] in ('+', '-'), 'Missing strand for '+t)
                if 'strand' in record:
                    require(record['strand'] == f[6] and record['chrom'] == f[0],
                            'Mixed chromosome/strand for '+t)
                record.update(strand=f[6], chrom=f[0])
                if f[2] in ('mRNA', 'transcript'):
                    record['name'] = a.get('Name', '')
                else:
                    segment = {'start': int(f[3]), 'end': int(f[4]),
                               'id': a.get('exon_id', a.get('ID', '')),
                               'phase': f[7]}
                    record['exons' if f[2] == 'exon' else 'cds'].append(segment)
    return tx, retained


def order_segments(segments, strand, label):
    result = sorted(segments, key=lambda s: s['start'], reverse=(strand == '-'))
    genomic = sorted(result, key=lambda s: s['start'])
    require(bool(result), 'No '+label+' features found')
    require(all(s['start'] <= s['end'] for s in result), 'Invalid '+label+' interval')
    require(all(a['end'] < b['start'] for a, b in zip(genomic, genomic[1:])),
            'Overlapping or duplicate '+label+' features')
    offset = 0
    for s in result:
        s['offset'] = offset
        s['length'] = s['end'] - s['start'] + 1
        offset += s['length']
    return result


def position_in(segments, position, strand):
    found = [s for s in segments if s['start'] <= position <= s['end']]
    require(len(found) == 1, 'Position %s does not map uniquely to an exon/CDS' % position)
    s = found[0]
    return s['offset'] + (s['end']-position if strand == '-' else position-s['start']) + 1


def genomic_from_offset(segments, position, strand):
    found = [s for s in segments if s['offset'] < position <= s['offset']+s['length']]
    require(len(found) == 1, 'Invalid CDS offset: '+str(position))
    s = found[0]
    delta = position-s['offset']-1
    return s['end']-delta if strand == '-' else s['start']+delta


def build_transcript(isoform, txid, annotation, rows, protein):
    require('strand' in annotation, 'Transcript absent from annotation: '+txid)
    strand = annotation['strand']
    require(not annotation.get('name') or annotation['name'] == isoform,
            'Annotation isoform name disagrees with TSV: '+txid)
    exons = order_segments(annotation['exons'], strand, isoform+' exon')
    cds = order_segments(annotation['cds'], strand, isoform+' CDS')
    total_cds = sum(s['length'] for s in cds)
    require(total_cds in (len(protein)*3, (len(protein)+1)*3),
            isoform+': CDS length does not match reference protein')
    # Internal phase bases finish a preceding codon: never trim each CDS block.
    for c in cds:
        require(c['phase'] in ('0', '1', '2'), isoform+': missing CDS phase')
        require(int(c['phase']) == (-c['offset']) % 3,
                isoform+': inconsistent CDS phase or incomplete 5-prime CDS')
        require(sum(e['start'] <= c['start'] and c['end'] <= e['end'] for e in exons) == 1,
                isoform+': CDS is not contained in exactly one exon')
    for i, e in enumerate(exons, 1):
        e['number'] = i
        e['rna_start'], e['rna_end'] = e['offset']+1, e['offset']+e['length']
        e['coding'] = []
        for c in cds:
            lo, hi = max(e['start'], c['start']), min(e['end'], c['end'])
            if lo <= hi:
                aa = position_in(exons, lo, strand)
                bb = position_in(exons, hi, strand)
                e['coding'].append([min(aa, bb)-e['rna_start'], max(aa, bb)-e['rna_start']+1])
    first_cds_genomic = genomic_from_offset(cds, 1, strand)
    last_cds_genomic = genomic_from_offset(cds, total_cds, strand)
    cds_rna_start = position_in(exons, first_cds_genomic, strand)
    cds_rna_end = position_in(exons, last_cds_genomic, strand)
    require(cds_rna_end-cds_rna_start+1 == total_cds,
            isoform+': discontinuous CDS in mature RNA')
    variants = []
    seen = set()
    for r in sorted(rows, key=lambda r: int(r['CDS_position'])):
        require(r['consequence'] == 'missense', isoform+': unexpected consequence')
        require(stable(r['transcript']) == txid, isoform+': transcript mismatch')
        require(re.sub(r'^chr', '', r['chrom']) == re.sub(r'^chr', '', annotation['chrom']),
                isoform+': chromosome mismatch')
        pos = int(r['position'])
        require(pos not in seen, isoform+': duplicate SNP position')
        seen.add(pos)
        cpos = position_in(cds, pos, strand)
        rpos = position_in(exons, pos, strand)
        ppos = (cpos-1)//3+1
        require(cpos == int(r['CDS_position']), isoform+': TSV/GFF CDS mismatch at '+str(pos))
        require(ppos == int(r['protein_position']), isoform+': TSV/GFF protein mismatch')
        require(rpos-cds_rna_start+1 == cpos, isoform+': mRNA/CDS offset mismatch')
        require(1 <= ppos <= len(protein) and protein[ppos-1] == r['reference_AA'],
                isoform+': reference residue mismatch')
        rc, ac = r['reference_codon'].upper(), r['alternate_codon'].upper()
        require(CODE.get(rc) == r['reference_AA'] and CODE.get(ac) == r['alternate_AA'],
                isoform+': codon translation mismatch')
        codon_offset = (cpos-1)%3
        ref, alt = r['genomic_REF'].upper(), r['genomic_ALT'].upper()
        require(len(ref) == len(alt) == 1 and ref in 'ACGT' and alt in 'ACGT', 'Not a SNP')
        cr, ca = (ref.translate(COMP), alt.translate(COMP)) if strand == '-' else (ref, alt)
        require(rc[codon_offset] == cr and ac == rc[:codon_offset]+ca+rc[codon_offset+1:],
                isoform+': coding/genomic strand mismatch')
        exon = next(e for e in exons if e['start'] <= pos <= e['end'])
        codon_genome = [genomic_from_offset(cds, (ppos-1)*3+j, strand) for j in (1, 2, 3)]
        codon_rna = [position_in(exons, p, strand) for p in codon_genome]
        variants.append(dict(chrom=r['chrom'], position=pos, genomic_ref=ref, genomic_alt=alt,
            exon=exon['number'], exon_id=exon['id'], exon_start=exon['start'], exon_end=exon['end'],
            rna_position=rpos, cds_position=cpos, protein_position=ppos,
            reference_aa=r['reference_AA'], alternate_aa=r['alternate_AA'],
            protein_change=r['reference_AA']+str(ppos)+r['alternate_AA'],
            reference_rna_codon=rc.replace('T', 'U'), alternate_rna_codon=ac.replace('T', 'U'),
            codon_offset=codon_offset, codon_genomic_positions=codon_genome,
            codon_rna_positions=codon_rna))
    return dict(isoform=isoform, transcript=txid, strand=strand, chrom=annotation['chrom'],
        exons=exons, cds=cds, rna_length=sum(e['length'] for e in exons),
        cds_length=total_cds, cds_rna_start=cds_rna_start, cds_rna_end=cds_rna_end,
        protein_length=len(protein), variants=variants)


HTML = r'''<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>LCORL · DNA to RNA to protein</title>
<style>
:root{color-scheme:light;font-family:Arial,Helvetica,sans-serif;color:#17293c;background:#f0f4f7}
*{box-sizing:border-box}body{margin:0}main{max-width:1440px;margin:auto;padding:26px 32px}
.eyebrow{font-size:12px;font-weight:700;letter-spacing:2px;color:#547083}h1{font-size:30px;margin:8px 0}
.sub{font-size:14px;color:#52677a;margin:0 0 20px;line-height:1.5}.panel{background:white;border:1px solid #d9e3eb;border-radius:16px;box-shadow:0 5px 22px #19364c08}
.controls{display:flex;gap:14px;align-items:end;flex-wrap:wrap;padding:18px 22px;border-bottom:1px solid #e4ebf0}
label{display:flex;flex-direction:column;gap:6px;font-size:12px;font-weight:700;color:#50657a}select,button{font:inherit;border:1px solid #cbd7e1;border-radius:8px;padding:10px 12px;background:#fff;color:#19364c}
button{font-weight:700;cursor:pointer}button.primary{background:#164d60;color:#fff;border-color:#164d60}select{min-width:190px}#variant{min-width:310px}
.progress{display:flex;align-items:center;gap:14px;padding:14px 22px}input[type=range]{flex:1;accent-color:#147e8c}#stage{font-size:13px;min-width:250px;color:#355369}
svg{display:block;width:100%;height:auto}svg text{font-family:Arial,Helvetica,sans-serif}.details{padding:16px 22px 20px;border-top:1px solid #e4ebf0;display:grid;grid-template-columns:repeat(3,1fr);gap:12px}.detail{background:#f5f8fb;padding:14px;border-radius:10px}.detail strong{display:block;font-size:12px;color:#587085;margin-bottom:7px}.detail span{font-size:15px;font-weight:700;line-height:1.5;overflow-wrap:anywhere}.notes{font-size:12px;line-height:1.6;color:#62778b;margin:16px 2px}
@media(max-width:750px){main{padding:14px}.details{grid-template-columns:1fr}#stage{min-width:120px}select,#variant{min-width:180px}h1{font-size:24px}}
</style></head><body><main>
<div class="eyebrow">LCORL · VARIANT COORDINATE TRACKING</div><h1>From genomic site to protein residue</h1>
<p class="sub">ARS-UCD1.2 · Ensembl 105 reference transcript models · Genomic → spliced mRNA → translated protein</p>
<div class="panel"><div class="controls">
<label>Reference isoform<select id="isoform"></select></label><label>Track a variant<select id="variant"></select></label>
<button class="primary" id="play">Play animation</button><button id="reset">Replay</button><button id="snapshot">Save SVG</button>
</div><div class="progress"><span id="stage"></span><input aria-label="Animation timeline" id="timeline" type="range" min="0" max="1000" value="0"><span id="clock">0.0 / 16 s</span></div>
<svg id="scene" xmlns="http://www.w3.org/2000/svg" font-family="Arial,Helvetica,sans-serif" viewBox="0 0 1280 815" role="img" aria-label="Continuous splicing and variant mapping animation"></svg>
<div class="details"><div class="detail"><strong>GENOMIC DNA · 1-BASED</strong><span id="dnaInfo"></span></div><div class="detail"><strong>SPLICED mRNA · 1-BASED</strong><span id="rnaInfo"></span></div><div class="detail"><strong>REFERENCE PROTEIN · 1-BASED</strong><span id="proteinInfo"></span></div></div></div>
<p class="notes">Teal: coding sequence (CDS). Pale blue: untranslated exon sequence. Grey: introns. Magenta: tracked SNP.<br>
Exon numbers refer to the selected transcript, in 5′→3′ order. The genomic/pre-mRNA layouts are schematic (exons enlarged, introns compressed); the mature mRNA and protein use their own linear coordinate scales. Exact coordinates are exported in the TSV tables.<br>
This depicts annotated reference splicing and a single candidate substitution at a time. It does not demonstrate tissue expression, variant-induced splice changes or a combined mutant haplotype. Protein graphics show sequence position, not a predicted fold. Codon/allele labels come from the archived variant-consequence table; GFF coordinates and reference protein residues are independently cross-checked.</p>
</main><script id="data" type="application/json">__DATA__</script><script>
'use strict';
const DATA=JSON.parse(document.getElementById('data').textContent), NS='http://www.w3.org/2000/svg';
const $=id=>document.getElementById(id), svg=$('scene'), iso=$('isoform'), variant=$('variant');
const C={ink:'#18354c',muted:'#627a8c',cds:'#128797',utr:'#c7deec',intron:'#a6b5c1',rose:'#c33378',line:'#dde7ee'};
let tr, v, progress=0, playing=false, previous=null, shapes=[], lines=[], dots=[], layout=[];
const left=95,right=1185,width=right-left, YG=120, YPRE=275,YRNA=425, YP=600;
function el(tag,attrs={},parent=svg){const n=document.createElementNS(NS,tag);Object.entries(attrs).forEach(([k,x])=>n.setAttribute(k,x));parent.appendChild(n);return n}
function txt(x,y,s,size=14,color=C.ink,anchor='start',parent=svg){const n=el('text',{x,y,fill:color,'font-size':size,'text-anchor':anchor},parent);n.textContent=s;return n}
function rect(x,y,w,h,fill,parent=svg){return el('rect',{x,y,width:Math.max(.1,w),height:h,rx:2,fill},parent)}
function line(x1,y1,x2,y2,color=C.line,parent=svg){return el('line',{x1,y1,x2,y2,stroke:color,'stroke-width':2},parent)}
function fade(node,a){node.setAttribute('opacity',Math.max(0,Math.min(1,a)))}
const clamp=x=>Math.max(0,Math.min(1,x)), smooth=x=>{x=clamp(x);return x*x*(3-2*x)}, mix=(a,b,t)=>a+(b-a)*t;
const rnaX=p=>left+(p-.5)/tr.rna_length*width;
const aaX=p=>left+(p-.5)/tr.protein_length*width;
function labelRow(y,n,title,sub){txt(32,y,n,13,'#14808e');txt(left,y,title,16);txt(left,y+22,sub,12,C.muted)}
function clear(){svg.replaceChildren();shapes=[];lines=[];dots=[]}
function rebuild(){
 clear();
 const neg=tr.strand==='-', exonSum=tr.exons.reduce((s,e)=>s+Math.sqrt(e.length),0), gap=tr.exons.length>1?Math.min(32,250/(tr.exons.length-1)):0;
 const available=width-gap*(tr.exons.length-1), minW=Math.min(28,available/tr.exons.length*.5), rest=available-minW*tr.exons.length;
 let x=left;
 layout=tr.exons.map(e=>{const w=minW+Math.sqrt(e.length)/exonSum*rest;const l={e,x,w,rx:left+e.offset/tr.rna_length*width,rw:e.length/tr.rna_length*width};x+=w+gap;return l});
 labelRow(35,'01','Genomic DNA',`Chromosome ${tr.chrom} · ${neg?'minus':'plus'} strand · transcription direction → (genomic coordinates ${neg?'decrease':'increase'})`);
 line(left,YG,right,YG,C.intron);
 txt(left,YG+52,neg?'Higher genomic coordinate':'Lower genomic coordinate',11,C.muted);
 txt(right,YG+52,neg?'Lower genomic coordinate':'Higher genomic coordinate',11,C.muted,'end');
 txt(left-30,YG+5,'5′',13,C.muted);txt(right+15,YG+5,'3′',13,C.muted);
 labelRow(215,'02','RNA processing',`Pre-mRNA retains exon/intron order; splicing joins the ${tr.exons.length} annotated exons`);
 labelRow(360,'03','Mature mRNA',`UTRs retained · ${tr.rna_length.toLocaleString()} nt · CDS starts at mRNA nucleotide ${tr.cds_rna_start.toLocaleString()}`);
 line(left,YRNA,right,YRNA);
 txt(left,YRNA+49,'1 nt · 5′',11,C.muted);txt(right,YRNA+49,`${tr.rna_length.toLocaleString()} nt · 3′`,11,C.muted,'end');
 labelRow(535,'04','Protein sequence',`${tr.isoform} · ${tr.protein_length.toLocaleString()} amino acids · sequence position only`);
 proteinGroup=el('g');rect(left,YP-12,width,24,'#d5e5ed',proteinGroup);
 txt(left,YP+39,'N-terminus · 1 aa',11,C.muted,'start',proteinGroup);txt(right,YP+39,`${tr.protein_length.toLocaleString()} aa · C-terminus`,11,C.muted,'end',proteinGroup);
 bridge=el('path',{fill:'none',stroke:C.rose,'stroke-width':2,'stroke-dasharray':'6 5'});
 layout.forEach((l,i)=>{
   let g=el('g');rect(l.x,YG-16,l.w,32,C.utr,g);
   l.e.coding.forEach(([a,b])=>rect(l.x+a/l.e.length*l.w,YG-16,(b-a)/l.e.length*l.w,32,C.cds,g));
   txt(l.x+l.w/2,YG+5,`E${l.e.number}`,11,'#18354c','middle',g);
   const title=el('title',{},g);title.textContent=`Exon ${l.e.number}: ${tr.chrom}:${l.e.start}-${l.e.end}; mRNA ${l.e.rna_start}-${l.e.rna_end}`;
   const moving=el('g'), box=rect(l.x,YPRE-16,l.w,32,C.utr,moving), coding=l.e.coding.map(()=>rect(0,0,1,32,C.cds,moving)), number=txt(0,0,`E${l.e.number}`,11,C.ink,'middle',moving);
   shapes.push({l,g:moving,box,coding,number});
   if(i<layout.length-1) lines.push({a:l,b:layout[i+1],node:line(0,0,0,0,C.intron)});
 });
 tr.variants.forEach(site=>{const l=layout[site.exon-1],fraction=(neg?l.e.end-site.position:site.position-l.e.start)+.5;
   const gx=l.x+fraction/l.e.length*l.w,selected=site.position===v.position;
   const fixed=line(gx,YG-22,gx,YG+22,selected?C.rose:'#e8a3c3');fixed.setAttribute('stroke-width',selected?3:1.5);
   const node=el('circle',{r:selected?6:3,fill:selected?C.rose:'#e8a3c3'});dots.push({site,l,gx,node,selected});
 });
 selectedGenomeX=dots.find(d=>d.selected).gx;
 const anchor=selectedGenomeX>right-250?'end':selectedGenomeX<left+250?'start':'middle';
 txt(selectedGenomeX,YG-39,`${tr.chrom}:${v.position.toLocaleString()} ${v.genomic_ref}>${v.genomic_alt}`,14,C.rose,anchor);
 selectedLabel=txt(0,0,`Exon ${v.exon} → mRNA nt ${v.rna_position.toLocaleString()}`,13,C.rose);
 aaDot=el('circle',{cx:aaX(v.protein_position),cy:YP,r:7,fill:C.rose},proteinGroup);
 txt(aaX(v.protein_position),YP-25,v.protein_change,15,C.rose,aaX(v.protein_position)>right-160?'end':'start',proteinGroup);
 rect(left,685,width,98,'#f1f6f9');txt(left+18,710,'CODON DETAIL · mature mRNA in the coding direction',12,C.muted);
 let cx=left+22;[...v.reference_rna_codon].forEach((base,j)=>{rect(cx+j*42,728,34,34,j===v.codon_offset?C.rose:'#e0eaf1');txt(cx+j*42+17,752,base,21,j===v.codon_offset?'white':C.ink,'middle')});
 txt(left+164,752,'→',22,C.muted);cx=left+198;[...v.alternate_rna_codon].forEach((base,j)=>{rect(cx+j*42,728,34,34,j===v.codon_offset?C.rose:'#e0eaf1');txt(cx+j*42+17,752,base,21,j===v.codon_offset?'white':C.ink,'middle')});
 txt(left+350,749,`${v.reference_aa} → ${v.alternate_aa} · residue ${v.protein_position}`,19,C.ink);
 txt(left+350,772,`CDS nt ${v.cds_position}; full mRNA nt ${v.rna_position} (includes 5′ UTR)`,12,C.muted);
 $('dnaInfo').textContent=`${tr.chrom}:${v.position.toLocaleString()} ${v.genomic_ref}>${v.genomic_alt} · exon ${v.exon}`;
 $('rnaInfo').textContent=`nt ${v.rna_position.toLocaleString()} of ${tr.rna_length.toLocaleString()} · CDS nt ${v.cds_position.toLocaleString()}`;
 $('proteinInfo').textContent=`${v.protein_change} · ${tr.transcript}`;
 render();
}
let proteinGroup,bridge,selectedGenomeX,selectedLabel,aaDot;
function render(){
 const t=smooth((progress-.12)/.5), translation=smooth((progress-.67)/.25);
 shapes.forEach(s=>{const x=mix(s.l.x,s.l.rx,t),w=mix(s.l.w,s.l.rw,t),y=mix(YPRE,YRNA,t);s.box.setAttribute('x',x);s.box.setAttribute('y',y-16);s.box.setAttribute('width',w);
   s.coding.forEach((n,j)=>{const [a,b]=s.l.e.coding[j];n.setAttribute('x',x+a/s.l.e.length*w);n.setAttribute('width',(b-a)/s.l.e.length*w);n.setAttribute('y',y-16)});
   s.number.setAttribute('x',x+w/2);s.number.setAttribute('y',y+5);fade(s.number,w<24?1-t:1);
 });
 lines.forEach(l=>{const a=shapes[layout.indexOf(l.a)],b=shapes[layout.indexOf(l.b)],y=mix(YPRE,YRNA,t);l.node.setAttribute('x1',mix(l.a.x+l.a.w,l.a.rx+l.a.rw,t));l.node.setAttribute('x2',mix(l.b.x,l.b.rx,t));l.node.setAttribute('y1',y);l.node.setAttribute('y2',y);fade(l.node,1-t)});
 dots.forEach(d=>{d.node.setAttribute('cx',mix(d.gx,rnaX(d.site.rna_position),t));d.node.setAttribute('cy',mix(YPRE,YRNA,t))});
 const mx=mix(selectedGenomeX,rnaX(v.rna_position),t),my=mix(YPRE,YRNA,t);
 selectedLabel.setAttribute('x',Math.min(right-5,Math.max(left+5,mx)));selectedLabel.setAttribute('y',my-31);
 selectedLabel.setAttribute('text-anchor',mx>right-310?'end':mx<left+310?'start':'middle');
 selectedLabel.textContent=t<.99?`Exon ${v.exon} · introns removed; exon retained`:`Exon ${v.exon} · mRNA nt ${v.rna_position.toLocaleString()}`;
 fade(proteinGroup,translation);
 const sx=rnaX(v.rna_position),ex=aaX(v.protein_position),ey=mix(YRNA+18,YP-19,translation);
 bridge.setAttribute('d',`M ${sx} ${YRNA+18} C ${sx} ${YRNA+80}, ${ex} ${YP-85}, ${mix(sx,ex,translation)} ${ey}`);fade(bridge,translation);
 $('stage').textContent=progress<.12?'Pre-mRNA · exon/intron structure':progress<.63?'Splicing · exons join continuously':progress<.67?'Mature mRNA · introns removed':'Translation · follow the SNP to its residue';
 $('timeline').value=Math.round(progress*1000);$('clock').textContent=`${(progress*16).toFixed(1)} / 16 s`;
}
function chooseIsoform(){const old=v?v.position:null;tr=DATA.transcripts[Number(iso.value)];variant.replaceChildren();tr.variants.forEach((site,i)=>{const o=document.createElement('option');o.value=i;o.textContent=`${tr.chrom}:${site.position} ${site.genomic_ref}>${site.genomic_alt} · ${site.protein_change}`;variant.appendChild(o)});const index=tr.variants.findIndex(x=>x.position===old);variant.value=index>=0?index:0;chooseVariant()}
function chooseVariant(){v=tr.variants[Number(variant.value)];progress=0;playing=false;previous=null;$('play').textContent='Play animation';rebuild()}
DATA.transcripts.forEach((t,i)=>{const o=document.createElement('option');o.value=i;o.textContent=`${t.isoform} · ${t.protein_length} aa`;iso.appendChild(o)});iso.value=Math.max(0,DATA.transcripts.findIndex(x=>x.isoform==='LCORL-206'));
iso.onchange=chooseIsoform;variant.onchange=chooseVariant;
$('play').onclick=()=>{if(progress>=1)progress=0;playing=!playing;previous=null;$('play').textContent=playing?'Pause':'Play animation';render()};
$('reset').onclick=()=>{progress=0;playing=true;previous=null;$('play').textContent='Pause';render()};
$('timeline').oninput=()=>{playing=false;previous=null;progress=Number($('timeline').value)/1000;$('play').textContent='Play animation';render()};
$('snapshot').onclick=()=>{const s=svg.cloneNode(true);s.setAttribute('xmlns',NS);const blob=new Blob([new XMLSerializer().serializeToString(s)],{type:'image/svg+xml'});const u=URL.createObjectURL(blob),a=document.createElement('a');a.href=u;a.download=`${tr.isoform}_${v.position}_mapping.svg`;a.click();setTimeout(()=>URL.revokeObjectURL(u),1000)};
function animate(now){if(playing){if(previous!==null)progress=Math.min(1,progress+Math.min(now-previous,100)/16000);previous=now;render();if(progress>=1){playing=false;$('play').textContent='Play animation'}}requestAnimationFrame(animate)}
chooseIsoform();requestAnimationFrame(animate);
</script></body></html>'''


def write_tsv(path, records):
    require(bool(records), 'No rows to write to '+str(path))
    with path.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(records[0]), delimiter='\t')
        writer.writeheader()
        writer.writerows(records)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source', type=Path, default=SOURCE,
                   help='Archived folder containing variant_consequences.long.tsv and reference_proteins.fasta')
    p.add_argument('--gff', type=Path, default=GFF)
    p.add_argument('--isoforms', default='202,203,204,205,206')
    p.add_argument('--expected-sites', type=int, default=11)
    p.add_argument('--out', type=Path, help='New output directory (must not exist)')
    args = p.parse_args()
    started = datetime.now(timezone.utc)
    out = args.out or PROJECT / 'results' / ('LCORL_splicing_animation_'+started.strftime('%Y%m%d_%H%M%S_%f'))
    require(not out.exists(), 'Output directory already exists: '+str(out))
    mapping = args.source / 'variant_consequences.long.tsv'
    proteins_path = args.source / 'reference_proteins.fasta'
    for f in (args.gff, mapping, proteins_path):
        require(f.is_file(), 'Missing input: '+str(f))
    isoforms = ['LCORL-'+s.strip().replace('LCORL-', '') for s in args.isoforms.split(',')]
    require(len(isoforms) == len(set(isoforms)), 'Duplicate isoform selection')
    proteins = fasta(proteins_path)
    with mapping.open(newline='') as stream:
        rows = list(csv.DictReader(stream, delimiter='\t'))
    txids = {}
    for name in isoforms:
        require(name in proteins, 'Missing reference protein '+name)
        selected = [r for r in rows if r['isoform'] == name]
        require(len(selected) == args.expected_sites, name+': unexpected number of variant rows')
        ids = {stable(r['transcript']) for r in selected}
        require(len(ids) == 1, name+': ambiguous transcript ID')
        txids[name] = ids.pop()
    print('Reading archived annotation:', args.gff, flush=True)
    annotation, retained = read_gff(args.gff, set(txids.values()))
    transcripts, checks, variant_rows, exon_rows = [], [], [], []
    for name in isoforms:
        selected = [r for r in rows if r['isoform'] == name]
        t = build_transcript(name, txids[name], annotation[txids[name]], selected, proteins[name])
        transcripts.append(t)
        print('%s: %d exons; mRNA %d nt; CDS %d nt; protein %d aa; %d sites PASS' % (
            name, len(t['exons']), t['rna_length'], t['cds_length'], t['protein_length'], len(t['variants'])), flush=True)
        checks.append(dict(isoform=name, transcript=t['transcript'], strand=t['strand'],
            exon_count=len(t['exons']), mRNA_nt=t['rna_length'], CDS_nt=t['cds_length'],
            protein_aa=t['protein_length'], variants_checked=len(t['variants']),
            status='PASS: GFF mapping, CDS phase, RNA/CDS offsets, protein residue, codon/alleles'))
        for v in t['variants']:
            variant_rows.append(dict(isoform=name, transcript=t['transcript'], strand=t['strand'],
                **{k:(','.join(map(str,value)) if isinstance(value,list) else value) for k,value in v.items()}))
        for e in t['exons']:
            exon_rows.append(dict(isoform=name, transcript=t['transcript'], chrom=t['chrom'],
                strand=t['strand'], exon_number_5to3=e['number'], exon_id=e['id'],
                genomic_start=e['start'], genomic_end=e['end'],
                mRNA_start=e['rna_start'], mRNA_end=e['rna_end'], exon_nt=e['length']))
    site_sets = [{(v['chrom'],v['position'],v['genomic_ref'],v['genomic_alt']) for v in t['variants']} for t in transcripts]
    require(all(s == site_sets[0] for s in site_sets), 'Isoforms do not contain the same genomic SNPs')
    data = dict(title='LCORL DNA / splicing / protein position mapping',
                coordinate_convention='1-based inclusive', transcripts=transcripts)
    out.mkdir(parents=True)
    write_tsv(out/'variant_DNA_exon_mRNA_CDS_protein.tsv', variant_rows)
    write_tsv(out/'exon_coordinates.tsv', exon_rows)
    write_tsv(out/'validation.tsv', checks)
    (out/'animation_data.json').write_text(json.dumps(data, indent=2)+'\n')
    embedded = json.dumps(data).replace('<', '\\u003c').replace('>', '\\u003e')
    (out/'LCORL_splicing_animation.html').write_text(HTML.replace('__DATA__', embedded), encoding='utf-8')
    (out/'input_LCORL_features.gff3').write_text('##gff-version 3\n'+'\n'.join(retained)+'\n')
    for source_file in (mapping, proteins_path, args.source/'reference_inputs.txt'):
        if source_file.is_file():
            shutil.copy2(source_file, out/source_file.name)
    script = Path(__file__).resolve()
    shutil.copy2(script, out/'build_lcorl_splicing.py')
    provenance = dict(created_utc=started.isoformat(), hostname=platform.node(),
        python=sys.version, command=sys.argv, inputs=[dict(path=str(f.resolve()),sha256=sha(f))
        for f in (args.gff,mapping,proteins_path,script)],
        limitations=[
            'GFF-based coordinate reconstruction; genomic FASTA bases are not re-extracted in this script.',
            'Codons and genomic alleles use the archived consequence table, cross-checked internally and against protein FASTA.',
            'Reference transcript illustration; no tissue expression or splice-impact inference.',
            'Genomic/pre-mRNA drawing is schematic; mature RNA and protein have independent linear scales.'])
    (out/'provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')
    log = '\n'.join(str(r) for r in checks)+'\nOUTPUT: '+str(out.resolve())+'\n'
    (out/'run_summary.txt').write_text(log)
    (out/'README.txt').write_text(
        'Open LCORL_splicing_animation.html locally. It needs no server or internet.\n'
        'Choose an isoform and SNP, press Play, or drag the timeline. Save SVG exports the current vector frame.\n'
        'The genomic DNA stays intact. The pre-mRNA exon blocks move together as introns disappear.\n'
        'Each tracked site is retained in its exon and mapped to mature mRNA, CDS and protein.\n'
        'mRNA includes UTRs; CDS coordinates do not. All coordinates are 1-based.\n'
        'On the negative strand, exon 1 is at the high-coordinate end and RNA alleles are complemented.\n'
        'Exon numbering is transcript-specific. Codons crossing splice junctions retain all three mapped genomic positions in the TSV.\n'
        'Exact exon coordinates are from the supplied GFF. Do not use the schematic widths as physical genomic distances.\n'
        'Reference splicing only: this is not evidence that any variant alters splicing.\n'
        'The build uses only the Python standard library and leaves the script, input subset, hashes, tables and logs here.\n')
    (out/'output_sha256.json').write_text(json.dumps({f.name:sha(f) for f in sorted(out.iterdir()) if f.is_file()},indent=2)+'\n')
    archive = out.with_suffix('.zip')
    require(not archive.exists(), 'Archive already exists: '+str(archive))
    with zipfile.ZipFile(archive, 'x', zipfile.ZIP_DEFLATED) as z:
        for f in sorted(out.iterdir()):
            z.write(f, arcname=out.name+'/'+f.name)
    print('\nBUILD PASSED.\nOUTPUT DIRECTORY: '+str(out.resolve())+'\nDOWNLOAD ZIP: '+str(archive.resolve()), flush=True)


if __name__ == '__main__':
    try:
        main()
    except (OSError, ValueError, KeyError, csv.Error) as error:
        print('STOP: '+str(error), file=sys.stderr)
        sys.exit(1)
