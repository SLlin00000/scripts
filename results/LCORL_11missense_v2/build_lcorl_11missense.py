#!/usr/bin/env python3
"""Build the 11-site LCORL DNA -> coding RNA -> protein animation on Bunya.

Python >=3.8; standard library only. Reads the existing SNP-consequence table
and validated reference proteins. No exon reconstruction or GFF input is used.
All positions are 1-based; RNA positions refer to CDS, excluding UTRs.
"""
import argparse
import csv
import hashlib
import json
import platform
import shutil
import sys
import zipfile
from datetime import datetime, timezone
from pathlib import Path

PROJECT = Path('/QRISdata/Q9774/lcorl_project')
SOURCE = PROJECT / 'results/LCORL_11SNP_isoform_comparison_20261025_231527'
COMP = str.maketrans('ACGT', 'TGCA')
CODE = dict(zip((a+b+c for a in 'TCAG' for b in 'TCAG' for c in 'TCAG'),
    'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'))


def require(ok, message):
    if not ok:
        raise ValueError(message)


def sha(path):
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1048576), b''):
            digest.update(block)
    return digest.hexdigest()


def read_proteins(path):
    sequences, name = {}, None
    for line in path.read_text().splitlines():
        if line.startswith('>'):
            name = line[1:].split('|')[0]
            require(name not in sequences, 'Duplicate FASTA entry '+name)
            sequences[name] = ''
        elif line.strip():
            require(name is not None, 'FASTA sequence before header')
            sequences[name] += line.strip().upper()
    return sequences


def validate(rows, proteins):
    transcripts, seen_sets = [], []
    for number in (206,):
        isoform = 'LCORL-'+str(number)
        require(isoform in proteins, 'Missing reference protein '+isoform)
        selected = sorted((r for r in rows if r['isoform'] == isoform),
                          key=lambda r: int(r['CDS_position']))
        require(len(selected) == 11, isoform+': expected 11 SNPs')
        ids = {r['transcript'] for r in selected}
        chroms = {r['chrom'] for r in selected}
        require(len(ids) == len(chroms) == 1, isoform+': inconsistent transcript/chromosome')
        require(ids == {'ENSBTAT00000066859'}, 'Unexpected LCORL-206 transcript ID')
        require(len(proteins[isoform]) == 1904, 'Expected the established 1904-aa LCORL-206 reference')
        sites, strands, keys = [], set(), set()
        for r in selected:
            require(r['consequence'] == 'missense', isoform+': unexpected consequence')
            dna, cds, aa = int(r['position']), int(r['CDS_position']), int(r['protein_position'])
            require(1 <= cds <= len(proteins[isoform])*3, isoform+': invalid CDS position')
            require((cds-1)//3+1 == aa, isoform+': CDS/protein position mismatch')
            require(proteins[isoform][aa-1] == r['reference_AA'], isoform+': reference residue mismatch')
            rc, ac = r['reference_codon'].upper(), r['alternate_codon'].upper()
            ref, alt = r['genomic_REF'].upper(), r['genomic_ALT'].upper()
            offset = (cds-1)%3
            require(len(ref) == len(alt) == 1 and ref in 'ACGT' and alt in 'ACGT', 'Non-SNP input')
            require(CODE.get(rc) == r['reference_AA'] and CODE.get(ac) == r['alternate_AA'], 'Codon translation mismatch')
            possible = []
            for strand in ('+', '-'):
                coding_ref = ref if strand == '+' else ref.translate(COMP)
                coding_alt = alt if strand == '+' else alt.translate(COMP)
                if rc[offset] == coding_ref and ac == rc[:offset]+coding_alt+rc[offset+1:]:
                    possible.append((strand, coding_ref, coding_alt))
            require(len(possible) == 1, 'Ambiguous genomic/coding allele mapping')
            strand, coding_ref, coding_alt = possible[0]
            strands.add(strand)
            key = (r['chrom'], dna, ref, alt)
            require(key not in keys, 'Duplicate genomic SNP')
            keys.add(key)
            sites.append(dict(chrom=r['chrom'], dna=dna, dna_ref=ref, dna_alt=alt,
                cds=cds, rna_ref=coding_ref.replace('T','U'), rna_alt=coding_alt.replace('T','U'),
                aa=aa, change=r['reference_AA']+str(aa)+r['alternate_AA'],
                rna_codon_ref=rc.replace('T','U'), rna_codon_alt=ac.replace('T','U')))
        require(len(strands) == 1, isoform+': inconsistent strand')
        seen_sets.append(keys)
        transcripts.append(dict(isoform=isoform, transcript=ids.pop(), chrom=chroms.pop(),
            protein_length=len(proteins[isoform]), strand=strands.pop(), sites=sites))
    require(all(s == seen_sets[0] for s in seen_sets), 'Isoforms contain different genomic SNP sets')
    return dict(coordinates='1-based; RNA=CDS position excluding UTRs', transcripts=transcripts)


HTML = r'''<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>LCORL-206 · 11 missense sites</title><style>
:root{font-family:-apple-system,BlinkMacSystemFont,"Segoe UI",Arial,sans-serif;color:#162e42;background:#edf3f7;color-scheme:light}*{box-sizing:border-box}[hidden]{display:none!important}body{margin:0}main{max-width:1640px;margin:auto;padding:24px}h1{font-size:30px;margin:0 0 8px}p{font-size:18px;line-height:1.5;margin:0 0 18px;color:#486176}.panel{background:#fff;border:1px solid #ccdae5;border-radius:14px;overflow:hidden}.toolbar{display:flex;align-items:end;gap:16px;flex-wrap:wrap;padding:20px;border-bottom:1px solid #dbe5ed}label{display:flex;flex-direction:column;gap:8px;font-size:17px;font-weight:650}select,button{font:inherit;font-size:18px;padding:11px 13px;border:1px solid #aec2d2;border-radius:8px;color:#19364c;background:#fff}select{min-width:225px}button{font-weight:650;cursor:pointer}button.primary{background:#164f63;color:#fff}button:disabled{opacity:.4;cursor:default}.timeline{padding:16px 20px;display:flex;gap:16px;align-items:center;font-size:18px}input{flex:1;min-width:160px;accent-color:#137d8c}#stage{min-width:270px}.scroll{overflow-x:auto}#scene{width:100%;min-width:1120px;display:block;background:#fff}svg text{font-family:Arial,Helvetica,sans-serif}#empty{padding:65px 24px;text-align:center;font-size:23px;color:#49677d}h2{font-size:23px;margin:0 0 7px}.tableHead{padding:22px 20px 10px}.tableHead p{margin:0;font-size:17px}table{width:100%;border-collapse:collapse;font-size:18px;white-space:nowrap}th{text-align:left;background:#edf4f8;color:#24485e;font-size:17px}th,td{padding:13px 15px;border-top:1px solid #dce6ed}td{font-variant-numeric:tabular-nums}tbody tr:nth-child(even){background:#f7fafc}tbody tr{cursor:default}tbody tr.highlight{background:#e8f2fc;outline:2px solid #28678e;outline-offset:-2px}.swatch{display:inline-block;width:12px;height:12px;border-radius:50%;margin-right:9px}.notes{font-size:16px;color:#496478;margin:18px 0 0;line-height:1.6}#selection{font-size:18px;color:#2c536e;margin-bottom:12px}.status{font-size:17px;color:#44647a;padding:12px 20px;background:#f3f8fb;border-top:1px solid #dce6ed}#downloadTable{margin-left:auto}@media(max-width:800px){main{padding:12px}.toolbar{padding:14px;gap:12px}#stage{min-width:140px}.timeline{font-size:16px;flex-wrap:wrap}h1{font-size:26px}}
</style></head><body><main><h1>LCORL-206 · all 11 missense sites</h1>
<p>Genomic DNA → RNA coding position → protein residue · ENSBTAT00000066859 · 1904 aa</p>
<div class="panel"><div class="toolbar">
<button class="primary" id="play">Play animation</button><button id="replay">Replay</button><button id="save">Save SVG</button>
</div><div id="viewer"><div class="timeline"><span id="stage">DNA → coding RNA</span><input id="time" type="range" min="0" max="1000" value="0" aria-label="Animation progress"><span id="clock">0 / 12 s</span></div>
<div class="scroll"><svg id="scene" xmlns="http://www.w3.org/2000/svg" font-family="Arial,Helvetica,sans-serif" viewBox="0 0 1440 800" role="img" aria-label="All eleven missense sites mapped through DNA, coding RNA and protein"></svg></div>
<div class="status" id="selection"></div><div class="tableHead"><h2>All 11 site mappings</h2><p>Each site keeps the same colour across DNA, coding RNA and protein.</p></div>
<div class="scroll"><table><thead><tr><th>Site</th><th>Genomic DNA</th><th>RNA · CDS nucleotide</th><th>RNA codon</th><th>Protein substitution</th></tr></thead><tbody id="rows"></tbody></table></div>
</div></div><p class="notes">RNA coordinates refer to the coding sequence (CDS), excluding UTRs; they are not full-mRNA coordinates. All coordinates are 1-based. Each axis shows the interval spanned by the 11 target sites, with an independent linear scale. Labels use leader lines to keep nearby sites readable. Reference sequences are shown; marking candidate substitutions does not create a combined mutant haplotype.</p>
</main><script id="data" type="application/json">__DATA__</script><script>
'use strict';
const D=JSON.parse(document.getElementById('data').textContent),$=id=>document.getElementById(id),svg=$('scene'),NS='http://www.w3.org/2000/svg';
const palette=['#9b2468','#1768aa','#a0460d','#5260ae','#08796b','#b12834','#765329','#5c3e9f','#51751a','#006b83','#8c3d62'];
const L=104,R=1336,W=R-L,Y=[125,380,635],labelY=[[188,232],[443,487],[698,742]];
let tr=D.transcripts[0],progress=0,playing=false,last=null,marks=[],paths=[],balls=[],groups=[];
const clamp=x=>Math.max(0,Math.min(1,x)),ease=x=>{x=clamp(x);return x*x*(3-2*x)},mix=(a,b,t)=>a+(b-a)*t;
function node(tag,a={},parent=svg){const n=document.createElementNS(NS,tag);Object.entries(a).forEach(([k,v])=>n.setAttribute(k,v));parent.appendChild(n);return n}
function text(x,y,value,size=22,fill='#183b50',anchor='start',parent=svg){const n=node('text',{x,y,'font-size':size,fill,'text-anchor':anchor,'font-family':'Arial,Helvetica,sans-serif'},parent);n.textContent=value;return n}
function opacity(n,value){n.setAttribute('opacity',value)}
function scale(key,value){const values=tr.sites.map(s=>s[key]),a=Math.min(...values),b=Math.max(...values);const f=(value-a)/(b-a||1);return L+(key==='dna'&&tr.strand==='-'?1-f:f)*W}
function labelPositions(items,gap){let xs=items.map(s=>s.x);for(let i=0;i<xs.length;i++)xs[i]=Math.max(xs[i],i?xs[i-1]+gap:125);if(xs.length&&xs.at(-1)>1315)xs[xs.length-1]=1315;for(let i=xs.length-2;i>=0;i--)xs[i]=Math.min(xs[i],xs[i+1]-gap);return xs}
function heading(y,title,subtitle){node('rect',{x:L-5,y:y-29,width:W+10,height:70,fill:'white'});text(L,y,title,29);text(L,y+33,subtitle,20,'#46687f')}
function renderDiagram(){
 svg.replaceChildren();marks=[];paths=[];balls=[];groups=[];
 const neg=tr.strand==='-';
 const leaders=node('g');
 tr.sites.forEach((s,i)=>{const xs=[scale('dna',s.dna),scale('cds',s.cds),scale('aa',s.aa)],color=palette[i];
   const first=node('path',{d:`M ${xs[0]} ${Y[0]+17} C ${xs[0]} ${Y[0]+82},${xs[1]} ${Y[1]-82},${xs[1]} ${Y[1]-17}`,stroke:color,'stroke-width':2,fill:'none'},leaders);
   const second=node('path',{d:`M ${xs[1]} ${Y[1]+17} C ${xs[1]} ${Y[1]+82},${xs[2]} ${Y[2]-82},${xs[2]} ${Y[2]-17}`,stroke:color,'stroke-width':2,fill:'none'},leaders);
   paths.push({first,second,i});balls.push({node:node('circle',{r:7,fill:color,stroke:'white','stroke-width':2}),xs,i});
 });
 for(let k=0;k<3;k++){
   const g=node('g');groups.push(g);node('rect',{x:L,y:Y[k]-14,width:W,height:28,rx:8,fill:'#e2edf4'},g);
   const items=tr.sites.map((s,i)=>({s,i,x:scale(['dna','cds','aa'][k],s[['dna','cds','aa'][k]])}));
   items.forEach(a=>marks.push({node:node('circle',{cx:a.x,cy:Y[k],r:6,fill:palette[a.i],stroke:'white','stroke-width':2,'data-layer':k,'data-site':a.i},g),i:a.i,k}));
   const labels=[];
   for(let lane=0;lane<2;lane++){
     const batch=items.filter(a=>a.i%2===lane),xs=labelPositions(batch,k===0?214:175);
     batch.forEach((a,j)=>{const x=xs[j],y=labelY[k][lane],color=palette[a.i];node('path',{d:`M ${a.x} ${Y[k]+15} L ${a.x} ${Y[k]+27} L ${x} ${y-24}`,stroke:color,'stroke-width':1.8,fill:'none'},g);labels.push({a,x,y,color})});
   }
   labels.forEach(({a,x,y,color})=>{
     const value=k===0?`${a.s.dna.toLocaleString()} ${a.s.dna_ref}>${a.s.dna_alt}`:k===1?`${a.s.cds.toLocaleString()} ${a.s.rna_ref}>${a.s.rna_alt}`:a.s.change;
     const w=k===0?211:k===1?178:128;
     node('rect',{x:x-w/2,y:y-25,width:w,height:34,rx:5,fill:'white',stroke:color,'stroke-width':1.4},g);
     text(x,y,value,k===0?21:23,'#102d43','middle',g).setAttribute('data-site-label',k+':'+a.i);
   });
 }
 // The moving markers stay above the tracks and label leaders.
 balls.forEach(b=>svg.appendChild(b.node));
 heading(38,'Genomic DNA',`Chromosome ${tr.chrom} · reported genomic alleles · coordinates ${neg?'decrease':'increase'} left to right`);
 heading(293,'RNA coding sequence',`CDS nucleotide positions in ${tr.isoform} · 5′ → 3′ · UTRs excluded`);
 heading(548,'Protein sequence',`${tr.isoform} · reference length ${tr.protein_length.toLocaleString()} aa · amino-acid positions`);
 text(L,790,'All 11 sites shown together · labels point to their exact positions',20,'#45667c');
 drawFrame();
}
function drawFrame(){
 if(!tr)return;const a=ease(progress*2),b=ease((progress-.5)*2);
 groups.forEach((g,k)=>opacity(g,k===0?1:k===1?.25+.75*a:.25+.75*b));
 paths.forEach(p=>{opacity(p.first,.20*a);opacity(p.second,.20*b)});
 balls.forEach(ball=>{const first=progress<.5,t=first?a:b,k=first?0:1;ball.node.setAttribute('cx',mix(ball.xs[k],ball.xs[k+1],t));ball.node.setAttribute('cy',mix(Y[k],Y[k+1],t))});
 $('time').value=Math.round(progress*1000);$('clock').textContent=`${(progress*12).toFixed(1)} / 12 s`;$('stage').textContent=progress<.5?'DNA → coding RNA':'Coding RNA → protein';
}
function initialise(){
 $('selection').textContent=`${tr.isoform} · ${tr.transcript} · ARS-UCD1.2 / Ensembl 105 · 11 of 11 sites shown`;
 $('rows').replaceChildren();
 tr.sites.forEach((s,i)=>{
   const r=document.createElement('tr');const values=[String(i+1),`${s.chrom}:${s.dna.toLocaleString()} ${s.dna_ref}>${s.dna_alt}`,`${s.cds.toLocaleString()} ${s.rna_ref}>${s.rna_alt}`,`${s.rna_codon_ref} → ${s.rna_codon_alt}`,s.change];
   values.forEach((value,j)=>{const td=document.createElement('td');if(j===0){const dot=document.createElement('span');dot.className='swatch';dot.style.background=palette[i];td.appendChild(dot)}const span=document.createElement('span');span.textContent=value;td.appendChild(span);r.appendChild(td)});
   $('rows').appendChild(r);
 });renderDiagram();
}
$('play').onclick=()=>{if(progress===1)progress=0;playing=!playing;last=null;$('play').textContent=playing?'Pause':'Play animation';drawFrame()};
$('replay').onclick=()=>{progress=0;playing=true;last=null;$('play').textContent='Pause';drawFrame()};
$('time').oninput=()=>{playing=false;last=null;progress=Number($('time').value)/1000;$('play').textContent='Play animation';drawFrame()};
$('save').onclick=()=>{const copy=svg.cloneNode(true);copy.setAttribute('xmlns',NS);const url=URL.createObjectURL(new Blob([new XMLSerializer().serializeToString(copy)],{type:'image/svg+xml'})),a=document.createElement('a');a.href=url;a.download=`${tr.isoform}_all_11_missense.svg`;a.click();setTimeout(()=>URL.revokeObjectURL(url),1000)};
function animate(now){if(playing&&tr){if(last!==null)progress=Math.min(1,progress+Math.min(100,now-last)/12000);last=now;drawFrame();if(progress===1){playing=false;$('play').textContent='Play animation'}}requestAnimationFrame(animate)}initialise();requestAnimationFrame(animate);
</script></body></html>'''


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, default=SOURCE)
    parser.add_argument('--out', type=Path, default=PROJECT/'results/LCORL_11missense_v2')
    args = parser.parse_args()
    out, source = args.out, args.source
    require(not out.exists() and not out.with_suffix('.zip').exists(), 'Output already exists: '+str(out))
    mapping, proteins = source/'variant_consequences.long.tsv', source/'reference_proteins.fasta'
    require(mapping.is_file() and proteins.is_file(), 'Missing archived mapping table or reference FASTA in '+str(source))
    with mapping.open(newline='') as stream:
        rows = list(csv.DictReader(stream, delimiter='\t'))
    data = validate(rows, read_proteins(proteins))
    out.mkdir(parents=True)
    (out/'LCORL_11missense.html').write_text(HTML.replace('__DATA__', json.dumps(data).replace('<','\\u003c')), encoding='utf-8')
    (out/'animation_data.json').write_text(json.dumps(data,indent=2)+'\n')
    flattened = [dict(isoform=t['isoform'],transcript=t['transcript'],**s) for t in data['transcripts'] for s in t['sites']]
    with (out/'LCORL206_11_site_mappings.tsv').open('w',newline='') as stream:
        writer = csv.DictWriter(stream,fieldnames=list(flattened[0]),delimiter='\t');writer.writeheader();writer.writerows(flattened)
    for f in (mapping,proteins,Path(__file__).resolve()):
        shutil.copy2(f,out/f.name)
    provenance=dict(created_utc=datetime.now(timezone.utc).isoformat(),hostname=platform.node(),python=sys.version,
        command=sys.argv,inputs=[dict(path=str(f.resolve()),sha256=sha(f)) for f in (mapping,proteins,Path(__file__))],
        scope='Display the existing 11 missense site mappings using the established 1904-aa LCORL-206 reference.',
        validation='11 mappings checked: CDS/protein offsets, reference residues, codon translations, genomic/coding alleles.',
        reference_choice='LCORL-206 is retained to keep reference sequence and residue numbering consistent with the earlier 1904-aa analysis; no claim of predominant expression or superior structural confidence.',
        coordinates='1-based. RNA=CDS nucleotide positions excluding UTRs.',
        note='No GFF or exon reconstruction. Existing variant consequences are not newly derived from a genome assembly.')
    (out/'provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')
    (out/'README.txt').write_text(
        'Open LCORL_11missense.html offline. Reference fixed to LCORL-206 / ENSBTAT00000066859 / 1904 aa.\n'
        'All 11 sites are labelled together on DNA, coding RNA and protein.\n'
        'RNA positions are CDS nucleotide positions, excluding UTRs. All positions are 1-based.\n'
        'Each track uses an independent linear scale over the first-to-last target interval.\n'
        'Play/replay or drag the timeline. Save SVG exports all 11 site labels on the current frame.\n'
        'The source script, input table, reference FASTA, hashes and all 11 mappings are retained here.\n')
    summary='\n'.join('%s: %d / 11 sites PASS' % (t['isoform'],len(t['sites'])) for t in data['transcripts'])+'\n'
    (out/'validation.txt').write_text(summary)
    (out/'output_sha256.json').write_text(json.dumps({f.name:sha(f) for f in sorted(out.iterdir())},indent=2)+'\n')
    archive=out.with_suffix('.zip')
    with zipfile.ZipFile(archive,'x',zipfile.ZIP_DEFLATED) as z:
        for f in sorted(out.iterdir()):
            z.write(f,arcname=out.name+'/'+f.name)
    print(summary+'BUILD PASSED\nOUTPUT: '+str(out.resolve())+'\nZIP: '+str(archive.resolve()),flush=True)


if __name__ == '__main__':
    try:
        main()
    except (OSError,ValueError,KeyError,IndexError,csv.Error) as error:
        print('STOP: '+str(error),file=sys.stderr)
        sys.exit(1)
