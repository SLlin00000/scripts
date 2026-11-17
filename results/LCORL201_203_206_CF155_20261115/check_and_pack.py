import json
from pathlib import Path
import sys
import zipfile

folder = Path(sys.argv[1]).resolve()
expected = int(sys.argv[2])
out = folder / 'predictions'
pdbs = sorted(out.glob('*unrelaxed_rank_*.pdb'))
scores = sorted(out.glob('*scores_rank_*.json'))
if len(pdbs) != 5 or len(scores) != 5:
    raise SystemExit('INCOMPLETE: expected five PDBs and five score JSONs')
for p in pdbs + scores:
    if p.stat().st_size == 0:
        raise SystemExit('EMPTY: ' + str(p))
for p in pdbs:
    ca = [s for s in p.read_text().splitlines()
          if s.startswith('ATOM  ') and s[12:16].strip() == 'CA']
    if len(ca) != expected:
        raise SystemExit('PDB length mismatch: ' + str(p))
for p in scores:
    data = json.loads(p.read_text())
    if len(data.get('plddt', [])) != expected:
        raise SystemExit('pLDDT length mismatch: ' + str(p))
(folder / 'output_check.txt').write_text('PASS: five nonempty PDBs and score files; lengths match input.\n')
archive = folder.parent / (folder.name + '_confidence_review.zip')
temporary = archive.with_suffix('.zip.partial')
with zipfile.ZipFile(temporary, 'w', compression=zipfile.ZIP_DEFLATED) as z:
    for name in ('reference_proteins.fasta', 'variant_consequences.long.tsv', 'check_and_pack.py'):
        p = folder.parent / name
        z.write(p, name)
    for p in sorted(folder.iterdir()):
        if p.is_file() and not p.name.startswith('slurm-'):
            z.write(p, str(p.relative_to(folder.parent)))
    for p in sorted(out.iterdir()):
        if p.is_file() and p.suffix in {'.pdb', '.json', '.png', '.a3m', '.txt', '.bibtex'}:
            z.write(p, str(p.relative_to(folder.parent)))
temporary.replace(archive)
print('Output check: PASS', flush=True)
print('Download archive:', archive, flush=True)
