#!/bin/bash -l
#SBATCH --job-name=LCORL205_CF
#SBATCH --partition=gpu_cuda
#SBATCH --gres=gpu:h100:1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=256G
#SBATCH --time=24:00:00
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err

set -eo pipefail
cd "$SLURM_SUBMIT_DIR"

trap 'rc=$?; printf "exit_code=%s\n" "$rc" > job_exit_status.txt' EXIT

module load localcolabfold/1.5.5

export OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK"
export PYTHONUNBUFFERED=1

date -Is
hostname
module list 2>&1
command -v colabfold_batch

nvidia-smi > gpu_info.txt
cat gpu_info.txt
free -h

colabfold_batch --help > colabfold_batch_help.txt 2>&1

# Stop immediately if JAX cannot use the allocated GPU.
python - <<'PY'
import jax
devices = jax.devices()
print("JAX devices:", devices)
if not any(d.platform == "gpu" for d in devices):
    raise SystemExit("STOP: JAX cannot access the allocated GPU.")
PY

colabfold_batch \
    LCORL205_REFERENCE.fasta \
    predictions \
    --model-type alphafold2_ptm \
    --msa-mode mmseqs2_uniref_env \
    --max-msa 512:1024 \
    --num-models 5 \
    --num-recycle 3 \
    --num-seeds 1 \
    --random-seed 0 \
    2>&1 | tee prediction_console.log

# Check actual output files; do not rely only on process exit status.
python - <<'PY'
from pathlib import Path

out = Path("predictions")
pdbs = sorted(out.glob("*unrelaxed*.pdb"))
scores = sorted(out.glob("*scores*.json"))

print("Saved model PDBs:", len(pdbs))
print("Saved score JSONs:", len(scores))
for p in pdbs:
    print(p)

if len(pdbs) < 5 or len(scores) < 5:
    raise SystemExit(
        "INCOMPLETE: fewer than five models or score files. "
        "Keep this directory and inspect prediction_console.log."
    )
PY

date -Is
echo "DONE: five model PDBs and score files found."
