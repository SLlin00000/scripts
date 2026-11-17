#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=04:00:00
set -eo pipefail
cd /scratch/user/s4922196/LCORL206_CF155_qim53006
trap 'rc=$?; printf "exit_code=%s\n" "$rc" > /QRISdata/Q9774/lcorl_project/results/LCORL201_203_206_CF155_20261115/LCORL206/job_exit_status.txt' EXIT
module load localcolabfold/1.5.5
export OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK"
export PYTHONUNBUFFERED=1
date -Is
hostname
module list 2>&1
nvidia-smi > /QRISdata/Q9774/lcorl_project/results/LCORL201_203_206_CF155_20261115/LCORL206/gpu_info.txt
python - <<'JAX_CHECK'
import jax
devices = jax.devices()
print("JAX devices:", devices, flush=True)
if not any(d.platform == "gpu" for d in devices):
    raise SystemExit("STOP: JAX cannot access the allocated GPU")
JAX_CHECK
colabfold_batch /QRISdata/Q9774/lcorl_project/results/LCORL201_203_206_CF155_20261115/LCORL206/LCORL206_REFERENCE.fasta /QRISdata/Q9774/lcorl_project/results/LCORL201_203_206_CF155_20261115/LCORL206/predictions --model-type alphafold2_ptm --msa-mode mmseqs2_uniref_env --max-msa 512:1024 --num-models 5 --num-recycle 3 --num-seeds 1 --random-seed 0 --rank plddt 2>&1 | tee /QRISdata/Q9774/lcorl_project/results/LCORL201_203_206_CF155_20261115/LCORL206/prediction_console.log
python /QRISdata/Q9774/lcorl_project/results/LCORL201_203_206_CF155_20261115/check_and_pack.py /QRISdata/Q9774/lcorl_project/results/LCORL201_203_206_CF155_20261115/LCORL206 1904
date -Is
