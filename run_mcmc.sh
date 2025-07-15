#!/bin/bash
#SBATCH --job-name=cs2
#SBATCH --time=4-00:00
#SBATCH --partition=cosmoobs
#SBATCH --output=./projects/cs2-project/logs/%x_%a_%A.out
#SBATCH --error=./projects/cs2-project/logs/%x_%a_%A.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=5
#SBATCH --exclusive
#SBATCH --mail-user=joao.reboucas@unesp.br
#SBATCH --mail-type=ALL

YAML=./projects/cs2-project/yamls/MCMC${SLURM_ARRAY_TASK_ID}.yaml

echo "Job started in `hostname`"

cd ~/cocoa/Cocoa
conda activate cocoa
source start_cocoa.sh

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OMP_PROC_BIND=close
export OMP_PLACES=cores

# NOTE: some libraries may use these variables to control threading
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export NUMEXPR_NUM_THREADS=${SLURM_CPUS_PER_TASK}

mpirun -n ${SLURM_NTASKS_PER_NODE} --mca btl tcp,self --bind-to core --map-by socket:PE=${OMP_NUM_THREADS} cobaya-run ${YAML} -r
