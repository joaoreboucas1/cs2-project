#!/bin/bash
#SBATCH --job-name=benchmark_atabey
#SBATCH --time=00:10:00
#SBATCH --partition=atabey-ext
#SBATCH --output=./projects/cs2-project/logs/%x_%a_%A.out
#SBATCH --error=./projects/cs2-project/logs/%x_%a_%A.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=8
#SBATCH --mail-user=joao.reboucas@unesp.br
#SBATCH --mail-type=ALL
#SBATCH --hint=nomultithread

YAML=./projects/cs2-project/yamls/TEST101.yaml

echo "Job started in `hostname` at `date`"

# NOTE: sometimes `source start_cocoa.sh` fails if many jobs are loaded simultaneously.
# This sleep prevents bugs from happening
# sleep $(( 10 + SLURM_ARRAY_TASK_ID*20 ))

cd ~/cocoa2/Cocoa
conda init bash
conda activate cocoa2
source start_cocoa.sh

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OMP_PROC_BIND=close
export OMP_PLACES=cores

# NOTE: some libraries may use these variables to control threading
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export NUMEXPR_NUM_THREADS=${SLURM_CPUS_PER_TASK}

mpirun -n ${SLURM_NTASKS_PER_NODE} --mca btl tcp,self --bind-to core --map-by socket:PE=${OMP_NUM_THREADS} cobaya-run ${YAML} -r

srun --mpi=pmix --cpu_bind=cores --distribution=block -n ${SLURM_NTASKS} cobaya-run ${YAML} -r

echo "Job ended at `date`"
