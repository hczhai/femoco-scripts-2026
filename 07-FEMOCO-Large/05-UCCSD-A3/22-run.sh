#!/bin/bash
#SBATCH --job-name F11-AC22
#SBATCH -o out/F11-AC22-%A_%a.out
#SBATCH --nodes=1
#SBATCH --qos=gen
#SBATCH --partition=icc
#SBATCH --constraint=icelake
#SBATCH --time=7-0:00:00
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=64
#SBATCH --no-requeue
#SBATCH --mem=1000G

source ~/.bashrc

module purge
module load intel-oneapi-mkl/2024.0.0
module load python/3.11.7
module load openmpi/4.1.6
module load intel-oneapi-compilers/2024.1.0
module load disBatch

. ~/pybase/bin/activate

which mpirun
which python3

XPWD=`pwd`

echo "===================================================="
echo "        Job ID is:         $SLURM_JOBID"
echo "        Job name is:       $SLURM_JOB_NAME"
echo "        Hostname is:       "`hostname`
echo "        This dir is:       $XPWD"
echo "        CPUs per Task is:  $SLURM_CPUS_PER_TASK"
echo "        Tasks per Node is: $SLURM_TASKS_PER_NODE"
echo "        Mem per CPU is:    $SLURM_MEM_PER_CPU"
echo "        Mem per GPU is:    $SLURM_MEM_PER_GPU"
echo "        Mem per node is:   $SLURM_MEM_PER_NODE"

export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export PYTHONUNBUFFERED=1
export LD_PRELOAD=

which orterun
which python3

python -u 22-uhf-nuc-density.py |& tee 22-uhf-nuc-density.out
