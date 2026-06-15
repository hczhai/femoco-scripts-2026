#!/bin/bash
#SBATCH --job-name M8-35
#SBATCH -o M8-35-%A_%a.out
#SBATCH --nodes=8
#SBATCH --qos=long
#SBATCH --time=14-0:00:00
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=56
#SBATCH --array=26-26
#SBATCH --partition=expansion
#SBATCH --reservation=changroup_himem
#SBATCH --no-requeue
#SBATCH --mem=1500G

source ~/.bashrc

module list
module purge
module load vim/9.0.0045-gcc-13.2.0-3naifpx
module load gcc/13.2.0-gcc-13.2.0-w55nxkl
module load openmpi/4.1.5-gcc-13.2.0-24q3ap2
module load intel-oneapi-mkl/2023.2.0-gcc-13.2.0-ohvyk7g

export PYTHONPATH=/home/hczhai/changroup/program/block2-new:/home/hczhai/changroup/program/block2-new/build:$PYTHONPATH

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
echo "        Tasks per Node is: $SLURM_MEM_PER_CPU"
echo "        Tasks per Node is: $SLURM_MEM_PER_GPU"
echo "        Tasks per Node is: $SLURM_MEM_PER_NODE"

export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export PYTHONUNBUFFERED=1
export LD_PRELOAD=

which orterun
which python3

IX=${SLURM_ARRAY_TASK_ID}

(orterun --map-by ppr:$SLURM_TASKS_PER_NODE:node:pe=$OMP_NUM_THREADS python3 -u 35-uhf-dmrg-12k-res-re.py ${IX}) |& tee 35-uhf-dmrg-12k-re-${IX}.out
