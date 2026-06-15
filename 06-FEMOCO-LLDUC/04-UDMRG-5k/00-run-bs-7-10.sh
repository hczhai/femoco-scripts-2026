#!/bin/bash
#SBATCH --job-name FM-02-00
#SBATCH -o out/FM-02-00-%A_%a.out
#SBATCH --nodes=2
#SBATCH --qos=preempt
#SBATCH --partition=preempt
#SBATCH --constraint=genoa
#SBATCH --time=7-0:00:00
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=96
#SBATCH --no-requeue
#SBATCH --array=2-4,6-11,16-21,25-27
#SBATCH --mem=1500G

source ~/.bashrc

module purge
module load intel-oneapi-mkl/2024.0.0
module load python/3.11.7
module load openmpi/4.1.6
module load intel-oneapi-compilers/2024.1.0

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

IX=${SLURM_ARRAY_TASK_ID}

(orterun --map-by ppr:$SLURM_TASKS_PER_NODE:node:pe=$OMP_NUM_THREADS python3 -u 00-fm-dmrg-uccsd.py ${IX}) |& tee 00-fm-dmrg-uccsd-${IX}.out
