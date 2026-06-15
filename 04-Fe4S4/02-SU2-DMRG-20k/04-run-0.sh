#!/bin/bash
#SBATCH --job-name F44-04
#SBATCH -o out/F44-04-%A_%a.out
#SBATCH --nodes=12
#SBATCH --qos=gen
#SBATCH --partition=icc
#SBATCH --constraint=genoa
#SBATCH --time=7-0:00:00
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=96
#SBATCH --no-requeue
#SBATCH --array=0-0
#SBATCH --exclude=worker[7022,7061,7085-7087,7096,7103,7114,7125,7139,7156,7176,7183,7203,7223,7227,7255,7276,7338,7360,7368,7375,7334,7329,7206,7213,7121,7105,7253,7298,7299,7030,7389,7095]
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

export OMPI_MCA_mpi_thread_multiple=1
export OMPI_MCA_mpi_async_progress=1
export OMPI_MCA_coll_tuned_use_dynamic_rules=true
export OMPI_MCA_coll_tuned_allgatherv_algorithm=3

which orterun
which python3

IX=${SLURM_ARRAY_TASK_ID}

(orterun --mca orte_keepalive_timeout 600 \
     --mca orte_base_help_aggregate 0 --map-by ppr:$SLURM_TASKS_PER_NODE:node:pe=$OMP_NUM_THREADS \
     numactl --interleave=all python3 -u 04-fe4-su2-dmrg-re.py ${IX}) |& tee 04-fe4-su2-dmrg-re-${IX}.out
