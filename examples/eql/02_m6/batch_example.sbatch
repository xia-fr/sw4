#!/bin/csh
#SBATCH -t 16:00:00
#SBATCH -N 160
#SBATCH --ntasks-per-node=40
#SBATCH --mem=200G
#SBATCH -A sasims
#SBATCH -p pbatch

#SBATCH --exclusive
#SBATCH -o slurm/eql_5hz_x50_y40_m6.out
#SBATCH -J eql

set CPUS = ` echo $SLURM_NNODES $SLURM_NTASKS_PER_NODE | awk '{ print $1 * $2 }' `
echo NODES: $SLURM_NNODES CPUS: $CPUS

set MEMperCPU = ` echo $SLURM_MEM_PER_NODE $SLURM_NTASKS_PER_NODE | awk '{ print $1 / $2 / 1000 }' `
echo MEM PER CPU: $MEMperCPU GB
echo `date`

set progs = /g/g11/xia8/sw4-3.0/optimize_quartz_mp
set scenario = earthquake_eql_5hz.sw4input
echo $scenario
srun -N $SLURM_NNODES -n $CPUS $progs/sw4 $scenario

cd ..

echo `date`
squeue -j $SLURM_JOBID
echo done successfully
