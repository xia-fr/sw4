#!/bin/csh
#SBATCH -t 00:30:00
#SBATCH -N 10
#SBATCH --ntasks-per-node=40
#SBATCH -A <ACCOUNT HERE>
#SBATCH -p pdebug

#SBATCH -o example_one_eql.out
#SBATCH -J eql

set CPUS = ` echo $SLURM_NNODES $SLURM_NTASKS_PER_NODE | awk '{ print $1 * $2 }' `
echo NODES: $SLURM_NNODES CPUS: $CPUS

set MEMperCPU = ` echo $SLURM_MEM_PER_NODE $SLURM_NTASKS_PER_NODE | awk '{ print $1 / $2 / 1000 }' `
echo MEM PER CPU: $MEMperCPU GB
echo `date`

#set progs = /usr/workspace/wave/ramesh/SW4/TOSS4_QUARTZ/sw4_fftw_bugfix/optimize_quartz_mp
set progs = /usr/workspace/wave/sw4-eql/optimize_quartz_mp
set scenario = example_one_eql.sw4input
echo $scenario
srun -N $SLURM_NNODES -n $CPUS $progs/sw4 $scenario

cd ..

echo `date`
squeue -j $SLURM_JOBID
echo done successfully