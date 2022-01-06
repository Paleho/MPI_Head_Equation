#!/bin/bash

## Give the Job a descriptive name
#PBS -N run_jacobi

## Output and error files
#PBS -o run.out
#PBS -e run.err

## Limit memory, runtime etc.
#PBS -l walltime=02:30:00

## How many nodes:processors_per_node should we get?
#PBS -l nodes=8:ppn=8

## Start
##echo "PBS_NODEFILE = $PBS_NODEFILE"
##cat $PBS_NODEFILE

## Run the job (use full paths to make sure we execute the correct thing)
## NOTE: Fix the path to show to your executable!

mkdir -p ${HOME}/tmp
export TMPDIR=${HOME}/tmp

cd /home/parallel/parlab09/a3/mpi/Jacobi/

module load openmpi/1.8.3

for i in 1 2 3
do
  mpirun -np 64 --map-by node --mca btl self,tcp jacobi_conv 1024 1024 8 8 >> conv.out
done

rm -r ${HOME}/tmp
