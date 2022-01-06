#!/bin/bash

## Give the Job a descriptive name
#PBS -N run_jacobi

## Output and error files
#PBS -o run.out
#PBS -e run.err

## Limit memory, runtime etc.
#PBS -l walltime=00:60:00

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
for size in 2048 #4096 6144
do
   for i in 1 2 3
   do
      mpirun -np 1 --map-by node --mca btl self,tcp jacobi ${size} ${size} 1 1 >> time1_${size}.out
      mpirun -np 2 --map-by node --mca btl self,tcp jacobi ${size} ${size} 2 1 >> time2_${size}.out
      mpirun -np 4 --map-by node --mca btl self,tcp jacobi ${size} ${size} 2 2 >> time4_${size}.out
      mpirun -np 8 --map-by node --mca btl self,tcp jacobi ${size} ${size} 4 2 >> time8_${size}.out
      mpirun -np 16 --map-by node --mca btl self,tcp jacobi ${size} ${size} 4 4 >> time16_${size}.out
      mpirun -np 32 --map-by node --mca btl self,tcp jacobi ${size} ${size} 4 8 >> time32_${size}.out
      mpirun -np 64 --map-by node --mca btl self,tcp jacobi ${size} ${size} 8 8 >> time64_${size}.out
   done
done

rm -r ${HOME}/tmp
