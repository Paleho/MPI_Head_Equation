#!/bin/bash

## Give the Job a descriptive name
#PBS -N make_script

## Output and error files
#PBS -o make.out
#PBS -e make.err

## Limit memory, runtime etc.
#PBS -l walltime=00:10:00

## How many nodes:processors_per_node should we get?
#PBS -l nodes=1:ppn=1

module load openmpi/1.8.3
cd /home/parallel/parlab09/a3/mpi/report/RedBlackSOR
make
