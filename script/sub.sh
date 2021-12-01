#!/bin/sh
#SBATCH -p G1Part_sce
#SBATCH -N 1
#SBATCH -n 56
#SBATCH -t 10
mpirun -np 56 MMC