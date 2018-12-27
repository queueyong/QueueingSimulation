#!/bin/bash
#SBATCH -J SL15ns6
#SBATCH -p bdwall
#SBATCH -A NEXTGENOPT
#SBATCH -N 1
#SBATCH -t 2:00:00

julia /blues/gpfs/home/choy/QS/experiments_distributed/task/instance15_ns25/main_scripts/main_6.jl
