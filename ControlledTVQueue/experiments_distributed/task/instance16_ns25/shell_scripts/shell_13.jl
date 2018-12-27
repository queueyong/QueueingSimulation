#!/bin/bash
#SBATCH -J SL16ns13
#SBATCH -p bdwall
#SBATCH -A NEXTGENOPT
#SBATCH -N 1
#SBATCH -t 1:00:00

julia /blues/gpfs/home/choy/QS/experiments_distributed/task/instance16_ns25/main_scripts/main_13.jl
