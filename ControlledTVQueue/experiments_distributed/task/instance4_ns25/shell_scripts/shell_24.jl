#!/bin/bash
#SBATCH -J SL4ns24
#SBATCH -p bdwall
#SBATCH -A NEXTGENOPT
#SBATCH -N 1
#SBATCH -t 1:00:00

julia /blues/gpfs/home/choy/QS/experiments_distributed/task/instance4_ns25/main_scripts/main_24.jl