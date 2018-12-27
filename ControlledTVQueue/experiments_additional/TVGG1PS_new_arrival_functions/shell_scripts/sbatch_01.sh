#!/bin/bash
#SBATCH -J SL01
#SBATCH -p knlall
#SBATCH -A NEXTGENOPT
#SBATCH -N 1
#SBATCH -t 20:00:00

julia /home/choy/QS/TVGG1PS_new_arrival_functions/main_functions/main_01.jl
