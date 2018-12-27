#!/bin/bash
#SBATCH -J SL05
#SBATCH -p bdwall
#SBATCH -A NEXTGENOPT
#SBATCH -N 1
#SBATCH -t 40:00:00


julia /home/choy/QS/TVGG1PS_new_arrival_functions/main_functions/main_05.jl
