#!/bin/bash                                                                                          
#SBATCH --partition=sixhour
#SBATCH --ntasks=1
#SBATCH --mem=2gb
#SBATCH --time=02:00:00
#SBATCH --array=1-200

R=$RANDOM
./precip.x $R
echo $R >> ../scratch/Jup_Sul/Output/Trials.dat
exit 0
