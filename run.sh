#MSUB -q sixhour
#MSUB -l nodes=1:ppn=1,walltime=06:00:00
#MSUB -t 1-200

R=$RANDOM
./precip.x $R
echo $R >> ./Output/Trials.dat
exit 0
