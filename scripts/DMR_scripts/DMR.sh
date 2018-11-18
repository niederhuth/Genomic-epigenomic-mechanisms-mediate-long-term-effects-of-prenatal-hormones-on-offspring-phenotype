#!/bin/bash -login
#PBS -l walltime=60:00:00 
#PBS -l nodes=1:ppn=6
#PBS -l mem=100gb
#PBS -N Find_DMRs

cd $PBS_O_WORKDIR
module load GSL/1.15

SAMPLES="FCH1 FCH2 FCH3 FCNT1 FCNT2 FCNT3 FTH1 FTH2 FTH3 FTNT1 FTNT2 FTNT3 MCH1 MCH2 MCH3 MCNT1 MCNT2 MCNT3 MTH1 MTH2 MTH3 MTNT1 MTNT2 MTNT3"
SAMPLE_CATEGORY=$(seq 1 $(echo $SAMPLES | wc -w) | tr '\n' ',' | sed s/,$//i | sed s/,/,\ /g)
CHROMS="1 10 11 12 13 14 15 17 18 19 1A 1B 2 20 21 22 23 24 25 26 27 28 3 4 4A 5 6 7 8 9 LG2 LG5 LGE22 Z Un 4_random 8_random Z_random 13_random 5_random 6_random 21_random 2_random 26_random 3_random 1_random 22_random 1A_random 7_random 10_random 23_random 18_random 25_random LGE22_random 9_random 15_random 2_random 20_random 11_random 4A_random 14_random 17_random 27_random 19_random 28_random 16_random 24_random 1B_random"
ALLC_FILES=$(for i in $SAMPLES; do awk -v a=$i 'BEGIN {print "../"a"/methylCseq/allc_"a".tsv.gz"}'; done | tr '\n' ' ')
echo $ALLC_FILES

for i in CNN
do
	echo "Calling $i DMRs"
	methylpy DMRfind --allc-files $ALLC_FILES \
        	         --samples $SAMPLES \
			 --output-prefix $i \
			 --chroms $CHROMS \
                	 --mc-type $i \
                	 --num-procs 6 \
			 --min-cov 1 \
                	 --dmr-max-dist 250 \
                	 --sig-cutoff 0.01 \
			 --num-sims 3000 \
                	 --min-tests 100 \
			 --min-num-dms 0 \
                	 --mc-max-dist 0 \
                	 --keep-temp-files False \
                	 --min-cluster 2 \
			 --seed -1 \
			 --sample-category $SAMPLE_CATEGORY 
			 #--resid-cutoff 0.01 #This is default value, do not change. This is only here to show what the value used in analyses is   
done
