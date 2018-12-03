#!/bin/bash --login
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=100GB
#SBATCH --job-name DMRfind
#SBATCH --output=job_reports/%x-%j.SLURMout

cd $PBS_O_WORKDIR

export TMPDIR=$PBS_O_WORKDIR
export TMP=$PBS_O_WORKDIR
export TEMP=$PBS_O_WORKDIR

#Define Variables

SAMPLES='FCH1 FCH2 FCH3 FCNT1 FCNT2 FCNT3 FTH1 FTH2 FTH3 FTNT1 FTNT2 FTNT3 MCH1 MCH2 MCH3 MCNT1 MCNT2 MCNT3 MTH1 MTH2 MTH3 MTNT1 MTNT2 MTNT3'

SAMPLE_CATEGORY=$(echo $(for i in $(seq 1 $(expr $(echo $SAMPLES | wc -w) / 1)); do printf "%0.s$i, " {1..1}; done) | sed s/,$//)

CHROMS='1 10 11 12 13 14 15 17 18 19 1A 1B 2 20 21 22 23 24 25 26 27 28 3 4 4A 5 6 7 8 9 LG2 LG5 LGE22 Z Un 4_random 8_random Z_random 13_random 5_random 6_random 21_random 2_random 26_random 3_random 1_random 22_random 1A_random 7_random 10_random 23_random 18_random 25_random LGE22_random 9_random 15_random 2_random 20_random 11_random 4A_random 14_random 17_random 27_random 19_random 28_random 16_random 24_random 1B_random'

ALLC_FILES=$(for i in $SAMPLES; do awk -v a=$i 'BEGIN {print "../"a"/methylCseq/allc_"a".tsv.gz"}'; done | tr '\n' ' ')

#Run DMRfind

for i in CG CH CN
do
	echo "Calling $i DMRs"
	methylpy DMRfind \
	--allc-files $ALLC_FILES \
	--samples $SAMPLES \
	--output-prefix $i \
	--chroms $CHROMS \
	--mc-type $i \
	--num-procs 20 \
	--min-cov 0 \
	--dmr-max-dist 250 \
	--sig-cutoff 0.01 \
	--num-sims 3000 \
	--min-tests 100 \
	--min-num-dms 0 \
	--mc-max-dist 10 \
	--keep-temp-files False \
	#--min-cluster 1 \
	#--seed -1 \ 
	#--sample-category $SAMPLE_CATEGORY 
	#--resid-cutoff 0.01   
done
