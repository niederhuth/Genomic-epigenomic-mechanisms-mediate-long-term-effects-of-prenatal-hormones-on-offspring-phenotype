#!/bin/bash --login
#SBATCH --time=3:59:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=25GB
#SBATCH --job-name reidentify-DMR
#SBATCH --output=job_reports/%x-%j.SLURMout

cd $PBS_O_WORKDIR

#Define variables

SAMPLES='FCH1 FCH2 FCH3 FCNT1 FCNT2 FCNT3 FTH1 FTH2 FTH3 FTNT1 FTNT2 FTNT3 MCH1 MCH2 MCH3 MCNT1 MCNT2 MCNT3 MTH1 MTH2 MTH3 MTNT1 MTNT2 MTNT3'

#SAMPLES="MCH1 MCH2 MCH3 MTH1 MTH2 MTH3"

SAMPLE_CATEGORY=$(echo $(for i in $(seq 1 $(expr $(echo $SAMPLES | wc -w) / 3)); do printf "%0.s$i, " {1..3}; done) | sed s/,$//)

echo $SAMPLE_CATEGORY

for i in CG
do
	methylpy reidentify-DMR \
		--input-rms-file "$i"_rms_results.tsv.gz \
		--output-file "$i"_recalled \
		--collapse-samples False \
		--sample-category $SAMPLE_CATEGORY \
		--min-cluster 2 \
		--sig-cutoff 0.01 \
		--dmr-max-dist 250 \
		--min-num-dms 0 \
		--num-sims 3000 \
		--min-tests 100 
		#--resid-cutoff 0.01
done

