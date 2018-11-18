SAMPLES="FCH1 FCH2 FCH3 FCNT1 FCNT2 FCNT3 FTH1 FTH2 FTH3 FTNT1 FTNT2 FTNT3 MCH1 MCH2 MCH3 MCNT1 MCNT2 MCNT3 MTH1 MTH2 MTH3 MTNT1 MTNT2 MTNT3"
SAMPLE_CATEGORY=$(echo $(for i in $(seq 1 $(expr $(echo $SAMPLES | wc -w) / 3)); do printf "%0.s$i, " {1..1}; done) | sed s/,$//)

methylpy reidentify-DMR \
--input-rms-file CGN_rms_results.tsv.gz \
--output-file test2 \
--collapse-samples True \
--sample-category $SAMPLE_CATEGORY \
--min-cluster 2 \
--sig-cutoff 0.01 \
--dmr-max-dist 1000 \
--min-num-dms 0 \
--num-sims 3000 \
--min-tests 100

