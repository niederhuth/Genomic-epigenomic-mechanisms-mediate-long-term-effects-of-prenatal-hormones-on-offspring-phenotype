#!/bin/bash

echo "Starting"
#Count reads per feature
echo "Running htseq-count"
export PYTHONPATH=${PYTHONPATH}:/usr/local/htseq/0.6.1p1/lib/python
time /usr/local/htseq/0.6.1p1/bin/htseq-count -f bam -s no -t exon -i gene_id -m union tophat/accepted_hits.bam ../ref/transcriptome/Gmax_v2.gff > gene_counts.tsv
echo "Done"
