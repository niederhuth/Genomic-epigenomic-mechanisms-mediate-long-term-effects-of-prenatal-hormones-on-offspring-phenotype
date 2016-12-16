#!/bin/bash

export PATH=/usr/local/samtools/0.1.18/:${PATH}
export LD_LIBRARY_PATH=/usr/local/boost/1.54.0/gcc447/lib:/usr/local/gcc/4.7.1/lib:/usr/local/gcc/4.7.1/lib64:${LD_LIBRARY_PATH}
export PATH=/usr/local/cufflinks/2.2.1/bin/:${PATH}
/usr/local/cufflinks/latest/bin/cuffquant -p 10 -b ../ref/bowtie2/Gmax_v2.fa -u --library-type fr-firststrand ../ref/transcriptome/Gmax_v2.gff tophat/accepted_hits.bam -o cuffquant
