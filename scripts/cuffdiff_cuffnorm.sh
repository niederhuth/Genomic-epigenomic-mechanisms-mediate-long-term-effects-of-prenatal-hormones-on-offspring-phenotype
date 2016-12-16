#!/bin/bash

echo "Starting"
export PATH=/usr/local/samtools/0.1.18/:${PATH}
export LD_LIBRARY_PATH=/usr/local/boost/1.54.0/gcc447/lib:/usr/local/gcc/4.7.1/lib:/usr/local/gcc/4.7.1/lib64:${LD_LIBRARY_PATH}
export PATH=/usr/local/cufflinks/2.2.1/bin/:${PATH}
#cuffdiff
echo "Running cuffdiff"
/usr/local/cufflinks/2.2.1/bin/cuffdiff -p 10 -L PI408105A-Control,PI408105A-Flooding,S99-2281-Control,S99-2281-Flooding --dispersion-method pooled --library-norm-method geometric --min-reps-for-js-test 3 -b ../ref/bowtie2/Gmax_v2.fa -u --library-type fr-firststrand -o cuffdiff ../ref/transcriptome/Gmax_v2.gff ../PI408105A_C1/cuffquant/abundances.cxb,../PI408105A_C2/cuffquant/abundances.cxb,../PI408105A_C3/cuffquant/abundances.cxb ../PI408105A_F1/cuffquant/abundances.cxb,../PI408105A_F2/cuffquant/abundances.cxb,../PI408105A_F3/cuffquant/abundances.cxb ../S99-2281_C1/cuffquant/abundances.cxb,../S99-2281_C2/cuffquant/abundances.cxb,../S99-2281_C3/cuffquant/abundances.cxb ../S99-2281_F1/cuffquant/abundances.cxb,../S99-2281_F2/cuffquant/abundances.cxb,../S99-2281_F3/cuffquant/abundances.cxb
#cuffnorm
echo "Running Cuffnorm"
/usr/local/cufflinks/2.2.1/bin/cuffnorm -p 10 -L PI408105A-Control,PI408105A-Flooding,S99-2281-Control,S99-2281-Flooding --library-norm-method geometric --library-type fr-firststrand -o cuffnorm ../ref/transcriptome/Gmax_v2.gff ../PI408105A_C1/cuffquant/abundances.cxb,../PI408105A_C2/cuffquant/abundances.cxb,../PI408105A_C3/cuffquant/abundances.cxb ../PI408105A_F1/cuffquant/abundances.cxb,../PI408105A_F2/cuffquant/abundances.cxb,../PI408105A_F3/cuffquant/abundances.cxb ../S99-2281_C1/cuffquant/abundances.cxb,../S99-2281_C2/cuffquant/abundances.cxb,../S99-2281_C3/cuffquant/abundances.cxb ../S99-2281_F1/cuffquant/abundances.cxb,../S99-2281_F2/cuffquant/abundances.cxb,../S99-2281_F3/cuffquant/abundances.cxb
