import os
import sys
import pandas as pd
import pybedtools as pbt

functionsfile = '../../../scripts/functions.py'
sys.path.append(os.path.dirname(os.path.expanduser(functionsfile)))

import functions

allc="allc/"+sys.argv[1]+"_allc_total.tsv"
fasta="../../ref/bowtie2/Tguttata_v3.2.4.fa"
genes_gff="../../ref/misc/goi_and_diff.gff"
genome_file="../../ref/bowtie2/Tguttata_v3.2.4.genome"
filter_chr=['MT','ChrL']
context=['CG','CH','CA','CT','CC']

if os.path.exists(genes_gff):
    print("subsetting data")
    functions.map2features(allc,genes_gff,genome_file,updown_stream=5000,
                           first_feature='gene',second_feature='gene',filter_chr=filter_chr)
else:
    print("No gene annotations found")

for i in ['c_tmp','f_tmp']:
    os.remove(i)
