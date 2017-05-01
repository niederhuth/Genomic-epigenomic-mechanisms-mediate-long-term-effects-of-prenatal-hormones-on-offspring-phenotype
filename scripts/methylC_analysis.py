import os
import sys
import pandas as pd
import pybedtools as pbt

functionsfile = '../../../scripts/functions.py'
sys.path.append(os.path.dirname(os.path.expanduser(functionsfile)))

import functions

allc="allc/"+sys.argv[1]+"_allc_total.tsv"
fasta="../../ref/bowtie2/Tguttata_v3.2.4.fa"
genome_file="../../ref/bowtie2/Tguttata_v3.2.4.genome"
genes_gff="../../ref/misc/Tguttata_v3.2.4.gff"
filter_chr=['MT','ChrL']
context=['CG','CH','CA','CT','CC']

print("Finding Total Weighted Methylation")
functions.weighted_mC(allc,output="results/total_weighted_mC.tsv",cutoff=0)

print("Getting per-site methylation levels")
functions.per_site_mC(allc,output_path='results/',context=context)

print("Analyzing Subcontext Methylation")
functions.subcontext_methylation(allc,fasta,context=context,output='results/subcontext_methylation.tsv',filter_chr=filter_chr)

if os.path.exists(genes_gff):
    print("Gene metaplot")
    functions.gene_metaplot(allc,genes_gff,genome_file,output="results/gene_metaplot.tsv",
                            ignoreStrand=False,windows=60,updown_stream=2000,cutoff=0,
                            first_feature='mRNA',second_feature='CDS',filter_chr=filter_chr)
else:
    print("No gene annotations found")
