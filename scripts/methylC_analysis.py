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
diff_exp="results/diff_exp.gff"
no_diff_exp="results/no_diff_exp.gff"
goi="../../ref/misc/gof.gff"
filter_chr=['MT','ChrL']
context=['CG','CH','CA','CT','CC']

print("Finding Total Weighted Methylation")
functions.weighted_mC(allc,output="results/total_weighted_mC.tsv",cutoff=0)

if os.path.exists(genes_gff):
    print("Gene metaplot")
    functions.gene_metaplot(allc,genes_gff,genome_file,output="results/gene_metaplot.tsv",
                            ignoreStrand=False,windows=60,updown_stream=2000,cutoff=0,
                            first_feature='mRNA',second_feature='exon',filter_chr=filter_chr,
                            remove_tmp=False)
else:
    print("No gene annotations found")

if os.path.exists(diff_exp):
    print("Differentially expressed gene metaplot")
    functions.feature_metaplot('CDS_allc.tmp',diff_exp,genome_file,output="results/diff_exp_metaplot.tsv",
                            ignoreStrand=False,windows=60,updown_stream=2000,cutoff=0,
                            filter_features='gene',filter_chr=filter_chr)
else:
    print("No gene annotations found")

if os.path.exists(no_diff_exp):
    print("Non-differentially expressed gene metaplot")
    functions.feature_metaplot('CDS_allc.tmp',no_diff_exp,genome_file,output="results/no_diff_exp_metaplot.tsv",
                            ignoreStrand=False,windows=60,updown_stream=2000,cutoff=0,
                            filter_features='gene',filter_chr=filter_chr)
else:
    print("No gene annotations found")

if os.path.exists(goi):
    print("Gene of interest metaplot")
    functions.feature_metaplot('CDS_allc.tmp',goi,genome_file,output="results/goi_metaplot.tsv",
                            ignoreStrand=False,windows=60,updown_stream=2000,cutoff=0,
                            filter_features='gene',filter_chr=filter_chr)
else:
    print("No gene annotations found")

if os.path.exists(genes_gff):
    print("Gene methylation levels")
    functions.feature_mC_levels('CDS_allc.tmp',genes_gff,output="results/gene_methylation_levels.tsv",
                                cutoff=0,filter_features='gene',filter_chr=filter_chr)
else:
    print("No gene annotations found")

#os.remove('CDS_allc.tmp')
