import os
import sys
import pandas as pd

functionsfile = '../../../scripts/functions.py'
sys.path.append(os.path.dirname(os.path.expanduser(functionsfile)))

import functions

#define variables
allc='annotation_filtered_allc.tmp'
annotations='../../ref/annotations/Tguttata.gff'
genome_file='../../ref/sequences/Tguttata.genome'
filter_chr=['ChrL','MT']
mc_type=['CG','CA','CT','CC','CH']
updown_stream=0
cutoff=0
feature='gene'
output='results/gene_exon_methylation.txt'

#get chromosome list
chrs = list(pd.read_table(genome_file,header=None,usecols=[0],dtype='str')[0])
chrs = list(set(chrs).difference(filter_chr))

#get gene exon methylation
print('getting gene exon methylation levels')
functions.gene_methylation(allc,annotations,genome_file,output=output,mc_type=mc_type,updown_stream=updown_stream,feature=feature,cutoff=cutoff,chrs=chrs)

