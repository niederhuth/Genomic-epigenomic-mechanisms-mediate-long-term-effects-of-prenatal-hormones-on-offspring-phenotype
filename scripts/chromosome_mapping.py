import os
import sys
import pandas as pd

functionsfile = '../../../scripts/functions.py'
sys.path.append(os.path.dirname(os.path.expanduser(functionsfile)))

import functions

#define variables
allc='allc_'+sys.argv[1]+'.tsv.gz'
genome_file='../../ref/sequences/Tguttata.genome'
filter_chr=['ChrL','MT']
mc_type=['CG','CA','CT','CC','CH']
window_size=100000
stepsize=50000
cutoff=0
output='results/chromosome_mapping.txt'

#get chromosome list
chrs = list(pd.read_table(genome_file,header=None,usecols=[0],dtype='str')[0])
chrs = list(set(chrs).difference(filter_chr))

#get gene metaplot data
functions.genome_window_methylation(allc,genome_file,output=output,mc_type=mc_type,window_size=window_size,stepsize=stepsize,cutoff=cutoff,chrs=chrs)

