import os
import sys
import pandas as pd

functionsfile = '../../../scripts/functions.py'
sys.path.append(os.path.dirname(os.path.expanduser(functionsfile)))

import functions

#define variables
allc='allc_'+sys.argv[1]+'.tsv.gz'
genome_file='../../ref/sequences/Tguttata.genome'
filter_chr=['ChrL','MT','Z','Z_random']
mc_type=['CG','CA','CT','CC','CH']
cutoff=0
output='results/total_weighted_methylation.txt'

#get chromosome list
chrs = list(pd.read_table(genome_file,header=None,usecols=[0],dtype='str')[0])
chrs = list(set(chrs).difference(filter_chr))

#get total weighted mC
functions.total_weighted_mC(allc, output=output, mc_type=mc_type, cutoff=cutoff, chrs=chrs)

#run for Z chromosome
output='results/chrZ_weighted_methylation.txt'
chrs=['Z','Z_random']
functions.total_weighted_mC(allc, output=output, mc_type=mc_type, cutoff=cutoff, chrs=chrs)

