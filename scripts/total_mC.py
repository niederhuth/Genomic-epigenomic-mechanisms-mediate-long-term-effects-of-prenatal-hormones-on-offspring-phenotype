import os
import sys

functionsfile = '../../../scripts/functions.py'
sys.path.append(os.path.dirname(os.path.expanduser(functionsfile)))

import functions
allc='allc_'+sys.argv[1]+'.tsv.gz'
#genome_file='../../ref/Tguttata.genome'
filter_chr=['ChrL','MT']
context=['CG','CA','CT','CC','CH']

functions.weighted_mC(allc, output="results/Total_methylation.tsv", cutoff=0)

