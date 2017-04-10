import os
import sys
import pandas as pd
import pybedtools as pbt

functionsfile = '../../../../scripts/functions.py'
sys.path.append(os.path.dirname(os.path.expanduser(functionsfile)))

import functions

dmr_file="CG_rms_results_collapsed.tsv"
filtered_dmr="CG_filtered.tsv"
genes_gff=""
repeats_gff=""

functions.filter_dmr(dmr_file,filtered_dmr,min_dms=5,min_mC_diff=0.1)
genes = pbt.BedTool(genes_gff)
repeats = pbt.BedTool(repeats_gff)
dmr = pd.read_table(filtered_dmr)
dmr2 = pbt.BedTool.from_dataframe(dmr[['chr', 'start', 'end']])
mapping = pbt.bedtool.BedTool.intersect(repeats,dmr2,wa=True,wb=True)
