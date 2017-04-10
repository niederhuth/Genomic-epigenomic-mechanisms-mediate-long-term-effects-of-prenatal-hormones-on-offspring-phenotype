#!/usr/local/python/2.7.2/bin/python
from methylpy.DMRfind import DMRfind
import sys

#Samples included in the comparison
samples = ["FCH1","FCH2","FCH3","FCNT1","FCNT2","FCNT3","FTH1","FTH2","FTH3",
"FTNT1","FTNT2","FTNT3","MCH1","MCH2","MCH3","MCNT1","MCNT2","MCNT3","MTH1",
"MTH2","MTH3","MTNT1","MTNT2","MTNT3"]

#Regions included in the DMR analysis
chrom ={"28"}
region_dict={}
for chr in chrom:
    region_dict[chr]=[0,5000000000000000]

#methylation type
mc_type=['CG']

#Number of processors
num_procs=20

path_to_allc="../allc"
output_prefix = "CG"

DMRfind(
mc_type = mc_type,
region_dict = region_dict,
samples = samples,
path_to_allc = path_to_allc,
num_procs = num_procs,
save_result=output_prefix,
min_cov=0,
keep_temp_files=False,
mc_max_dist=0,
dmr_max_dist=1000,
resid_cutoff=.01,
sig_cutoff=.01,
num_sims=3000,
num_sig_tests=100,
seed=-1,
min_num_dms=0,
collapse_samples=False,
sample_category=False,
min_cluster=0,
use_mc_status=False,
max_iterations=1000,
convergence_diff=1)

#End of script
