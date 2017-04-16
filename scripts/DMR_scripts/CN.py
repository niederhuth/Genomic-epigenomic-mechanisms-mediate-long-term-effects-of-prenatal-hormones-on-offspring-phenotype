#!/usr/local/python/2.7.2/bin/python
from methylpy.DMRfind import DMRfind
import sys

#Samples included in the comparison
samples = ["FCH1","FCH2","FCH3","FCNT1","FCNT2","FCNT3","FTH1","FTH2","FTH3",
"FTNT1","FTNT2","FTNT3","MCH1","MCH2","MCH3","MCNT1","MCNT2","MCNT3","MTH1",
"MTH2","MTH3","MTNT1","MTNT2","MTNT3"]

#Regions included in the DMR analysis
chrom ={"1","10","11","12","13","14","15","16","17","18","19","1A","1B","2","20","21",
"22","23","24","25","26","27","28","3","4","4A","5","6","7","8","9","LG2","LG5",
"LGE22","Z","Un","4random","8random","Zrandom","13random","5random","6random",
"21random","2random","26random","3random","1random","22random","1Arandom","7random",
"10random","23random","18random","25random","LGE22random","9random","15random",
"12random","20random","11random","4Arandom","14random","17random","27random",
"19random","28random","16random","24random","1Brandom"}
region_dict={}
for chr in chrom:
    region_dict[chr]=[0,5000000000000000]

#methylation type
mc_type=['CN']

#Number of processors
num_procs=20

path_to_allc="../allc"
output_prefix = "CN"

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
