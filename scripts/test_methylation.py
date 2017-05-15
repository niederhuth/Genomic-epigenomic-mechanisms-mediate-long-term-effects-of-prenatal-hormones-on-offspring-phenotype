import os
import sys
import pybedtools as pbt
import pandas as pd
import itertools
import numpy as np
from scipy.stats.stats import pearsonr
from collections import Counter
from Bio import SeqIO

#interpret sequence context, taken from methylpy.utils
def expand_nucleotide_code(mc_type=["C"]):
    iub_dict = {"N":["A","C","G","T"],"H":["A","C","T"],"C":["C"],"G":["G"],"T":["T"],"A":["A"]}
    for type in mc_type[:]:
        type += "N" * (3 - len(type))
        mc_type.extend(["".join(i) for i in itertools.product(*[iub_dict[nuc] for nuc in type])])
    if "C" in mc_type:
        mc_type.extend(["CG", "CHG", "CHH","CNN"])
    if "CG" in mc_type:
        mc_type.extend(["CGN"])
    if "CH" in mc_type:
        mc_type.extend(["CHN"])
    if "CA" in mc_type:
        mc_type.extend(["CAN"])
    if "CT" in mc_type:
        mc_type.extend(["CTN"])
    if "CC" in mc_type:
        mc_type.extend(["CCN"])
    return mc_type

#filter allc file based on sequence context
def filter_context(allc,context=["C"]):
    a = pd.read_table(allc)
    a = a[a.mc_class.isin(expand_nucleotide_code(context))]
    return a

def allc2bed(allc,context=["C"],bed=True):
    a = filter_context(allc,context)
    a['pos2'] = a.pos
    a['name'] = a.index
    a['score'] = "."
    a = a[['chr','pos','pos2','name','score','strand','mc_class','mc_count','total','methylated']]
    if bed is True:
        a = pbt.BedTool.from_dataframe(a)
    return a

#simple function for filtering gff files based on feature (gene, exon, mRNA, etc)
def feat_filter(x,feature):
    if feature:
        return x[2] == feature
    else:
        return x

#simple function for filtering gff files based on strand
def strand_filter(x,strand):
    return x.strand == strand

#simple function for filtering gff files based on chromosome
def chr_filter(x,chr):
    return x.chrom not in chr

mC_bed = allc2bed(allc)
bed = pbt.BedTool(genes_gff).filter(feat_filter,"mRNA").filter(chr_filter,filter_chr)
cds_bed = pbt.BedTool(genes_gff).filter(feat_filter,"exon").filter(chr_filter,filter_chr)
mapping = pbt.bedtool.BedTool.intersect(mC_bed,cds_bed,wa=True)
m = pd.read_table(mapping.fn, header=None, usecols = [0,1,5,6,7,8,9])
m.columns = ['chr','pos','strand','mc_class','mc_count','total','methylated']
m = m.drop_duplicates()
m.to_csv('CDS_allc.tmp', sep='\t', index=False)
mC_bed2 = allc2bed('CDS_allc.tmp')
mapping2 = pbt.bedtool.BedTool.intersect(mC_bed2,bed,wa=True,wb=True)
m = pd.read_table(mapping2.fn, header=None)
name="none"
t = 0
m = 0
for c in m.itertuples():
    if c[18] == name:
        
    else:
        print("false")
        name = c[18]


c[6]
c[7]
c[8]
