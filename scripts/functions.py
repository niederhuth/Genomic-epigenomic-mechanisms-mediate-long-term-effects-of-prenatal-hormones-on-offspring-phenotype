import pandas as pd
import itertools
import numpy as np

#interpret sequence context, taken from methylpy utilities
def expand_nucleotide_code(mc_type=["C"]):
	iub_dict = {"N":["A","C","G","T","N"],
	"H":["A","C","T","H"],
	"D":["A","G","T","D"],
	"B":["C","G","T","B"],
	"A":["A","C","G","A"],
	"R":["A","G","R"],
	"Y":["C","T","Y"],
	"K":["G","T","K"],
	"M":["A","C","M"],
	"S":["G","C","S"],
	"W":["A","T","W"],
	"C":["C"],
	"G":["G"],
	"T":["T"],
	"A":["A"]}

	mc_class = list(mc_type) # copy
	if "C" in mc_type:
		mc_class.extend(["CGN", "CHG", "CHH","CNN"])
	elif "CG" in mc_type:
		mc_class.extend(["CGN"])

	mc_class_final = []
	for motif in mc_class:
		mc_class_final.extend(["".join(i) for i in itertools.product(*[iub_dict[nuc] for nuc in motif])])
	return(set(mc_class_final))

#Collect mC data for a context
def get_mC_data(a,mc_type=['C'],cutoff=0):
	b = expand_nucleotide_code(mc_type)
	d1 = d2 = d3 = d4 = 0
	for c in a.itertuples():
		if c[4] in b:
			if int(c[6]) >= int(cutoff):
				d1 = d1 + 1
				d2 = d2 + int(c[7])
				d3 = d3 + int(c[6])
				d4 = d4 + int(c[5])
	e = [mc_type,d1,d2,d3,d4]
	return e

#Collect total methylation data for genome or sets of chromosomes
def total_weighted_mC(allc,output=(),mc_type=['CG','CHG','CHH'],cutoff=0,chrs=[]):
	a =  pd.read_table(allc,names=['chr','pos','strand','mc_class','mc_count','total','methylated'],dtype={'chr':str,'pos':int,'strand':str,'mc_class':str,'mc_count':int,'total':int,'methylated':int})
	if chrs:
		a = a[a.chr.isin(chrs)]
	b = pd.DataFrame(columns=['Context','Total_sites','Methylated_sites','Total_reads','Methylated_reads','Weighted_mC'])
	for c in mc_type:
		d = get_mC_data(a,mc_type=c,cutoff=cutoff)
		b = b.append({'Context':print(d[0]),'Total_sites':d[1],'Methylated_sites':d[2],'Total_reads':d[3],'Methylated_reads':d[4],'Weighted_mC':(np.float64(d[4])/np.float64(d[3]))}, ignore_index=True)
	if output:
		b.to_csv(output, sep='\t', index=False)
	else:
		return b
