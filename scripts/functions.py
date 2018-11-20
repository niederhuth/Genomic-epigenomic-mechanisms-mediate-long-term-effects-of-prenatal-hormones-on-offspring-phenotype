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
	#expand nucleotide list for a given context
	b = expand_nucleotide_code(mc_type)
	d1 = d2 = d3 = d4 = 0
	#iterate over each line
	for c in a.itertuples():
		#check if line is correct context
		if c[4] in b:
			#check if meets minimum cutoff for read depth
			if int(c[6]) >= int(cutoff):
				#add up number of sites
				d1 = d1 + 1
				#add up number of sites called methylated by methylpy
				d2 = d2 + int(c[7])
				#add up total reads covering a site
				d3 = d3 + int(c[6])
				#add up total methylated reads covering a site
				d4 = d4 + int(c[5])
	#create list
	e = [mc_type,d1,d2,d3,d4]
	#return that list
	return e

#Collect total methylation data for genome or sets of chromosomes
def total_weighted_mC(allc,output=(),mc_type=['CG','CHG','CHH'],cutoff=0,chrs=[]):
	#read allc file
	a =  pd.read_table(allc,names=['chr','pos','strand','mc_class','mc_count','total','methylated'],dtype={'chr':str,'pos':int,'strand':str,'mc_class':str,'mc_count':int,'total':int,'methylated':int})
	#filter chromosome sequences
	if chrs:
		a = a[a.chr.isin(chrs)]
	#create data frame
	columns=['Context','Total_sites','Methylated_sites','Total_reads','Methylated_reads','Weighted_mC']
	b = pd.DataFrame(columns=columns)
	#iterate over each mC type and run get_mC_data
	for c in mc_type:
		d=[]
		d = get_mC_data(a,mc_type=c,cutoff=cutoff)
		#calculate weighted methylation
		d = d + [np.float64(d[4])/np.float64(d[3])]
		b = b.append(pd.DataFrame(d,columns=columns), ignore_index=True)
	#output results
	if output:
		b.to_csv(output, sep='\t', index=False)
	else:
		return b
