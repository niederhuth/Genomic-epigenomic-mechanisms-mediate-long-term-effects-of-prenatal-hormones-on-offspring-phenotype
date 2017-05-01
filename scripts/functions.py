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

#generic function for calculating methylation levels in windows
def window_methylation_levels(m,cutoff=0,filter=0.5,nuc_bed=()):
    a = pd.DataFrame(columns=['window','mCG','mCH'])
    name = "none"
    CG = mCG = CH = mCH = 0
    if nuc_bed:
        nuc = pd.read_table(nuc_bed.fn, usecols = [3,7,8])
        m = pd.merge(m,nuc,left_on=13,right_on='4_usercol')
    for c in m.itertuples():
      if name == "none":
        name = c[5]
        if nuc_bed:
            GC = int(c[7]) + int(c[8])
        if int(c[3]) >= int(cutoff):
          if c[1].startswith("CN"):
            continue
          elif c[1].startswith("CG"):
            CG = CG + int(c[3])
            mCG = mCG + int(c[2])
          else:
            CH = CH + int(c[3])
            mCH = mHH + int(c[2])
      elif c[5] != name:
        if nuc_bed:
            if ((CG + CH)/GC) >= filter:
                a = a.append({'window':str(name), 'mCG':(np.float64(mCG)/np.float64(CG)), 'mCH':(np.float64(mCH)/np.float64(CH)), ignore_index=True)
        else:
            a = a.append({'window':str(name), 'mCG':(np.float64(mCG)/np.float64(CG)), 'mCH':(np.float64(mCH)/np.float64(CH)), ignore_index=True)
        name = c[5]
        if nuc_bed:
            GC = int(c[7]) + int(c[8])
        CG = mCG = CH = mCH = 0
        if int(c[3]) >= int(cutoff):
            if c[1].startswith("CN"):
              continue
            elif c[1].startswith("CG"):
              CG = CG + int(c[3])
              mCG = mCG + int(c[2])
            else:
              CH = CH + int(c[3])
              mCH = mCH + int(c[2])
      elif c[5] == name:
        if int(c[3]) >= int(cutoff):
          if c[1].startswith("CN"):
            continue
          elif c[1].startswith("CG"):
            CG = CG + int(c[3])
            mCG = mCG + int(c[2])
          else:
            CH = CH + int(c[3])
            mCH = mCH + int(c[2])
    if nuc_bed:
        if ((CG + CH)/GC) >= filter:
            a = a.append({'window':str(name), 'mCG':(np.float64(mCG)/np.float64(CG)), 'mCH':(np.float64(mCH)/np.float64(CH)), ignore_index=True)
    else:
        a = a.append({'window':str(name), 'mCG':(np.float64(mCG)/np.float64(CG)), 'mCH':(np.float64(mCH)/np.float64(CH)), ignore_index=True)
    return a

#get total weighted methylation
def weighted_mC(allc, output=(), cutoff=0):
    CG = mCG = CH = mCH = CN = mCN = 0
    with open(allc) as f:
        next(f)
        for l in f:
            c = l.split('\t')
            if int(c[5]) >= int(cutoff):
                if c[3].startswith("CN"):
                    CN = CN + int(c[5])
                    mCN = mCN + int(c[4])
                elif c[3].startsswith("CG"):
                    CG = CG + int(c[5])
                    mCG = mCG + int(c[4])
                else:
                    CH = CH + int(c[5])
                    mCH = mCH + int(c[4])
    a = pd.DataFrame(columns=['Context','Total','Methylated','Weighted_mC'])
    a = a.append({'Context': 'mCG', 'Total': CG, 'Methylated': mCG, 'Weighted_mC': (np.float64(mCG)/np.float64(CG))}, ignore_index=True)
    a = a.append({'Context': 'mCH', 'Total': CH, 'Methylated': mCH, 'Weighted_mC': (np.float64(mCH)/np.float64(CH))}, ignore_index=True)
    a = a.append({'Context': 'mCN', 'Total': CN, 'Methylated': mCN, 'Weighted_mC': (np.float64(mCN)/np.float64(CN))}, ignore_index=True)
    if output:
        a.to_csv(output, sep='\t', index=False)
    else:
        return a

#plot methylation level for features
def feature_metaplot(allc,features,genome_file,output=(),ignoreStrand=False,windows=60,updown_stream=2000,cutoff=0,filter_features=(),filter_chr=[]):
    counter = 1
    mC_bed = allc2bed(allc)
    a = pbt.BedTool(features)
    if ignoreStrand:
        p_bed = a.filter(feat_filter,filter_features).filter(chr_filter,filter_chr).saveas('p_tmp')
    else:
        p_bed = a.filter(strand_filter,strand='+').filter(feat_filter,filter_features).filter(chr_filter,filter_chr).saveas('p_tmp')
        n_bed = a.filter(strand_filter,strand='-').filter(feat_filter,filter_features).filter(chr_filter,filter_chr).saveas('n_tmp')
    CG = mCG = CHG = mCHG = CHH = mCHH = 0
    metaplot = pd.DataFrame(columns=['Bin','mCG','mCH'])
    for y in ['u','f','d']:
        if y == 'u':
            if ignoreStrand:
                pf_bed = pbt.bedtool.BedTool.flank(p_bed,g=genome_file,l=updown_stream,r=0,s=True)
            else:
                pf_bed = pbt.bedtool.BedTool.flank(p_bed,g=genome_file,l=updown_stream,r=0,s=True)
                nf_bed = pbt.bedtool.BedTool.flank(n_bed,g=genome_file,l=updown_stream,r=0,s=True)
                pw_bed = pbt.bedtool.BedTool.window_maker(pf_bed,b=pf_bed,n=int(windows/3),i='srcwinnum')
                nw_bed = pbt.bedtool.BedTool.window_maker(nf_bed,b=nf_bed,n=int(windows/3),i='srcwinnum',reverse=True)
        elif y == 'f':
            if ignoreStrand:
                w_bed = pbt.bedtool.BedTool.window_maker(p_bed,b=p_bed,n=int(windows/3),i='srcwinnum')
            else:
                pw_bed = pbt.bedtool.BedTool.window_maker(p_bed,b=p_bed,n=int(windows/3),i='srcwinnum')
                nw_bed = pbt.bedtool.BedTool.window_maker(n_bed,b=n_bed,n=int(windows/3),i='srcwinnum',reverse=True)
        elif y == 'd':
            if ignoreStrand:
                pf_bed = pbt.bedtool.BedTool.flank(p_bed,g=genome_file,l=0,r=updown_stream,s=True)
            else:
                pf_bed = pbt.bedtool.BedTool.flank(p_bed,g=genome_file,l=0,r=updown_stream,s=True)
                nf_bed = pbt.bedtool.BedTool.flank(n_bed,g=genome_file,l=0,r=updown_stream,s=True)
                pw_bed = pbt.bedtool.BedTool.window_maker(pf_bed,b=pf_bed,n=int(windows/3),i='srcwinnum')
                nw_bed = pbt.bedtool.BedTool.window_maker(nf_bed,b=nf_bed,n=int(windows/3),i='srcwinnum',reverse=True)
        if ignoreStrand:
            w_bed = pbt.bedtool.BedTool.window_maker(pf_bed,b=pf_bed,n=int(windows/3),i='srcwinnum')
        else:
            w_bed = pw_bed.cat(nw_bed, postmerge=False)
        mapping = pbt.bedtool.BedTool.intersect(mC_bed,w_bed,wa=True,wb=True)
        m = pd.read_table(mapping.fn, header=None, usecols = [13,6,7,8])
        for x in list(range(1,int(windows/3)+1)):
            for c in m.itertuples():
                if c[4].endswith("_"+str(x)):
                    if int(c[3]) >= int(cutoff):
                        if c[1].startswith("CN"):
                            continue
                        elif c[1].startswith("CG"):
                            CG = CG + int(c[3])
                            mCG = mCG + int(c[2])
                        else:
                            CH = CH + int(c[3])
                            mCH = mCH + int(c[2])
            metaplot = metaplot.append({'Bin': counter,'mCG': (np.float64(mCG)/np.float64(CG)),
                                        'mCH': (np.float64(mCH)/np.float64(CH)), ignore_index=True)
            counter = counter + 1
            CG = mCG = CH = mCH = 0
    os.remove('p_tmp')
    os.remove('n_tmp')
    if output:
        metaplot.to_csv(output, sep='\t', index=False)
    else:
        return metaplot

#plot methylation levels for genes
def gene_metaplot(allc,features,genome_file,output=(),ignoreStrand=False,windows=60,updown_stream=2000,cutoff=0,first_feature=(),second_feature=(),filter_chr=[]):
    mC_bed = allc2bed(allc)
    bed = pbt.BedTool(unique_genes).filter(feat_filter,first_feature).filter(chr_filter,filter_chr)
    flank_bed = pbt.bedtool.BedTool.flank(bed,g=genome_file,l=2000,r=2000,s=True).saveas('f_tmp')
    cds_bed = bed.filter(feat_filter,second_feature).filter(chr_filter,filter_chr).saveas('c_tmp')
    bed = cds_bed.cat(flank_bed, postmerge=False)
    mapping = pbt.bedtool.BedTool.intersect(mC_bed,bed,wa=True)
    m = pd.read_table(mapping.fn, header=None, usecols = [0,1,5,6,7,8,9])
    m.columns = ['chr','pos','strand','mc_class','mc_count','total','methylated']
    m = m.drop_duplicates()
    m.to_csv('CDS_allc.tmp', sep='\t', index=False)
    feature_metaplot('CDS_allc.tmp',features,genome_file,output,ignoreStrand,
                     windows,updown_stream,cutoff,first_feature,filter_chr)
    for i in ['CDS_allc.tmp','c_tmp','f_tmp']:
        os.remove(i)

#output per-site methylation levels for mCs in each specified context
def per_site_mC(allc,output_path,context=['CG','CH']):
    for i in context:
        a = filter_context(allc,[i])
        a = a[a['methylated'] == 1]
        a['mc_level'] = a['mc_count']/a['total']
        a.to_csv(output_path+i+'_site_methylation.txt', sep='\t', index=False)

# count the subcontexts in fasta
def count_subcontext_fasta(fasta,context=['CG','CH'],output=(),filter_chr=[]):
    df = pd.DataFrame(columns=['context','total_bases'])
    for c in context:
        count = 0
        for i in expand_nucleotide_code([c]):
            for sequence in SeqIO.parse(fasta, "fasta"):
                if sequence.name not in filter_chr:
                    count = count + sequence.seq.count(i) + sequence.seq.reverse_complement().count(i)
        df = df.append({'context': c, 'total_bases': count}, ignore_index=True)
    if output:
        df.to_csv(output, sep='\t', index=False)
    else:
        return df

#Count up subcontext methylation
def count_subcontext_allc(allc,context=['CG','CH'],output=(),cutoff=3,filter_chr=[]):
    a = pd.read_table(allc)
    df = pd.DataFrame(columns=['context','total_passed','methylated','weighted_mC'])
    for c in context:
        tp = mC = t = m = 0
        for i in expand_nucleotide_code([c]):
            i_table = a[(a['mc_class']==i) & (~a['chr'].isin(filter_chr))]
            t = t + i_table['total'].values.sum()
            m = m + i_table['mc_count'].values.sum()
            i_table = i_table[i_table['total']>=cutoff]
            tp = tp + len(i_table.index)
            mC = mC + len(i_table[i_table['methylated']==1].index)
        df = df.append({'context': c, 'total_passed': tp, 'methylated': mC,
                        "weighted_mC": np.float64(m)/np.float64(t)}, ignore_index=True)
    if output:
        df.to_csv(output, sep='\t', index=False)
    else:
        return df

#Get subcontexts of genome and methylation
def subcontext_methylation(allc,fasta,context=['CG','CH'],output=(),cutoff=3,filter_chr=[]):
    a = count_subcontext_fasta(fasta,context,output=(),filter_chr=filter_chr)
    b = count_subcontext_allc(allc,context,output=(),cutoff=3,filter_chr=filter_chr)
    df = pd.merge(a,b,on='context')
    if output:
        df.to_csv(output, sep='\t', index=False)
    else:
        return df

#filter DMR output from methylpy
#use filter_dmr( collapsed DMR file from methylpy, output file for filtered data, minimum number DMS (default = 5)
#minimum difference in methylation between highest and lowest methylated sample (default = 0.1 or 10% difference)
def filter_dmr(dmr_file,output=(),min_dms=5,min_mC_diff=0.1):
    a=pd.read_table(dmr_file)
    list=[]
    for c in a.itertuples():
        if c[4] >= min_dms:
            if max(c[7:])-min(c[7:]) >= min_mC_diff:
                list.append(c[0])
    a=a.ix[list]
    if output:
        a.to_csv(output, sep='\t', index=False)
    else:
        return a
