#!/usr/local/apps/python/2.7.8/bin/python
import sys
import getopt
from methylpy.call_mc import run_methylation_pipeline_pe


def main(argv):

	#Sample name
	sample = argv[0]

	#Path to fastq files
	read1_files = [argv[1]]
	read2_files = [argv[2]]

	print(read1_files)
	print(read2_files)

	#Denoting corresponding library
	libraries = ["libA"]

	#Genome (take col0_t9 as example)
	f_ref = argv[3] + "_f"
	r_ref = argv[3] + "_r"
	ref_fasta = argv[3] + ".fasta"

	#Number of processors
	num_procs = int(argv[4])

	#Sort memory
	sort_mem = argv[5] + "G"

	#Adapter
	adapter_seq_R1 = "AGATCGGAAGAGCACACGTCTGAAC"
	adapter_seq_R2 = "AGATCGGAAGAGCGTCGTGTAGGGA"

	#Control "chrC:" or "chrL:"
	m_control = argv[6] + ":"

	run_methylation_pipeline_pe(read1_files=read1_files,
								read2_files=read2_files,
								libraries=libraries,
								sample=sample,
								forward_reference=f_ref,
								reverse_reference=r_ref,
								reference_fasta=ref_fasta,
								unmethylated_control=m_control,
                                quality_version="1.8",
								path_to_samtools="",
								path_to_bowtie="",
                                bowtie_options=["-X 1000", "-k 2", "--no-mixed", "--no-discordant"],
                                num_procs=1,
								trim_reads=True,
								path_to_cutadapt="",
                                adapter_seq_R1 = adapter_seq_R1,
                                adapter_seq_R2 = adapter_seq_R2,
                                max_adapter_removal=None,
                                overlap_length=None,
								zero_cap=None,
								error_rate=None,
								min_qual_score=10,
                                min_read_len=30,
								sig_cutoff=.01,
								min_cov=3,
								binom_test=True,
								bh=True,
                                keep_temp_files=False,
								num_reads=-1,
								save_space=True,
                                bowtie2=True,
								sort_mem=sort_mem,
								path_to_output="",
								in_mem=False,
								min_base_quality=1,
                                path_to_MarkDuplicates=False)

if __name__ == "__main__":
   main(sys.argv[1:])
