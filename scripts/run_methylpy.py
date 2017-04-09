#!/usr/local/apps/python/2.7.8/bin/python
import sys
import getopt
from methylpy.call_mc import run_methylation_pipeline


def main(argv):

	#Sample name
	sample = argv[0]

	#Path to fastq files
	files = [argv[1]]

	#Denoting corresponding library
	libraries = ["libA"]

	#Genome (take col0_t9 as example)
	f_ref = argv[2] + "_f"
	r_ref = argv[2] + "_r"
	ref_fasta = argv[2] + ".fa"

	#Number of processors
	num_procs = int(argv[3])

	#Sort memory
	sort_mem = argv[4] + "G"

	#Adapter
	adapter_seq = 'AGATCGGAAGAGCTCGTATGCC'

	#Control "chrC:" or "chrL:"
	m_control = argv[5] + ":"

	run_methylation_pipeline(files=files,
						 	 libraries=libraries,
						 	 sample=sample,
						 	 forward_reference=f_ref,
						 	 reverse_reference=r_ref,
						 	 reference_fasta=ref_fasta,
						 	 unmethylated_control=m_control,
                    	 	 quality_version="1.8",
                    	 	 path_to_samtools="",
                    	 	 path_to_bowtie="",
                         	 bowtie_options=["-S","-k 1","-m 1","--chunkmbs 3072","--best","--strata","-o 4","-e 80","-l 20","-n 0"],
                        	 num_procs=num_procs,
                         	 trim_reads=True,
                         	 path_to_cutadapt="",
                         	 adapter_seq=adapter_seq,
                         	 max_adapter_removal=None,
                         	 overlap_length=None,
                         	 zero_cap=None,
                         	 error_rate=None,
                         	 min_qual_score=10,
                         	 min_read_len=30,
                         	 sig_cutoff=0.01,
                         	 min_cov=3,
                        	 binom_test=True,
                         	 bh=True,
                         	 keep_temp_files=False,
                         	 num_reads=-1,
                         	 save_space=True,
                         	 bowtie2=False,
                         	 sort_mem=sort_mem,
                         	 path_to_output="",
                         	 in_mem=False,
                        	 min_base_quality=1,
                         	 path_to_MarkDuplicates=False
                         	 )


if __name__ == "__main__":
   main(sys.argv[1:])
