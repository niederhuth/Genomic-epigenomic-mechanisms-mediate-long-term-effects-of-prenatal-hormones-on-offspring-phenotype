---
Daily status of the analysis
---
<<<<<<< HEAD
#### April 13, 2017
1. Run test DMR analysis -> failed
2. MTH3 and MTNT1 for some reason failed during methylCseq mapping
3. Restart methylCseq mapping for MTH3 and MTNT1
=======
#### To do
1. DMR analysis
  1. Filter DMRs
  2. Cluster Samples based on DMRs
    1. Does it cluster by treatment, sex, or testosterone level?
    2. Associate with differentially expressed genes
  3. Gene metaplots
    1. Compare differentially expressed genes
  4. etc
2. RNAseq
  1. Quality control metrics
  2. Differential expression
  3. GO term enrichment

#### April 18, 2017
1. Make scripts for preparing browser files
2. filter DMRs and examine in browser
3. Add functions.py
4. Add scripts for DMR analysis

#### April 16, 2017
1. Run cuffdiff_cuffnorm.sh
2. Work on differential expression

#### April 15, 2017
1. methylCseq mapping completed
2. Run prepare_DMR.sh
3. Fix DMR scripts (missing parentheses in chromosomes)
4. Start DMR analysis
5. Run cuffquant.sh
6. Make HTseq-count.sh script
7. Run HTseq-count.sh

#### April 14, 2017
1. Fixed initial trimming and reverse complementation of methylCseq reads
2. Rerun methylCseq mapping

#### April 10, 2017
1. Added scripts for DMR analysis
2. methylCseq mapping still running

#### April 9, 2017
1. Fixed scripts for methylCseq and RNAseq mapping
2. Started methylCseq mapping
3. Started RNAseq mapping

#### Dec 16, 2016
1. Edited scripts
2. Remade all reference files
3. Restarted trim_align_rnaseq on all samples
4. Restarted methylpy on all samples
