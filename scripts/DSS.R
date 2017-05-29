setwd("Analysis/git/zebra-finch-brains/data/")


library(DSS)
require(bsseq)

MCH1 <- read.table("MCH1/methylCseq/MCH1_dss_format.tsv",header=T)
MCH2 <- read.table("MCH2/methylCseq/MCH2_dss_format.tsv",header=T)
MCH3 <- read.table("MCH3/methylCseq/MCH3_dss_format.tsv",header=T)
MTH1 <- read.table("MTH1/methylCseq/MTH1_dss_format.tsv",header=T)
MTH2 <- read.table("MTH2/methylCseq/MTH2_dss_format.tsv",header=T)
MTH3 <- read.table("MTH3/methylCseq/MTH3_dss_format.tsv",header=T)

BSobj <- makeBSseqData(list(MCH1,MCH2,MCH3,MTH1,MTH2,MTH3),
                       c("C1","C2","C3","T1","T2","T3"))

dmlTest <- DMLtest(BSobj, group1=c("C1","C2","C3"), group2=c("T1","T2","T3"))
dmlTest.sm <- DMLtest(BSobj, group1=c("C1","C2","C3"), group2=c("T1","T2","T3"),smoothing=TRUE)

dmls <- callDML(dmlTest, p.threshold=0.001)
dmls.sm <- callDML(dmlTest.sm, p.threshold=0.001)

dmrs <- callDMR(dmlTest, p.threshold=0.01)
dmrs.sm <- callDMR(dmlTest.sm, p.threshold=0.01)
