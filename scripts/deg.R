#Define functions
makeResultsTable <- function(x,conditionA,conditionB,filter=FALSE){
    require(DESeq2)
    bml <- sapply(levels(x$condition),function(lvl) rowMeans(counts(x,normalized=TRUE)[,x$condition == lvl]))
    bml <- as.data.frame(bml)
    y <- results(x,contrast=c("condition",conditionA,conditionB),independentFiltering=filter)
    y <- data.frame(id=row.names(y),
                    sampleA=c(conditionA),sampleB=c(conditionB),
                    baseMeanA=bml[,conditionA],baseMeanB=bml[,conditionB],
                    log2FC=y$log2FoldChange,pval=y$pvalue,padj=y$padj)
    row.names(y) <- c(1:nrow(y))
    return(y)
}

##
setwd("rnaseq/")

library(readr)
library(tximport)
library(DESeq2)


tx2gene <- read.csv("../ref/misc/tx2gene.csv",header=T)
files <- c("../FCH1/rnaseq/salmon_gibbs/quant.sf","../FCH2/rnaseq/salmon_gibbs/quant.sf",
"../FCH3/rnaseq/salmon_gibbs/quant.sf","../FTH1/rnaseq/salmon_gibbs/quant.sf",
"../FTH2/rnaseq/salmon_gibbs/quant.sf","../FTH3/rnaseq/salmon_gibbs/quant.sf",
"../FCNT1/rnaseq/salmon_gibbs/quant.sf","../FCNT2/rnaseq/salmon_gibbs/quant.sf",
"../FCNT3/rnaseq/salmon_gibbs/quant.sf","../FTNT1/rnaseq/salmon_gibbs/quant.sf",
"../FTNT2/rnaseq/salmon_gibbs/quant.sf","../FTNT3/rnaseq/salmon_gibbs/quant.sf",
"../MCH1/rnaseq/salmon_gibbs/quant.sf","../MCH2/rnaseq/salmon_gibbs/quant.sf",
"../MCH3/rnaseq/salmon_gibbs/quant.sf","../MTH1/rnaseq/salmon_gibbs/quant.sf",
"../MTH2/rnaseq/salmon_gibbs/quant.sf","../MTH3/rnaseq/salmon_gibbs/quant.sf",
"../MCNT1/rnaseq/salmon_gibbs/quant.sf","../MCNT2/rnaseq/salmon_gibbs/quant.sf",
"../MCNT3/rnaseq/salmon_gibbs/quant.sf","../MTNT1/rnaseq/salmon_gibbs/quant.sf",
"../MTNT2/rnaseq/salmon_gibbs/quant.sf","../MTNT3/rnaseq/salmon_gibbs/quant.sf")
names(files) <- c("FCH1","FCH2","FCH3","FTH1","FTH2","FTH3","FCNT1","FCNT2",
"FCNT3","FTNT1","FTNT2","FTNT3","MCH1","MCH2","MCH3","MTH1","MTH2","MTH3",
"MCNT1","MCNT2","MCNT3","MTNT1","MTNT2","MTNT3")

txi <- tximport(files, type = "salmon", tx2gene = tx2gene, reader = read_tsv)
sampleTable <- data.frame(condition = factor(rep(c("FCH", "FTH","FCNT","FTNT","MCH","MTH","MCNT","MTNT"), each = 3)))
rownames(sampleTable) <- colnames(txi$counts)
dds <- DESeqDataSetFromTximport(txi, sampleTable, ~condition)
dds <- dds[ rowSums(counts(dds)) > 1, ]
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds,fitType="local")

txi <- tximport(files[1:6], type = "salmon", tx2gene = tx2gene, reader = read_tsv)
sampleTable <- data.frame(condition = factor(rep(c("FCH", "FTH"), each = 3)))
rownames(sampleTable) <- colnames(txi$counts)
FH <- DESeqDataSetFromTximport(txi, sampleTable, ~condition)
FH <- FH[ rowSums(counts(FH)) > 1, ]
FH <- estimateSizeFactors(FH)
FH <- estimateDispersions(FH,fitType="local")
FH <- nbinomWaldTest(FH)
resFH <- makeResultsTable(FH,"FCH","FTH",filter=FALSE)
FHcounts <- counts(FH, normalized=TRUE)

txi <- tximport(files[7:12], type = "salmon", tx2gene = tx2gene, reader = read_tsv)
sampleTable <- data.frame(condition = factor(rep(c("FCNT", "FTNT"), each = 3)))
rownames(sampleTable) <- colnames(txi$counts)
FNT <- DESeqDataSetFromTximport(txi, sampleTable, ~condition)
FNT <- FNT[ rowSums(counts(FNT)) > 1, ]
FNT <- estimateSizeFactors(FNT)
FNT <- estimateDispersions(FNT,fitType="local")
FNT <- nbinomWaldTest(FNT)
resFNT <- makeResultsTable(FNT,"FCNT","FTNT",filter=FALSE)
FNTcounts <- counts(FNT, normalized=TRUE)

txi <- tximport(files[13:18], type = "salmon", tx2gene = tx2gene, reader = read_tsv)
sampleTable <- data.frame(condition = factor(rep(c("MCH", "MTH"), each = 3)))
rownames(sampleTable) <- colnames(txi$counts)
MH <- DESeqDataSetFromTximport(txi, sampleTable, ~condition)
MH <- MH[ rowSums(counts(MH)) > 1, ]
MH <- estimateSizeFactors(MH)
MH <- estimateDispersions(MH,fitType="local")
MH <- nbinomWaldTest(MH)
resMH <- makeResultsTable(MH,"MCH","MTH",filter=FALSE)
MHcounts <- counts(MH, normalized=TRUE)

txi <- tximport(files[19:24], type = "salmon", tx2gene = tx2gene, reader = read_tsv)
sampleTable <- data.frame(condition = factor(rep(c("MCNT", "MTNT"), each = 3)))
rownames(sampleTable) <- colnames(txi$counts)
MNT <- DESeqDataSetFromTximport(txi, sampleTable, ~condition)
MNT <- MNT[ rowSums(counts(MNT)) > 1, ]
MNT <- estimateSizeFactors(MNT)
MNT <- estimateDispersions(MNT,fitType="local")
MNT <- nbinomWaldTest(MNT)
resMNT <- makeResultsTable(MNT,"MCNT","MTNT",filter=FALSE)
MNTcounts <- counts(MNT, normalized=TRUE)

total <- as.data.frame(rbind(resFH,resFNT,resMH,resMNT))
total$padj <- p.adjust(total$pval,method="BH")
sig <- na.omit(total[total$padj <= 0.05,])

write.table(total,"../../figures_tables/deseq.tsv",quote=F,sep="\t",row.names=F)
write.table(sig,"../../figures_tables/sig_genes.tsv",quote=F,sep="\t",row.names=F)

path <- c("../../figures_tables/sig_genes/")
ifelse(!dir.exists(path),
dir.create(path), FALSE)
for(gene in gof$gene){
    x <- plotCounts(dds, gene, intgroup = "condition", normalized = TRUE,
                    transform = FALSE, returnData = TRUE)
    p <- ggplot(x, aes(x = condition, y = count)) + geom_jitter(aes(color=condition)) +
         theme(panel.background=element_blank(),
               axis.line=element_line(color="black"),
               axis.text=element_text(color="black"),
               axis.title=element_text(color="black",face="bold"),
               legend.position="none")
    ggsave(paste(path,gene,"_counts.pdf",sep=""), p, width=5, height=4)
}
