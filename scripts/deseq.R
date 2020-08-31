#Define functions
makeResultsTable <- function(x,conditionA,conditionB,filter=FALSE){
  require(DESeq2)
  bml <- sapply(levels(dds$condition),function(lvl) rowMeans(counts(dds,normalized=TRUE)[,dds$condition == lvl]))
  bml <- as.data.frame(bml)
  y <- results(x,contrast=c("condition",conditionA,conditionB),independentFiltering=filter)
  y <- data.frame(id=gsub(pattern = "gene:", replacement = "", row.names(y)),
                  sampleA=c(conditionA),sampleB=c(conditionB),
                  baseMeanA=bml[,conditionA],baseMeanB=bml[,conditionB],
                  log2FC=y$log2FoldChange,pval=y$pvalue,padj=y$padj)
  row.names(y) <- c(1:nrow(y))
  return(y)
}

sampleHeatMap <- function(x){
  require(pheatmap)
  require("RColorBrewer")
  sampleDists <- dist(t(assay(x)))
  sampleDistMatrix <- as.matrix(sampleDists)
  rownames(sampleDistMatrix) <- colnames(x)
  colnames(sampleDistMatrix) <- NULL
  colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
  pheatmap(sampleDistMatrix,
           clustering_distance_rows=sampleDists,
           clustering_distance_cols=sampleDists,
           col=colors)
}

geneHeatMap <- function(dds,geneList){
  require(pheatmap)
  select <- select <- row.names(counts(dds,normalized=TRUE)) %in% geneList
  nt <- normTransform(dds) # defaults to log2(x+1)
  log2.norm.counts <- assay(nt)[select,]
  COL <- as.data.frame(colData(dds)[,c("condition")])
  pheatmap(log2.norm.counts, cluster_rows=TRUE, show_rownames=TRUE,
           cluster_cols=FALSE, annotation_col=df)
}
sampleNames <- c("FCH1","FCH2","FCH3","FCNT1","FCNT2","FCNT3","FTH1","FTH2","FTH3","FTNT1","FTNT2","FTNT3","MCH1","MCH2","MCH3","MCNT1","MCNT2","MCNT3","MTH1","MTH2","MTH3","MTNT1","MTNT2","MTNT3")
sampleFiles <- c("FCH1.counts.tab","FCH2.counts.tab","FCH3.counts.tab","FCNT1.counts.tab","FCNT2.counts.tab","FCNT3.counts.tab","FTH1.counts.tab","FTH2.counts.tab","FTH3.counts.tab","FTNT1.counts.tab","FTNT2.counts.tab","FTNT3.counts.tab","MCH1.counts.tab","MCH2.counts.tab","MCH3.counts.tab","MCNT1.counts.tab","MCNT2.counts.tab","MCNT3.counts.tab","MTH1.counts.tab","MTH2.counts.tab","MTH3.counts.tab","MTNT1.counts.tab","MTNT2.counts.tab","MTNT3.counts.tab")
sampleGroup <- c("FCH","FCH","FCH","FCNT","FCNT","FCNT","FTH","FTH","FTH","FTNT","FTNT","FTNT","MCH","MCH","MCH","MCNT","MCNT","MCNT","MTH","MTH","MTH","MTNT","MTNT","MTNT")
sampleCondition <- c("control","control","control","control","control","control","treatment","treatment","treatment","treatment","treatment","treatment","control","control","control","control","control","control","treatment","treatment","treatment","treatment","treatment","treatment")
sampleSex <- c("female","female","female","female","female","female","female","female","female","female","female","female","male","male","male","male","male","male","male","male","male","male","male","male")
sampleTissue <- c("hypothalamus","hypothalamus","hypothalamus","nucleus_taenia","nucleus_taenia","nucleus_taenia","hypothalamus","hypothalamus","hypothalamus","nucleus_taenia","nucleus_taenia","nucleus_taenia","hypothalamus","hypothalamus","hypothalamus","nucleus_taenia","nucleus_taenia","nucleus_taenia","hypothalamus","hypothalamus","hypothalamus","nucleus_taenia","nucleus_taenia","nucleus_taenia")
sampleTable <- data.frame(sampleName = sampleNames, fileName = sampleFiles, group = sampleGroup, condition = sampleGroup, sex = sampleSex, tissue = sampleTissue)
dds <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, directory = 'htseq/', design= ~ group)
dds <- dds[ rowSums(counts(dds)) > 10, ]

dds <- DESeq(dds)
res <- results(dds)
summary(res)

vsd <- vst(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE)
head(assay(vsd), 3)
                                       
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds,fitType="local")
dds <- nbinomWaldTest(dds)
res <- results(dds)
summary(res)

makeResultsTable(dds,"MCH","MTH",filter=FALSE)
