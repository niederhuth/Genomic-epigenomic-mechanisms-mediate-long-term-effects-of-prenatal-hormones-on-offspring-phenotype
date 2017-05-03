topGO <- function(genelist,goTerms,nodeSize,filename,writeData=FALSE){
    require(topGO)
    require(GO.db)
    path <- c("../../figures_tables/goTerms/")
    ifelse(!dir.exists(path),
    dir.create(path), FALSE)
    BP <- new("topGOdata",description="Biological Process",ontology="BP",
              allGenes=genelist,annot=annFUN.gene2GO,nodeSize=nodeSize,gene2GO=goTerms)
    MF <- new("topGOdata",description="Molecular Function",ontology="MF",
              allGenes=genelist,annot=annFUN.gene2GO,nodeSize=nodeSize,gene2GO=goTerms)
    CC <- new("topGOdata",description="Cellular Compartment",ontology="CC",
              allGenes=genelist,annot=annFUN.gene2GO,nodeSize=nodeSize,gene2GO=goTerms)
    FisherBP <- runTest(BP,algorithm="parentchild",statistic="fisher")
    FisherMF <- runTest(MF,algorithm="parentchild",statistic="fisher")
    FisherCC <- runTest(CC,algorithm="parentchild",statistic="fisher")
    BPgenTable <- GenTable(BP,Fisher=FisherBP,ranksOf="Fisher",topNodes=length(score(FisherBP)))
    MFgenTable <- GenTable(MF,Fisher=FisherMF,ranksOf="Fisher",topNodes=length(score(FisherMF)))
    CCgenTable <- GenTable(CC,Fisher=FisherCC,ranksOf="Fisher",topNodes=length(score(FisherCC)))
    write.csv(BPgenTable[BPgenTable$Fisher < 0.01,],paste(path,filename,"_BP.csv",sep=""),row.names=FALSE,quote=FALSE)
    write.csv(MFgenTable[MFgenTable$Fisher < 0.01,],paste(path,filename,"_MF.csv",sep=""),row.names=FALSE,quote=FALSE)
    write.csv(CCgenTable[CCgenTable$Fisher < 0.01,],paste(path,filename,"_CC.csv",sep=""),row.names=FALSE,quote=FALSE)
    pdf(paste(path,filename,"_BP_top5_nodes.pdf",sep=""),width=6,height=6,paper='special')
    showSigOfNodes(BP, score(FisherBP), firstSigNodes = 5, useInfo = 'all')
    dev.off()
    pdf(paste(path,filename,"_MF_top5_nodes.pdf",sep=""),width=6,height=6,paper='special')
    showSigOfNodes(MF, score(FisherMF), firstSigNodes = 5, useInfo = 'all')
    dev.off()
    pdf(paste(path,filename,"_CC_top5_nodes.pdf",sep=""),width=6,height=6,paper='special')
    showSigOfNodes(CC, score(FisherCC), firstSigNodes = 5, useInfo = 'all')
    dev.off()
    if(writeData){
      return(list(BP=BPgenTable,MF=MFgenTable,CC=CCgenTable))
    }
}


#run
setwd("rnaseq/")
library(cummeRbund)

#Redo FDR, prepare samples
cuff <- readCufflinks("cuffdiff_all2/")
df <- read.table("cuffdiff_all2/gene_exp.diff",header=T,sep="\t")
df$test_id <- gsub("gene:","",df$test_id)
df$gene_id <- gsub("gene:","",df$gene_id)
df <- df[df$sample_1 == "FCH" & df$sample_2 == "FTH" | df$sample_1 == "MCH" & df$sample_2 == "MTH" | df$sample_1 == "FCNT" & df$sample_2 == "FTNT" | df$sample_1 == "MCNT" & df$sample_2 == "MTNT",]
df <- df[df$status == "OK",]
df$q_value <- p.adjust(df$p_value,method="BH")
df$significant <- ifelse(df$q_value <= 0.05, "yes", "no")
sig <- df[df$q_value <= 0.05,]

FH <- sig[sig$sample_1 == "FCH" & sig$sample_2 == "FTH",]
FNT <- sig[sig$sample_1 == "FCNT" & sig$sample_2 == "FTNT",]
MH <- sig[sig$sample_1 == "MCH" & sig$sample_2 == "MTH",]
MNT <- sig[sig$sample_1 == "MCNT" & sig$sample_2 == "MTNT",]

#write.table(total,"../../figures_table/.tsv",quote=F,sep="\t",row.names=F)
write.table(sig,"../../figures_tables/sig_genes.tsv",quote=F,sep="\t",row.names=F)

#venn
d <- data.frame(gene_id=unique(sig$gene_id))
d <- data.frame(gene_id=d$gene_id,FH=ifelse(d$gene_id %in% sig[sig$sample_1 == "FCH" & sig$sample_2 == "FTH",]$gene_id, 1, 0),
  FNT=ifelse(d$gene_id %in% sig[sig$sample_1 == "FCNT" & sig$sample_2 == "FTNT",]$gene_id, 1, 0),
  MH=ifelse(d$gene_id %in% sig[sig$sample_1 == "MCH" & sig$sample_2 == "MTH",]$gene_id, 1, 0),
  MNT=ifelse(d$gene_id %in% sig[sig$sample_1 == "MCNT" & sig$sample_2 == "MTNT",]$gene_id, 1, 0)
)
library(VennDiagram)
pdf("../../figures_tables/comparisonVennDiagram.pdf",width=6,height=6,paper='special')
draw.quad.venn( area1 = nrow(subset(d, FH == 1)),
  area2 = nrow(subset(d, MH == 1)),
  area3 = nrow(subset(d, FNT == 1)),
  area4 = nrow(subset(d, MNT == 1)),
  n12 = nrow(subset(d, FH == 1 & MH == 1)),
  n13 = nrow(subset(d, FH == 1 & FNT == 1)),
  n14 = nrow(subset(d, FH == 1 & MNT == 1)),
  n23 = nrow(subset(d, MH == 1 & FNT == 1)),
  n24 = nrow(subset(d, MH == 1 & MNT == 1)),
  n34 = nrow(subset(d, FNT == 1 & MNT == 1)),
  n123 = nrow(subset(d, FH == 1 & MH == 1 & FNT == 1 )),
  n134 = nrow(subset(d, FH == 1 & FNT == 1 & MNT == 1)),
  n124 = nrow(subset(d, FH == 1 & MH == 1 & MNT == 1)),
  n234 = nrow(subset(d, MH == 1 & FNT == 1 & MNT == 1)),
  n1234 = nrow(subset(d, FH == 1 & MH == 1 & FNT == 1 & MNT == 1)),
  category = c("Female Hypothalamus", "Male Hypothalamus",
  "Female Nucleus taeniae", "Male Nucleus taeniae"),
  lty = "blank", fill = c("skyblue", "pink1", "mediumorchid", "orange"))
dev.off()

#GO terms
library(topGO)
library(GO.db)
goTerms <- readMappings(file="../../misc/topGO.txt")

FHgo <- factor(as.integer(names(goTerms) %in% FH$gene_id ))
names(FHgo) <- names(goTerms)
FNTgo <- factor(as.integer(names(goTerms) %in% FNT$gene_id ))
names(FNTgo) <- names(goTerms)
MHgo <- factor(as.integer(names(goTerms) %in% MH$gene_id ))
names(MHgo) <- names(goTerms)
MNTgo <- factor(as.integer(names(goTerms) %in% MNT$gene_id ))
names(MNTgo) <- names(goTerms)

FHgotest <- topGO(FHgo,goTerms,nodeSize=0,"FH",writeData=TRUE)
FNTgotest <- topGO(FNTgo,goTerms,nodeSize=0,"FNT",writeData=TRUE)
MHgotest <- topGO(MHgo,goTerms,nodeSize=0,"MH",writeData=TRUE)
MNTgotest <- topGO(MNTgo,goTerms,nodeSize=5,"MNT",writeData=TRUE)

#gof
path <- c("../../figures_tables/genes_of_interest/")
ifelse(!dir.exists(path),
dir.create(path), FALSE)
gof <- read.csv("../../misc/gof.txt")
fpkm <- read.table("cuffdiff_all2/fpkm.tsv",header=T,sep="\t",row.names=1)
row.names(fpkm) <- gsub("gene:","",row.names(fpkm))
for(gene in gof$gene){
    x <- data.frame(t(fpkm[row.names(fpkm) == gene,]),
    condition = rep(c("FCH","FTH","FCNT","FTNT","MCH","MTH","MCNT","MTNT"),each = 3),
    order = rep(c(1,2,3,4,5,6,7,8), each = 3))
    p <- ggplot(x, aes(x = reorder(x[,2],x[,3]), y = x[,1])) +
                geom_point(aes(color=condition), size = 4) +
                theme(panel.background=element_blank(),
                axis.line=element_line(color="black"),
                axis.text=element_text(color="black"),
                axis.title=element_text(color="black",face="bold"),
                legend.position="none") + xlab("Sample") + ylab("FPKM")
    ggsave(paste(path,gene,"_fpkm.pdf",sep=""), p, width=5, height=4)
}
