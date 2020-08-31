setwd('methylation/')
library(ggplot2)
library(reshape2)

path <- c("../../figures_tables/methylation/genes_of_interest/")
ifelse(!dir.exists(path),
dir.create(path), FALSE)
df <- read.table('goi_methylation_levels.tsv',header=T,sep='\t')
df2 <- melt(df)
df2$condition = rep(c("MCH","MTH","MCNT","MTNT","FCH","FTH","FCNT","FTNT"),each = 216)
df2$order = rep(c(5,6,7,8,1,2,3,4), each = 216)
df2$mC_type = rep(c('CG','CH'),each=36)
for(gene in df$Gene){
    #CG
    x <- df2[df2$Gene==gene & df2$mC_type=='CG',]
    p <- ggplot(x, aes(x = reorder(condition,order), y = value)) +
                geom_point(aes(color=condition), size = 4) +
                theme(panel.background=element_blank(),
                axis.line=element_line(color="black"),
                axis.text=element_text(color="black"),
                axis.title=element_text(color="black",face="bold"),
                legend.position="none") + xlab("Sample") + ylab("Methylation level") +
                scale_y_continuous(limits=c(0,1),expand=c(0,0),
                breaks=c(0.25,0.5,0.75,1),labels=c("25%","50%","75%","100%")) +
    ggsave(paste(path,gene,"_CG_methylation_level.pdf",sep=""), p, width=5, height=4)
}
for(gene in df$Gene){
    x <- df2[df2$Gene==gene & df2$mC_type=='CH',]
    p <- ggplot(x, aes(x = reorder(condition,order), y = value)) +
                geom_point(aes(color=condition), size = 4) +
                theme(panel.background=element_blank(),
                axis.line=element_line(color="black"),
                axis.text=element_text(color="black"),
                axis.title=element_text(color="black",face="bold"),
                legend.position="none") + xlab("Sample") + ylab("Methylation level") +
                scale_y_continuous(limits=c(0,0.1),expand=c(0,0),
                breaks=c(0.025,0.05,0.075,0.1),labels=c("2.5%","5.0%","7.5%","10.0%"))
    ggsave(paste(path,gene,"_CH_methylation_level.pdf",sep=""), p, width=5, height=4)
}

