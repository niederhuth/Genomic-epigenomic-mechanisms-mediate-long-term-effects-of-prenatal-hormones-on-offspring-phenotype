#!/usr/bin/env Rscript

library(ggplot2)
library(reshape2)

samples <- read.csv('../../../misc/samples.csv',header=T)

#
df <- matrix(,nrow=5,ncol=0)
for(i in names <- samples$Name){
  input=paste(i,"total_weighted_methylation.txt",sep="_")
  df <- cbind(df,read.table(input,header=T,sep="\t")[6])
}
colnames(df) <- names
row.names(df) <- c('CG','CA','CT','CC','CH')
df <- melt(df)
df$context <- c('CG','CA','CT','CC','CH')

#gene metaplots
df <- matrix(,nrow=60,ncol=0)
for(i in names <- samples[samples$Sex=='female' & samples$Tissue=='Nucleus taeniae',]$Name){
  input=paste(i,"gene_metaplot.txt",sep="_")
  df <- cbind(df,read.table(input,header=T,sep="\t")[2])
}
colnames(df) <- names
df <- melt(df)
df$Bin <- 1:60
print(names)

ggplot(df, aes(x=Bin, y=value, color=variable)) + geom_line(size=0.8) +
              theme(panel.background=element_blank(), panel.grid=element_blank(),
                    axis.text.y=element_text(color="black"), axis.text.x=element_text(color="black"),
                    axis.ticks=element_line(color="black"), axis.title=element_text(color="black"),
                    legend.position="none", axis.line=element_line(color="black")) +
              ylab("Percent methylation") + 
              xlab( "" ) +
              scale_y_continuous(limits=c(0,1), expand=c(0,0), breaks=c(0.05,0.1,0.15,0.2),labels=c("5%","10%","15%","20%")) +
              geom_vline(xintercept=20, linetype="longdash", color="grey55") +
              geom_vline(xintercept=40, linetype="longdash", color="grey55") +
              scale_x_continuous(labels=c("-2000","TSS'","TTS'","+2000"), breaks=c(1, 20, 40, 60)) +
              scale_color_manual("",values=c("black","black","black","tomato3","tomato3","tomato3"))

