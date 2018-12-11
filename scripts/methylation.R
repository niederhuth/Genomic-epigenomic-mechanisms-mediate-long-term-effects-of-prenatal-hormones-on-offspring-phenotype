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
#CG group metaplots
for(x in c('male','female')){
  for(y in c('Nucleus taeniae','Hypothalamus')){
    
    df <- matrix(,nrow=60,ncol=0)
    for(i in names <- samples[samples$Sex==x & samples$Tissue==y,]$Name){
      input=paste(i,"gene_metaplot.txt",sep="_")
      df <- cbind(df,read.table(input,header=T,sep="\t")[2])
    }
    colnames(df) <- names
    df <- melt(df)
    df$Bin <- 1:60
    print(names)
    
    plot <- ggplot(df, aes(x=Bin, y=value, color=variable)) + geom_line(size=0.8) +
      theme(panel.background=element_blank(), panel.grid=element_blank(),
            axis.text.y=element_text(color="black"), axis.text.x=element_text(color="black"),
            axis.ticks=element_line(color="black"), axis.title=element_text(color="black"),
            legend.position="right", axis.line=element_line(color="black")) +
      ylab("Percent methylation") + 
      xlab( "" ) +
      scale_y_continuous(limits=c(0,1), expand=c(0,0), breaks=c(0,0.25,0.75,1),labels=c("25%","50%","75%","100%")) +
      geom_vline(xintercept=20, linetype="longdash", color="grey55") +
      geom_vline(xintercept=40, linetype="longdash", color="grey55") +
      scale_x_continuous(labels=c("-2000","TSS'","TTS'","+2000"), breaks=c(1, 20, 40, 60)) +
      scale_color_manual("",values=c("steelblue4","steelblue2","blue","goldenrod4","goldenrod2","gold"))
    
    ggsave(paste(x,y,"CG_gene_metaplot.pdf",sep="_"),plot,path="group_metaplots/")
  }
}  

#CH group metaplots
for(x in c('male','female')){
  for(y in c('Nucleus taeniae','Hypothalamus')){
    
    df <- matrix(,nrow=60,ncol=0)
    for(i in names <- samples[samples$Sex==x & samples$Tissue==y,]$Name){
      input=paste(i,"gene_metaplot.txt",sep="_")
      df <- cbind(df,read.table(input,header=T,sep="\t")[6])
    }
    colnames(df) <- names
    df <- melt(df)
    df$Bin <- 1:60
    print(names)
    
    plot <- ggplot(df, aes(x=Bin, y=value, color=variable)) + geom_line(size=0.8) +
      theme(panel.background=element_blank(), panel.grid=element_blank(),
            axis.text.y=element_text(color="black"), axis.text.x=element_text(color="black"),
            axis.ticks=element_line(color="black"), axis.title=element_text(color="black"),
            legend.position="right", axis.line=element_line(color="black")) +
      ylab("Percent methylation") + 
      xlab( "" ) +
      scale_y_continuous(limits=c(0,1), expand=c(0,0), breaks=c(0,0.25,0.75,1),labels=c("25%","50%","75%","100%")) +
      geom_vline(xintercept=20, linetype="longdash", color="grey55") +
      geom_vline(xintercept=40, linetype="longdash", color="grey55") +
      scale_y_continuous(limits=c(0,0.10), expand=c(0,0), breaks=c(0.025,0.05,0.075,0.10),labels=c("2.5%","5%","7.5%","10%")) +     
      scale_color_manual("",values=c("steelblue4","steelblue2","blue","goldenrod4","goldenrod2","gold"))
    
    ggsave(paste(x,y,"CH_gene_metaplot.pdf",sep="_"),plot,path="group_metaplots/")
  }
}  


#sample metaplots
for(i in names <- samples$Name){
  input=paste(i,"gene_metaplot.txt",sep="_")
  df <- read.table(input,header=T,sep="\t")[c(2,6)]
  df <- melt(df)
  df$Bin <- 1:60
  df$variable <- gsub("_Weighted_mC","",df$variable)
  
  plot <- ggplot(df, aes(x=Bin, y=value, color=variable)) + geom_line(size=0.8) +
    theme(panel.background=element_blank(), panel.grid=element_blank(),
          axis.text.y=element_text(color="black"), axis.text.x=element_text(color="black"),
          axis.ticks=element_line(color="black"), axis.title=element_text(color="black"),
          legend.position="right", axis.line=element_line(color="black")) +
    ylab("Weighted methylation") + 
    xlab( "" ) +
    scale_y_continuous(limits=c(0,1), expand=c(0,0), breaks=c(0.25,0.5,0.75,1),labels=c("25%","50%","75%","100%")) +
    geom_vline(xintercept=20, linetype="longdash", color="grey55") +
    geom_vline(xintercept=40, linetype="longdash", color="grey55") +
    scale_x_continuous(labels=c("-2000","TSS","TTS","+2000"), breaks=c(1, 20, 40, 60)) +
    scale_color_manual("",values=c("black","tomato3"))
  
  ggsave(paste(i,"_gene_metaplot.pdf",sep="_"),plot,path="sample_metaplots/")
}

#CH sample metaplots
for(i in names <- samples$Name){
  input=paste(i,"gene_metaplot.txt",sep="_")
  df <- read.table(input,header=T,sep="\t")[3:6]
  df <- melt(df)
  df$Bin <- 1:60
  df$order <- c(rep(1,60),rep(2,60),rep(3,60),rep(4,60))
  df$variable <- gsub("_Weighted_mC","",df$variable)
  
  plot <- ggplot(df, aes(x=Bin, y=value, color=reorder(variable,order))) + geom_line(size=0.8) +
    theme(panel.background=element_blank(), panel.grid=element_blank(),
          axis.text.y=element_text(color="black"), axis.text.x=element_text(color="black"),
          axis.ticks=element_line(color="black"), axis.title=element_text(color="black"),
          legend.position="right", axis.line=element_line(color="black")) +
    ylab("Weighted methylation") + 
    xlab( "" ) +
    scale_y_continuous(limits=c(0,0.10), expand=c(0,0), breaks=c(0.025,0.05,0.075,0.10),labels=c("2.5%","5%","7.5%","10%")) +
    geom_vline(xintercept=20, linetype="longdash", color="grey55") +
    geom_vline(xintercept=40, linetype="longdash", color="grey55") +
    scale_x_continuous(labels=c("-2000","TSS","TTS","+2000"), breaks=c(1, 20, 40, 60)) +
    scale_color_manual("",values=c("tomato3","orange3","magenta3","royalblue3"))
  
  ggsave(paste(i,"CH_gene_metaplot.pdf",sep="_"),plot,path="sample_metaplots/")
}

