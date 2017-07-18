#!/usr/bin/env Rscript
library(ggplot2)

args = commandArgs(trailingOnly=TRUE)
#test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("Sample name must be given", call.=FALSE)
}

#define plotting function for features
plot_features <- function(df){
  require(ggplot2)
  ggplot(df, aes(x=Bin)) + geom_line(aes(y=mCG), color="dodgerblue4", size=0.8) +
         geom_line(aes(y=mCH), color="hotpink4", size=0.8) +
         theme(panel.background=element_blank(), panel.grid=element_blank(),
         axis.text.y=element_text(color="black"), axis.text.x=element_text(color="black"),
         axis.ticks=element_line(color="black"), axis.title=element_text(color="black"),
         legend.position="none", axis.line=element_line(color="black")) +
         ylab("Percent methylation") + xlab( "" ) +
         scale_y_continuous(limits=c(0,1), expand=c(0,0),
         breaks=c(0.25,0.5,0.75,1),labels=c("25%","50%","75%","100%")) +
         geom_vline(xintercept=20, linetype="longdash", color="grey55") +
         geom_vline(xintercept=40, linetype="longdash", color="grey55")
}

#plot gene metaplots
if(file.exists("results/gene_metaplot.tsv")){
  print("Making gene metaplot")
  df <- read.table("results/gene_metaplot.tsv",header=T,sep="\t")
  plot <- plot_features(df) +
          scale_x_continuous(labels=c("-2000","TSS'","TTS'","+2000"), breaks=c(1, 20, 40, 60))
  filename=paste("figures_tables/", args[1], "_gene_metaplot.pdf", sep="")
  ggsave(filename=filename, plot, height=4, width=4, useDingbats=F)
  rm(df,plot)
}

#plot diff exp gene metaplots
if(file.exists("results/diff_exp_metaplot.tsv")){
  print("Making diff exp gene metaplot")
  df <- read.table("results/diff_exp_metaplot.tsv",header=T,sep="\t")
  plot <- plot_features(df) +
          scale_x_continuous(labels=c("-2000","TSS'","TTS'","+2000"), breaks=c(1, 20, 40, 60))
  filename=paste("figures_tables/", args[1], "_diff_exp_metaplot.pdf", sep="")
  ggsave(filename=filename, plot, height=4, width=4, useDingbats=F)
  rm(df,plot)
}

#plot no diff exp gene metaplots
if(file.exists("results/no_diff_exp_metaplot.tsv")){
  print("Making no diff exp gene metaplot")
  df <- read.table("results/no_diff_exp_metaplot.tsv",header=T,sep="\t")
  plot <- plot_features(df) +
          scale_x_continuous(labels=c("-2000","TSS'","TTS'","+2000"), breaks=c(1, 20, 40, 60))
  filename=paste("figures_tables/", args[1], "_no_diff_exp_metaplot.pdf", sep="")
  ggsave(filename=filename, plot, height=4, width=4, useDingbats=F)
  rm(df,plot)
}

#plot goi gene metaplots
if(file.exists("results/goi_metaplot.tsv")){
  print("Making goi gene metaplot")
  df <- read.table("results/goi_metaplot.tsv",header=T,sep="\t")
  plot <- plot_features(df) +
          scale_x_continuous(labels=c("-2000","TSS'","TTS'","+2000"), breaks=c(1, 20, 40, 60))
  filename=paste("figures_tables/", args[1], "_goi_metaplot.pdf", sep="")
  ggsave(filename=filename, plot, height=4, width=4, useDingbats=F)
  rm(df,plot)
}

if(file.exists("results/up_metaplot.tsv")){
  if(file.exists("results/down_metaplot.tsv")){
    print("Making combined metaplot")
    up <- read.table("results/up_metaplot.tsv",header=T,sep="\t")
    down <- read.table("results/down_metaplot.tsv",header=T,sep="\t")
    no <- read.table("results/no_diff_exp_metaplot.tsv",header=T,sep="\t")
    df <- data.frame(Bin=up$Bin,up_CG=up$mCG,up_CH=up$mCH,down_CG=down$mCG,
                     down_CH=down$mCH,no_CG=no$mCG,no_CH=no$mCH)
    plot <- ggplot(df, aes(x=Bin)) + geom_line(aes(y=up_CG), color="hotpink4", size=0.8) +
                geom_line(aes(y=down_CG), color="dodgerblue4", size=0.8) +
                geom_line(aes(y=no_CG), color="grey54", size=0.8) +
                theme(panel.background=element_blank(), panel.grid=element_blank(),
                axis.text.y=element_text(color="black"), axis.text.x=element_text(color="black"),
                axis.ticks=element_line(color="black"), axis.title=element_text(color="black"),
                legend.position="none", axis.line=element_line(color="black")) +
                ylab("Percent methylation") + xlab( "" ) +
                scale_y_continuous(limits=c(0,1), expand=c(0,0),
                breaks=c(0.25,0.5,0.75,1),labels=c("25%","50%","75%","100%")) +
                geom_vline(xintercept=20, linetype="longdash", color="grey55") +
                geom_vline(xintercept=40, linetype="longdash", color="grey55")
    filename=paste("figures_tables/", args[1], "_combined_CG_metaplot.pdf", sep="")
    ggsave(filename=filename, plot, height=4, width=4, useDingbats=F)
    plot <- ggplot(df, aes(x=Bin)) + geom_line(aes(y=up_CH), color="hotpink4", size=0.8) +
                geom_line(aes(y=down_CH), color="dodgerblue4", size=0.8) +
                geom_line(aes(y=no_CH), color="grey54", size=0.8) +
                theme(panel.background=element_blank(), panel.grid=element_blank(),
                axis.text.y=element_text(color="black"), axis.text.x=element_text(color="black"),
                axis.ticks=element_line(color="black"), axis.title=element_text(color="black"),
                legend.position="none", axis.line=element_line(color="black")) +
                ylab("Percent methylation") + xlab( "" ) +
                scale_y_continuous(limits=c(0,0.2), expand=c(0,0),
                breaks=c(0.05,0.1,0.15,0.2),labels=c("5.0%","10%","15%","20%")) +
                geom_vline(xintercept=20, linetype="longdash", color="grey55") +
                geom_vline(xintercept=40, linetype="longdash", color="grey55")
    filename=paste("figures_tables/", args[1], "_combined_CH_metaplot.pdf", sep="")
    ggsave(filename=filename, plot, height=4, width=4, useDingbats=F)
  }
}
