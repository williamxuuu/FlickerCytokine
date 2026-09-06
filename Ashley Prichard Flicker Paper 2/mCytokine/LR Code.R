setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


library('xlsx')
library('fastcluster')
library("gplots")
library("dendextend")
library(reshape)
library(matrixStats)
library("readxl")
library("heatmap3")
library("compare")
library("grid")
library("gridGraphics")
#library("DESeq2")
library("ggplot2")
library("scatr")
library("ggpubr")
#library("GSVA")
library("RColorBrewer")
library(boot)
library(caret)
library(hrbrthemes)
library("ClassDiscovery")
#library("tidyverse")
#library("plsRglm")
library("plsdof")
library("gplots")
library("lme4")
library("foreign")
#source("summarySE.R")
library("corrplot")
library("Hmisc")
library("PerformanceAnalytics")
library("dplyr")
library("caret")
library("ggpubr")
dataCytInALL=read_excel("BCA_Analyzed_PrichardSamples_11425_V2.xlsx", sheet="WT", range="B2:AH12")
dataCytBackground=as.matrix(dataCytInALL[dataCytInALL$Sample=="Background0",2:dim(dataCytInALL)[2]])
dataCytIn=dataCytInALL[1:(dim(dataCytInALL)[1]-2),2:dim(dataCytInALL)[2]]
dataCyt=as.matrix(dataCytIn)
for (i in 1:dim(dataCyt)[2])
{
  dataCyt[,i]=dataCyt[,i]-dataCytBackground[,i]
}
dataCyt[dataCyt<0]=0

concentration=dataCytInALL$Sample[1:(dim(dataCytInALL)[1]-2)]
df_var=as.numeric(concentration)
x_regress=df_var
x_label="Protein Concentration"

for(i in 1:ncol(dataCyt)){
  
  y_regress=dataCyt[,i] 
  y_label=colnames(dataCytIn)[i]
  
  
  dataReg=data.frame(y_regress,x_regress)
  linmodel=lm(y_regress~x_regress,dataReg)
  # RSQ=summary(linmodel)$r.squared
  # f <- summary(linmodel)$fstatistic
  # pval <- pf(f[1],f[2],f[3],lower.tail=F)
  # attributes(pval) <- NULL
  
  outputFolder<-"./WT/"
  
  #Plot regression
  
  
  png(paste(outputFolder,"/",y_label,' vs ',x_label,".png",sep= ""), width=4,height=4,units="in",res=600,pointsize =10)
  print({
    ggplot(dataReg, aes(x=x_regress, y=y_regress)) +
      geom_smooth(method='',se=FALSE,color='red',linetype="longdash")+
      geom_point(size=4 , color='royalblue')+
      # xlim(0.95*min(x_regress),1.05*max(x_regress))+
      # ylim(0.95*min(y_regress),1.05*max(y_regress))+
      theme(panel.background = element_rect(fill = 'white'),
            text = element_text(size=12),
            axis.text.x = element_text(color = "black",
                                       size = 14),
            axis.text.y = element_text(color = "black",
                                       size = 14),
            plot.title = element_text(hjust = 0.5,face='bold'),
            panel.border=element_blank(),
            axis.line=element_line(color='black'),
            plot.margin=margin(10,20,10,10)
      )+
      xlab(x_label)+
      ylab(y_label)+
     
      ggtitle(paste(y_label," vs ",x_label,sep=""))
  })
  dev.off()
  
}
