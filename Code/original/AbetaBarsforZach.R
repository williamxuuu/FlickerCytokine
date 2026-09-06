# PRELIMINARIES ################################

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

source("WoodLabFunctions.R")
source("summarySE.R")
if (!require("pacman")) install.packages("pacman")
if (!require("BiocManager")) install.packages("BiocManager")
library("readxl")
library("tidyverse")
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
library("corrplot")
library("Hmisc")
library("PerformanceAnalytics")
library("dplyr")
library("caret")
library("ggpubr")
library("ggpmisc")
library("RColorBrewer")

library("SAMtool")
source("summarySE.R")



p_load(readxl,tidyverse,ropls,ggpubr)

filterCyt = FALSE
breakBarColors=c(-200,seq(-1.5, 1.5, 0.01),200)
barColors = colorpanel(length(breakBarColors)-1, "blue", "white", "red2")

dataCytInALL=read_excel("030625_mCytokineZach_copy.xlsx", sheet="Female Cytokine 5xFAD", range="D1:AV15")
dataCytBackground=as.matrix(dataCytInALL[dataCytInALL$Sample=="Background0",14:dim(dataCytInALL)[2]])
dataCytIn=dataCytInALL[1:(dim(dataCytInALL)[1]-2),14:dim(dataCytInALL)[2]]
dataCyt=as.matrix(dataCytIn)

#Convert cytokines to matrix and subtract background
dataCyt=as.matrix(dataCytIn)
for (i in 1:dim(dataCyt)[2])
{
  dataCyt[,i]=dataCyt[,i]-mean((dataCytBackground[,i]))
} 
dataCyt[dataCyt<0]=0

#### CODE BELOW IS HOW I NORMALIZE BASED ON NO FLICKER GROUP 

#This makes the matrix into a dataframe
dataCytWithGroup=as.data.frame(dataCyt)
#This adds a column with the groups to the dataframe
dataCytWithGroup['Flicker'] <- dataCytInALL$Flicker[1:((nrow)(dataCytInALL)-2)]



# Susbet Condition is what we are normalizing to 
subset_condition <- dataCytWithGroup$Flicker == "NoFlicker"

# Function to calculate statistics based on a subset
calculate_subset_stats <- function(data, condition) {
  # Get the subset data
  subset_data <- data[condition, ]
  
  # Get only numeric columns (excluding grouping variables)
  numeric_cols <- sapply(subset_data, is.numeric)
  
  # Calculate stats for each numeric column
  stats_df <- data.frame(
    column = names(subset_data)[numeric_cols],
    mean = sapply(subset_data[, numeric_cols], mean),
    sd = sapply(subset_data[, numeric_cols], sd),
    stringsAsFactors = FALSE
  )
  
  return(stats_df)
}

# Calculate stats for the subset
subset_stats <- calculate_subset_stats(dataCytWithGroup, subset_condition)
print(subset_stats)

# Function to apply normalization using the subset statistics to the ENTIRE dataset
normalize_with_subset_stats <- function(data, stats) {
  result <- data  # Create a copy
  
  for(col in stats$column) {
    # Get the mean and sd from our subset statistics
    col_mean <- stats$mean[stats$column == col]
    col_sd <- stats$sd[stats$column == col]
    
    # Apply normalization to the entire column
    result[[col]] <- (data[[col]] - col_mean) / col_sd
  }
  
  return(result)
}

# Apply the normalization using statistics from only group "A"
normalized_df <- normalize_with_subset_stats(dataCytWithGroup, subset_stats)

#Get rid of the group marker...
print(normalized_df)
dataCytZ <- normalized_df[, -ncol(normalized_df)]

#### THE END OF THE CODE FOR NORMLAIZING BASED ON NO FLICKER 

rownames(dataCytZ)=dataCytInALL$Sample[1:(dim(dataCytInALL)[1]-2)]
#dataCytZ=(apply(t(dataCyt),1,scale)); #z-score the data
dataCytZ[is.na(dataCytZ)]=0
rownames(dataCytZ)=rownames(dataCyt)



Condition=dataCytInALL$Condition[1:(dim(dataCytInALL)[1]-2)]
#Sex=dataCytInALL$Sex[1:(dim(dataCytInALL)[1]-2)]
#dose=dataCytInALL$dose[1:(dim(dataCytInALL)[1]-2)]
#genotype=dataCytInALL$genotype[1:(dim(dataCytInALL)[1]-2)]
ab42=dataCytInALL$Abeta42[1:(dim(dataCytInALL)[1]-2)]
ab40=dataCytInALL$Abeta40[1:(dim(dataCytInALL)[1]-2)]
ratio=dataCytInALL$AbetaRatio[1:(dim(dataCytInALL)[1]-2)]
#Sex_genotype=dataCytInALL$"Sex_genotype"[1:(dim(dataCytInALL)[1]-2)]
#Sample=dataCytInALL$Sample[1:(dim(dataCytInALL)[1]-2)]
#Treatment=dataCytInALL$Treatment[1:(dim(dataCytInALL)[1]-2)]
#Sort based on xxxx
tpSort=sort(Condition, decreasing=FALSE, index.return=TRUE)
# ix based on top to bottom 
indSort=tpSort$ix
dataCyt=dataCyt[indSort,]
dataCytZ=dataCytZ[indSort,]
#Define a group
Condition=Condition[indSort]
# Sex=Sex[indSort]
# dose=dose[indSort]
# genotype=genotype[indSort]
ab42=ab42[indSort]
ab40=ab40[indSort]
ratio=ratio[indSort]
# Sex_genotype=Sex_genotype[indSort]
# Sample=Sample[indSort]
# Treatment=Treatment[indSort]
#Setup color bars
# So instead of sex, genotype, do GROUPS instead of Condition... and assign colors to each one of them 
colorBarCondition=matrix(nrow = length(Condition),1)
colorBarCondition[Condition==1]="#E298D6"
colorBarCondition[Condition==2]="#1A1513"
colorBarCondition[Condition==3]="#76A33E"
#colorBarCondition[Condition==4]="darkgoldenrod4"
#colorBarCondition[Condition==5]="mediumorchid4"
#colorBarCondition[Condition==6]="chartreuse4"

#colorBarCondition[Condition==7]="black"
#colorBarCondition[Condition==8]="orangered3"
#colorBarCondition[Condition==9]="forestgreen"
# colorBarCondition[Condition==10]="darkturquoise"
# colorBarCondition[Condition==11]="royalblue"
# colorBarCondition[Condition==12]="maroon2"
# 
# 
# colorBarSex=matrix(nrow = length(Sex),1)
# colorBarSex[Sex=="female"]="hotpink"
# colorBarSex[Sex=="male"]="deepskyblue2"
# 
# colorBardose=matrix(nrow = length(dose),1)
# colorBardose[dose=="4"]="darkred"
# colorBardose[dose=="1"]="darkgreen"
# colorBardose[dose=="0"]="SNOW4"
# 
# colorBarGenotype=matrix(nrow = length(genotype),1)
# colorBarGenotype[genotype=="FAD"]="lightslateblue"
# colorBarGenotype[genotype=="WT"]="olivedrab2"


breakBarColorsAbeta42=c(-1,seq(0, 1200, 0.1))*10 #Outside numbers clip outliers. This is for zscoring.
barColorsAbeta42Vec = c("black",colorpanel(length(breakBarColorsAbeta42), "white", "gold4"))
indMatchColorone=match(round(ab42*10),round(breakBarColorsAbeta42))
barColorsAb42=barColorsAbeta42Vec[indMatchColorone]

breakBarColorsAbeta40=c(-1,seq(0, 15, 0.01))*10 #Outside numbers clip outliers. This is for zscoring.
barColorsAbeta40Vec = c("black",colorpanel(length(breakBarColorsAbeta40), "white", "orange"))
indMatchColortwo=match(round(ab40*10),round(breakBarColorsAbeta40))
barColorsAb40=barColorsAbeta40Vec[indMatchColortwo]


breakBarColorsAbetaratio=c(-1,seq(0, 100, 0.01))*10 #Outside numbers clip outliers. This is for zscoring.
barColorsAbetaratioVec = c("black",colorpanel(length(breakBarColorsAbetaratio), "white", "gray4"))
indMatchColorthree=match(round(ratio*10),round(breakBarColorsAbetaratio))
barColorsAbratio=barColorsAbetaratioVec[indMatchColorthree]


# I will just have group and Abeta..... 
colorBarCombined=cbind(c(colorBarCondition),c(barColorsAb42),c(barColorsAb40),c(barColorsAbratio))
colnames(colorBarCombined)=c("Condition","Ab42","Ab40","ratio")


#Cluster the cytokine data by columns
hc= hclust(dist(t(dataCytZ),method = "euclidean"), method="average")

pdf("NORAMLIZEDFemale5xFADOnlyWithAbeta.pdf", width=8,height=8,pointsize = 6)
heatmap3(dataCytZ,main="Female 5xFAD Cytokine Profile",col=barColors, breaks=breakBarColors,
         Rowv=NA, Colv=as.dendrogram(hc),  scale="none", na.rm = TRUE,
         labRow = rownames(dataCytZ), labCol = colnames(dataCytZ),
         RowSideColors=colorBarCombined)


dev.off()


dataCytOut=data.frame(Sample,Condition,Treatment,Sex_genotype,dataCyt)


write.xlsx(dataCytOut, file="ignore.xlsx", sheetName="sheet1", col.names=TRUE, row.names=FALSE, append=FALSE)



#----Make bar plots Temporal



dataCyt_groupVars=data.frame(dataCyt,dose,genotype)
Sex_genotype=as.factor(Sex_genotype)
dose=as.factor(dose)
Sex=as.factor(Sex)
genotype=as.factor(genotype)

for (ind in 1:44)
{
  dataSum=summarySE(data.frame(dataCyt_groupVars), measurevar=colnames(dataCyt_groupVars)[ind], groupvars = c("dose","genotype"))
  
  png(paste("Analyte_",ind,".png",sep= ""), width=6,height=4,units="in",res=600, pointsize = 20)
  print({
    ggplot(dataSum, aes(x=dose, y=mean, color=genotype)) + 
      geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2) +
      geom_line(size=1) +
      geom_point()+
      #geom_jitter(aes(colour=genotype),width=0.25)+ #Can remove
      theme(panel.background = element_rect(fill = 'white', colour = 'black'),
            text = element_text(size=20),
            axis.text.x = element_text(face = "bold", color = "black",
                                       size = 16),
            axis.text.y = element_text(face = "bold", color = "black",
                                       size = 16),
            plot.title = element_text(hjust = 0.5))+
     # scale_color_manual(values=c("deeppink3","orchid3","deepskyblue3","darkgreen"))+
      xlab("LPS dose (Day)")+
      ylab("Flourescent Int (a.u.)")+
      labs(color="genotype")+
      ggtitle(colnames(dataCyt_groupVars)[ind])
    
  })
  dev.off()
}


######
#PCA analysis :


# All Samples PCA ####

variabletochange=dataCytInALL$'Treatment'[1:(dim(dataCytInALL)[1]-2)]
variabletochange=as.character(variabletochange)
variabletochange=variabletochange[indSort]
pheno<-as.factor(variabletochange)



PCA_folder = "New 041121"
#change here
df_samples<-dataCyt
df_PCA <- df_samples


df_PCA[,] <- apply(df_PCA[,],2,"scale")
#df_PCA[is.na(df_PCA)] <- 0
PCA_mat <- df_PCA[,] %>% as.matrix()

oplsOut <- opls(PCA_mat,predI=2)
rotateOut <- rotate.opls(oplsOut,degrees=0)


png(paste0("./Figures/PCA/",PCA_folder,"/PCA.png"),height=3,width=6,units="in",res=600)
scores.plot(rotateOut$T1,rotateOut$T2,pheno,"condition",analysis="PCA")

dev.off()


png(paste0("./Figures/PCA/",PCA_folder,"/p1.png"),height=3,width=6,units="in",res=600)
scores.plot(rotateOut$T1,rotateOut$T2,pheno,"condition",analysis="PCA")+
  xlab(paste0("Scores on PC1"))+
  ylab(paste0("Scores on PC2"))
#scale_color_manual(values=c("deeppink3","orchid3","deepskyblue3","darkgreen"))
dev.off()



png(paste0("./Figures/PCA/",PCA_folder,"/PCA_PC1.png"),height=3,width=6,units="in",res=600)
loadings.chart(rotateOut$P1,"PC1",cex=0.5)
dev.off()

sp <- scores.plot(rotateOut$T1,rotateOut$T2,pheno,"Condition",analysis="PCA")
x1 <- loadings.chart2(rotateOut$P1,"PC1",x_size = 9)
x2 <- loadings.chart2(rotateOut$P2,"PC2",x_size = 9)
#X3 <-loadings.chart2(rotateOut$P3,"PC3",x_size = 9)
#X4 <-loadings.chart2(rotateOut$P4,"PC4",x_size = 9)

ggarrange(sp,ggarrange(x1,x2,now=2,ncol=1),
          ncol=2,nrow=1)

png(paste0("./Figures/PCA/",PCA_folder,"/PCA_PC2.png"),height=3,width=6,units="in",res=600)
loadings.chart(rotateOut$P2,"PC2",cex=0.5)
dev.off()


png(paste0("./Figures/PCA/",PCA_folder,"/PCA_PC1_dot.png"),height=5,width=8,units="in",res=600)
ggpar(
  dotplot.mean.sem(x=pheno,y=rotateOut$T1,x.str="",y.str="PC1",title="PC1"),
  legend="none"
)+
  scale_color_manual(values=c("deeppink3","orchid3","deepskyblue3","darkgreen"))
dev.off()

png(paste0("./Figures/PCA/",PCA_folder,"/PCA_PC2_dot.png"),height=5,width=8,units="in",res=600)
ggpar(
  dotplot.mean.sem(x=pheno,y=rotateOut$T2,x.str="",y.str="PC2",title="PC2"),
  legend="none"
)+scale_color_manual(values=c("deeppink3","orchid3","deepskyblue3","darkgreen"))
dev.off()



#No Lines

png(paste0("./Figures/PCA/",PCA_folder,"/PCA_PC2withNoLine_dot.png"),height=5,width=8,units="in",res=600)
ggpar(
  errorplot.meanSE.compare(x=pheno,y=rotateOut$T2,colors=pheno,title="PC2 no line",
                           p.display = "p.format",
                           my_comparisons = list(c(1,2),c(2,3),c(1,3)))
)+scale_color_manual(values=c("deeppink3","orchid3","deepskyblue3","darkgreen"))

dev.off()

png(paste0("./Figures/PCA/",PCA_folder,"/PCA_PC1withNoLine_dot.png"),height=5,width=8,units="in",res=600)
ggpar(
  errorplot.meanSE.compare(x=pheno,y=rotateOut$T1,colors=pheno,title="PC1 no line",
                           p.display = "p.format",
                           my_comparisons = list(c(1,2),c(2,3),c(1,3)))
)+scale_color_manual(values=c("deeppink3","orchid3","deepskyblue3","darkgreen"))

dev.off()

png(paste0("./Figures/PCA/",PCA_folder,"/PCA_PC2withNoLine_dot.png"),height=5,width=8,units="in",res=600)
ggpar(
  errorplot.meanSE.compare(x=pheno,y=rotateOut$T2,colors=pheno,title="PC2 no line",
                           p.display = "p.format",
                           my_comparisons = list(c(1,2),c(2,3),c(1,3)))
)+scale_color_manual(values=c("deeppink3","orchid3","deepskyblue3","darkgreen"))

dev.off()


#----Make bar plots Temporal



dataCyt_groupVars=data.frame(dataCyt,dose,Sex_genotype)
Sex_genotype=as.factor(Sex_genotype)
dose=as.factor(dose)
Sex=as.factor(Sex)

for (ind in 1:44)
{
  dataSum=summarySE(data.frame(dataCyt_groupVars), measurevar=colnames(dataCyt_groupVars)[ind], groupvars = c("dose","Sex_genotype"))
  
  png(paste("Analyte_",ind,".png",sep= ""), width=6,height=4,units="in",res=600, pointsize = 20)
  print({
    ggplot(dataSum, aes(x=dose, y=mean, color=Sex_genotype)) + 
      geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2) +
      geom_line(size=1) +
      geom_point()+
      geom_jitter(aes(colour=Sex_genotype),width=0.25)+ #Can remove
      theme(panel.background = element_rect(fill = 'white', colour = 'black'),
            text = element_text(size=20),
            axis.text.x = element_text(face = "bold", color = "black",
                                       size = 16),
            axis.text.y = element_text(face = "bold", color = "black",
                                       size = 16),
            plot.title = element_text(hjust = 0.5))+
      scale_color_manual(values=c("deeppink3","orchid3","deepskyblue3","darkgreen"))+
      xlab("LPS dose (Day)")+
      ylab("Flourescent Int (a.u.)")+
      labs(color="Sex_genotype")+
      ggtitle(colnames(dataCyt_groupVars)[ind])
    
  })
  dev.off()
}

