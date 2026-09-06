
#Sara Bitarafan March 2025#
#---Ashely Prichard---Flicker Papper 2----Analaysis--#

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("WoodLabFunctions.R")
source("summarySE.R")
if (!require("pacman")) install.packages("pacman")
if (!require("BiocManager")) install.packages("BiocManager")

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("WGCNA")
# 
# BiocManager::install("htmlwidgets")
#if (!require("WGCNA")) install.packages("WGCNA")
# Load contributed packages with pacman
pacman::p_load(pacman,rio,tidyverse,readxl,matrixStats,ggpubr,heatmap3,Rtsne,limma,
               gplots,GSVA,ropls,RColorBrewer,ggsci, EnhancedVolcano,DESeq2,readxl,xlsx,ggvenn,grid,forcats,viridis,org.Mm.eg.db,ggbeeswarm,WGCNA,impute,preprocessCore) # install/load these packages
# pacman: for package handling
# rio: for smart import/export functions
# tidyverse: for handling tibble data frames and data wrangling
# readxl: for excel importing
# matrixStats: for matrix operations
# ggpubr: for publication-ready figures
# heatmap3: our favorite heatmap function
# Rtsne: t-distributed stochastic neighbor embedding
# gplots: more graphical things for figures
# GSVA: gene set variation analysis, Bioconductor
# ropls: for PCA, PLSR, PLSDA, etc., Bioconductor
# RColorBrewer: for color palettes
# ggsci: scientific color palettes
# Enhanced Volcano: volcano plots, Bioconductor
# DESeq2: for normalizing RNAseq read counts, Bioconductor

mypal <- pal_npg()(10)   # a palette from ggsci 


breakBarColors=c(-200,seq(-1.5, 1.5, 0.01),200)
barColors = colorpanel(length(breakBarColors)-1, "blue", "white", "red2")

dfInALL=read_excel("Metadata_mCytokine_AmyloidBeta_MicroglialMorphology_SBversionv2.xlsx", sheet="Simplified", range="A2:BP47")
dfBackground=as.matrix(dfInALL[dfInALL$Sample=="Background0",14:dim(dfInALL)[2]])
dfIn=dfInALL[1:(dim(dfInALL)[1]-2),14:dim(dfInALL)[2]]
df=as.matrix(dfIn)

#Convert cytokines to matrix and subtract background
df=as.matrix(dfIn)
for (i in 1:dim(df)[2])
{
  df[,i]=df[,i]-mean((dfBackground[,i]))
} 
df[df<0]=0


Group=dfInALL$`Group`[1:(dim(dfInALL)[1]-2)]
Genotype=dfInALL$Genotype[1:(dim(dfInALL)[1]-2)]
Sex=dfInALL$Sex[1:(dim(dfInALL)[1]-2)]
Age=dfInALL$Age[1:(dim(dfInALL)[1]-2)]
Flicker=dfInALL$Flicker[1:(dim(dfInALL)[1]-2)]
DOB=dfInALL$DOB[1:(dim(dfInALL)[1]-2)]
Lit=dfInALL$litter[1:(dim(dfInALL)[1]-2)]

Exp_data<-dfInALL[-c(44,45), c(1:13)]

Exp_data<- Exp_data %>%
  mutate(Condition = paste0(Genotype," ",Sex , " ", Flicker))

Exp_data<- Exp_data %>%
  mutate(SexGenotype = paste0(Genotype," ",Sex))


dataCyt=df[,-c(1:24)]
rownames(dataCyt)<-dfInALL$`Sample Number`[1:(dim(dfInALL)[1]-2)]
dataAbeta=df[, c(19:23)]
rownames(dataAbeta)<-dfInALL$`Sample Number`[1:(dim(dfInALL)[1]-2)]

dataMorph=df[, c(1:18)]
rownames(dataMorph)<-dfInALL$`Sample Number`[1:(dim(dfInALL)[1]-2)]


Df_CytZ=(apply(t(dataCyt),1,scale)); #z-score the data
Df_CytZ[is.na(Df_CytZ)]=0
rownames(Df_CytZ)=rownames(dataCyt)






# ---------- Outlier detection method -------------------

output_folder=paste0("./OutlierDetection")

#5xFAD Group
sortby <- which(Exp_data$Genotype =="5XFAD")

Df_Cyt_Sub <-  dataCyt[sortby ,] 
Df_Abeta_Sub <-dataAbeta[sortby ,]
Exp_data_sub <- Exp_data[sortby ,]

indSort=sort(Exp_data_sub$Group, index.return=TRUE)$ix
Group=Group[indSort]
Df_Cyt_Sub=Df_Cyt_Sub[indSort,]
Df_Cyt_Sub=Df_Cyt_Sub[, -c(1:8)]
Df_Abeta_Sub=Df_Abeta_Sub[indSort,]
Exp_data_sub=Exp_data_sub[indSort,]


#Z scoring data

Df_Cyt_SubZ=(apply(t(Df_Cyt_Sub),1,scale)); #z-score the data
Df_Cyt_SubZ[is.na(Df_Cyt_SubZ)]=0
rownames(Df_Cyt_SubZ)=rownames(Df_Cyt_Sub)


df_PCA <- Df_Abeta_Sub
PCA_mat <- df_PCA[,] %>% as.matrix()

oplsOut <- opls(PCA_mat,predI=2)
rotateOut <- rotate.opls(oplsOut,degrees=0)
scores_data=data.frame(rotateOut$T1,rotateOut$T2)
# 

p<- scores.plot(rotateOut$T1,rotateOut$T2,factor(Exp_data_sub$Flicker),"Flicker")
#ggpar(p,palette=c("#5C81BD","#3E2A47","#E8A941"))
print(p)




df_plot = data.frame(T1=rotateOut$T1,T2=rotateOut$T2)

Groupss= df_plot[,c("T1" , "T2")]
rownames(Groupss)=rownames(df_plot)
Group.center= colMeans(Groupss)
Group.cov= cov(Groupss)
rad=qchisq(p=0.95, df=ncol(Groupss))
rad=sqrt(rad)
ellipse <- car::ellipse(center = Group.center , shape = Group.cov , radius = rad ,
                        segments = 150 , draw = FALSE)

ellipse <- as.data.frame(ellipse)
colnames(ellipse) <- colnames(Groupss)



png(paste0(output_folder,"OutlierDetection 5xFAD.png"),height=4,width=4,units="in",res=600)
figure <- ggplot(Groupss , aes(x = T1 , y = T2)) +
  geom_point(size = 2) +
  geom_polygon(data = ellipse , fill = "orange" , color = "orange" , alpha = 0.5)+
  geom_point(aes(Group.center[1] , Group.center[2]) , size = 5 , color = "blue") +
  geom_text( aes(label = row.names(Groupss)) , hjust = 1 , vjust = -1.5 ,size = 2.5 ) +
  ylab("PC2 scores") + xlab("PC1 scores")+ggtitle("PCA 5xFAD ")
dev.off()
figure


# Finding distances
distances <- mahalanobis(x = Groupss , center = Group.center , cov = Group.cov)

# Cutoff value for distances from Chi-Sqaure Dist. 
# with p = 0.95 df = 2 which in ncol
cutoff <- qchisq(p = 0.95 , df = ncol(Groupss))

## Display observation whose distance greater than cutoff value
Groupss[distances > cutoff ,]


df_PCA <- Df_CytZ
PCA_mat <- df_PCA[,] %>% as.matrix()

oplsOut <- opls(PCA_mat,predI=2)
rotateOut <- rotate.opls(oplsOut,degrees=0)
scores_data=data.frame(rotateOut$T1,rotateOut$T2)
# 


pdfPath <- paste("./Aug25/ScorePlot_opls_all_Litter.pdf", sep = '')
pdf(pdfPath, width = 5, height = 4, pointsize = 14)
p<- scores.plot(rotateOut$T1,rotateOut$T2,factor(Exp_data$litter),"Litter")
#ggpar(p,palette=c("#5C81BD","#3E2A47","#E8A941"))
print(p)
dev.off()




#----PCA---
df_PCA <- Df_CytZ
PCA_mat <- df_PCA[,] %>% as.matrix()

oplsOut <- opls(PCA_mat,predI=2)
rotateOut <- rotate.opls(oplsOut,degrees=0)
scores_data=data.frame(rotateOut$T1,rotateOut$T2,Exp_data$Condition)
Lit=factor(Exp_data$litter)

oplsVar = oplsOut@modelDF$R2X * 100



df_plot = data.frame(T1=rotateOut$T1,T2=rotateOut$T2,Condition=Exp_data$Condition)


statsOut_T1 = summarySE(df_plot,measurevar="T1", groupvars="Condition")
statsOut_T2 = summarySE(df_plot,measurevar="T2", groupvars="Condition")


pdf(paste0("./Aug25/PC1_All.pdf"),height=4,width=6,useDingbats = FALSE)
ggplot(data=df_plot, aes( y = T1, x = Condition, fill = Condition)) +
  geom_violin()+
  geom_point()+
  
 # stat_summary(fun = mean, geom = "point", shape = 18, size = 3, color = "black", position = position_dodge(width = 0.9)) +
  #stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.2, color = "black", position = position_dodge(width = 0.9)) +
  
  geom_errorbar(data=statsOut_T1, aes(x=Condition , ymin=mean-se, ymax=mean +se),
  width=0.4, colour="black", alpha=1, size=1) +
 # geom_crossbar(data=df_plot, aes(x=Condition ,ymin = T1, ymax = T1),
 # size=0.5,col="red", width = .5)+
  theme(panel.background = element_rect(fill = 'white', colour = 'black'),
        text = element_text(size=20), 
        axis.text.x = element_text(face = "plain", color = "black",
                                   size = 10,angle=45,hjust=1),
        axis.text.y = element_text(face = "plain", color = "black",
                                   size = 10),
        legend.position = "none", plot.title = element_text(hjust = 0.5))+
  xlab("")+
  ylab(paste0("Scores on PC1 (",oplsVar[1],"%)"))
  #ylim(c(-5, 5))+
 # scale_fill_manual(values = c("#5C81BD","#3E2A47","#E8A941"))+
 # stat_compare_means(method="anova",p.adjust.methods = "bonferroni",label="p.format",tip.length = 0.03)
#coord_flip()
dev.off()



pdf(paste0("PC2_all.pdf"),height=4,width=6,useDingbats = FALSE)
ggplot(data=df_plot, aes( y = T2, x = Condition, fill = Condition)) +
  geom_point()+
  
 # stat_summary(fun = mean, geom = "point", shape = 18, size = 3, color = "black", position = position_dodge(width = 0.9)) +
#  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.2, color = "black", position = position_dodge(width = 0.9)) +
  
  #geom_errorbar(data=statsOut_T1, aes(x=Condition, ymin=mean-se, ymax=mean +se),
  # width=0.4, colour="black", alpha=1, size=1) +
  #geom_crossbar(data=statsOut_T1, aes(x=Condition,ymin = T1, ymax = T1),
  #  size=0.5,col="red", width = .5)+
  theme(panel.background = element_rect(fill = 'white', colour = 'black'),
        text = element_text(size=20), 
        axis.text.x = element_text(face = "plain", color = "black",
                                   size = 10,angle=45,hjust=1),
        axis.text.y = element_text(face = "plain", color = "black",
                                   size = 10),
        legend.position = "none", plot.title = element_text(hjust = 0.5))+
  xlab("")+
  ylab(paste0("Scores on PC2 (",oplsVar[2],"%)"))+
 # ylim(c(-5, 5))+
  #scale_fill_manual(values = "#E298D6","#1A1513","#76A33E")+
 # scale_fill_manual(values = c("#5C81BD","#3E2A47","#E8A941"))+
 # stat_compare_means(method="anova",p.adjust.methods = "bonferroni",label="p.format",tip.length = 0.03)
#coord_flip()
dev.off()



pdfPath <- paste("./ScorePlot_PCA_All.pdf", sep = '')
pdf(pdfPath, width = 10, height = 8, pointsize = 14)

png("ScorePlot_PCA_All_litter.png", width=16,height=16,units="in",res=600, pointsize = 20)
p <- scores.plot(rotateOut$T1,rotateOut$T2,Exp_data$litter,"litter",analysis = "PCA")
#ggpar(p,palette=c("#5C81BD","#3E2A47","#E8A941"))

print(p)
dev.off()



GroupColorBarF <- case_when(
  Exp_data$Condition == "WT F NoFlicker" ~ "#E8AA5B",
  Exp_data$Condition == "WT F 20Hz" ~ "#6BACDB",
  Exp_data$Condition == "WT F 40Hz" ~ "#D27FC2",
  Exp_data$Condition == "5XFAD F NoFlicker" ~  "#B76945",
  Exp_data$Condition == "5XFAD F 20Hz"~ "#4E5491",
  Exp_data$Condition == "5XFAD F 40Hz" ~ "#863775"
)


GroupColorBarM <- case_when(
  Exp_data$Condition == "WT M NoFlicker" ~ "#E8AA5B",
  Exp_data$Condition == "WT M 20Hz" ~ "#6BACDB",
  Exp_data$Condition == "WT M 40Hz" ~ "#D27FC2",
  Exp_data$Condition == "5XFAD M NoFlicker" ~  "#B76945",
  Exp_data$Condition == "5XFAD M 20Hz"~ "#4E5491",
  Exp_data$Condition == "5XFAD M 40Hz" ~ "#863775"
)





GroupColorBarM <- case_when(
 Exp_data$Condition == "WT M NoFlicker" ~ "#E8AA5B",
 Exp_data$Condition == "WT M 20Hz" ~ "#6BACDB",
  Exp_data$Condition == "WT M 40Hz" ~ "#D27FC2",
  Exp_data$Condition == "5XFAD M NoFlicker" ~  "#B76945",
  Exp_data$Condition == "5XFAD M 20Hz"~ "#4E5491",
  Exp_data$Condition == "5XFAD M 40Hz" ~ "#863775"
)





#---- Select a Genotype or Group (5xFAD vs WT)---#

sortby <- which(Exp_data$Genotype =="WT")

Df_Cyt_Sub <-  dataCyt[sortby ,] 
Df_Abeta_Sub <-dataAbeta[sortby ,]
Exp_data_sub <- Exp_data[sortby ,]
Df_Morph_sub <- dataMorph [sortby ,]

# RemovingOutlierSamples (samples 2 (cher4) ,5(dolly1) ,8(xena3) are outliers)
# Df_Cyt_Sub2<- Df_Cyt_Sub[-c(5), ]
# Df_Abeta_Sub2 <-Df_Abeta_Sub[-c(5),]
# Exp_data_sub2 <- Exp_data_sub[-c(5),]

Flicker=Exp_data_sub$Flicker
 indSort=sort(Flicker, index.return=TRUE)$ix
 Flicker=Flicker[indSort]
 Df_Cyt_Sub2=Df_Cyt_Sub[indSort ,]
 Df_Abeta_Sub2=Df_Abeta_Sub[indSort ,]
 Df_Morph_Sub2=Df_Morph_sub[indSort ,]
Exp_data_sub2=Exp_data_sub[indSort ,]






#Z scoring data

Df_Cyt_SubZ=(apply(t(Df_Cyt_Sub2),1,scale)); #z-score the data
Df_Cyt_SubZ[is.na(Df_Cyt_SubZ)]=0
rownames(Df_Cyt_SubZ)=rownames(Df_Cyt_Sub2)



Litterpal<-brewer.pal(8,"Paired")
LitterColorBar <- case_when(
  Exp_data_sub2$litter == "a" ~ Litterpal[1],
  Exp_data_sub2$litter == "b" ~ Litterpal[2],
  Exp_data_sub2$litter == "c" ~ Litterpal[3],
  Exp_data_sub2$litter == "d" ~ Litterpal[4],
  Exp_data_sub2$litter == "e" ~ Litterpal[5],
  Exp_data_sub2$litter == "f" ~ Litterpal[6],
  Exp_data_sub2$litter == "g" ~ Litterpal[7],
)


# FlickerColorBar <- case_when(
#   Exp_data_sub$Flicker == "20Hz" ~ "black",
#   Exp_data_sub$Flicker == "40Hz" ~ "#76A33E",
#   Exp_data_sub$Flicker == "NoFlicker" ~ "#E298D6"
# )


GroupColorBarF <- case_when(
  Exp_data_sub2$Flicker == "NoFlicker" ~ "#E8AA5B",
  Exp_data_sub2$Flicker == "20Hz" ~ "#6BACDB",
  Exp_data_sub2$Flicker == "40Hz" ~ "#D27FC2"
#  Exp_data_sub2$Condition == "5XFAD F NoFlicker" ~  "#B76945",
 # Exp_data_sub2$Condition == "5XFAD F 20Hz"~ "#4E5491",
 # Exp_data_sub2$Condition == "5XFAD F 40Hz" ~ "#863775"
)

GroupColorBarF2 <- case_when(
  Exp_data_sub2$Flicker == "NoFlicker" ~ "#B76945",
  Exp_data_sub2$Flicker == "20Hz" ~ "#4E5491",
  Exp_data_sub2$Flicker == "40Hz" ~  "#863775"
)

GroupColorBarSex <- case_when(
  Exp_data_sub2$Sex == "F" ~ "pink",
  Exp_data_sub2$Sex == "M" ~ "blue"
)


colorBarCombined=cbind(c(GroupColorBarF),c(GroupColorBarSex))
colnames(colorBarCombined)=c("Flicker","Sex")


#Cluster the cytokine data by columns
hc= hclust(dist(t(Df_Cyt_SubZ),method = "euclidean"), method="average")

##Make Heatmap
pdf("./Aug25/WT/Cytokine profile WT sep25.pdf", width=12,height=10,pointsize = 6)
#png("Cytokine profile Female 5xFAD Only.png", width=16,height=16,units="in",res=600, pointsize = 20)
heatmap3(Df_Cyt_SubZ,main="Cytokine Profile",col=barColors, breaks=breakBarColors,legendfun=function()showLegend(),
         Rowv=NA, Colv=as.dendrogram(hc),  scale="none", na.rm = TRUE,
         labRow = rownames(Df_Cyt_SubZ), labCol = colnames(Df_Cyt_SubZ),
         RowSideColors=colorBarCombined)


dev.off()



#Amyloid beta side bar


ab42<-Df_Abeta_Sub2[,2]
ab40<-Df_Abeta_Sub2[,1]
ratio<-Df_Abeta_Sub2[,3]


breakBarColorsAbeta42=c(-1,seq(0, 800, 0.01))*10 #Outside numbers clip outliers. This is for zscoring.
barColorsAbeta42Vec = c("black",colorpanel(length(breakBarColorsAbeta42), "white", "gold4"))
indMatchColorone=match(round(ab42*10),round(breakBarColorsAbeta42))
barColorsAb42=barColorsAbeta42Vec[indMatchColorone]

breakBarColorsAbeta40=c(-1,seq(0, 12, 0.01))*10 #Outside numbers clip outliers. This is for zscoring.
barColorsAbeta40Vec = c("black",colorpanel(length(breakBarColorsAbeta40), "white", "orange2"))
indMatchColortwo=match(round(ab40*10),round(breakBarColorsAbeta40))
barColorsAb40=barColorsAbeta40Vec[indMatchColortwo]


breakBarColorsAbetaratio=c(-1,seq(0, 180, 0.01))*10 #Outside numbers clip outliers. This is for zscoring.
barColorsAbetaratioVec = c("black",colorpanel(length(breakBarColorsAbetaratio), "white", "gold3"))
indMatchColorthree=match(round(ratio*10),round(breakBarColorsAbetaratio))
barColorsAbratio=barColorsAbetaratioVec[indMatchColorthree]

colorBarCombined=cbind(c(GroupColorBarF2),c(barColorsAb42),c(barColorsAb40))
colnames(colorBarCombined)=c("Group","Ab42","Ab40")



##Make Heatmap
pdf("./Aug25/FAD/Cytokine profile  5xFAD with Side Bar for Abeta Sep25 v3.pdf", width=12,height=10,pointsize = 6)
#png("Cytokine profile Female 5xFAD Only_with side bar.png", width=16,height=16,units="in",res=600, pointsize = 20)
heatmap3(Df_Cyt_SubZ,main="Female 5xFAD Cytokine Profile",col=barColors, breaks=breakBarColors,legendfun=function()showLegend(),
         Rowv=NA, Colv=as.dendrogram(hc),  scale="none", na.rm = TRUE,
         labRow = rownames(Df_Cyt_SubZ), labCol = colnames(Df_Cyt_SubZ),
         RowSideColors=colorBarCombined)


dev.off()



Condition=Exp_data_sub2$Flicker

#datagroupVars=data.frame(Df_Cyt_Sub2, Df_Abeta_Sub2, Condition)
datagroupVars=data.frame(Df_Cyt_Sub2, Condition)
datagroupVars$Condition <- factor(datagroupVars$Condition, levels = c("NoFlicker", "20Hz", "40Hz"))  

for (i in 1:3) {
  pdf(paste0("./Aug25/WT/BarPlots_WT_SEP25_TTEST_check", colnames(datagroupVars)[i], ".pdf"), height = 4, width = 5)
  # png(paste0("./Aug25/FAD/Female/BarPlots_Female_40HzandCTRL_",colnames(datagroupVars)[i],".png"),height=5,width=6,units="in",res=600)
  
  print({
    ggplot(datagroupVars, aes(x = Condition, y = datagroupVars[, i], fill = Condition)) +
      
      stat_summary(fun = mean, geom = "bar", position = position_dodge(width = 0.75),
                   width = 0.6, color = "black", alpha=0.4) +
      
      stat_summary(fun.data = mean_se, geom = "errorbar",
                   position = position_dodge(width = 0.75), width = 0.2, stroke=1, color = "black") +
      
      # geom_jitter(width = 0.2, size = 2, stroke=1 , aes(color = as.factor(Condition))) +
      # geom_jitter(width = 0.2, size = 3.5, shape = 21, stroke =  1, color = "black", aes(fill = Condition)) +
      geom_beeswarm(aes(fill = Condition),size =4 ,shape=21, stroke = 1, cex = 2,groupOnX = TRUE)+
      theme(panel.background = element_rect(fill = 'white', colour = 'black'),
            text = element_text(size = 10),
            axis.text.x = element_text(face = "bold", color = "black",
                                       size = 8, angle = 45, hjust = 1),
            axis.text.y = element_text(face = "bold", color = "black", size = 16),
            plot.title = element_text(hjust = 0.5)) +
      #scale_color_manual(values=c("#B76945","#4E5491","#863775"))+
     # scale_fill_manual(values=c("#B76945","#4E5491","#863775"))+
      scale_color_manual(values=c("#E8AA5B","#6BACDB","#D27FC2"))+
      scale_fill_manual(values=c("#E8AA5B","#6BACDB","#D27FC2"))+
      xlab("") + ylab("")+
      labs(fill = "Condition", color = "Condition") +
      ggtitle(colnames(datagroupVars)[i]) +
      stat_compare_means(method="t.test",comparison=list(c(1,3),c(1,2),c(2,3)), label="p.format",tip.length = 0.03,paired = FALSE) 
  })
  dev.off()
}



#----Make color bars----
color.bar <- function(lut, min, max=-min, nticks=11, ticks=seq(min, max, len=nticks), title='') {
  scale = (length(lut)-1)/(max-min)
  
  #dev.new(width=1.75, height=5)
  plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title)
  axis(2, ticks, las=1)
  for (i in 1:(length(lut)-1)) {
    y = (i-1)/scale + min
    rect(0,y,10,y+1/scale, col=lut[i], border=NA)
  }
}


pdf("barColorsAbeta42Vec_ColorBar.pdf", width=2,height=8, pointsize = 14)
color.bar(barColorsAbeta42Vec,min=0, max=800, nticks=3)
dev.off()







df_PCA <- Df_Cyt_SubZ
PCA_mat <- df_PCA[,] %>% as.matrix()
Var=as.factor(Exp_data_sub2$Flicker)
#df_var=as.numeric(log2(PathVar))

#na_indices <- which(is.na(df_var), arr.ind = TRUE)



oplsOut <- opls(PCA_mat,Var, predI=2)
rotateOut <- rotate.opls(oplsOut,degrees=0)
scores_data=data.frame(rotateOut$T1,rotateOut$T2)
#shooka

# GroupColorBarF <- case_when(
#   Exp_data_sub2$Flicker == "NoFlicker" ~  "#B76945",
#   Exp_data_sub2$Flicker == "20Hz"~ "#4E5491",
#   Exp_data_sub2$Flicker == "40Hz" ~ "#863775"
# )

GroupColorBarF <- case_when(
  Exp_data_sub2$Flicker == "NoFlicker" ~ "#E8AA5B",
  Exp_data_sub2$Flicker == "20Hz" ~ "#6BACDB",
  Exp_data_sub2$Flicker == "40Hz" ~ "#D27FC2"
  #  Exp_data_sub2$Condition == "5XFAD F NoFlicker" ~  "#B76945",
  # Exp_data_sub2$Condition == "5XFAD F 20Hz"~ "#4E5491",
  # Exp_data_sub2$Condition == "5XFAD F 40Hz" ~ "#863775"
)


pdfPath <- paste("./Aug25/WT/ScorePlot_PLSDA_WT.pdf", sep = '')
pdf(pdfPath, width = 5, height = 4, pointsize = 14)
p<- scores.plot(rotateOut$T1,rotateOut$T2,Exp_data_sub2$Flicker,"Flicker",analysis = "PLSDA")
ggpar(p,palette=c("#6BACDB","#D27FC2","#E8AA5B"))
#print(p)
dev.off()
# 
# # PathVar=Df_Abeta_Sub[,1]
# # df_var=as.numeric(PathVar)
# # 
# # oplsOut <- opls(PCA_mat,df_var,predI=2)
# # rotateOut <- rotate.opls(oplsOut,degrees=0)
# # scores_data=data.frame(rotateOut$T1,rotateOut$T2)
# 
# pdfPath <- paste("./.pdf", sep = '')
# pdf(pdfPath, width = 5, height = 4, pointsize = 14)
# scores.plot(rotateOut$T1,rotateOut$T2,df_var,"df_var")
# print(p)
# dev.off()
# 
# 
# #Partial Least Squares Regression PLSR
# df_PCA <- Df_Cyt_SubZ
# PCA_mat <- as.matrix(df_PCA) 
# PathVar <- Df_Abeta_Sub[, 1]  
# df_var <- as.numeric(PathVar) 
# 
# 
# plsr_model <- plsr(df_var ~ PCA_mat, ncomp = 2, validation = "CV")
# 
# 
# scores_data <- data.frame(T1 = plsr_model$scores[, 1], T2 = plsr_model$scores[, 2])
# 
# 
# 
# pdfPath <- paste("./FiguresOut/OPLS and PLSR/Scoreplot_PLSR_GFAP", i, ".pdf", sep = '')
# pdf(pdfPath, width = 5, height = 4, pointsize = 14)
# scores.plot(scores_data$T1,scores_data$T2,PathVar,"PathVar")
# ;
# dev.off()
# 
# 
# 
# pdfPath <- paste("./FiguresOut/OPLS and PLSR/PLSR_GFAP_LV1", i, ".pdf", sep = '')
# pdf(pdfPath, width = 5, height = 4, pointsize = 14)
# loadings.chart(rotateOut$P1,"LV1");
# dev.off()
# 
# pdfPath <- paste("./FiguresOut/OPLS and PLSR/PLSR_GFAP_LV2", i, ".pdf", sep = '')
# pdf(pdfPath, width = 5, height = 4, pointsize = 14)
# loadings.chart(rotateOut$P2,"LV2");
# dev.off()
# 
# 
# 
# 
# score_df <- data.frame(
#   PC1 = rotateOut$T1,
#   PC2 = rotateOut$T2,
#   Flicker = factor(Exp_data_sub$Flicker, levels = c("NoFlicker", "20Hz", "40Hz"))
# )
# 
# 
# flicker_colors <- c("NoFlicker" = "#E298D6", "20Hz" = "#1A1513", "40Hz" = "#76A33E")
# 
# p <- ggplot(score_df, aes(x = PC1, y = PC2, color = Flicker)) +
#   geom_point(size = 3, alpha = 0.8) +
#   scale_color_manual(values = flicker_colors) +
#   geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
#   geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
#   theme_minimal(base_size = 14) +
#   labs(x = "PC1", y = "PC2", title = "Score Plot") +
#   theme(
#     plot.title = element_text(hjust = 0.5, face = "bold"),
#     legend.title = element_blank(),
#     panel.border = element_rect(color = "black", fill = NA, size = 0.8),
#     panel.grid = element_blank()
#   )
# 
# pdfPath <- "./ScorePlot_OPLS_5xFAD_Female.pdf"
# pdf(pdfPath, width = 5, height = 4, pointsize = 14)
# #png(paste0("./ScorePlot_OPLS_5xFAD_Fem.png"),height=4,width=5,units="in",res=600)
# 
# print(p)
# dev.off()







pdfPath <- paste("./Aug25/WT/PC1_WT_Cytokines_PLSDA.pdf", sep = '')
pdf(pdfPath, width = 5, height = 4, pointsize = 6)

#png(paste0("./PC1_5xFAD_Cytokines_Female.png"),height=10,width=12,units="in",res=600)
loadings.chart(rotateOut$P1,"LV1");

dev.off()

pdfPath <- paste("./Aug25/WT/PC2_WT_Cytokines_PLSDA.pdf", sep = '')
pdf(pdfPath, width = 5, height = 4, pointsize = 6)
#png(paste0("./PC2_5xFAD_Cytokines_Female.png"),height=10,width=12,units="in",res=600)
loadings.chart(rotateOut$P2,"LV2");
dev.off()



# 
# df_PCA <- Df_Cyt_SubZ
# PCA_mat <- df_PCA[,] %>% as.matrix()
# 
# oplsOut <- opls(PCA_mat,predI=2)
# rotateOut <- rotate.opls(oplsOut,degrees=0)
# scores_data=data.frame(rotateOut$T1,rotateOut$T2,Exp_data_sub$Flicker)
# Flicker=factor(Exp_data_sub$Flicker)
# 
# oplsVar = oplsOut@modelDF$R2X * 100
# 




df_plot = data.frame(T1=rotateOut$T1,T2=rotateOut$T2,Condition=Exp_data_sub2$Flicker)


statsOut_T1 = summarySE(df_plot,measurevar="T1", groupvars="Flicker")
statsOut_T2 = summarySE(df_plot,measurevar="T2", groupvars="Flicker")


export(data.frame(statsOut_T1),file="./Aug25/WT/PC1 Stats_WT_PLSDA.csv")
export(data.frame(statsOut_T2),file="./Aug25/WT/PC2 Stats_WT_PLSDA.csv")

df_plot$Condition <- factor(df_plot$Condition, levels = c("NoFlicker", "20Hz","40Hz"))  


pdf(paste0("./Aug25/WT/PC1_Violin_WT_PLSDA.pdf"),height=4,width=6,useDingbats = FALSE)
ggplot(data=df_plot, aes( y = T1, x = Condition, fill = Condition)) +
  geom_violin()+
  
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3, color = "black", position = position_dodge(width = 0.9)) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.2, color = "black", position = position_dodge(width = 0.9)) +
  
  #geom_errorbar(data=statsOut_T1, aes(x=Condition, ymin=mean-se, ymax=mean +se),
  # width=0.4, colour="black", alpha=1, size=1) +
  #geom_crossbar(data=statsOut_T1, aes(x=Condition,ymin = T1, ymax = T1),
  #  size=0.5,col="red", width = .5)+
  theme(panel.background = element_rect(fill = 'white', colour = 'black'),
        text = element_text(size=20), 
        axis.text.x = element_text(face = "plain", color = "black",
                                   size = 10,angle=45,hjust=1),
        axis.text.y = element_text(face = "plain", color = "black",
                                   size = 10),
        legend.position = "none", plot.title = element_text(hjust = 0.5))+
  xlab("")+
  ylab(paste0("Scores on LV1 (",oplsVar[1],"%)"))+
  #ylim(c(-5, 5))+
  scale_fill_manual(values = c("#E8AA5B","#6BACDB","#D27FC2"))+
  stat_compare_means(method="t.test",comparisons = list(c(1,2),c(1,3),c(2,3)),label="p.format",tip.length = 0.03,paired = FALSE)
#coord_flip()
dev.off()




pdf(paste0("./Aug25/WT/PC2_Violin_WT_PLSDA.pdf"),height=4,width=6,useDingbats = FALSE)
ggplot(data=df_plot, aes( y = T2, x = Condition, fill = Condition)) +
  geom_violin()+
  
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3, color = "black", position = position_dodge(width = 0.9)) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.2, color = "black", position = position_dodge(width = 0.9)) +
  
  #geom_errorbar(data=statsOut_T1, aes(x=Condition, ymin=mean-se, ymax=mean +se),
  # width=0.4, colour="black", alpha=1, size=1) +
  #geom_crossbar(data=statsOut_T1, aes(x=Condition,ymin = T1, ymax = T1),
  #  size=0.5,col="red", width = .5)+
  theme(panel.background = element_rect(fill = 'white', colour = 'black'),
        text = element_text(size=20), 
        axis.text.x = element_text(face = "plain", color = "black",
                                   size = 10,angle=45,hjust=1),
        axis.text.y = element_text(face = "plain", color = "black",
                                   size = 10),
        legend.position = "none", plot.title = element_text(hjust = 0.5))+
  xlab("")+
  ylab(paste0("Scores on LV2 (",oplsVar[2],"%)"))+
  ylim(c(-5, 5))+
  #scale_fill_manual(values = "#E298D6","#1A1513","#76A33E")+
  #scale_fill_manual(values = c("#B76945","#4E5491","#863775"))+
  scale_fill_manual(values = c("#E8AA5B","#6BACDB","#D27FC2"))+
  stat_compare_means(method="t.test",comparisons = list(c(1,2),c(1,3),c(2,3)),label="p.format",tip.length = 0.03,paired = FALSE)
#coord_flip()
dev.off()


#--------------Dont Need----------#
# pdfPath <- paste("./mCytokine/PCA/PCA new/ScorePlot_PCA_Female 5xFAD.pdf", sep = '')
# pdf(pdfPath, width = 10, height = 8, pointsize = 14)
# p <- scores.plot(rotateOut$T1,rotateOut$T2,Exp_data_sub$Flicker,"Flicker",analysis = "PCA")
# ggpar(p,palette=c("#1A1513","#76A33E","#E298D6"))
# 
# print(p)
# dev.off()
# 
# #blah
# 
# pdfPath <- paste("./mCytokine/PCA/PCA new/PC1_PCA_Female 5xFAD.pdf", sep = '')
# pdf(pdfPath, width = 5, height = 4, pointsize = 6)
# 
# #png(paste0("./mCytokine/PCA/PCA new/PC1_PCA_Female 5xFAD.png"),height=10,width=12,units="in",res=600)
# loadings.chart(rotateOut$P1,"PC1");
# 
# dev.off()
# 
# pdfPath <- paste("./mCytokine/PCA/PCA new/PC2_PCA_Female 5xFAD.pdf", sep = '')
# pdf(pdfPath, width = 5, height = 4, pointsize = 6)
# #png(paste0("./PC2_PCA_All_Group1.png"),height=10,width=12,units="in",res=600)
# loadings.chart(rotateOut$P2,"PC2");
# dev.off()




#Correlation of PC1 vs Amyloid beta

ind40=which(Exp_data_sub$Flicker=="40Hz")
indCTRL=which(Exp_data_sub$Flicker=="NoFlicker")

df_var=as.numeric(rotateOut$T1)
x_regress=df_var
x_regress=x_regress[indCTRL]
x_label="PC1"


for(i in 1:ncol(Df_Abeta_Sub)){
  y_regress=Df_Abeta_Sub[indCTRL,i] 
  y_label=colnames(Df_Abeta_Sub)[i]
  dataReg=data.frame(y_regress,x_regress)
  linmodel=lm(y_regress~x_regress,dataReg)
  RSQ=summary(linmodel)$r.squared
  f <- summary(linmodel)$fstatistic
  pval <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(pval) <- NULL
  
  outputFolder<-"mCytokine/PCA/PCA new"
  
  #Plot regression
  Flicker=as.factor(Exp_data_sub$Flicker[indCTRL])
  
  #pdfPath <- paste("./mCytokine/PCA/PCA new/PC1_Regression_Abeta.pdf", sep = '')
  pdf(paste("mCytokine/PCA/PCA new","/",y_label,' vs ',x_label,"_Ctrl.pdf",sep= ""), width = 5, height = 4, pointsize = 6)
  
  # pdf("mCytokine/PCA/PCA new/40Hz.pdf",width=10,height=10,pointsize = 10)
  
  #png(paste(outputFolder,"/",y_label,' vs ',x_label,"_40Hz.png",sep= ""), width=6,height=4,units="in",res=600,pointsize = 10)
  print({
    ggplot(dataReg, aes(x=x_regress, y=y_regress)) +
      geom_smooth(method='lm',se=FALSE,color="#76A33E",linetype="longdash")+
      # geom_rug(aes(color =condition))+
      geom_point(aes(color="#76A33E"),size=4)+
      #xlim(0.95*min(x_regress),1.05*max(x_regress))+
      #ylim(0.95*min(y_regress),1.05*max(y_regress))+
      theme(panel.background = element_rect(fill = 'white'),
            text = element_text(size=16),
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
      #annotate("text",size=6, x=0.8*max(x_regress)+0.2*min(x_regress), y=0.1*min(y_regress)+.9*max(y_regress), label=bquote(R^2==.(round(RSQ,digits=2))))+
      #annotate("text",size=6, x=0.8*max(x_regress)+0.2*min(x_regress), y=0*min(y_regress)+1*max(y_regress), label= bquote(p==.(round(pval,digits=4))))+
      #annotate("text",size=6, x=0.8*max(x_regress)+0.2*min(x_regress), y=0.98*min(y_regress)+0.02*max(y_regress), label=bquote(R^2==.(round(RSQ,digits=2))))+
      #annotate("text",size=6, x=0.8*max(x_regress)+0.2*min(x_regress), y=0.9*min(y_regress)+.1*max(y_regress), label= bquote(p==.(round(pval,digits=4))))+
      ggtitle(paste(y_label," vs ",x_label,sep=""))+
      ggpubr::stat_cor(method="pearson")
  })
  
  # + ggpar(p,palette=c("#76A33E"))
  dev.off()
  
}



#Correlation of PC1 vs Amyloid beta

# ind40=which(Exp_data_sub2$Flicker=="40Hz")
# indCTRL=which(Exp_data_sub2$Flicker=="NoFlicker")
# 
# df_var=as.numeric(Df_Cyt_Sub2[,5])
# x_regress=df_var
# x_regress=x_regress[ind40]
# x_label="Ab40"


df_var=as.numeric(Df_Morph_Sub2[, 18])
x_regress=df_var
x_regress=x_regress
x_label="median_Astrocyte Cell Volume" 

for(i in 1:ncol(Df_Cyt_Sub2)){
  y_regress=Df_Cyt_Sub2[,i] 
  y_label=colnames(Df_Cyt_Sub2)[i]
  dataReg=data.frame(y_regress,x_regress)
  linmodel=lm(y_regress~x_regress,dataReg)
  RSQ=summary(linmodel)$r.squared
  f <- summary(linmodel)$fstatistic
  pval <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(pval) <- NULL
  
  outputFolder<-"./CorrelationFigures/WT/Astrocyte/Median Ast Cell Volume"
  
  #Plot regression
  Flicker=as.factor(Exp_data_sub2$Flicker)
  
  #pdf(paste(outputFolder,"/",y_label,' vs ',x_label,"_Nov25_median_Astrocyte Cell Volume_vs_Cytokines_.pdf",sep= ""), width = 7, height = 4, pointsize = 6)
png(paste(outputFolder,"/",y_label,' vs ',x_label,"_Nov25_median_Astrocyte Cell Volume _vs_Cytokines_.png",sep= ""), width=6,height=4,units="in",res=600,pointsize = 10)
  print(
    ggplot(dataReg, aes(x=x_regress, y=y_regress)) +
      geom_smooth(method='lm',se=FALSE,color="black",linetype="longdash")+
      # geom_rug(aes(color =condition))+
      geom_point(aes(fill=Exp_data_sub2$Flicker),size =4 ,shape=21, stroke = 1, cex = 2)+
      #xlim(0.95*min(x_regress),1.05*max(x_regress))+
      #ylim(0.95*min(y_regress),1.05*max(y_regress))+
      theme(panel.background = element_rect(fill = 'white'),
            text = element_text(size=16),
            axis.text.x = element_text(color = "black",
                                       size = 14),
            axis.text.y = element_text(color = "black",
                                       size = 14),
            plot.title = element_text(hjust = 0.5,face='bold'),
            panel.border=element_blank(),
            axis.line=element_line(color='black'),
            plot.margin=margin(10,20,10,10)
      )+
      #scale_color_manual(values=c("#4E5491","#863775","#B76945"))+
     # scale_fill_manual(values=c("#4E5491","#863775","#B76945"))+
      scale_color_manual(values=c("#6BACDB","#D27FC2","#E8AA5B"))+
      scale_fill_manual(values = c("#6BACDB","#D27FC2","#E8AA5B"))+
      xlab(x_label)+
      ylab(y_label)+
      #annotate("text",size=6, x=0.8*max(x_regress)+0.2*min(x_regress), y=0.1*min(y_regress)+.9*max(y_regress), label=bquote(R^2==.(round(RSQ,digits=2))))+
      #annotate("text",size=6, x=0.8*max(x_regress)+0.2*min(x_regress), y=0*min(y_regress)+1*max(y_regress), label= bquote(p==.(round(pval,digits=4))))+
      #annotate("text",size=6, x=0.8*max(x_regress)+0.2*min(x_regress), y=0.98*min(y_regress)+0.02*max(y_regress), label=bquote(R^2==.(round(RSQ,digits=2))))+
      #annotate("text",size=6, x=0.8*max(x_regress)+0.2*min(x_regress), y=0.9*min(y_regress)+.1*max(y_regress), label= bquote(p==.(round(pval,digits=4))))+
      ggtitle(paste(y_label," vs ",x_label,sep=""))+
      ggpubr::stat_cor(method="pearson")
  )
  
  # + ggpar(p,palette=c("#76A33E"))
  dev.off()
  
}



for(i in 1:ncol(Df_Cyt_Sub2)){
  y_label=colnames(Df_Cyt_Sub2)[i]
  y_regress <- Df_Cyt_Sub2[, i]
  x_regress=x_regress
  dataReg=data.frame(x_regress, y_regress)
  dataReg <- drop_na(dataReg)
  
  cor_test <- cor.test(dataReg$x_regress, dataReg$y_regress, method = "spearman")
  cor_label <- paste0("rho == ", round(cor_test$estimate, 2), 
                      "~~~p == ", formatC(cor_test$p.value, format = "e", digits = 2))
  
  poly_model <- lm(y_regress ~ poly(x_regress, 2), data = dataReg)

  pdf(paste(outputFolder,"/",y_label,' vs ',x_label,"_Sep25_AB40vsCyt_v3.pdf",sep= ""), width = 7, height = 4, pointsize = 6)
  print(
    ggplot(dataReg, aes(x = x_regress, y = y_regress)) +
      geom_point(aes(fill=Exp_data_sub2$Flicker),size =5 ,shape=21, stroke = 1, cex = 2)+
      geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "longdash") +
      
      annotate("text", x = Inf, y = Inf, label = cor_label, hjust = 1.1, vjust = 1.5, parse = TRUE, size = 5) +
      
      theme(panel.background = element_rect(fill = 'white'),
            text = element_text(size = 16),
            axis.text.x = element_text(color = "black", size = 16),
            axis.text.y = element_text(color = "black", size = 16),
            plot.title = element_text(hjust = 0.5, face = 'bold'),
            axis.line = element_line(color = 'black'),
            plot.margin = margin(10, 20, 10, 20),
            strip.text.x = element_text(size = 16),
            panel.border = element_rect(colour = "black", fill = NA, size = 1.5)
      ) +
      scale_color_manual(values=c("#4E5491","#863775","#B76945"))+
      scale_fill_manual(values=c("#4E5491","#863775","#B76945"))+
      xlab("Ab40") +
      ylab(colnames(Df_Cyt_Sub2)[i]) +
      ylim(0, NA)
  )
  dev.off()
}




new_meanSEM3 <- function(x,y,colors=x,x_str="",y_str="",title=y_str,axis.text.x.size=20,axis.text.y.size=20,axis.text.x.angle=90,title.size=20,text.size=20){
  dataPlot <- data.frame(factor(x),y,factor(colors))
  colnames(dataPlot)=c("x_dot","y_dot","colors_dot")
  plotOut <- ggplot(dataPlot,aes(x=x_dot,y=y_dot)) +
    geom_boxplot(aes(color=colors))+
    geom_point(aes(colour=colors,fill=colors), size=6)+
    geom_beeswarm(aes(colour=colors,fill = colors,groupOnX = TRUE))+
    #stat_summary(geom="errorbar",color="black",size=1.5,aes(width=0.3),)+
    #stat_summary(fun = "mean", geom="point", size=4,show.legend = F, color=colors) +
    theme(panel.background = element_rect(fill = 'transparent'),
          panel.grid.major=element_blank(),
          axis.line=element_line(color="black"),
          text = element_text(size=text.size),
          legend.key = element_blank(),
          axis.text.x = element_text(color = "black",
                                     size = axis.text.x.size,angle= axis.text.x.angle),
          axis.text.y = element_text( color = "black",
                                      size = axis.text.y.size),
          plot.title=element_text(hjust=0.5,face = "bold",size=title.size))+
    xlab(x_str)+
    ylab(y_str)+
    #  ylim(c(-2, 2))+
    ggtitle(paste(title))
  return(plotOut)
}

new_meanSEM4 <- function(x, y, colors = x, 
                         x_str = "", y_str = "", title = y_str,
                         axis.text.x.size = 20, axis.text.y.size = 20, 
                         axis.text.x.angle = 90, title.size = 20, text.size = 20) {
  
  dataPlot <- data.frame(x_dot = factor(x), 
                         y_dot = y, 
                         colors_dot = factor(colors))
  
  plotOut <- ggplot(dataPlot, aes(x = x_dot, y = y_dot)) +
    
    geom_bar(aes(color = colors_dot), size = 1.5, outlier.shape = NA) +
    
    geom_beeswarm(aes(color = colors_dot, fill = colors_dot), 
                  groupOnX = TRUE, size = 6, shape = 21, stroke = 0.8) +
    
    
    # geom_point(aes(color = colors_dot, fill = colors_dot), 
    #  size = 4.5, shape = 21, stroke = 0.8, alpha = 0.7, position = position_jitter(width = 0.1)) +
    
    # mean + error bar
    #stat_summary(fun = mean, geom = "point", shape = 23, size = 8, fill = "white", color = "black") +
    # stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3, color = "black", size = 1.2) +
    
    scale_color_manual(values = c("NoFlicker" = "#E298D6", "20Hz" = "#1A1513", "40Hz" = "#76A33E")) +
    scale_fill_manual(values = c("NoFlicker" = "#E298D6", "20Hz" = "#1A1513", "40Hz" = "#76A33E")) +
    
    theme_bw(base_size = text.size) +
    theme(
      panel.border = element_rect(color = "black", fill = NA, size = 1.2),
      axis.line = element_line(color = "black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(color = "black", size = axis.text.x.size, angle = axis.text.x.angle),
      axis.text.y = element_text(color = "black", size = axis.text.y.size),
      plot.title = element_text(hjust = 0.5, face = "bold", size = title.size),
      legend.position = "none"  # Turn off legend unless you want to keep it
    ) +
    xlab(x_str) +
    ylab(y_str) +
    ggtitle(title)
  
  return(plotOut)
}

#----Make bar  plots for each cytokine


dataCyt_groupVars=data.frame(Df_Cyt_Sub, Exp_data_sub)


# 
# for(i in 1:33){
#   plotOut <- new_meanSEM3(x=factor(Exp_data_sub$Flicker, levels = c("NoFlicker", "20Hz", "40Hz")),
#                           y=Df_Cyt_Sub[,i],
#                           colors = factor(Exp_data_sub$Flicker, levels = c("NoFlicker", "20Hz", "40Hz")),
#                           title = colnames(Df_Cyt_Sub)[i],y_str = "Flourescent Int (a.u.)",
#                           axis.text.x.size = 12,axis.text.y.size = 12,title.size = 15,text.size = 15)+
#     
#   #  scale_color_manual(values=c("#3D92CF","#B385BC","#3F8A3F","#66D8FF","#A579FF","#3FFF9F"))+
#    # scale_fill_manual(values=c("#3D92CF","#B385BC","#3F8A3F","#66D8FF","#A579FF","#3FFF9F"))+
#      stat_compare_means(method="anova",p.adjust.method = "bonferroni",label="p.format",tip.length = 0.03)
#  #  pdf(paste0("./mCytokine/DotBox plots All/Zscored/Female ",colnames(dataCyt)[i],".pdf"),height=5,width=6,useDingbats = FALSE,useKerning = FALSE)
# png(paste0("./mCytokine/DotBox plots/NonZscored/Female/5xFAD/", colnames(dataCytZ)[i],".png"),height=10,width=12,units="in",res=600)
#   
#   print(plotOut)
#   dev.off()
# }
# 









for (i in 1:33) {
  flicker_levels <- c("NoFlicker", "20Hz", "40Hz")
  condition_ordered <- factor(Exp_data_sub$Flicker, levels = flicker_levels)
  
  # if (length(unique(condition_ordered)) < 2) {
  #   message(paste("Skipping", colnames(Df_Cyt_Sub)[i], "- insufficient group levels"))
  #   next
  # }
  
  plotOut <- new_meanSEM4(
    x = condition_ordered,
    y = Df_Cyt_Sub[, i],
    colors = condition_ordered,
    title = colnames(Df_Cyt_Sub)[i],
    y_str = "Fluorescent Int (a.u.)",
    axis.text.x.size = 12,
    axis.text.y.size = 12,
    title.size = 15,
    text.size = 15
  ) +
  #  stat_compare_means(method = "t-test",comparisons = c(2,3) , label = "p.format", tip.length = 0.03)
  
  png(
  filename = paste0("./Aug25/FAD/Female/BarPlot_", colnames(Df_Cyt_Sub)[i], ".png"),
   height = 10, width = 12, units = "in", res = 600
   )
 # pdf(paste0("./mCytokine/DotBox plots/NonZscored/Female/5xFAD/ ",colnames(Df_Cyt_Sub)[i],".pdf"),height=5,width=6,useDingbats = FALSE,useKerning = FALSE)
  
  
  print(plotOut)
  dev.off()
}
# for (i in 1:37)
# {
#   png(paste0("./mCytokine/DotBox plots/NonZscored/Female/5xFAD/",colnames(dataCyt)[i],".png"),height=5,width=6,units="in",res=600)
#   # pdf(paste0("./Individual Plots/",colnames(dataCyt)[i],".pdf"),height=4,width=5)
#   print({
#     ggplot(datagroupVars, aes(x=Order, y=dataCyt[,i], color=Order)) +
#       #geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2) +
#       geom_boxplot(aes(color=as.factor(Order)))+
#       #geom_line(size=1) +
#       geom_point(size=4, aes(color=as.factor(Order)))+
#       theme(panel.background = element_rect(fill = 'white', colour = 'black'),
#             text = element_text(size=10),
#             axis.text.x = element_text(face = "bold", color = "black",
#                                        size = 8),
#             axis.text.y = element_text(face = "bold", color = "black",
#                                        size = 16),
#             plot.title = element_text(hjust = 0.5))+
#       # scale_color_manual(values=c("purple3","#00A087FF"))+
#       #scale_fill_manual(values=c("purple3","#00A087FF"))+
#       stat_compare_means(method="anova",label="p.format",tip.length = 0.03) +
#       xlab("")+
#       ylab("Fluorescent Int (a.u.)")+
#       labs(color="Treatment")+
#       # ylim(c(-2, 2)) +
#       ggtitle(colnames(datagroupVars)[i])
#     # geom_text()+
#     # annotate("text", label = p)
#   })
#   dev.off()
# }
# 
# 
# 



for (i in 1:4) {
  flicker_levels <- c("NoFlicker", "20Hz", "40Hz")
  condition_ordered <- factor(Exp_data_sub$Flicker, levels = flicker_levels)
  
  if (length(unique(condition_ordered)) < 2) {
    message(paste("Skipping", colnames(Df_Abeta_Sub)[i], "- insufficient group levels"))
    next
  }
  
  plotOut <- new_meanSEM4(
    x = condition_ordered,
    y = Df_Abeta_Sub[, i],
    colors = condition_ordered,
    title = colnames(Df_Abeta_Sub)[i],
    y_str = "Fluorescent Int (a.u.)",
    axis.text.x.size = 12,
    axis.text.y.size = 12,
    title.size = 15,
    text.size = 15
  ) +
    stat_compare_means(method = "anova", label = "p.format", tip.length = 0.03)
  
  #png(
  #filename = paste0("./mCytokine/DotBox plots/NonZscored/Female/5xFAD/", colnames(dataCytZ)[i], ".png"),
  # height = 10, width = 12, units = "in", res = 600
  # )
  pdf(paste0("./mCytokine/DotBox plots/NonZscored/Female/5xFAD/Abeta_ ",colnames(Df_Abeta_Sub)[i],".pdf"),height=5,width=6,useDingbats = FALSE,useKerning = FALSE)
  
  
  print(plotOut)
  dev.off()
}


#----Correlation analysis ---- Correaltion Heatmap-----------#
library(reshape2)

dfMorph<- Df_Cyt_Sub2 [, c(1:11)]

cor_matrix <- cor(Df_Cyt_Sub2,dfMorph ,use = "pairwise.complete.obs", method = "pearson")
nSamples <- nrow(Df_Cyt_Sub)
Pvalue <- corPvalueStudent(cor_matrix, nSamples)

cor_df <- melt(cor_matrix)
colnames(cor_df) <- c("Analyte", "morphology", "Correlation")

pval_df <- melt(Pvalue)
pval_df$padj <- p.adjust(pval_df$value, method = "fdr")
colnames(pval_df) <- c("Analyte", "morphology", "Pvalue","fdr pval")

df <- merge(cor_df, pval_df, by = c("Analyte", "morphology"))

df$Significance <- cut(df$Pvalue,
                       breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                       labels = c("***", "**", "*", ""),
                       right = FALSE)


df$Label <- paste0(df$Significance)

#-------or--------------
df$Label <- paste0("R=", round(df$Correlation, 2), 
                   #"\n", df$Significance,
                   "\nPvalue=", signif(df$Pvalue, 2))

#-------------
df$padj <- p.adjust(pval_df$Pvalue, method = "fdr")

row_dist <- dist(cor_matrix)
row_clust <- hclust(row_dist)
row_order <- row_clust$labels[row_clust$order]

col_dist <- dist(t(cor_matrix))
col_clust <- hclust(col_dist)
col_order <- col_clust$labels[col_clust$order]

df$Analyte <- factor(df$Analyte, levels = row_order)
df$morphology <- factor(df$morphology, levels = col_order)

#df_sub=df %>% filter(!(Cytokine == c("IL-6",
                                    # "IL-4","IL-2","IL-1B","MCP-1","KC","IL-5","IP-10","IFNy","MIP-1a","IL-9","TNFa","IL-13","IL-1a","IL-17"))) 


# --- PLOT ---
pdf("./Aug25/PearsonCorrealtion_Morphology_vs_Cytokines_with Pvalues_sep25.pdf", width = 20, height = 16, pointsize = 10)
ggplot(df, aes(x = morphology, y = Analyte, fill = Correlation)) +
  geom_tile(color = "white") +
  geom_text(aes(label = Label), size = 2) +
  scale_fill_gradient2(
    low = "blue", mid = "white", high = "red", midpoint = 0,
    limit = c(-1, 1), space = "Lab", name = "pearson R"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text = element_text(size = 10),
    axis.title = element_blank()
  )

dev.off()



write_csv(df[, c("Analyte", "morphology", "Correlation", "Pvalue", "padj", "Significance")],
          "./PearsonCorrelation_Morphology vs Cytokines_CorrelationMatrix_Pvals_FDR_SummaryTable_sep25.csv")




#---------------Old analysis----------------#

sortby <- which(Exp_data$SexGenotype =="WT M")

Df_Cyt_Sub <-  dataCyt[sortby ,] 
Df_Abeta_Sub <-dataAbeta[sortby ,]
Exp_data_sub <- Exp_data[sortby ,]



#Z scoring data

Df_Cyt_SubZ=(apply(t(Df_Cyt_Sub),1,scale)); #z-score the data
Df_Cyt_SubZ[is.na(Df_Cyt_SubZ)]=0
rownames(Df_Cyt_SubZ)=rownames(Df_Cyt_Sub)




FlickerColorBar <- case_when(
  Exp_data_sub$Flicker == "20Hz" ~ "#5C81BD",
  Exp_data_sub$Flicker == "40Hz" ~ "#3E2A47",
  Exp_data_sub$Flicker == "NoFlicker" ~ "#E8A941"
)



colorBarCombined=cbind(c(FlickerColorBar),c(FlickerColorBar))
colnames(colorBarCombined)=c("Flicker","Flicker")


#Cluster the cytokine data by columns
hc= hclust(dist(t(Df_Cyt_SubZ),method = "euclidean"), method="average")

##Make Heatmap
pdf("Male WT/Cytokine profile Male WT Only.pdf", width=12,height=10,pointsize = 6)
#png("Male WT/Cytokine profile Male WT Only.png", width=16,height=16,units="in",res=600, pointsize = 20)
heatmap3(Df_Cyt_SubZ,main="Male WT Cytokine Profile",col=barColors, breaks=breakBarColors,legendfun=function()showLegend(),
         Rowv=NA, Colv=as.dendrogram(hc),  scale="none", na.rm = TRUE,
         labRow = rownames(Df_Cyt_SubZ), labCol = colnames(Df_Cyt_SubZ),
         RowSideColors=colorBarCombined)


dev.off()


df_PCA <- Df_Cyt_SubZ
PCA_mat <- df_PCA[,] %>% as.matrix()

oplsOut <- opls(PCA_mat,predI=2)
rotateOut <- rotate.opls(oplsOut,degrees=0)
scores_data=data.frame(rotateOut$T1,rotateOut$T2)
# 


pdfPath <- paste("./Male WT/ScorePlot_OPLS_WT_Male.pdf", sep = '')
pdf(pdfPath, width = 5, height = 4, pointsize = 14)
p<- scores.plot(rotateOut$T1,rotateOut$T2,factor(Exp_data_sub$Flicker),"Flicker")
ggpar(p,palette=c("#5C81BD","#3E2A47","#E8A941"))
print(p)
dev.off()

# PathVar=Df_Abeta_Sub[,1]
# df_var=as.numeric(PathVar)
# 
# oplsOut <- opls(PCA_mat,df_var,predI=2)
# rotateOut <- rotate.opls(oplsOut,degrees=0)
# scores_data=data.frame(rotateOut$T1,rotateOut$T2)
# 







pdfPath <- paste("./Male WT/PC1_WT_Cytokines_Male.pdf", sep = '')
pdf(pdfPath, width = 5, height = 4, pointsize = 6)

#png(paste0("./Male WT/PC1_WT_Cytokines_Male.png"),height=10,width=12,units="in",res=600)
loadings.chart(rotateOut$P1,"PC1");

dev.off()

pdfPath <- paste("./Male WT/PC2_WT_Cytokines_Female.pdf", sep = '')
pdf(pdfPath, width = 5, height = 4, pointsize = 6)
#png(paste0("./Male WT/PC_WT_Cytokines_Male.png"),height=10,width=12,units="in",res=600)
loadings.chart(rotateOut$P2,"PC2");
dev.off()







#PCA Analaysis

#----PCA---
df_PCA <- Df_Cyt_SubZ
PCA_mat <- df_PCA[,] %>% as.matrix()

oplsOut <- opls(PCA_mat,predI=2)
rotateOut <- rotate.opls(oplsOut,degrees=0)
scores_data=data.frame(rotateOut$T1,rotateOut$T2,Exp_data_sub$Flicker)
Flicker=factor(Exp_data_sub$Flicker)

oplsVar = oplsOut@modelDF$R2X * 100



df_plot = data.frame(T1=rotateOut$T1,T2=rotateOut$T2,Condition=Exp_data_sub$Flicker)


statsOut_T1 = summarySE(df_plot,measurevar="T1", groupvars="Flicker")
statsOut_T2 = summarySE(df_plot,measurevar="T2", groupvars="Flicker")


export(data.frame(statsOut_T1),file="./Male WT/PC1 Stats_Male wt.csv")
export(data.frame(statsOut_T2),file="./Male WT/PC2 Stats Male wt.csv")


pdf(paste0("./Male WT/PC1_Violin_Male wt.pdf"),height=4,width=6,useDingbats = FALSE)
ggplot(data=df_plot, aes( y = T1, x = Condition, fill = Condition)) +
  geom_violin()+
  
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3, color = "black", position = position_dodge(width = 0.9)) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.2, color = "black", position = position_dodge(width = 0.9)) +
  
  #geom_errorbar(data=statsOut_T1, aes(x=Condition, ymin=mean-se, ymax=mean +se),
  # width=0.4, colour="black", alpha=1, size=1) +
  #geom_crossbar(data=statsOut_T1, aes(x=Condition,ymin = T1, ymax = T1),
  #  size=0.5,col="red", width = .5)+
  theme(panel.background = element_rect(fill = 'white', colour = 'black'),
        text = element_text(size=20), 
        axis.text.x = element_text(face = "plain", color = "black",
                                   size = 10,angle=45,hjust=1),
        axis.text.y = element_text(face = "plain", color = "black",
                                   size = 10),
        legend.position = "none", plot.title = element_text(hjust = 0.5))+
  xlab("")+
  ylab(paste0("Scores on PC1 (",oplsVar[1],"%)"))+
  #ylim(c(-5, 5))+
  scale_fill_manual(values = c("#5C81BD","#3E2A47","#E8A941"))+
  stat_compare_means(method="anova",p.adjust.methods = "bonferroni",label="p.format",tip.length = 0.03)
#coord_flip()
dev.off()



pdf(paste0("Male WT/PC2_Violin_Male WT.pdf"),height=4,width=6,useDingbats = FALSE)
ggplot(data=df_plot, aes( y = T2, x = Condition, fill = Condition)) +
  geom_violin()+
  
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3, color = "black", position = position_dodge(width = 0.9)) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.2, color = "black", position = position_dodge(width = 0.9)) +
  
  #geom_errorbar(data=statsOut_T1, aes(x=Condition, ymin=mean-se, ymax=mean +se),
  # width=0.4, colour="black", alpha=1, size=1) +
  #geom_crossbar(data=statsOut_T1, aes(x=Condition,ymin = T1, ymax = T1),
  #  size=0.5,col="red", width = .5)+
  theme(panel.background = element_rect(fill = 'white', colour = 'black'),
        text = element_text(size=20), 
        axis.text.x = element_text(face = "plain", color = "black",
                                   size = 10,angle=45,hjust=1),
        axis.text.y = element_text(face = "plain", color = "black",
                                   size = 10),
        legend.position = "none", plot.title = element_text(hjust = 0.5))+
  xlab("")+
  ylab(paste0("Scores on PC2 (",oplsVar[2],"%)"))+
  ylim(c(-5, 5))+
  #scale_fill_manual(values = "#E298D6","#1A1513","#76A33E")+
  scale_fill_manual(values = c("#5C81BD","#3E2A47","#E8A941"))+
  stat_compare_means(method="anova",p.adjust.methods = "bonferroni",label="p.format",tip.length = 0.03)
#coord_flip()
dev.off()



pdfPath <- paste("./Male WT/ScorePlot_PCA_Male WT.pdf", sep = '')
pdf(pdfPath, width = 10, height = 8, pointsize = 14)
p <- scores.plot(rotateOut$T1,rotateOut$T2,Exp_data_sub$Flicker,"Flicker",analysis = "PCA")
ggpar(p,palette=c("#5C81BD","#3E2A47","#E8A941"))

print(p)
dev.off()


pdfPath <- paste("./Male WT/PC1_PCA_Male WT.pdf", sep = '')
pdf(pdfPath, width = 5, height = 4, pointsize = 6)

#png(paste0("./Male WT/PC1_PCA_Male WT.png"),height=10,width=12,units="in",res=600)
loadings.chart(rotateOut$P1,"PC1");

dev.off()

pdfPath <- paste("./Male WT/PC2_PCA_Male WT.pdf", sep = '')
pdf(pdfPath, width = 5, height = 4, pointsize = 6)
#png(paste0("./PC2_PCA_All_Group1.png"),height=10,width=12,units="in",res=600)
loadings.chart(rotateOut$P2,"PC2");
dev.off()





new_meanSEM3 <- function(x,y,colors=x,x_str="",y_str="",title=y_str,axis.text.x.size=20,axis.text.y.size=20,axis.text.x.angle=90,title.size=20,text.size=20){
  dataPlot <- data.frame(factor(x),y,factor(colors))
  colnames(dataPlot)=c("x_dot","y_dot","colors_dot")
  plotOut <- ggplot(dataPlot,aes(x=x_dot,y=y_dot)) +
    geom_boxplot(aes(color=colors))+
    geom_point(aes(colour=colors,fill=colors), size=6)+
    geom_beeswarm(aes(colour=colors,fill = colors,groupOnX = TRUE))+
    #stat_summary(geom="errorbar",color="black",size=1.5,aes(width=0.3),)+
    #stat_summary(fun = "mean", geom="point", size=4,show.legend = F, color=colors) +
    theme(panel.background = element_rect(fill = 'transparent'),
          panel.grid.major=element_blank(),
          axis.line=element_line(color="black"),
          text = element_text(size=text.size),
          legend.key = element_blank(),
          axis.text.x = element_text(color = "black",
                                     size = axis.text.x.size,angle= axis.text.x.angle),
          axis.text.y = element_text( color = "black",
                                      size = axis.text.y.size),
          plot.title=element_text(hjust=0.5,face = "bold",size=title.size))+
    xlab(x_str)+
    ylab(y_str)+
    #  ylim(c(-2, 2))+
    ggtitle(paste(title))
  return(plotOut)
}

new_meanSEM4 <- function(x, y, colors = x, 
                         x_str = "", y_str = "", title = y_str,
                         axis.text.x.size = 20, axis.text.y.size = 20, 
                         axis.text.x.angle = 90, title.size = 20, text.size = 20) {
  
  dataPlot <- data.frame(x_dot = factor(x), 
                         y_dot = y, 
                         colors_dot = factor(colors))
  
  plotOut <- ggplot(dataPlot, aes(x = x_dot, y = y_dot)) +
    
    geom_boxplot(aes(color = colors_dot), size = 1.5, outlier.shape = NA) +
    
    geom_beeswarm(aes(color = colors_dot, fill = colors_dot), 
                  groupOnX = TRUE, size = 6, shape = 21, stroke = 0.8) +

    scale_color_manual(values = c("NoFlicker" = "#E8A941", "20Hz" = "#5C81BD", "40Hz" = "#3E2A47")) +
    scale_fill_manual(values = c("NoFlicker" = "#E8A941", "20Hz" = "#5C81BD", "40Hz" = "#3E2A47")) +
    
    theme_bw(base_size = text.size) +
    theme(
      panel.border = element_rect(color = "black", fill = NA, size = 1.2),
      axis.line = element_line(color = "black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(color = "black", size = axis.text.x.size, angle = axis.text.x.angle),
      axis.text.y = element_text(color = "black", size = axis.text.y.size),
      plot.title = element_text(hjust = 0.5, face = "bold", size = title.size),
      legend.position = "none"  # Turn off legend unless you want to keep it
    ) +
    xlab(x_str) +
    ylab(y_str) +
    ggtitle(title)
  
  return(plotOut)
}

#----Make bar  plots for each cytokine

dataCyt_groupVars=data.frame(Df_Cyt_Sub, Exp_data_sub)


for (i in 1:33) {
  flicker_levels <- c("NoFlicker", "20Hz", "40Hz")
  condition_ordered <- factor(Exp_data_sub$Flicker, levels = flicker_levels)
  
  if (length(unique(condition_ordered)) < 2) {
    message(paste("Skipping", colnames(Df_Cyt_Sub)[i], "- insufficient group levels"))
    next
  }
  
  plotOut <- new_meanSEM4(
    x = condition_ordered,
    y = Df_Cyt_Sub[, i],
    colors = condition_ordered,
    title = colnames(Df_Cyt_Sub)[i],
    y_str = "Fluorescent Int (a.u.)",
    axis.text.x.size = 12,
    axis.text.y.size = 12,
    title.size = 15,
    text.size = 15
  ) +
    stat_compare_means(method = "anova", label = "p.format", tip.length = 0.03)
  
  #png(
  #filename = paste0("./Male WT/Boxplots/", colnames(dataCytZ)[i], ".png"),
  # height = 10, width = 12, units = "in", res = 600
  # )
  pdf(paste0("./Male WT/Boxplots/ ",colnames(Df_Cyt_Sub)[i],".pdf"),height=5,width=6,useDingbats = FALSE,useKerning = FALSE)
  
  
  print(plotOut)
  dev.off()
}






#----Female WT---#
sortby <- which(Exp_data$SexGenotype =="WT F")

Df_Cyt_Sub <-  dataCyt[sortby ,] 
Df_Abeta_Sub <-dataAbeta[sortby ,]
Exp_data_sub <- Exp_data[sortby ,]



#Z scoring data

Df_Cyt_SubZ=(apply(t(Df_Cyt_Sub),1,scale)); #z-score the data
Df_Cyt_SubZ[is.na(Df_Cyt_SubZ)]=0
rownames(Df_Cyt_SubZ)=rownames(Df_Cyt_Sub)




FlickerColorBar <- case_when(
  Exp_data_sub$Flicker == "20Hz" ~ "#5C81BD",
  Exp_data_sub$Flicker == "40Hz" ~ "#3E2A47",
  Exp_data_sub$Flicker == "NoFlicker" ~ "#E8A941"
)



colorBarCombined=cbind(c(FlickerColorBar),c(FlickerColorBar))
colnames(colorBarCombined)=c("Flicker","Flicker")


#Cluster the cytokine data by columns
hc= hclust(dist(t(Df_Cyt_SubZ),method = "euclidean"), method="average")

##Make Heatmap
#pdf("Female WT/Cytokine profile Female WT Only.pdf", width=12,height=10,pointsize = 6)
png("Female WT/Cytokine profile Female WT Only.png", width=16,height=16,units="in",res=600, pointsize = 20)
heatmap3(Df_Cyt_SubZ,main="Female WT Cytokine Profile",col=barColors, breaks=breakBarColors,legendfun=function()showLegend(),
         Rowv=NA, Colv=as.dendrogram(hc),  scale="none", na.rm = TRUE,
         labRow = rownames(Df_Cyt_SubZ), labCol = colnames(Df_Cyt_SubZ),
         RowSideColors=colorBarCombined)


dev.off()


df_PCA <- Df_Cyt_SubZ
PCA_mat <- df_PCA[,] %>% as.matrix()

oplsOut <- opls(PCA_mat,predI=2)
rotateOut <- rotate.opls(oplsOut,degrees=0)
scores_data=data.frame(rotateOut$T1,rotateOut$T2)
# 


pdfPath <- paste("./Female WT/ScorePlot_OPLS_WT_Female.pdf", sep = '')
pdf(pdfPath, width = 5, height = 4, pointsize = 14)
p<- scores.plot(rotateOut$T1,rotateOut$T2,factor(Exp_data_sub$Flicker),"Flicker")
ggpar(p,palette=c("#5C81BD","#3E2A47","#E8A941"))
print(p)
dev.off()

# PathVar=Df_Abeta_Sub[,1]
# df_var=as.numeric(PathVar)
# 
# oplsOut <- opls(PCA_mat,df_var,predI=2)
# rotateOut <- rotate.opls(oplsOut,degrees=0)
# scores_data=data.frame(rotateOut$T1,rotateOut$T2)
# 




#Partial Least Squares Regression PLSR
df_PCA <- Df_Cyt_SubZ
PCA_mat <- as.matrix(df_PCA) 
PathVar <- Df_Abeta_Sub[, 1]  
df_var <- as.numeric(PathVar) 


plsr_model <- plsr(df_var ~ PCA_mat, ncomp = 2, validation = "CV")


scores_data <- data.frame(T1 = plsr_model$scores[, 1], T2 = plsr_model$scores[, 2])



pdfPath <- paste("./FiguresOut/OPLS and PLSR/Scoreplot_PLSR_GFAP", i, ".pdf", sep = '')
pdf(pdfPath, width = 5, height = 4, pointsize = 14)
scores.plot(scores_data$T1,scores_data$T2,PathVar,"PathVar")
;
dev.off()



pdfPath <- paste("./FiguresOut/OPLS and PLSR/PLSR_GFAP_LV1", i, ".pdf", sep = '')
pdf(pdfPath, width = 5, height = 4, pointsize = 14)
loadings.chart(rotateOut$P1,"LV1");
dev.off()

pdfPath <- paste("./FiguresOut/OPLS and PLSR/PLSR_GFAP_LV2", i, ".pdf", sep = '')
pdf(pdfPath, width = 5, height = 4, pointsize = 14)
loadings.chart(rotateOut$P2,"LV2");
dev.off()




score_df <- data.frame(
  PC1 = rotateOut$T1,
  PC2 = rotateOut$T2,
  Flicker = factor(Exp_data_sub$Flicker, levels = c("NoFlicker", "20Hz", "40Hz"))
)


flicker_colors <- c("NoFlicker" = "#E298D6", "20Hz" = "#1A1513", "40Hz" = "#76A33E")

p <- ggplot(score_df, aes(x = PC1, y = PC2, color = Flicker)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = flicker_colors) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  theme_minimal(base_size = 14) +
  labs(x = "PC1", y = "PC2", title = "Score Plot") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.title = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.8),
    panel.grid = element_blank()
  )

pdfPath <- "./ScorePlot_OPLS_5xFAD_Female.pdf"
pdf(pdfPath, width = 5, height = 4, pointsize = 14)
#png(paste0("./ScorePlot_OPLS_5xFAD_Fem.png"),height=4,width=5,units="in",res=600)

print(p)
dev.off()







pdfPath <- paste("./Female WT/PC1_WT_Cytokines_Female.pdf", sep = '')
pdf(pdfPath, width = 5, height = 4, pointsize = 6)

#png(paste0("./Female WT/PC1_WT_Cytokines_Female.png"),height=10,width=12,units="in",res=600)
loadings.chart(rotateOut$P1,"PC1");

dev.off()

pdfPath <- paste("./Female WT/PC2_WT_Cytokines_Female.pdf", sep = '')
pdf(pdfPath, width = 5, height = 4, pointsize = 6)
#png(paste0("./Female WT/PC_WT_Cytokines_Female.png"),height=10,width=12,units="in",res=600)
loadings.chart(rotateOut$P2,"PC2");
dev.off()







#PCA Analaysis

#----PCA---
df_PCA <- Df_Cyt_SubZ
PCA_mat <- df_PCA[,] %>% as.matrix()

oplsOut <- opls(PCA_mat,predI=2)
rotateOut <- rotate.opls(oplsOut,degrees=0)
scores_data=data.frame(rotateOut$T1,rotateOut$T2,Exp_data_sub$Flicker)
Flicker=factor(Exp_data_sub$Flicker)

oplsVar = oplsOut@modelDF$R2X * 100





df_plot = data.frame(T1=rotateOut$T1,T2=rotateOut$T2,Condition=Exp_data_sub$Flicker)


statsOut_T1 = summarySE(df_plot,measurevar="T1", groupvars="Flicker")
statsOut_T2 = summarySE(df_plot,measurevar="T2", groupvars="Flicker")


export(data.frame(statsOut_T1),file="./Female WT/PC1 Stats_Female wt.csv")
export(data.frame(statsOut_T2),file="./Female WT/PC2 Stats Female wt.csv")


pdf(paste0("./Female WT/PC1_Violin_Female wt.pdf"),height=4,width=6,useDingbats = FALSE)
ggplot(data=df_plot, aes( y = T1, x = Condition, fill = Condition)) +
  geom_violin()+
  
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3, color = "black", position = position_dodge(width = 0.9)) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.2, color = "black", position = position_dodge(width = 0.9)) +
  
  #geom_errorbar(data=statsOut_T1, aes(x=Condition, ymin=mean-se, ymax=mean +se),
  # width=0.4, colour="black", alpha=1, size=1) +
  #geom_crossbar(data=statsOut_T1, aes(x=Condition,ymin = T1, ymax = T1),
  #  size=0.5,col="red", width = .5)+
  theme(panel.background = element_rect(fill = 'white', colour = 'black'),
        text = element_text(size=20), 
        axis.text.x = element_text(face = "plain", color = "black",
                                   size = 10,angle=45,hjust=1),
        axis.text.y = element_text(face = "plain", color = "black",
                                   size = 10),
        legend.position = "none", plot.title = element_text(hjust = 0.5))+
  xlab("")+
  ylab(paste0("Scores on PC1 (",oplsVar[1],"%)"))+
  #ylim(c(-5, 5))+
  scale_fill_manual(values = c("#5C81BD","#3E2A47","#E8A941"))+
  stat_compare_means(method="anova",p.adjust.methods = "bonferroni",label="p.format",tip.length = 0.03)
#coord_flip()
dev.off()



pdf(paste0("Female WT/PC2_Violin_Female WT.pdf"),height=4,width=6,useDingbats = FALSE)
ggplot(data=df_plot, aes( y = T2, x = Condition, fill = Condition)) +
  geom_violin()+
  
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3, color = "black", position = position_dodge(width = 0.9)) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.2, color = "black", position = position_dodge(width = 0.9)) +
  
  #geom_errorbar(data=statsOut_T1, aes(x=Condition, ymin=mean-se, ymax=mean +se),
  # width=0.4, colour="black", alpha=1, size=1) +
  #geom_crossbar(data=statsOut_T1, aes(x=Condition,ymin = T1, ymax = T1),
  #  size=0.5,col="red", width = .5)+
  theme(panel.background = element_rect(fill = 'white', colour = 'black'),
        text = element_text(size=20), 
        axis.text.x = element_text(face = "plain", color = "black",
                                   size = 10,angle=45,hjust=1),
        axis.text.y = element_text(face = "plain", color = "black",
                                   size = 10),
        legend.position = "none", plot.title = element_text(hjust = 0.5))+
  xlab("")+
  ylab(paste0("Scores on PC2 (",oplsVar[2],"%)"))+
  ylim(c(-5, 5))+
  #scale_fill_manual(values = "#E298D6","#1A1513","#76A33E")+
  scale_fill_manual(values = c("#5C81BD","#3E2A47","#E8A941"))+
  stat_compare_means(method="anova",p.adjust.methods = "bonferroni",label="p.format",tip.length = 0.03)
#coord_flip()
dev.off()



pdfPath <- paste("./Female WT/ScorePlot_PCA_Female WT.pdf", sep = '')
pdf(pdfPath, width = 10, height = 8, pointsize = 14)
p <- scores.plot(rotateOut$T1,rotateOut$T2,Exp_data_sub$Flicker,"Flicker",analysis = "PCA")
ggpar(p,palette=c("#5C81BD","#3E2A47","#E8A941"))

print(p)
dev.off()

#blah

pdfPath <- paste("./Female WT/PC1_PCA_Female WT.pdf", sep = '')
pdf(pdfPath, width = 5, height = 4, pointsize = 6)

#png(paste0("./Female WT/PC1_PCA_Female WT.png"),height=10,width=12,units="in",res=600)
loadings.chart(rotateOut$P1,"PC1");

dev.off()

pdfPath <- paste("./Female WT/PC2_PCA_Female WT.pdf", sep = '')
pdf(pdfPath, width = 5, height = 4, pointsize = 6)
#png(paste0("./PC2_PCA_All_Group1.png"),height=10,width=12,units="in",res=600)
loadings.chart(rotateOut$P2,"PC2");
dev.off()



#----Make bar  plots for each cytokine

dataCyt_groupVars=data.frame(Df_Cyt_Sub, Exp_data_sub)


for (i in 1:33) {
  flicker_levels <- c("NoFlicker", "20Hz", "40Hz")
  condition_ordered <- factor(Exp_data_sub$Flicker, levels = flicker_levels)
  
  if (length(unique(condition_ordered)) < 2) {
    message(paste("Skipping", colnames(Df_Cyt_Sub)[i], "- insufficient group levels"))
    next
  }
  
  plotOut <- new_meanSEM4(
    x = condition_ordered,
    y = Df_Cyt_Sub[, i],
    colors = condition_ordered,
    title = colnames(Df_Cyt_Sub)[i],
    y_str = "Fluorescent Int (a.u.)",
    axis.text.x.size = 12,
    axis.text.y.size = 12,
    title.size = 15,
    text.size = 15
  ) +
    stat_compare_means(method = "anova", label = "p.format", tip.length = 0.03)
  
  #png(
  #filename = paste0("./Male WT/Boxplots/", colnames(dataCytZ)[i], ".png"),
  # height = 10, width = 12, units = "in", res = 600
  # )
  pdf(paste0("./Female WT/Boxplots/ ",colnames(Df_Cyt_Sub)[i],".pdf"),height=5,width=6,useDingbats = FALSE,useKerning = FALSE)
  
  
  print(plotOut)
  dev.off()
}




# 
# #Correlation of PC1 vs Amyloid beta
# 
# ind40=which(Exp_data_sub$Flicker=="40Hz")
# indCTRL=which(Exp_data_sub$Flicker=="NoFlicker")
# 
# df_var=as.numeric(rotateOut$T1)
# x_regress=df_var
# x_regress=x_regress[indCTRL]
# x_label="PC1"
# 
# 
# for(i in 1:ncol(Df_Abeta_Sub)){
#   y_regress=Df_Abeta_Sub[indCTRL,i] 
#   y_label=colnames(Df_Abeta_Sub)[i]
#   dataReg=data.frame(y_regress,x_regress)
#   linmodel=lm(y_regress~x_regress,dataReg)
#   RSQ=summary(linmodel)$r.squared
#   f <- summary(linmodel)$fstatistic
#   pval <- pf(f[1],f[2],f[3],lower.tail=F)
#   attributes(pval) <- NULL
#   
#   outputFolder<-"mCytokine/PCA/PCA new"
#   
#   #Plot regression
#   Flicker=as.factor(Exp_data_sub$Flicker[indCTRL])
#   
#   #pdfPath <- paste("./mCytokine/PCA/PCA new/PC1_Regression_Abeta.pdf", sep = '')
#   pdf(paste("mCytokine/PCA/PCA new","/",y_label,' vs ',x_label,"_Ctrl.pdf",sep= ""), width = 5, height = 4, pointsize = 6)
#   
#   # pdf("mCytokine/PCA/PCA new/40Hz.pdf",width=10,height=10,pointsize = 10)
#   
#   #png(paste(outputFolder,"/",y_label,' vs ',x_label,"_40Hz.png",sep= ""), width=6,height=4,units="in",res=600,pointsize = 10)
#   print({
#     ggplot(dataReg, aes(x=x_regress, y=y_regress)) +
#       geom_smooth(method='lm',se=FALSE,color="#76A33E",linetype="longdash")+
#       # geom_rug(aes(color =condition))+
#       geom_point(aes(color="#76A33E"),size=4)+
#       #xlim(0.95*min(x_regress),1.05*max(x_regress))+
#       #ylim(0.95*min(y_regress),1.05*max(y_regress))+
#       theme(panel.background = element_rect(fill = 'white'),
#             text = element_text(size=16),
#             axis.text.x = element_text(color = "black",
#                                        size = 14),
#             axis.text.y = element_text(color = "black",
#                                        size = 14),
#             plot.title = element_text(hjust = 0.5,face='bold'),
#             panel.border=element_blank(),
#             axis.line=element_line(color='black'),
#             plot.margin=margin(10,20,10,10)
#       )+
#       xlab(x_label)+
#       ylab(y_label)+
#       #annotate("text",size=6, x=0.8*max(x_regress)+0.2*min(x_regress), y=0.1*min(y_regress)+.9*max(y_regress), label=bquote(R^2==.(round(RSQ,digits=2))))+
#       #annotate("text",size=6, x=0.8*max(x_regress)+0.2*min(x_regress), y=0*min(y_regress)+1*max(y_regress), label= bquote(p==.(round(pval,digits=4))))+
#       #annotate("text",size=6, x=0.8*max(x_regress)+0.2*min(x_regress), y=0.98*min(y_regress)+0.02*max(y_regress), label=bquote(R^2==.(round(RSQ,digits=2))))+
#       #annotate("text",size=6, x=0.8*max(x_regress)+0.2*min(x_regress), y=0.9*min(y_regress)+.1*max(y_regress), label= bquote(p==.(round(pval,digits=4))))+
#       ggtitle(paste(y_label," vs ",x_label,sep=""))+
#       ggpubr::stat_cor(method="pearson")
#   })
#   
#   # + ggpar(p,palette=c("#76A33E"))
#   dev.off()
#   
# }
# 
# 
# 
# #Correlation of PC1 vs Amyloid beta
# 
# ind40=which(Exp_data_sub$Flicker=="40Hz")
# indCTRL=which(Exp_data_sub$Flicker=="NoFlicker")
# 
# df_var=as.numeric(Df_Abeta_Sub[,1])
# x_regress=df_var
# x_regress=x_regress[ind40]
# x_label="Ab40"
# 
# 
# for(i in 1:ncol(Df_Cyt_Sub)){
#   y_regress=Df_Cyt_Sub[ind40,i] 
#   y_label=colnames(Df_Cyt_Sub)[i]
#   dataReg=data.frame(y_regress,x_regress)
#   linmodel=lm(y_regress~x_regress,dataReg)
#   RSQ=summary(linmodel)$r.squared
#   f <- summary(linmodel)$fstatistic
#   pval <- pf(f[1],f[2],f[3],lower.tail=F)
#   attributes(pval) <- NULL
#   
#   outputFolder<-"Abeta vs Cytokines"
#   
#   #Plot regression
#   Flicker=as.factor(Exp_data_sub$Flicker[ind40])
#   
#   pdf(paste("Abeta vs Cytokines","/",y_label,' vs ',x_label,"_40Hz_ABETA40.pdf",sep= ""), width = 5, height = 4, pointsize = 6)
#   #png(paste(outputFolder,"/",y_label,' vs ',x_label,"40Hz.png",sep= ""), width=6,height=4,units="in",res=600,pointsize = 10)
#   print({
#     ggplot(dataReg, aes(x=x_regress, y=y_regress)) +
#       geom_smooth(method='lm',se=FALSE,color="#76A33E",linetype="longdash")+
#       # geom_rug(aes(color =condition))+
#       geom_point(color="#76A33E",size=4)+
#       #xlim(0.95*min(x_regress),1.05*max(x_regress))+
#       #ylim(0.95*min(y_regress),1.05*max(y_regress))+
#       theme(panel.background = element_rect(fill = 'white'),
#             text = element_text(size=16),
#             axis.text.x = element_text(color = "black",
#                                        size = 14),
#             axis.text.y = element_text(color = "black",
#                                        size = 14),
#             plot.title = element_text(hjust = 0.5,face='bold'),
#             panel.border=element_blank(),
#             axis.line=element_line(color='black'),
#             plot.margin=margin(10,20,10,10)
#       )+
#       xlab(x_label)+
#       ylab(y_label)+
#       #annotate("text",size=6, x=0.8*max(x_regress)+0.2*min(x_regress), y=0.1*min(y_regress)+.9*max(y_regress), label=bquote(R^2==.(round(RSQ,digits=2))))+
#       #annotate("text",size=6, x=0.8*max(x_regress)+0.2*min(x_regress), y=0*min(y_regress)+1*max(y_regress), label= bquote(p==.(round(pval,digits=4))))+
#       #annotate("text",size=6, x=0.8*max(x_regress)+0.2*min(x_regress), y=0.98*min(y_regress)+0.02*max(y_regress), label=bquote(R^2==.(round(RSQ,digits=2))))+
#       #annotate("text",size=6, x=0.8*max(x_regress)+0.2*min(x_regress), y=0.9*min(y_regress)+.1*max(y_regress), label= bquote(p==.(round(pval,digits=4))))+
#       ggtitle(paste(y_label," vs ",x_label,sep=""))+
#       ggpubr::stat_cor(method="pearson")
#   })
#   
#   # + ggpar(p,palette=c("#76A33E"))
#   dev.off()
#   
# }
# 
# 





#Amyloid beta side bar


ab42<-Df_Abeta_Sub[,2]
ab40<-Df_Abeta_Sub[,1]
ratio<-Df_Abeta_Sub[,3]


breakBarColorsAbeta42=c(-1,seq(0, 1200, 0.1))*10 #Outside numbers clip outliers. This is for zscoring.
barColorsAbeta42Vec = c("black",colorpanel(length(breakBarColorsAbeta42), "white", "gold4"))
indMatchColorone=match(round(ab42*10),round(breakBarColorsAbeta42))
barColorsAb42=barColorsAbeta42Vec[indMatchColorone]

breakBarColorsAbeta40=c(-1,seq(0, 10, 0.01))*10 #Outside numbers clip outliers. This is for zscoring.
barColorsAbeta40Vec = c("black",colorpanel(length(breakBarColorsAbeta40), "white", "orange"))
indMatchColortwo=match(round(ab40*10),round(breakBarColorsAbeta40))
barColorsAb40=barColorsAbeta40Vec[indMatchColortwo]


breakBarColorsAbetaratio=c(-1,seq(0, 180, 0.01))*10 #Outside numbers clip outliers. This is for zscoring.
barColorsAbetaratioVec = c("black",colorpanel(length(breakBarColorsAbetaratio), "white", "gray4"))
indMatchColorthree=match(round(ratio*10),round(breakBarColorsAbetaratio))
barColorsAbratio=barColorsAbetaratioVec[indMatchColorthree]

colorBarCombined=cbind(c(FlickerColorBar),c(barColorsAb42),c(barColorsAb40),c(barColorsAbratio))
colnames(colorBarCombined)=c("Flicker","Ab42","Ab40","ratio")



##Make Heatmap
pdf("Cytokine profile Female 5xFAD Only with Side Bar for Abeta.pdf", width=12,height=10,pointsize = 6)
#png("Cytokine profile Female 5xFAD Only_with side bar.png", width=16,height=16,units="in",res=600, pointsize = 20)
heatmap3(Df_Cyt_SubZ,main="Female 5xFAD Cytokine Profile",col=barColors, breaks=breakBarColors,legendfun=function()showLegend(),
         Rowv=NA, Colv=as.dendrogram(hc),  scale="none", na.rm = TRUE,
         labRow = rownames(Df_Cyt_SubZ), labCol = colnames(Df_Cyt_SubZ),
         RowSideColors=colorBarCombined)


dev.off()




new_meanSEM3 <- function(x,y,colors=x,x_str="",y_str="",title=y_str,axis.text.x.size=20,axis.text.y.size=20,axis.text.x.angle=90,title.size=20,text.size=20){
  dataPlot <- data.frame(factor(x),y,factor(colors))
  colnames(dataPlot)=c("x_dot","y_dot","colors_dot")
  plotOut <- ggplot(dataPlot,aes(x=x_dot,y=y_dot)) +
    geom_boxplot(aes(color=colors))+
    geom_point(aes(colour=colors,fill=colors), size=6)+
    geom_beeswarm(aes(colour=colors,fill = colors,groupOnX = TRUE))+
    #stat_summary(geom="errorbar",color="black",size=1.5,aes(width=0.3),)+
    #stat_summary(fun = "mean", geom="point", size=4,show.legend = F, color=colors) +
    theme(panel.background = element_rect(fill = 'transparent'),
          panel.grid.major=element_blank(),
          axis.line=element_line(color="black"),
          text = element_text(size=text.size),
          legend.key = element_blank(),
          axis.text.x = element_text(color = "black",
                                     size = axis.text.x.size,angle= axis.text.x.angle),
          axis.text.y = element_text( color = "black",
                                      size = axis.text.y.size),
          plot.title=element_text(hjust=0.5,face = "bold",size=title.size))+
    xlab(x_str)+
    ylab(y_str)+
    #  ylim(c(-2, 2))+
    ggtitle(paste(title))
  return(plotOut)
}

new_meanSEM4 <- function(x, y, colors = x, 
                         x_str = "", y_str = "", title = y_str,
                         axis.text.x.size = 20, axis.text.y.size = 20, 
                         axis.text.x.angle = 90, title.size = 20, text.size = 20) {
  
  dataPlot <- data.frame(x_dot = factor(x), 
                         y_dot = y, 
                         colors_dot = factor(colors))
  
  plotOut <- ggplot(dataPlot, aes(x = x_dot, y = y_dot)) +
    
    geom_boxplot(aes(color = colors_dot), size = 1.5, outlier.shape = NA) +
    
    geom_beeswarm(aes(color = colors_dot, fill = colors_dot), 
                  groupOnX = TRUE, size = 6, shape = 21, stroke = 0.8) +
    
    
    # geom_point(aes(color = colors_dot, fill = colors_dot), 
    #  size = 4.5, shape = 21, stroke = 0.8, alpha = 0.7, position = position_jitter(width = 0.1)) +
    
    # mean + error bar
    #stat_summary(fun = mean, geom = "point", shape = 23, size = 8, fill = "white", color = "black") +
    # stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3, color = "black", size = 1.2) +
    
    scale_color_manual(values = c("NoFlicker" = "#E298D6", "20Hz" = "#1A1513", "40Hz" = "#76A33E")) +
    scale_fill_manual(values = c("NoFlicker" = "#E298D6", "20Hz" = "#1A1513", "40Hz" = "#76A33E")) +
    
    theme_bw(base_size = text.size) +
    theme(
      panel.border = element_rect(color = "black", fill = NA, size = 1.2),
      axis.line = element_line(color = "black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(color = "black", size = axis.text.x.size, angle = axis.text.x.angle),
      axis.text.y = element_text(color = "black", size = axis.text.y.size),
      plot.title = element_text(hjust = 0.5, face = "bold", size = title.size),
      legend.position = "none"  # Turn off legend unless you want to keep it
    ) +
    xlab(x_str) +
    ylab(y_str) +
    ggtitle(title)
  
  return(plotOut)
}

#----Make bar  plots for each cytokine


dataCyt_groupVars=data.frame(Df_Cyt_Sub, Exp_data_sub)


# 
# for(i in 1:33){
#   plotOut <- new_meanSEM3(x=factor(Exp_data_sub$Flicker, levels = c("NoFlicker", "20Hz", "40Hz")),
#                           y=Df_Cyt_Sub[,i],
#                           colors = factor(Exp_data_sub$Flicker, levels = c("NoFlicker", "20Hz", "40Hz")),
#                           title = colnames(Df_Cyt_Sub)[i],y_str = "Flourescent Int (a.u.)",
#                           axis.text.x.size = 12,axis.text.y.size = 12,title.size = 15,text.size = 15)+
#     
#   #  scale_color_manual(values=c("#3D92CF","#B385BC","#3F8A3F","#66D8FF","#A579FF","#3FFF9F"))+
#    # scale_fill_manual(values=c("#3D92CF","#B385BC","#3F8A3F","#66D8FF","#A579FF","#3FFF9F"))+
#      stat_compare_means(method="anova",p.adjust.method = "bonferroni",label="p.format",tip.length = 0.03)
#  #  pdf(paste0("./mCytokine/DotBox plots All/Zscored/Female ",colnames(dataCyt)[i],".pdf"),height=5,width=6,useDingbats = FALSE,useKerning = FALSE)
# png(paste0("./mCytokine/DotBox plots/NonZscored/Female/5xFAD/", colnames(dataCytZ)[i],".png"),height=10,width=12,units="in",res=600)
#   
#   print(plotOut)
#   dev.off()
# }
# 


for (i in 1:33) {
  flicker_levels <- c("NoFlicker", "20Hz", "40Hz")
  condition_ordered <- factor(Exp_data_sub$Flicker, levels = flicker_levels)
  
  if (length(unique(condition_ordered)) < 2) {
    message(paste("Skipping", colnames(Df_Cyt_Sub)[i], "- insufficient group levels"))
    next
  }
  
  plotOut <- new_meanSEM4(
    x = condition_ordered,
    y = Df_Cyt_Sub[, i],
    colors = condition_ordered,
    title = colnames(Df_Cyt_Sub)[i],
    y_str = "Fluorescent Int (a.u.)",
    axis.text.x.size = 12,
    axis.text.y.size = 12,
    title.size = 15,
    text.size = 15
  ) +
    stat_compare_means(method = "anova", label = "p.format", tip.length = 0.03)
  
  #png(
  #filename = paste0("./mCytokine/DotBox plots/NonZscored/Female/5xFAD/", colnames(dataCytZ)[i], ".png"),
  # height = 10, width = 12, units = "in", res = 600
  # )
  pdf(paste0("./mCytokine/DotBox plots/NonZscored/Female/5xFAD/ ",colnames(Df_Cyt_Sub)[i],".pdf"),height=5,width=6,useDingbats = FALSE,useKerning = FALSE)
  
  
  print(plotOut)
  dev.off()
}
# for (i in 1:37)
# {
#   png(paste0("./mCytokine/DotBox plots/NonZscored/Female/5xFAD/",colnames(dataCyt)[i],".png"),height=5,width=6,units="in",res=600)
#   # pdf(paste0("./Individual Plots/",colnames(dataCyt)[i],".pdf"),height=4,width=5)
#   print({
#     ggplot(datagroupVars, aes(x=Order, y=dataCyt[,i], color=Order)) +
#       #geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2) +
#       geom_boxplot(aes(color=as.factor(Order)))+
#       #geom_line(size=1) +
#       geom_point(size=4, aes(color=as.factor(Order)))+
#       theme(panel.background = element_rect(fill = 'white', colour = 'black'),
#             text = element_text(size=10),
#             axis.text.x = element_text(face = "bold", color = "black",
#                                        size = 8),
#             axis.text.y = element_text(face = "bold", color = "black",
#                                        size = 16),
#             plot.title = element_text(hjust = 0.5))+
#       # scale_color_manual(values=c("purple3","#00A087FF"))+
#       #scale_fill_manual(values=c("purple3","#00A087FF"))+
#       stat_compare_means(method="anova",label="p.format",tip.length = 0.03) +
#       xlab("")+
#       ylab("Fluorescent Int (a.u.)")+
#       labs(color="Treatment")+
#       # ylim(c(-2, 2)) +
#       ggtitle(colnames(datagroupVars)[i])
#     # geom_text()+
#     # annotate("text", label = p)
#   })
#   dev.off()
# }
# 
# 
# 



for (i in 1:4) {
  flicker_levels <- c("NoFlicker", "20Hz", "40Hz")
  condition_ordered <- factor(Exp_data_sub$Flicker, levels = flicker_levels)
  
  if (length(unique(condition_ordered)) < 2) {
    message(paste("Skipping", colnames(Df_Abeta_Sub)[i], "- insufficient group levels"))
    next
  }
  
  plotOut <- new_meanSEM4(
    x = condition_ordered,
    y = Df_Abeta_Sub[, i],
    colors = condition_ordered,
    title = colnames(Df_Abeta_Sub)[i],
    y_str = "Fluorescent Int (a.u.)",
    axis.text.x.size = 12,
    axis.text.y.size = 12,
    title.size = 15,
    text.size = 15
  ) +
    stat_compare_means(method = "anova", label = "p.format", tip.length = 0.03)
  
  #png(
  #filename = paste0("./mCytokine/DotBox plots/NonZscored/Female/5xFAD/", colnames(dataCytZ)[i], ".png"),
  # height = 10, width = 12, units = "in", res = 600
  # )
  pdf(paste0("./mCytokine/DotBox plots/NonZscored/Female/5xFAD/Abeta_ ",colnames(Df_Abeta_Sub)[i],".pdf"),height=5,width=6,useDingbats = FALSE,useKerning = FALSE)
  
  
  print(plotOut)
  dev.off()
}





cor_matrix <- cor(Df_Cyt_Sub,Df_Abeta_Sub ,use = "pairwise.complete.obs", method = "spearman")
nSamples <- nrow(Df_Cyt_Sub)
Pvalue <- corPvalueStudent(cor_matrix, nSamples)

cor_df <- melt(cor_matrix)
colnames(cor_df) <- c("Cytokine", "Amyloid", "Correlation")

pval_df <- melt(Pvalue)
pval_df$padj <- p.adjust(pval_df$Pvalue, method = "fdr")
colnames(pval_df) <- c("Cytokine", "Amyloid", "Pvalue")

df <- merge(cor_df, pval_df, by = c("Cytokine", "Amyloid"))

df$Significance <- cut(df$Pvalue,
                       breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                       labels = c("***", "**", "*", ""),
                       right = FALSE)


df$Label <- paste0(df$Significance)

#-------or--------------
df$Label <- paste0("R=", round(df$Correlation, 2), 
                   #"\n", df$Significance,
                   "\nPvalue=", signif(df$Pvalue, 2))

#-------------
df$padj <- p.adjust(pval_df$Pvalue, method = "fdr")

row_dist <- dist(cor_matrix)
row_clust <- hclust(row_dist)
row_order <- row_clust$labels[row_clust$order]

col_dist <- dist(t(cor_matrix))
col_clust <- hclust(col_dist)
col_order <- col_clust$labels[col_clust$order]

df$Cytokine <- factor(df$Cytokine, levels = row_order)
df$Amyloid <- factor(df$Amyloid, levels = col_order)

# --- PLOT ---
pdf("./mCytokine/SpermanCorrealtion_Abeta_vs_Cytokines_with Pvalues.pdf", width = 20, height = 16, pointsize = 10)
ggplot(df, aes(x = Amyloid, y = Cytokine, fill = Correlation)) +
  geom_tile(color = "white") +
  geom_text(aes(label = Label), size = 2) +
  scale_fill_gradient2(
    low = "blue", mid = "white", high = "red", midpoint = 0,
    limit = c(-1, 1), space = "Lab", name = "Spearman R"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text = element_text(size = 10),
    axis.title = element_blank()
  )

dev.off()



write_csv(df[, c("Cytokine", "Amyloid", "Correlation", "Pvalue", "padj", "Significance")],
          "./mCytokine/SpearmanCorrelation_Abeta vs Cytokines_CorrelationMatrix_Pvals_FDR_SummaryTable.csv")





#----Male 5xFAD---#
sortby <- which(Exp_data$SexGenotype =="5XFAD M")

Df_Cyt_Sub <-  dataCyt[sortby ,] 
Df_Abeta_Sub <-dataAbeta[sortby ,]
Exp_data_sub <- Exp_data[sortby ,]



#Z scoring data

Df_Cyt_SubZ=(apply(t(Df_Cyt_Sub),1,scale)); #z-score the data
Df_Cyt_SubZ[is.na(Df_Cyt_SubZ)]=0
rownames(Df_Cyt_SubZ)=rownames(Df_Cyt_Sub)





FlickerColorBar <- case_when(
  Exp_data_sub$Flicker == "20Hz" ~ "#1A1513",
  Exp_data_sub$Flicker == "40Hz" ~ "#76A33E",
  Exp_data_sub$Flicker == "NoFlicker" ~ "#E298D6"
)


colorBarCombined=cbind(c(FlickerColorBar),c(FlickerColorBar))
colnames(colorBarCombined)=c("Flicker","Flicker")


#Cluster the cytokine data by columns
hc= hclust(dist(t(Df_Cyt_SubZ),method = "euclidean"), method="average")

##Make Heatmap
pdf("Cytokine profile Male 5xFAD Only.pdf", width=12,height=10,pointsize = 6)
#png("Cytokine profile Male 5xFAD Only.png", width=16,height=16,units="in",res=600, pointsize = 20)
heatmap3(Df_Cyt_SubZ,main="Female 5xFAD Cytokine Profile",col=barColors, breaks=breakBarColors,legendfun=function()showLegend(),
         Rowv=NA, Colv=as.dendrogram(hc),  scale="none", na.rm = TRUE,
         labRow = rownames(Df_Cyt_SubZ), labCol = colnames(Df_Cyt_SubZ),
         RowSideColors=colorBarCombined)


dev.off()




















#----Female 5xFAD---#
sortby <- which(Exp_data$Genotype =="5XFAD")
                  
Df_Cyt_Sub <-  dataCyt[sortby ,] 
Df_Abeta_Sub <-dataAbeta[sortby ,]
Exp_data_sub <- Exp_data[sortby ,]



#Z scoring data

Df_Cyt_SubZ=(apply(t(Df_Cyt_Sub),1,scale)); #z-score the data
Df_Cyt_SubZ[is.na(Df_Cyt_SubZ)]=0
rownames(Df_Cyt_SubZ)=rownames(Df_Cyt_Sub)



Litterpal<-brewer.pal(8,"Paired")
LitterColorBar <- case_when(
  Exp_data_sub$litter == "a" ~ Litterpal[1],
  Exp_data_sub$litter == "b" ~ Litterpal[2],
  Exp_data_sub$litter == "c" ~ Litterpal[3],
  Exp_data_sub$litter == "d" ~ Litterpal[4],
  Exp_data_sub$litter == "e" ~ Litterpal[5]
)


FlickerColorBar <- case_when(
  Exp_data_sub$Flicker == "20Hz" ~ "black",
  Exp_data_sub$Flicker == "40Hz" ~ "#76A33E",
  Exp_data_sub$Flicker == "NoFlicker" ~ "#E298D6"
)


colorBarCombined=cbind(c(FlickerColorBar),c(FlickerColorBar))
colnames(colorBarCombined)=c("Flicker","Flicker")


#Cluster the cytokine data by columns
hc= hclust(dist(t(Df_Cyt_SubZ),method = "euclidean"), method="average")

##Make Heatmap
pdf("Cytokine profile Female 5xFAD Only.pdf", width=12,height=10,pointsize = 6)
#png("Cytokine profile Female 5xFAD Only.png", width=16,height=16,units="in",res=600, pointsize = 20)
heatmap3(Df_Cyt_SubZ,main="Female 5xFAD Cytokine Profile",col=barColors, breaks=breakBarColors,legendfun=function()showLegend(),
         Rowv=NA, Colv=as.dendrogram(hc),  scale="none", na.rm = TRUE,
         labRow = rownames(Df_Cyt_SubZ), labCol = colnames(Df_Cyt_SubZ),
         RowSideColors=colorBarCombined)


dev.off()






#Amyloid beta side bar


ab42<-Df_Abeta_Sub[,2]
ab40<-Df_Abeta_Sub[,1]
ratio<-Df_Abeta_Sub[,3]


breakBarColorsAbeta42=c(-1,seq(0, 1400, 0.01))*10 #Outside numbers clip outliers. This is for zscoring.
barColorsAbeta42Vec = c("black",colorpanel(length(breakBarColorsAbeta42), "white", "gold4"))
indMatchColorone=match(round(ab42*10),round(breakBarColorsAbeta42))
barColorsAb42=barColorsAbeta42Vec[indMatchColorone]

breakBarColorsAbeta40=c(-1,seq(0, 12, 0.01))*10 #Outside numbers clip outliers. This is for zscoring.
barColorsAbeta40Vec = c("black",colorpanel(length(breakBarColorsAbeta40), "white", "orange2"))
indMatchColortwo=match(round(ab40*10),round(breakBarColorsAbeta40))
barColorsAb40=barColorsAbeta40Vec[indMatchColortwo]


breakBarColorsAbetaratio=c(-1,seq(0, 180, 0.01))*10 #Outside numbers clip outliers. This is for zscoring.
barColorsAbetaratioVec = c("black",colorpanel(length(breakBarColorsAbetaratio), "white", "gold2"))
indMatchColorthree=match(round(ratio*10),round(breakBarColorsAbetaratio))
barColorsAbratio=barColorsAbetaratioVec[indMatchColorthree]

colorBarCombined=cbind(c(FlickerColorBar),c(LitterColorBar),c(barColorsAb42),c(barColorsAb40),c(barColorsAbratio))
colnames(colorBarCombined)=c("Flicker","Litter","Ab42","Ab40","ratio")



##Make Heatmap
pdf("Cytokine profile Female 5xFAD Only with Side Bar for Abeta.pdf", width=12,height=10,pointsize = 6)
#png("Cytokine profile Female 5xFAD Only_with side bar.png", width=16,height=16,units="in",res=600, pointsize = 20)
heatmap3(Df_Cyt_SubZ,main="Female 5xFAD Cytokine Profile",col=barColors, breaks=breakBarColors,legendfun=function()showLegend(),
         Rowv=NA, Colv=as.dendrogram(hc),  scale="none", na.rm = TRUE,
         labRow = rownames(Df_Cyt_SubZ), labCol = colnames(Df_Cyt_SubZ),
         RowSideColors=colorBarCombined)


dev.off()






df_PCA <- Df_Cyt_SubZ
PCA_mat <- df_PCA[,] %>% as.matrix()
#PathVar=pathData[, 1]
#df_var=as.numeric(log2(PathVar))

#na_indices <- which(is.na(df_var), arr.ind = TRUE)



oplsOut <- opls(PCA_mat,predI=2)
rotateOut <- rotate.opls(oplsOut,degrees=0)
scores_data=data.frame(rotateOut$T1,rotateOut$T2)
# 


pdfPath <- paste("./ScorePlot_OPLS_5xFAD_Female_v2.pdf", sep = '')
pdf(pdfPath, width = 5, height = 4, pointsize = 14)
p<- scores.plot(rotateOut$T1,rotateOut$T2,factor(Exp_data_sub$Flicker),"Flicker")

print(p)
dev.off()

# PathVar=Df_Abeta_Sub[,1]
# df_var=as.numeric(PathVar)
# 
# oplsOut <- opls(PCA_mat,df_var,predI=2)
# rotateOut <- rotate.opls(oplsOut,degrees=0)
# scores_data=data.frame(rotateOut$T1,rotateOut$T2)
# 


pdfPath <- paste("./.pdf", sep = '')
pdf(pdfPath, width = 5, height = 4, pointsize = 14)
scores.plot(rotateOut$T1,rotateOut$T2,df_var,"df_var")
print(p)
dev.off()


#Partial Least Squares Regression PLSR
df_PCA <- Df_Cyt_SubZ
PCA_mat <- as.matrix(df_PCA) 
PathVar <- Df_Abeta_Sub[, 1]  
df_var <- as.numeric(PathVar) 


plsr_model <- plsr(df_var ~ PCA_mat, ncomp = 2, validation = "CV")


scores_data <- data.frame(T1 = plsr_model$scores[, 1], T2 = plsr_model$scores[, 2])



pdfPath <- paste("./FiguresOut/OPLS and PLSR/Scoreplot_PLSR_GFAP", i, ".pdf", sep = '')
pdf(pdfPath, width = 5, height = 4, pointsize = 14)
scores.plot(scores_data$T1,scores_data$T2,PathVar,"PathVar")
;
dev.off()



pdfPath <- paste("./FiguresOut/OPLS and PLSR/PLSR_GFAP_LV1", i, ".pdf", sep = '')
pdf(pdfPath, width = 5, height = 4, pointsize = 14)
loadings.chart(rotateOut$P1,"LV1");
dev.off()

pdfPath <- paste("./FiguresOut/OPLS and PLSR/PLSR_GFAP_LV2", i, ".pdf", sep = '')
pdf(pdfPath, width = 5, height = 4, pointsize = 14)
loadings.chart(rotateOut$P2,"LV2");
dev.off()




score_df <- data.frame(
  PC1 = rotateOut$T1,
  PC2 = rotateOut$T2,
  Flicker = factor(Exp_data_sub$Flicker, levels = c("NoFlicker", "20Hz", "40Hz"))
)


flicker_colors <- c("NoFlicker" = "#E298D6", "20Hz" = "#1A1513", "40Hz" = "#76A33E")

p <- ggplot(score_df, aes(x = PC1, y = PC2, color = Flicker)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = flicker_colors) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  theme_minimal(base_size = 14) +
  labs(x = "PC1", y = "PC2", title = "Score Plot") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.title = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.8),
    panel.grid = element_blank()
  )

pdfPath <- "./ScorePlot_OPLS_5xFAD_Female.pdf"
pdf(pdfPath, width = 5, height = 4, pointsize = 14)
#png(paste0("./ScorePlot_OPLS_5xFAD_Fem.png"),height=4,width=5,units="in",res=600)

print(p)
dev.off()







pdfPath <- paste("./PC1_5xFAD_Cytokines_Female.pdf", sep = '')
pdf(pdfPath, width = 5, height = 4, pointsize = 6)

png(paste0("./PC1_5xFAD_Cytokines_Female.png"),height=10,width=12,units="in",res=600)
loadings.chart(rotateOut$P1,"PC1");

dev.off()

pdfPath <- paste("./PC2_5xFAD_Cytokines_Female.pdf", sep = '')
pdf(pdfPath, width = 5, height = 4, pointsize = 6)
png(paste0("./PC2_5xFAD_Cytokines_Female.png"),height=10,width=12,units="in",res=600)
loadings.chart(rotateOut$P2,"PC2");
dev.off()







#PCA Analaysis

#----PCA---
df_PCA <- Df_Cyt_SubZ
PCA_mat <- df_PCA[,] %>% as.matrix()

oplsOut <- opls(PCA_mat,predI=2)
rotateOut <- rotate.opls(oplsOut,degrees=0)
scores_data=data.frame(rotateOut$T1,rotateOut$T2,Exp_data_sub$Flicker)
Flicker=factor(Exp_data_sub$Flicker)

oplsVar = oplsOut@modelDF$R2X * 100





df_plot = data.frame(T1=rotateOut$T1,T2=rotateOut$T2,Condition=Exp_data_sub$Flicker)


statsOut_T1 = summarySE(df_plot,measurevar="T1", groupvars="Flicker")
statsOut_T2 = summarySE(df_plot,measurevar="T2", groupvars="Flicker")


export(data.frame(statsOut_T1),file="mCytokine/PCA/PCA new/PC1 Stats_Female 5xFAD.csv")
export(data.frame(statsOut_T2),file="mCytokine/PCA/PCA new/PC2 Stats Female 5xFAD.csv")


pdf(paste0("mCytokine/PCA/PCA new/PC1_Violin_Female 5xFAD.pdf"),height=4,width=6,useDingbats = FALSE)
ggplot(data=df_plot, aes( y = T1, x = Condition, fill = Condition)) +
  geom_violin()+
  
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3, color = "black", position = position_dodge(width = 0.9)) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.2, color = "black", position = position_dodge(width = 0.9)) +
  
  #geom_errorbar(data=statsOut_T1, aes(x=Condition, ymin=mean-se, ymax=mean +se),
  # width=0.4, colour="black", alpha=1, size=1) +
  #geom_crossbar(data=statsOut_T1, aes(x=Condition,ymin = T1, ymax = T1),
  #  size=0.5,col="red", width = .5)+
  theme(panel.background = element_rect(fill = 'white', colour = 'black'),
        text = element_text(size=20), 
        axis.text.x = element_text(face = "plain", color = "black",
                                   size = 10,angle=45,hjust=1),
        axis.text.y = element_text(face = "plain", color = "black",
                                   size = 10),
        legend.position = "none", plot.title = element_text(hjust = 0.5))+
  xlab("")+
  ylab(paste0("Scores on PC1 (",oplsVar[1],"%)"))+
  #ylim(c(-5, 5))+
  scale_fill_manual(values = c("#1A1513","#76A33E","#E298D6"))+
  stat_compare_means(method="anova",p.adjust.methods = "bonferroni",label="p.format",tip.length = 0.03)
#coord_flip()
dev.off()




pdf(paste0("mCytokine/PCA/PCA new/PC2_Violin_Female 5xFAD.pdf"),height=4,width=6,useDingbats = FALSE)
ggplot(data=df_plot, aes( y = T2, x = Condition, fill = Condition)) +
  geom_violin()+
  
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3, color = "black", position = position_dodge(width = 0.9)) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.2, color = "black", position = position_dodge(width = 0.9)) +
  
  #geom_errorbar(data=statsOut_T1, aes(x=Condition, ymin=mean-se, ymax=mean +se),
  # width=0.4, colour="black", alpha=1, size=1) +
  #geom_crossbar(data=statsOut_T1, aes(x=Condition,ymin = T1, ymax = T1),
  #  size=0.5,col="red", width = .5)+
  theme(panel.background = element_rect(fill = 'white', colour = 'black'),
        text = element_text(size=20), 
        axis.text.x = element_text(face = "plain", color = "black",
                                   size = 10,angle=45,hjust=1),
        axis.text.y = element_text(face = "plain", color = "black",
                                   size = 10),
        legend.position = "none", plot.title = element_text(hjust = 0.5))+
  xlab("")+
  ylab(paste0("Scores on PC2 (",oplsVar[2],"%)"))+
  ylim(c(-5, 5))+
  #scale_fill_manual(values = "#E298D6","#1A1513","#76A33E")+
  scale_fill_manual(values = c("#1A1513","#76A33E","#E298D6"))+
  stat_compare_means(method="anova",p.adjust.methods = "bonferroni",label="p.format",tip.length = 0.03)
#coord_flip()
dev.off()



pdfPath <- paste("./mCytokine/PCA/PCA new/ScorePlot_PCA_Female 5xFAD.pdf", sep = '')
pdf(pdfPath, width = 10, height = 8, pointsize = 14)
p <- scores.plot(rotateOut$T1,rotateOut$T2,Exp_data_sub$Flicker,"Flicker",analysis = "PCA")
ggpar(p,palette=c("#1A1513","#76A33E","#E298D6"))

print(p)
dev.off()

#blah

pdfPath <- paste("./mCytokine/PCA/PCA new/PC1_PCA_Female 5xFAD.pdf", sep = '')
pdf(pdfPath, width = 5, height = 4, pointsize = 6)

#png(paste0("./mCytokine/PCA/PCA new/PC1_PCA_Female 5xFAD.png"),height=10,width=12,units="in",res=600)
loadings.chart(rotateOut$P1,"PC1");

dev.off()

pdfPath <- paste("./mCytokine/PCA/PCA new/PC2_PCA_Female 5xFAD.pdf", sep = '')
pdf(pdfPath, width = 5, height = 4, pointsize = 6)
#png(paste0("./PC2_PCA_All_Group1.png"),height=10,width=12,units="in",res=600)
loadings.chart(rotateOut$P2,"PC2");
dev.off()


#Correlation of PC1 vs Amyloid beta

ind40=which(Exp_data_sub$Flicker=="40Hz")
indCTRL=which(Exp_data_sub$Flicker=="NoFlicker")

df_var=as.numeric(rotateOut$T1)
x_regress=df_var
x_regress=x_regress[indCTRL]
x_label="PC1"


for(i in 1:ncol(Df_Abeta_Sub)){
  y_regress=Df_Abeta_Sub[indCTRL,i] 
  y_label=colnames(Df_Abeta_Sub)[i]
  dataReg=data.frame(y_regress,x_regress)
  linmodel=lm(y_regress~x_regress,dataReg)
  RSQ=summary(linmodel)$r.squared
  f <- summary(linmodel)$fstatistic
  pval <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(pval) <- NULL
  
  outputFolder<-"mCytokine/PCA/PCA new"
  
  #Plot regression
  Flicker=as.factor(Exp_data_sub$Flicker[indCTRL])
  
  #pdfPath <- paste("./mCytokine/PCA/PCA new/PC1_Regression_Abeta.pdf", sep = '')
  pdf(paste("mCytokine/PCA/PCA new","/",y_label,' vs ',x_label,"_Ctrl.pdf",sep= ""), width = 5, height = 4, pointsize = 6)
  
 # pdf("mCytokine/PCA/PCA new/40Hz.pdf",width=10,height=10,pointsize = 10)
  
 #png(paste(outputFolder,"/",y_label,' vs ',x_label,"_40Hz.png",sep= ""), width=6,height=4,units="in",res=600,pointsize = 10)
  print({
  ggplot(dataReg, aes(x=x_regress, y=y_regress)) +
      geom_smooth(method='lm',se=FALSE,color="#76A33E",linetype="longdash")+
      # geom_rug(aes(color =condition))+
      geom_point(aes(color="#76A33E"),size=4)+
      #xlim(0.95*min(x_regress),1.05*max(x_regress))+
      #ylim(0.95*min(y_regress),1.05*max(y_regress))+
      theme(panel.background = element_rect(fill = 'white'),
            text = element_text(size=16),
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
      #annotate("text",size=6, x=0.8*max(x_regress)+0.2*min(x_regress), y=0.1*min(y_regress)+.9*max(y_regress), label=bquote(R^2==.(round(RSQ,digits=2))))+
      #annotate("text",size=6, x=0.8*max(x_regress)+0.2*min(x_regress), y=0*min(y_regress)+1*max(y_regress), label= bquote(p==.(round(pval,digits=4))))+
      #annotate("text",size=6, x=0.8*max(x_regress)+0.2*min(x_regress), y=0.98*min(y_regress)+0.02*max(y_regress), label=bquote(R^2==.(round(RSQ,digits=2))))+
      #annotate("text",size=6, x=0.8*max(x_regress)+0.2*min(x_regress), y=0.9*min(y_regress)+.1*max(y_regress), label= bquote(p==.(round(pval,digits=4))))+
      ggtitle(paste(y_label," vs ",x_label,sep=""))+
      ggpubr::stat_cor(method="pearson")
  })
  
# + ggpar(p,palette=c("#76A33E"))
  dev.off()
  
}



#Correlation of PC1 vs Amyloid beta

ind40=which(Exp_data_sub$Flicker=="40Hz")
indCTRL=which(Exp_data_sub$Flicker=="NoFlicker")

df_var=as.numeric(Df_Abeta_Sub[,1])
x_regress=df_var
x_regress=x_regress[ind40]
x_label="Ab40"


for(i in 1:ncol(Df_Cyt_Sub)){
  y_regress=Df_Cyt_Sub[ind40,i] 
  y_label=colnames(Df_Cyt_Sub)[i]
  dataReg=data.frame(y_regress,x_regress)
  linmodel=lm(y_regress~x_regress,dataReg)
  RSQ=summary(linmodel)$r.squared
  f <- summary(linmodel)$fstatistic
  pval <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(pval) <- NULL
  
  outputFolder<-"Abeta vs Cytokines"
  
  #Plot regression
  Flicker=as.factor(Exp_data_sub$Flicker[ind40])
  
  pdf(paste("Abeta vs Cytokines","/",y_label,' vs ',x_label,"_40Hz_ABETA40.pdf",sep= ""), width = 5, height = 4, pointsize = 6)
  #png(paste(outputFolder,"/",y_label,' vs ',x_label,"40Hz.png",sep= ""), width=6,height=4,units="in",res=600,pointsize = 10)
  print({
    ggplot(dataReg, aes(x=x_regress, y=y_regress)) +
      geom_smooth(method='lm',se=FALSE,color="#76A33E",linetype="longdash")+
      # geom_rug(aes(color =condition))+
      geom_point(color="#76A33E",size=4)+
      #xlim(0.95*min(x_regress),1.05*max(x_regress))+
      #ylim(0.95*min(y_regress),1.05*max(y_regress))+
      theme(panel.background = element_rect(fill = 'white'),
            text = element_text(size=16),
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
      #annotate("text",size=6, x=0.8*max(x_regress)+0.2*min(x_regress), y=0.1*min(y_regress)+.9*max(y_regress), label=bquote(R^2==.(round(RSQ,digits=2))))+
      #annotate("text",size=6, x=0.8*max(x_regress)+0.2*min(x_regress), y=0*min(y_regress)+1*max(y_regress), label= bquote(p==.(round(pval,digits=4))))+
      #annotate("text",size=6, x=0.8*max(x_regress)+0.2*min(x_regress), y=0.98*min(y_regress)+0.02*max(y_regress), label=bquote(R^2==.(round(RSQ,digits=2))))+
      #annotate("text",size=6, x=0.8*max(x_regress)+0.2*min(x_regress), y=0.9*min(y_regress)+.1*max(y_regress), label= bquote(p==.(round(pval,digits=4))))+
      ggtitle(paste(y_label," vs ",x_label,sep=""))+
      ggpubr::stat_cor(method="pearson")
  })
  
  # + ggpar(p,palette=c("#76A33E"))
  dev.off()
  
}





#Amyloid beta side bar


ab42<-Df_Abeta_Sub[,2]
ab40<-Df_Abeta_Sub[,1]
ratio<-Df_Abeta_Sub[,3]


breakBarColorsAbeta42=c(-1,seq(0, 1200, 0.1))*10 #Outside numbers clip outliers. This is for zscoring.
barColorsAbeta42Vec = c("black",colorpanel(length(breakBarColorsAbeta42), "white", "gold4"))
indMatchColorone=match(round(ab42*10),round(breakBarColorsAbeta42))
barColorsAb42=barColorsAbeta42Vec[indMatchColorone]

breakBarColorsAbeta40=c(-1,seq(0, 10, 0.01))*10 #Outside numbers clip outliers. This is for zscoring.
barColorsAbeta40Vec = c("black",colorpanel(length(breakBarColorsAbeta40), "white", "orange"))
indMatchColortwo=match(round(ab40*10),round(breakBarColorsAbeta40))
barColorsAb40=barColorsAbeta40Vec[indMatchColortwo]


breakBarColorsAbetaratio=c(-1,seq(0, 180, 0.01))*10 #Outside numbers clip outliers. This is for zscoring.
barColorsAbetaratioVec = c("black",colorpanel(length(breakBarColorsAbetaratio), "white", "gray4"))
indMatchColorthree=match(round(ratio*10),round(breakBarColorsAbetaratio))
barColorsAbratio=barColorsAbetaratioVec[indMatchColorthree]

colorBarCombined=cbind(c(FlickerColorBar),c(barColorsAb42),c(barColorsAb40),c(barColorsAbratio))
colnames(colorBarCombined)=c("Flicker","Ab42","Ab40","ratio")



##Make Heatmap
pdf("Cytokine profile Female 5xFAD Only with Side Bar for Abeta.pdf", width=12,height=10,pointsize = 6)
#png("Cytokine profile Female 5xFAD Only_with side bar.png", width=16,height=16,units="in",res=600, pointsize = 20)
heatmap3(Df_Cyt_SubZ,main="Female 5xFAD Cytokine Profile",col=barColors, breaks=breakBarColors,legendfun=function()showLegend(),
         Rowv=NA, Colv=as.dendrogram(hc),  scale="none", na.rm = TRUE,
         labRow = rownames(Df_Cyt_SubZ), labCol = colnames(Df_Cyt_SubZ),
         RowSideColors=colorBarCombined)


dev.off()




new_meanSEM3 <- function(x,y,colors=x,x_str="",y_str="",title=y_str,axis.text.x.size=20,axis.text.y.size=20,axis.text.x.angle=90,title.size=20,text.size=20){
  dataPlot <- data.frame(factor(x),y,factor(colors))
  colnames(dataPlot)=c("x_dot","y_dot","colors_dot")
  plotOut <- ggplot(dataPlot,aes(x=x_dot,y=y_dot)) +
    geom_boxplot(aes(color=colors))+
    geom_point(aes(colour=colors,fill=colors), size=6)+
    geom_beeswarm(aes(colour=colors,fill = colors,groupOnX = TRUE))+
    #stat_summary(geom="errorbar",color="black",size=1.5,aes(width=0.3),)+
    #stat_summary(fun = "mean", geom="point", size=4,show.legend = F, color=colors) +
    theme(panel.background = element_rect(fill = 'transparent'),
          panel.grid.major=element_blank(),
          axis.line=element_line(color="black"),
          text = element_text(size=text.size),
          legend.key = element_blank(),
          axis.text.x = element_text(color = "black",
                                     size = axis.text.x.size,angle= axis.text.x.angle),
          axis.text.y = element_text( color = "black",
                                      size = axis.text.y.size),
          plot.title=element_text(hjust=0.5,face = "bold",size=title.size))+
    xlab(x_str)+
    ylab(y_str)+
  #  ylim(c(-2, 2))+
    ggtitle(paste(title))
  return(plotOut)
}

new_meanSEM4 <- function(x, y, colors = x, 
                         x_str = "", y_str = "", title = y_str,
                         axis.text.x.size = 20, axis.text.y.size = 20, 
                         axis.text.x.angle = 90, title.size = 20, text.size = 20) {
  
  dataPlot <- data.frame(x_dot = factor(x), 
                         y_dot = y, 
                         colors_dot = factor(colors))
  
  plotOut <- ggplot(dataPlot, aes(x = x_dot, y = y_dot)) +
    
    geom_boxplot(aes(color = colors_dot), size = 1.5, outlier.shape = NA) +
    
    geom_beeswarm(aes(color = colors_dot, fill = colors_dot), 
                  groupOnX = TRUE, size = 6, shape = 21, stroke = 0.8) +
    
 
   # geom_point(aes(color = colors_dot, fill = colors_dot), 
             #  size = 4.5, shape = 21, stroke = 0.8, alpha = 0.7, position = position_jitter(width = 0.1)) +
    
   # mean + error bar
     #stat_summary(fun = mean, geom = "point", shape = 23, size = 8, fill = "white", color = "black") +
   # stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3, color = "black", size = 1.2) +
    
    scale_color_manual(values = c("NoFlicker" = "#E298D6", "20Hz" = "#1A1513", "40Hz" = "#76A33E")) +
    scale_fill_manual(values = c("NoFlicker" = "#E298D6", "20Hz" = "#1A1513", "40Hz" = "#76A33E")) +
    
    theme_bw(base_size = text.size) +
    theme(
      panel.border = element_rect(color = "black", fill = NA, size = 1.2),
      axis.line = element_line(color = "black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(color = "black", size = axis.text.x.size, angle = axis.text.x.angle),
      axis.text.y = element_text(color = "black", size = axis.text.y.size),
      plot.title = element_text(hjust = 0.5, face = "bold", size = title.size),
      legend.position = "none"  # Turn off legend unless you want to keep it
    ) +
    xlab(x_str) +
    ylab(y_str) +
    ggtitle(title)
  
  return(plotOut)
}

#----Make bar  plots for each cytokine


dataCyt_groupVars=data.frame(Df_Cyt_Sub, Exp_data_sub)


# 
# for(i in 1:33){
#   plotOut <- new_meanSEM3(x=factor(Exp_data_sub$Flicker, levels = c("NoFlicker", "20Hz", "40Hz")),
#                           y=Df_Cyt_Sub[,i],
#                           colors = factor(Exp_data_sub$Flicker, levels = c("NoFlicker", "20Hz", "40Hz")),
#                           title = colnames(Df_Cyt_Sub)[i],y_str = "Flourescent Int (a.u.)",
#                           axis.text.x.size = 12,axis.text.y.size = 12,title.size = 15,text.size = 15)+
#     
#   #  scale_color_manual(values=c("#3D92CF","#B385BC","#3F8A3F","#66D8FF","#A579FF","#3FFF9F"))+
#    # scale_fill_manual(values=c("#3D92CF","#B385BC","#3F8A3F","#66D8FF","#A579FF","#3FFF9F"))+
#      stat_compare_means(method="anova",p.adjust.method = "bonferroni",label="p.format",tip.length = 0.03)
#  #  pdf(paste0("./mCytokine/DotBox plots All/Zscored/Female ",colnames(dataCyt)[i],".pdf"),height=5,width=6,useDingbats = FALSE,useKerning = FALSE)
# png(paste0("./mCytokine/DotBox plots/NonZscored/Female/5xFAD/", colnames(dataCytZ)[i],".png"),height=10,width=12,units="in",res=600)
#   
#   print(plotOut)
#   dev.off()
# }
# 


for (i in 1:33) {
  flicker_levels <- c("NoFlicker", "20Hz", "40Hz")
  condition_ordered <- factor(Exp_data_sub$Flicker, levels = flicker_levels)
  
  if (length(unique(condition_ordered)) < 2) {
    message(paste("Skipping", colnames(Df_Cyt_Sub)[i], "- insufficient group levels"))
    next
  }
  
  plotOut <- new_meanSEM4(
    x = condition_ordered,
    y = Df_Cyt_Sub[, i],
    colors = condition_ordered,
    title = colnames(Df_Cyt_Sub)[i],
    y_str = "Fluorescent Int (a.u.)",
    axis.text.x.size = 12,
    axis.text.y.size = 12,
    title.size = 15,
    text.size = 15
  ) +
    stat_compare_means(method = "anova", label = "p.format", tip.length = 0.03)
  
  #png(
    #filename = paste0("./mCytokine/DotBox plots/NonZscored/Female/5xFAD/", colnames(dataCytZ)[i], ".png"),
   # height = 10, width = 12, units = "in", res = 600
 # )
 pdf(paste0("./mCytokine/DotBox plots/NonZscored/Female/5xFAD/ ",colnames(Df_Cyt_Sub)[i],".pdf"),height=5,width=6,useDingbats = FALSE,useKerning = FALSE)
  
  
  print(plotOut)
  dev.off()
}
# for (i in 1:37)
# {
#   png(paste0("./mCytokine/DotBox plots/NonZscored/Female/5xFAD/",colnames(dataCyt)[i],".png"),height=5,width=6,units="in",res=600)
#   # pdf(paste0("./Individual Plots/",colnames(dataCyt)[i],".pdf"),height=4,width=5)
#   print({
#     ggplot(datagroupVars, aes(x=Order, y=dataCyt[,i], color=Order)) +
#       #geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2) +
#       geom_boxplot(aes(color=as.factor(Order)))+
#       #geom_line(size=1) +
#       geom_point(size=4, aes(color=as.factor(Order)))+
#       theme(panel.background = element_rect(fill = 'white', colour = 'black'),
#             text = element_text(size=10),
#             axis.text.x = element_text(face = "bold", color = "black",
#                                        size = 8),
#             axis.text.y = element_text(face = "bold", color = "black",
#                                        size = 16),
#             plot.title = element_text(hjust = 0.5))+
#       # scale_color_manual(values=c("purple3","#00A087FF"))+
#       #scale_fill_manual(values=c("purple3","#00A087FF"))+
#       stat_compare_means(method="anova",label="p.format",tip.length = 0.03) +
#       xlab("")+
#       ylab("Fluorescent Int (a.u.)")+
#       labs(color="Treatment")+
#       # ylim(c(-2, 2)) +
#       ggtitle(colnames(datagroupVars)[i])
#     # geom_text()+
#     # annotate("text", label = p)
#   })
#   dev.off()
# }
# 
# 
# 



for (i in 1:4) {
  flicker_levels <- c("NoFlicker", "20Hz", "40Hz")
  condition_ordered <- factor(Exp_data_sub$Flicker, levels = flicker_levels)
  
  if (length(unique(condition_ordered)) < 2) {
    message(paste("Skipping", colnames(Df_Abeta_Sub)[i], "- insufficient group levels"))
    next
  }
  
  plotOut <- new_meanSEM4(
    x = condition_ordered,
    y = Df_Abeta_Sub[, i],
    colors = condition_ordered,
    title = colnames(Df_Abeta_Sub)[i],
    y_str = "Fluorescent Int (a.u.)",
    axis.text.x.size = 12,
    axis.text.y.size = 12,
    title.size = 15,
    text.size = 15
  ) +
    stat_compare_means(method = "anova", label = "p.format", tip.length = 0.03)
  
  #png(
  #filename = paste0("./mCytokine/DotBox plots/NonZscored/Female/5xFAD/", colnames(dataCytZ)[i], ".png"),
  # height = 10, width = 12, units = "in", res = 600
  # )
  pdf(paste0("./mCytokine/DotBox plots/NonZscored/Female/5xFAD/Abeta_ ",colnames(Df_Abeta_Sub)[i],".pdf"),height=5,width=6,useDingbats = FALSE,useKerning = FALSE)
  
  
  print(plotOut)
  dev.off()
}





cor_matrix <- cor(Df_Cyt_Sub,Df_Abeta_Sub ,use = "pairwise.complete.obs", method = "spearman")
nSamples <- nrow(Df_Cyt_Sub)
Pvalue <- corPvalueStudent(cor_matrix, nSamples)

cor_df <- melt(cor_matrix)
colnames(cor_df) <- c("Cytokine", "Amyloid", "Correlation")

pval_df <- melt(Pvalue)
pval_df$padj <- p.adjust(pval_df$value, method = "fdr")
colnames(pval_df) <- c("Cytokine", "Amyloid", "Pvalue","fdr pval")

df <- merge(cor_df, pval_df, by = c("Cytokine", "Amyloid"))

df$Significance <- cut(df$Pvalue,
                       breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                       labels = c("***", "**", "*", ""),
                       right = FALSE)


df$Label <- paste0(df$Significance)

#-------or--------------
df$Label <- paste0("R=", round(df$Correlation, 2), 
                   #"\n", df$Significance,
                   "\nPvalue=", signif(df$Pvalue, 2))

#-------------
df$padj <- p.adjust(pval_df$Pvalue, method = "fdr")

row_dist <- dist(cor_matrix)
row_clust <- hclust(row_dist)
row_order <- row_clust$labels[row_clust$order]

col_dist <- dist(t(cor_matrix))
col_clust <- hclust(col_dist)
col_order <- col_clust$labels[col_clust$order]

df$Cytokine <- factor(df$Cytokine, levels = row_order)
df$Amyloid <- factor(df$Amyloid, levels = col_order)

df_sub=df %>% filter(!(Cytokine == c("IL-6",
                                  "IL-4","IL-2","IL-1B","MCP-1","KC","IL-5","IP-10","IFNy","MIP-1a","IL-9","TNFa","IL-13","IL-1a","IL-17"))) 


# --- PLOT ---
pdf("./mCytokine/SpermanCorrealtion_Abeta_vs_Cytokines_with Pvalues.pdf", width = 20, height = 16, pointsize = 10)
ggplot(df, aes(x = Amyloid, y = Cytokine, fill = Correlation)) +
  geom_tile(color = "white") +
  geom_text(aes(label = Label), size = 2) +
  scale_fill_gradient2(
    low = "blue", mid = "white", high = "red", midpoint = 0,
    limit = c(-1, 1), space = "Lab", name = "Spearman R"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text = element_text(size = 10),
    axis.title = element_blank()
  )

dev.off()



write_csv(df[, c("Cytokine", "Amyloid", "Correlation", "Pvalue", "padj", "Significance")],
          "./mCytokine/SpearmanCorrelation_Abeta vs Cytokines_CorrelationMatrix_Pvals_FDR_SummaryTable.csv")





#----Male 5xFAD---#
sortby <- which(Exp_data$SexGenotype =="5XFAD M")

Df_Cyt_Sub <-  dataCyt[sortby ,] 
Df_Abeta_Sub <-dataAbeta[sortby ,]
Exp_data_sub <- Exp_data[sortby ,]



#Z scoring data

Df_Cyt_SubZ=(apply(t(Df_Cyt_Sub),1,scale)); #z-score the data
Df_Cyt_SubZ[is.na(Df_Cyt_SubZ)]=0
rownames(Df_Cyt_SubZ)=rownames(Df_Cyt_Sub)





FlickerColorBar <- case_when(
  Exp_data_sub$Flicker == "20Hz" ~ "#1A1513",
  Exp_data_sub$Flicker == "40Hz" ~ "#76A33E",
  Exp_data_sub$Flicker == "NoFlicker" ~ "#E298D6"
)


colorBarCombined=cbind(c(FlickerColorBar),c(FlickerColorBar))
colnames(colorBarCombined)=c("Flicker","Flicker")


#Cluster the cytokine data by columns
hc= hclust(dist(t(Df_Cyt_SubZ),method = "euclidean"), method="average")

##Make Heatmap
pdf("Cytokine profile Male 5xFAD Only.pdf", width=12,height=10,pointsize = 6)
#png("Cytokine profile Male 5xFAD Only.png", width=16,height=16,units="in",res=600, pointsize = 20)
heatmap3(Df_Cyt_SubZ,main="Female 5xFAD Cytokine Profile",col=barColors, breaks=breakBarColors,legendfun=function()showLegend(),
         Rowv=NA, Colv=as.dendrogram(hc),  scale="none", na.rm = TRUE,
         labRow = rownames(Df_Cyt_SubZ), labCol = colnames(Df_Cyt_SubZ),
         RowSideColors=colorBarCombined)


dev.off()






#Amyloid beta side bar


ab42<-Df_Abeta_Sub[,2]
ab40<-Df_Abeta_Sub[,1]
ratio<-Df_Abeta_Sub[,3]


breakBarColorsAbeta42=c(-1,seq(0, 1200, 0.1))*10 #Outside numbers clip outliers. This is for zscoring.
barColorsAbeta42Vec = c("black",colorpanel(length(breakBarColorsAbeta42), "white", "gold4"))
indMatchColorone=match(round(ab42*10),round(breakBarColorsAbeta42))
barColorsAb42=barColorsAbeta42Vec[indMatchColorone]

breakBarColorsAbeta40=c(-1,seq(0, 10, 0.01))*10 #Outside numbers clip outliers. This is for zscoring.
barColorsAbeta40Vec = c("black",colorpanel(length(breakBarColorsAbeta40), "white", "orange"))
indMatchColortwo=match(round(ab40*10),round(breakBarColorsAbeta40))
barColorsAb40=barColorsAbeta40Vec[indMatchColortwo]


breakBarColorsAbetaratio=c(-1,seq(0, 180, 0.01))*10 #Outside numbers clip outliers. This is for zscoring.
barColorsAbetaratioVec = c("black",colorpanel(length(breakBarColorsAbetaratio), "white", "gray4"))
indMatchColorthree=match(round(ratio*10),round(breakBarColorsAbetaratio))
barColorsAbratio=barColorsAbetaratioVec[indMatchColorthree]

colorBarCombined=cbind(c(FlickerColorBar),c(barColorsAb42),c(barColorsAb40),c(barColorsAbratio))
colnames(colorBarCombined)=c("Flicker","Ab42","Ab40","ratio")



##Make Heatmap
pdf("Cytokine profile Male 5xFAD Only with Side Bar for Abeta.pdf", width=12,height=10,pointsize = 6)
#png("Cytokine profile Male 5xFAD Only_with side bar.png", width=16,height=16,units="in",res=600, pointsize = 20)
heatmap3(Df_Cyt_SubZ,main="Male 5xFAD Cytokine Profile",col=barColors, breaks=breakBarColors,legendfun=function()showLegend(),
         Rowv=NA, Colv=as.dendrogram(hc),  scale="none", na.rm = TRUE,
         labRow = rownames(Df_Cyt_SubZ), labCol = colnames(Df_Cyt_SubZ),
         RowSideColors=colorBarCombined)


dev.off()




new_meanSEM4 <- function(x, y, colors = x, 
                         x_str = "", y_str = "", title = y_str,
                         axis.text.x.size = 20, axis.text.y.size = 20, 
                         axis.text.x.angle = 90, title.size = 20, text.size = 20) {
  
  dataPlot <- data.frame(x_dot = factor(x), 
                         y_dot = y, 
                         colors_dot = factor(colors))
  
  plotOut <- ggplot(dataPlot, aes(x = x_dot, y = y_dot)) +
    
    geom_boxplot(aes(color = colors_dot), size = 1.5, outlier.shape = NA) +
    
    geom_beeswarm(aes(color = colors_dot, fill = colors_dot), 
                  groupOnX = TRUE, size = 6, shape = 21, stroke = 0.8) +
    
    
    # geom_point(aes(color = colors_dot, fill = colors_dot), 
    #  size = 4.5, shape = 21, stroke = 0.8, alpha = 0.7, position = position_jitter(width = 0.1)) +
    
    # mean + error bar
    #stat_summary(fun = mean, geom = "point", shape = 23, size = 8, fill = "white", color = "black") +
    # stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3, color = "black", size = 1.2) +
    
    scale_color_manual(values = c("NoFlicker" = "#E298D6", "20Hz" = "#1A1513", "40Hz" = "#76A33E")) +
    scale_fill_manual(values = c("NoFlicker" = "#E298D6", "20Hz" = "#1A1513", "40Hz" = "#76A33E")) +
    
    theme_bw(base_size = text.size) +
    theme(
      panel.border = element_rect(color = "black", fill = NA, size = 1.2),
      axis.line = element_line(color = "black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(color = "black", size = axis.text.x.size, angle = axis.text.x.angle),
      axis.text.y = element_text(color = "black", size = axis.text.y.size),
      plot.title = element_text(hjust = 0.5, face = "bold", size = title.size),
      legend.position = "none"  # Turn off legend unless you want to keep it
    ) +
    xlab(x_str) +
    ylab(y_str) +
    ggtitle(title)
  
  return(plotOut)
}

#----Make bar  plots for each cytokine


dataCyt_groupVars=data.frame(Df_Cyt_Sub, Exp_data_sub)





for (i in 1:33) {
  flicker_levels <- c("NoFlicker", "20Hz", "40Hz")
  condition_ordered <- factor(Exp_data_sub$Flicker, levels = flicker_levels)
  
  if (length(unique(condition_ordered)) < 2) {
    message(paste("Skipping", colnames(Df_Cyt_Sub)[i], "- insufficient group levels"))
    next
  }
  
  plotOut <- new_meanSEM4(
    x = condition_ordered,
    y = Df_Cyt_Sub[, i],
    colors = condition_ordered,
    title = colnames(Df_Cyt_Sub)[i],
    y_str = "Fluorescent Int (a.u.)",
    axis.text.x.size = 12,
    axis.text.y.size = 12,
    title.size = 15,
    text.size = 15
  ) +
    stat_compare_means(method = "anova", label = "p.format", tip.length = 0.03)
  
  #png(
  #filename = paste0("./mCytokine/DotBox plots/NonZscored/Male/5xFAD/", colnames(dataCytZ)[i], ".png"),
  # height = 10, width = 12, units = "in", res = 600
  # )
  pdf(paste0("./mCytokine/DotBox plots/NonZscored/Male/5xFAD/",colnames(Df_Cyt_Sub)[i],".pdf"),height=5,width=6,useDingbats = FALSE,useKerning = FALSE)
  
  
  print(plotOut)
  dev.off()
}






for (i in 1:4) {
  flicker_levels <- c("NoFlicker", "20Hz", "40Hz")
  condition_ordered <- factor(Exp_data_sub$Flicker, levels = flicker_levels)
  
  if (length(unique(condition_ordered)) < 2) {
    message(paste("Skipping", colnames(Df_Abeta_Sub)[i], "- insufficient group levels"))
    next
  }
  
  plotOut <- new_meanSEM4(
    x = condition_ordered,
    y = Df_Abeta_Sub[, i],
    colors = condition_ordered,
    title = colnames(Df_Abeta_Sub)[i],
    y_str = "Fluorescent Int (a.u.)",
    axis.text.x.size = 12,
    axis.text.y.size = 12,
    title.size = 15,
    text.size = 15
  ) +
    stat_compare_means(method = "anova", label = "p.format", tip.length = 0.03)
  
  #png(
  #filename = paste0("./mCytokine/DotBox plots/NonZscored/Female/5xFAD/", colnames(dataCytZ)[i], ".png"),
  # height = 10, width = 12, units = "in", res = 600
  # )
  pdf(paste0("./mCytokine/DotBox plots/NonZscored/Male/5xFAD/Abeta_ ",colnames(Df_Abeta_Sub)[i],".pdf"),height=5,width=6,useDingbats = FALSE,useKerning = FALSE)
  
  
  print(plotOut)
  dev.off()
}





cor_matrix <- cor(Df_Cyt_Sub,Df_Abeta_Sub ,use = "pairwise.complete.obs", method = "spearman")
nSamples <- nrow(Df_Cyt_Sub)
Pvalue <- corPvalueStudent(cor_matrix, nSamples)

cor_df <- melt(cor_matrix)
colnames(cor_df) <- c("Cytokine", "Amyloid", "Correlation")

pval_df <- melt(Pvalue)
colnames(pval_df) <- c("Cytokine", "Amyloid", "Pvalue")

df <- merge(cor_df, pval_df, by = c("Cytokine", "Amyloid"))

df$Significance <- cut(df$Pvalue,
                       breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                       labels = c("***", "**", "*", ""),
                       right = FALSE)


df$Label <- paste0(df$Significance)

#-------or--------------
df$Label <- paste0("R=", round(df$Correlation, 2), 
                   #"\n", df$Significance,
                   "\nPvalue=", signif(df$Pvalue, 2))

#-------------
df$padj <- p.adjust(pval_df$Pvalue, method = "fdr")

row_dist <- dist(cor_matrix)
row_clust <- hclust(row_dist)
row_order <- row_clust$labels[row_clust$order]

col_dist <- dist(t(cor_matrix))
col_clust <- hclust(col_dist)
col_order <- col_clust$labels[col_clust$order]

df$Cytokine <- factor(df$Cytokine, levels = row_order)
df$Amyloid <- factor(df$Amyloid, levels = col_order)

# --- PLOT ---
pdf("./mCytokine/SpermanCorrealtion_Abeta_vs_Cytokines_Male_5xFAD.pdf", width = 20, height = 16, pointsize = 10)
ggplot(df, aes(x = Amyloid, y = Cytokine, fill = Correlation)) +
  geom_tile(color = "white") +
  geom_text(aes(label = Label), size = 2) +
  scale_fill_gradient2(
    low = "blue", mid = "white", high = "red", midpoint = 0,
    limit = c(-1, 1), space = "Lab", name = "Spearman R"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text = element_text(size = 10),
    axis.title = element_blank()
  )

dev.off()



write_csv(df[, c("Cytokine", "Amyloid", "Correlation", "Pvalue", "padj", "Significance")],
          "./mCytokine/SpearmanCorrelation_Abeta vs Cytokines_CorrelationMatrix_Pvals_FDR_SummaryTable_Male_5xFAD.csv")











df_PCA <- Df_Cyt_SubZ
PCA_mat <- df_PCA[,] %>% as.matrix()
#PathVar=pathData[, 1]
#df_var=as.numeric(log2(PathVar))

#na_indices <- which(is.na(df_var), arr.ind = TRUE)



oplsOut <- opls(PCA_mat,predI=2)
rotateOut <- rotate.opls(oplsOut,degrees=0)
scores_data=data.frame(rotateOut$T1,rotateOut$T2)
# 


pdfPath <- paste("./ScorePlot_OPLS_5xFAD_Male.pdf", sep = '')
pdf(pdfPath, width = 5, height = 4, pointsize = 14)
scores.plot(rotateOut$T1,rotateOut$T2,factor(Exp_data_sub$Flicker))
print(p)
dev.off()



score_df <- data.frame(
  PC1 = rotateOut$T1,
  PC2 = rotateOut$T2,
  Flicker = factor(Exp_data_sub$Flicker, levels = c("NoFlicker", "20Hz", "40Hz"))
)


flicker_colors <- c("NoFlicker" = "#E298D6", "20Hz" = "#1A1513", "40Hz" = "#76A33E")

p <- ggplot(score_df, aes(x = PC1, y = PC2, color = Flicker)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = flicker_colors) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  theme_minimal(base_size = 14) +
  labs(x = "PC1", y = "PC2", title = "Score Plot") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.title = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.8),
    panel.grid = element_blank()
  )

pdfPath <- "./ScorePlot_OPLS_5xFAD_Male.pdf"
pdf(pdfPath, width = 5, height = 4, pointsize = 14)
png(paste0("./ScorePlot_OPLS_5xFAD_Male.png"),height=4,width=5,units="in",res=600)

print(p)
dev.off()







pdfPath <- paste("./PC1_5xFAD_Cytokines_Male.pdf", sep = '')
pdf(pdfPath, width = 5, height = 4, pointsize = 6)

png(paste0("./PC1_5xFAD_Cytokines_Male.png"),height=10,width=12,units="in",res=600)
loadings.chart(rotateOut$P1,"PC1");

dev.off()

pdfPath <- paste("./PC2_5xFAD_Cytokines_Male.pdf", sep = '')
pdf(pdfPath, width = 5, height = 4, pointsize = 6)
#png(paste0("./PC2_5xFAD_Cytokines_Male.png"),height=10,width=12,units="in",res=600)
loadings.chart(rotateOut$P2,"PC2");
dev.off()





































for(i in 1:33){
  plotOut <- new_meanSEM3(x=Condition,y=dataCytZ[,i],colors = as.factor(Condition),
                          title = colnames(dataCytZ)[i],y_str = "Flourescent Int (a.u.)",
                          axis.text.x.size = 12,axis.text.y.size = 12,title.size = 15,text.size = 15)+
    
    scale_color_manual(values=c("#657D7D","#7F0707","#327777","#BFBFBF","#ED899D","#7EFFFF"))+
    scale_fill_manual(values=c("#657D7D","#7F0707","#327777","#BFBFBF","#ED899D","#7EFFFF"))+
    stat_compare_means(method="anova",p.adjust.method = "bonferroni",label="p.format",tip.length = 0.03)
  pdf(paste0("./mCytokine/DotBox plots All/Zscored/Male ",colnames(dataCyt)[i],".pdf"),height=5,width=6,useDingbats = FALSE,useKerning = FALSE)
  #  png(paste0("./mCytokine/DotBox plots All/Zscored/Male  ", colnames(dataCytZ)[i],".png"),height=10,width=12,units="in",res=600)
  
  print(plotOut)
  dev.off()
}


#All Groups

for(i in 1:33){
  plotOut <- new_meanSEM3(x=Condition,y=dataCytZ[,i],colors = as.factor(Condition),
                          title = colnames(dataCytZ)[i],y_str = "Flourescent Int (a.u.)",
                          axis.text.x.size = 12,axis.text.y.size = 12,title.size = 15,text.size = 15)+
    
    scale_color_manual(values=c("#3D92CF","#B385BC","#3F8A3F","#657D7D","#7F0707","#327777","#66D8FF","#A579FF","#3FFF9F","#BFBFBF","#ED899D","#7EFFFF"))+
    scale_fill_manual(values=c("#3D92CF","#B385BC","#3F8A3F","#657D7D","#7F0707","#327777","#66D8FF","#A579FF","#3FFF9F","#BFBFBF","#ED899D","#7EFFFF"))+
    stat_compare_means(method="anova",p.adjust.method = "bonferroni",label="p.format",tip.length = 0.03)
  pdf(paste0("./mCytokine/DotBox plots All/Zscored/Female ",colnames(dataCyt)[i],".pdf"),height=5,width=6,useDingbats = FALSE,useKerning = FALSE)
  #  png(paste0("./mCytokine/DotBox plots All/Zscored/Female  ", colnames(dataCytZ)[i],".png"),height=10,width=12,units="in",res=600)
  
  print(plotOut)
  dev.off()
}





#---- PCA Analysis----#

pheno<-as.factor(Condition)
df_PCA<-dataCyt

df_PCA[,] <- apply(df_PCA[,],2,"scale")
#df_PCA[is.na(df_PCA)] <- 0
PCA_mat <- df_PCA[,] %>% as.matrix()

oplsOut <- opls(PCA_mat,predI=2)
rotateOut <- rotate.opls(oplsOut,degrees=0)

pdf(paste0("./mCytokine/PCA/PCA Female All Analytes.pdf"),height=8,width=12,useDingbats = FALSE)
#png(paste0("./mCytokine/PCA/PCA Female All Analytes.png"),height=3,width=6,units="in",res=600)
scores.plot(rotateOut$T1,rotateOut$T2,pheno,"Condition",analysis="PCA")+  geom_point(size=10) 

dev.off()


pdf(paste0("./mCytokine/PCA/PCA_PC1_Female All analyte.pdf"),height=16,width=16,useDingbats = FALSE)
#png(paste0("./mCytokine/PCA/PCA_PC1_Female All analyte.png"),height=3,width=6,units="in",res=600)
loadings.chart(rotateOut$P1,"PC1")
dev.off()

pdf(paste0("./mCytokine/PCA/PCA_PC2_Female All analyte.pdf"),height=16,width=16,useDingbats = FALSE)
#png(paste0("./mCytokine/PCA/PCA_PC2_Female All analyte.png"),height=3,width=6,units="in",res=600)
loadings.chart(rotateOut$P2,"PC2")
dev.off()

#pdf(paste0("./mCytokine/PCA/PCA_PC2 Score Dotplots_female All analyte.pdf"),height=16,width=16,useDingbats = FALSE)
png(paste0("./mCytokine/PCA/PCA_PC2 Score Dotplots_female All analyte.png"),height=5,width=8,units="in",res=600)
ggpar(
  errorplot.meanSE.compare(x=pheno,y=rotateOut$T2,colors=pheno,title="PC2 ",
                           p.display = "p.format"))
  + scale_color_manual(values=c("#657D7D","#7F0707","#327777","#BFBFBF","#ED899D","#7EFFFF"))
dev.off()


#pdf(paste0("./mCytokine/PCA/PCA_PC1 Score Dotplots_female All analyte.pdf"),height=16,width=16,useDingbats = FALSE)
png(paste0("./mCytokine/PCA/PCA_PC1 Score Dotplots_female All analyte.png"),height=5,width=8,units="in",res=600)
ggpar(
  errorplot.meanSE.compare(x=pheno,y=rotateOut$T1,colors=pheno,title="PC1 ",
                           p.display = "p.format"))
+ scale_color_manual(values=c("#657D7D","#7F0707","#327777","#BFBFBF","#ED899D","#7EFFFF"))
dev.off()








