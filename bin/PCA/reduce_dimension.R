## library packages
library(pcaMethods)
library(data.table)
library(ggplot2)

## input
percent <- 0.7
cols <- c("DMSO"="#2872c1", "1uM"="#eac015")

setwd("/path/to/data")
CG_df <- fread("nanoCool_5aza.GCH_genebody_NA_Per0.7.mtx",header=T)

## PCA
pos_df <- CG_df[,1:3]
CG_df <- CG_df[,4:ncol(CG_df)] %>% as.data.frame()
cellsAll <- colnames(CG_df)
use_Pos <- apply(CG_df,1,function(x){sum(is.na(x))<percent*ncol(CG_df)})
use_Cell <- apply(CG_df,2,function(x){sum(is.na(x))<nrow(CG_df)})

use_GC_Merge <- CG_df[use_Pos,use_Cell]
use_GC_pos <- pos_df[use_Pos,]
use_GC_Merge <- as.matrix(use_GC_Merge)

pc <- pca(use_GC_Merge,nPcs=5,method="ppca")
imputed <- completeObs(pc)
pca<-prcomp(t(imputed))

pca_df <- as.data.frame(pca$x[,1:5])
colnames(pca_df) <- paste("PC",1:5,sep="")

## optional analysis
# summary(pca)
# screeplot(pca, npcs = 10, type = "lines")

## plot the result
pca_df$Sample <- rownames(pca_df)
stat_df <- read.table("nanoCool_5aza.stat_info.txt", header = T, stringsAsFactors = F)
pca_df <- merge(pca_df, stat_df[stat_df$QC=="Pass",c("Sample", "Time", "Dose", "WCG_meth")], by="Sample")
pca_df$Group <- paste(pca_df$Time, pca_df$Dose, sep = "_")


pca_df %>% filter(Dose!="untreated") %>%
  ggplot(aes(x=PC1, y=PC2, color=Dose))+
  stat_ellipse(aes(fill=Dose),
               type ="norm", geom ="polygon",
               alpha=0.2,color=NA)+
  geom_point()+
  facet_wrap(.~Time, nrow = 1)+
  scale_color_manual(values = cols)+
  scale_fill_manual(values = cols)+
  theme_classic(base_size = 15)+
  theme(legend.position = "right")
  