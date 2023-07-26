source("/mnt/e/Project/nanoCool/Report/Manuscript/CellResearch/rebuttal_R/cellline_config.R")

library(dplyr)
library(ggplot2)
library(data.table)
library(RColorBrewer)

setwd(wkdir)

percent <- 0.3
min_ncell <- 5

#1) load raw data
df <- fread("Correlation/TGS_NGS_K562_HFF1.merged.1000bp_1X_CpG.bed", header = T)
df <- as.data.frame(df)
pos_df <- df[,1:2]
df <- df[,grep("K562|P14", colnames(df))]

#2) original filtering and correlation
use_Pos <- apply(df[,3:ncol(df)],1,function(x){sum(is.na(x))<percent*c(ncol(df)-2)})
cor_mt <- cor(df[use_Pos, 3:ncol(df)], use="pairwise.complete.obs", method="pearson")

#3) new filtering and correlation
ngs_idx <- grep("K562", colnames(df))
tgs_ids <- grep("P14", colnames(df))

use_Pos <- apply(
  df[,3:ncol(df)], 1,
  function(x){
    sum(is.na(x[ngs_idx])) < c(ncol(df[,ngs_idx]) - min_ncell) &
      sum(is.na(x[tgs_ids])) < c(ncol(df[,tgs_ids]) - min_ncell)
  }
)

mean_df <- data.frame(
  chr=pos_df[use_Pos, 1], pos=pos_df[use_Pos, 2],
  ngs_meth=rowMeans(df[use_Pos,grep("K562", colnames(df))], na.rm = T),
  tgs_meth=rowMeans(df[use_Pos,grep("P14", colnames(df))], na.rm = T)
)

#4) plot cor
plot_cor_mtx <- function(exclude_black=F,save_plot=F){
  # Set color palette for 2D heatmap
  rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
  r <- rf(32)
  
  # pal <- wesanderson::wes_palette("Zissou1", 100, type = "continuous")
  
  plot_df <- mean_df %>% filter(!is.na(ngs_meth) & !is.na(tgs_meth)) 
  # %>% filter(!(ngs_meth < 0.25 & tgs_meth > 0.5)) 
  
  plot_gr <- plot_df %>% mutate(end=pos+1e3-1)
  plot_gr <- makeGRangesFromDataFrame(
    plot_gr, start.field = "pos", keep.extra.columns = T
  )
  
  if(exclude_black){
    # exclude_gr <- genomation::readBed("hg38_rmsk_Simple_repeat.bed.gz")
    over_black <- overlapsAny(plot_gr, Signac::blacklist_hg38_unified)
    # over_black <- overlapsAny(plot_gr, exclude_gr)
    plot_gr <- plot_gr[!over_black,]
    print("Regions overlapping blacklist were excluded!")
  }
  
  plot_df <- as.data.frame(plot_gr)
  
  c <- cor(plot_df$ngs_meth, plot_df$tgs_meth, "pairwise.complete.obs")
  m <- lm(tgs_meth ~ ngs_meth, plot_df)
  
  title <- sprintf(
    "N = %d r = %.2f \ny = %.2f + %.2f * x (R2 = %.2f)", 
    nrow(plot_df), c, coef(m)[1], coef(m)[2], summary(m)$r.squared
  )
  
  p <- plot_df %>% 
    ggplot(aes(x=ngs_meth, y=tgs_meth)) +
    # geom_density_2d()+
    geom_bin2d(bins=50) +
    geom_smooth(method = "lm", se=FALSE, color="black", formula = y~x, linetype="longdash")+ 
    # scale_fill_gradientn(colours = pal) + 
    scale_fill_gradientn(colors=r, trans="log10") +
    ylab("scNanoCOOL-seq methylation level") +
    xlab("scCOOL-seq methylation level") +
    coord_equal(ratio = 1)+
    # geom_abline(slope = 1, intercept = 0.25, linetype="longdash")+
    theme_bw(base_size=15) +
    ggtitle(title)
  
  if(save_plot){
    plot_time <- strsplit(as.character(Sys.time()), " ")[[1]][1]
    ggsave(
      plot = p, width = 10, height = 8, 
      filename = paste0(plot_prefix, "/cool_K562.sc_meth_cor",plot_time,".pdf")
    )} else{print(p)}
}

plot_cor_mtx(exclude_black = T, save_plot = T)

