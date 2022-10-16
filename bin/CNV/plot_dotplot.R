## library packages
library(dplyr)
library(ggplot2)
library(data.table)

## meta info
stat_df <- read.table("/path/to/stat_info.txt",header = T, stringsAsFactors = F)
cellsPass <- stat_df %>%
  filter(WCG_genebody_profile=="Pass" &GCH_TSS_profile=="Pass" & Conversion_ratio >=0.98

## input
indir <- "/path/to/input"
filePattern <- ".clean_bismark_mm2.sort.rmdup.bam_ratio.txt"
fileNames <- list.files(indir,pattern = filePattern)

samples <- gsub(filePattern,"",ratio_files)
ratio_files <- fileNames[samples %in% cellsPass$Sample]

todo_index <- 1:length(ratio_files)
plot_dotplot <- function(x){
  fileNames <- paste0(indir, ratio_files[x])
  df <- fread(fileNames, header = T, stringsAsFactors = F)
  df$Chromosome <- paste0("chr", df$Chromosome)
  df <- df %>% arrange(
    factor(df$Chromosome, levels = c(paste0("chr", c(1:22, "X", "Y")))),Start
  )
  df$Pos <- 1:nrow(df)
  df$Cell <- gsub(filePattern,"", ratio_files[i])
}

cnv_res <- lapply(todo_index, plot_dotplot)
cnv_res <- do.call("rbind", cnv_res)

## prepare coordinate of chromosomes
tmp_df <- cnv_res[cnv_res$Cell==unique(cnv_res$Cell)[1],]
tmp_df$Chromosome <-factor(tmp_df$Chromosome, levels = paste0("chr", c(1:22, "X", "Y")))

chr_window_num <- as.data.frame(table(tmp_df$Chromosome))
colnames(chr_window_num)<-c("Chr","Number")
chr_window_num$position <-rep(0,24)
for(i in 2:24){
  chr_window_num[i,2] <-chr_window_num[i,2]+chr_window_num[i-1,2]
}
chr_window_num[1,3] <- floor(chr_window_num[1,2]/2)
for(i in 2:24){
  chr_window_num[i,3]<-floor(chr_window_num[i-1,2]+0.5*(chr_window_num[i,2]-chr_window_num[i-1,2]))
}
chr_window_num  <- chr_window_num %>% arrange(Chr)

## plot, e.g. K562
cnv_res$CNV <- "Normal"
cnv_res$CNV[cnv_res$CopyNumber > 3] <- "Gain"
cnv_res$CNV[cnv_res$CopyNumber < 3] <- "Loss"

select_cell <- cnv_res %>% 
  group_by(Cell) %>%
  filter(CopyNumber >0) %>%
  mutate(norm_ratio=Ratio/CopyNumber) %>%
  summarise(CV=sd(norm_ratio, na.rm = T)/mean(norm_ratio, na.rm = T)) %>%
  arrange(CV) %>%
  filter(Cell %in% paste0("cool_P13_", c(15,58,85)))

for (i in seq(1,60,4)) {
  p <- cnv_res %>% filter(cnv_res$Cell %in% select_cell$Cell[i:(i+3)]) %>%
    ggplot(aes(x=Pos))+
    geom_hline(yintercept=c(c(0.5,1,2,2.5)*2),alpha=0.5,colour="grey",linetype="longdash")+
    # geom_hline(yintercept=c(c(0.5,1.5)*2),alpha=0.5,colour="grey",linetype="longdash")+
    geom_hline(yintercept=3,alpha=0.4,colour="lightgrey",size=0.9)+
    # geom_point(aes(y=Ratio*2), color="darkgrey",size=1.2)+
    # geom_point(aes(y=round(CopyNumber/3*2), color=CNV),shape=15)+
    geom_point(aes(y=Ratio*3), color="grey",size=.5)+
    geom_point(aes(y=CopyNumber, color=CNV), shape=15, size=1)+
    geom_vline(xintercept = chr_window_num[1:23,2],colour="black", linetype = "longdash")+
    # scale_color_manual(values=c("1"="#E41A1C","2"="#377EB8"))+
    scale_color_manual(values=c("Gain"="#E41A1C","Loss"="#377EB8", "Normal"="darkgrey"))+
    scale_x_discrete(breaks = chr_window_num$position,labels=chr_window_num$Chr,
                     expand = c(0,0),guide = guide_axis(n.dodge=2))+
    facet_grid(Cell~.)+
    theme_bw()+theme(
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      panel.background = element_blank(), 
      axis.line = element_line(colour = "black"),
      strip.background = element_blank(),
      strip.text.y = element_text(angle = 0, size = 10),
      panel.border = element_rect(size = 1),
      axis.text = element_text(size = 10),
      legend.position = "none",
      axis.title = element_blank(),
      axis.ticks.x = element_blank()
    )+
    ylim(0,6)
  print(p)
  ggsave(paste0("Plot/cellline.K562_P13.5M.", i, "-", i+4, "dotplot.pdf"), width = 16, height = 4)
}