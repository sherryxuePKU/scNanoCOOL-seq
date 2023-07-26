################## presetting ##################
source(config_dir)

library(dplyr)
library(ggplot2)
library(stringr)
library(DescTools)

setwd(wkdir)

################## functions ##################
#1) read bam_ratio.txt file
read_ratio <- function(x) {
  chr_order <- c(paste0("chr", c(1:22, "X", "Y")))
  
  fname <- ratio_files[x]
  print(fname)
  df <- read.table(fname, header = T, stringsAsFactors = F)
  df$Chromosome <- paste0("chr", df$Chromosome)
  df <- df %>%
    arrange(factor(df$Chromosome, levels = chr_order),
            Start)
  df$Pos <- 1:nrow(df)
  tmp <- rev(str_split(ratio_files[x], "/")[[1]])[1]
  df$Cell <- gsub(file_pattern, "", tmp)
  return(df)
}
#2) make annotation for the x-axis of CNV dotplot
make_x_anno <- function(cnv_res){
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
  return(chr_window_num)
}
#3) group cells for plot
group_sc <- function(select_cell, step=4){
  
  sc_order_df <- data.frame(
    start=seq(1,nrow(select_cell),step), 
    end=seq(1,nrow(select_cell),step)+(step-1)
  )
  sc <- apply(
    sc_order_df, 1, function(x){
      lst <-c(select_cell$Cell[x[1]:x[2]])
      lst <- lst[!is.na(lst)]
      return(lst)
    }, simplify = FALSE)
}
#4) plot single CNV dotplot
CNV_dotplot <- function(cells){
  
  chr_window_num <- make_x_anno(cnv_res)
  
  p <- cnv_res %>% 
    # filter(cnv_res$Cell %in% select_cell$Cell[i:(i+3)]) %>%
    filter(cnv_res$Cell %in% cells) %>% 
    left_join(cv_cell %>% dplyr::select(Cell, cell_cv), by="Cell") %>% 
    ggplot(aes(x=Pos))+
    geom_hline(
      yintercept=c(c(0.5,1,2,2.5)*2),
      alpha=0.5,colour="grey",linetype="longdash")+
    geom_hline(yintercept=yint,alpha=0.4,colour="lightgrey",size=0.9)+
    geom_point(aes(y=Ratio), color="grey",size=.5)+
    geom_point(aes(y=CopyNumber, color=CNV), shape=15, size=1)+
    geom_vline(
      xintercept = chr_window_num[1:23,2],
      colour="black", linetype = "longdash")+
    scale_color_manual(values=cnv_color)+
    scale_x_discrete(
      breaks = chr_window_num$position,
      limits = chr_window_num$position,
      labels = gsub("chr", "", chr_window_num$Chr),
      expand = c(0,0), guide = guide_axis(n.dodge=2))+
    facet_grid(cell_cv~.) + cnv_theme + ylim(0,6)
  
  return(p)
}

################## metadata ##################
stat_df <- read.table(stat_dir,header = T, stringsAsFactors = F)
cellsPass <- stat_df$Sample[stat_df$QC=="Pass"]

file_pattern <- ".clean_bismark_mm2.sort.rmdup.bam_ratio.txt"

################## K562 (single cell) ################## 
## FREEC: ploidy = 2 & breakPointThreshold = .6 & control (HFF1 pseudobulk)

#1) load raw input
ratio_files <- list.files(
  path = "CNV/1M/Bed", full.names = T, pattern = ".bam_ratio.txt"
)
ratio_files <- ratio_files[grep("P13", ratio_files)]
ratio_files <- ratio_files[
  ratio_files %like any% paste0("%",cellsPass, file_pattern,"%")]

cnv_res <- lapply(1:length(ratio_files), read_ratio)
cnv_res <- do.call("rbind", cnv_res)

#2) copy number
cnv_res <- cnv_res %>% 
  mutate(
    CNV = case_when(
      CopyNumber > 3 ~ "Gain", CopyNumber < 3 ~ "Loss",
      CopyNumber == 3 | is.na(CopyNumber) ~ "Normal"
    )
  )

#3) select cells with top low CV
cv_cell <- cnv_res %>% group_by(Cell) %>% 
  filter(CopyNumber > 0) %>%
  # filter(Ratio > 0) %>%
  mutate(norm_ratio=Ratio/CopyNumber) %>% 
  summarise(CV=sd(norm_ratio, na.rm = T)/mean(norm_ratio, na.rm = T)) %>% 
  mutate(cell_cv = paste0(Cell, "\nCV=", round(CV, 3))) %>%
  arrange(CV)

select_cell <- cv_cell %>%
  filter(Cell %in% paste0("cool_P13_", c(15,58,85)))

#5) dotplot
cnv_res <- cnv_res %>% mutate(Ratio=3*Ratio)

sc <- group_sc(select_cell, step = 4)

n <- 1; yint = 3

lapply(sc[1], function(i){
  p <- CNV_dotplot(i)
  
  plot_time <- strsplit(as.character(Sys.time()), " ")[[1]][1]
  
  fname <- paste(
    plot_prefix,"/cellline.K562_P13.1M", 
    n, plot_time,"pdf", sep = ".")
  ggsave(
    plot = p, filename = fname, 
    width = 16, height = 4, units = "in"
  )
  n <- n + 1
})

################## HFF1 (single cell) ################## 
## FREEC: ploidy = 2 & breakPointThreshold = .6 & control (HFF1 pseudobulk)

#1) load raw input

ratio_files <- list.files(
  path = "CNV/1M/Bed", full.names = T, 
  pattern = ".CHlt80per.bam_ratio.txt"
)
ratio_files <- ratio_files[grep("P11", ratio_files)]
ratio_files <- ratio_files[
  ratio_files %like any% paste0("%",cellsPass, file_pattern,"%")]

cnv_res <- lapply(1:length(ratio_files), read_ratio)
cnv_res <- do.call("rbind", cnv_res)

#2) smoothing inferred copy number
cnv_res <- cnv_res %>% 
  mutate(
    CopyNumber = round(
      forecast::ma(CopyNumber, order = 50, centre = F)),
    CNV = case_when(
      CopyNumber > 2 ~ "Gain", 
      CopyNumber < 2 ~ "Loss",
      CopyNumber == 2 | is.na(CopyNumber) ~ "Normal"
    )
)

#3) select cells with no CNV
select_cell <- cnv_res %>%
  filter(Chromosome %in% paste0("chr", 1:22)) %>%
  group_by(Cell, CNV) %>% summarise(count=n()) %>% 
  group_by(Cell) %>% filter(count[CNV=="Normal"]==sum(count))

#4) calc the CV (sd/mean)
cv_cell <- cnv_res %>% group_by(Cell) %>% 
  filter(CopyNumber > 0) %>%
  # filter(Ratio > 0) %>%
  mutate(norm_ratio=Ratio/CopyNumber) %>% 
  summarise(CV=sd(norm_ratio, na.rm = T)/mean(norm_ratio, na.rm = T)) %>% 
  mutate(cell_cv = paste0(Cell, "\nCV=", round(CV, 3))) %>%
  arrange(CV)

#5) dotplot

cnv_res <- cnv_res %>% mutate(Ratio=2*Ratio)

sc <- group_sc(select_cell, step = 4)

n <- 1; yint = 2

lapply(sc[1], function(i){
  p <- CNV_dotplot(i)
  
  plot_time <- strsplit(as.character(Sys.time()), " ")[[1]][1]
  
  fname <- paste(
    plot_prefix,"/cellline.HFF1_P11.1M", 
    n, plot_time,"pdf", sep = ".")
  ggsave(
    plot = p, filename = fname, 
    width = 16, height = 4, units = "in"
  )
  n <- n + 1
})


################## K562 (pseudobulk) ##################
## FREEC: ploidy = 3 & breakPointThreshold = .8  & control (HFF1 pseudobulk | WGBS_GM12878)

#1) load raw input

ratio_files <- list.files(
  path = "CNV/1M/Bed", full.names = T, 
  pattern = ".CHlt80per.bam_ratio.txt"
)
ratio_files <- ratio_files[grep("5aza|ENCFF963XLT", ratio_files)]

cnv_res <- lapply(1:length(ratio_files), read_ratio)
cnv_res <- do.call("rbind", cnv_res)

#2) copy number
cnv_res <- cnv_res %>% 
  mutate(
    CNV = case_when(
      CopyNumber > 3 ~ "Gain", CopyNumber < 3 ~ "Loss",
      CopyNumber == 3 | is.na(CopyNumber) ~ "Normal"
    )
  )

#3) calc CV
cv_cell <- cnv_res %>% group_by(Cell) %>% 
  filter(CopyNumber > 0 & Ratio > 0) %>%
  mutate(norm_ratio=Ratio/CopyNumber) %>% 
  summarise(CV=sd(norm_ratio, na.rm = T)/mean(norm_ratio, na.rm = T)) %>% 
  mutate(cell_cv = paste0(Cell, "\nCV=", round(CV, 3))) %>% 
  arrange(CV)

#4) plot and output
chr_window_num <- make_x_anno(cnv_res)
cnv_res <- cnv_res %>% mutate(Ratio=3*Ratio)
yint <- 3
p <- CNV_dotplot(c("ENCFF963XLT", "nanoCool_5aza.0"))
print(p)

plot_time <- strsplit(as.character(Sys.time()), " ")[[1]][1]
ggsave(
  filename = paste0(plot_prefix,"/K562.pseudobulk_ENCODE.CNV.",plot_time,".pdf"),
  width = 16, height = 4, units = "in"
)
