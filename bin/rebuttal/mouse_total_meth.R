source("/mnt/e/Project/nanoCool/Report/Manuscript/CellResearch/configure.R")

library(data.table)
library(stringr)
library(ggplot2)

## profile, genebody
plot_pseudobulk_profile <- function(C_motif,save_plot=F){
  # C_motif <- "WCG"
  switch(C_motif,
         "WCG" = {y_lab = "DNA Methylation Level (%)"},
         "GCH" = {y_lab = "Chromatin accessibility (%)"}
  )
  
  # base_dir <- "/mnt/e/Project/nanoCool/Data/mouse_embryo/DNA/"
  indir <- paste0(wkdir, "/DNA/pseudobulk/",C_motif,"_Genebody")
  fileNames <- list.files(indir, full.names = T)
  
  for(i in 1:length(fileNames)){
    methy <- fread(fileNames[i], header = T)
    methy <- data.frame(methy)
    methy$data <- methy$data * 100
    
    if(i == 1)
      com.ave.methy <- methy else
        com.ave.methy <- rbind(com.ave.methy, methy)
  }
  
  com.ave.methy$Sample <- str_split(com.ave.methy$Sample, "[.]", simplify = T)[,1]
  colnames(com.ave.methy) <- c("methylation", "bin", "sample")
  # com.ave.methy$sample <- gsub("mouse_", "", com.ave.methy$sample)
  com.ave.methy$sample <- plyr::mapvalues(
    com.ave.methy$sample, 
    from = paste0("mouse_", c('Epi', "eTE", "ICM", "lTE", "PrE")), 
    to =  c('Epi', "Early_TE", "ICM", "Late_TE", "PrE")
  )
  
  p <- ggplot(
    data = com.ave.methy, 
    aes(x = bin, y = methylation, group=sample, color=sample)) + 
    geom_line(position = position_dodge(0.2), linewidth = 1.5) +
    scale_color_manual(values = celltype_col) + ylab(y_lab) + 
    # ylab("DNA Methylation Level (%)") +
    # ylab("Chromatin accessibility (%)") +
    geom_vline(aes(xintercept = 50), linetype = "dashed", linewidth = 1) +
    geom_vline(aes(xintercept = 250), linetype = "dashed", linewidth = 1)+
    scale_x_continuous(
      breaks = c(0,50,250,300),labels = c("-15kb", "TSS", "TES", "+15kb")) + 
    # facet_wrap(. ~ sample, scales = "free_y", ncol = 6) +
    # ggh4x::force_panelsizes(rows = unit(1, "in"), cols = unit(2, "in"), TRUE) +
    theme_classic(base_size = 15) + theme(axis.title.x = element_blank())
  # egg::theme_presentation() + theme(axis.title.x = element_blank())
  
  ratio <- diff(range(p$data$bin))/(diff(range(p$data$methylation))*1.25)
  p <- p + coord_equal(ratio = ratio)
  # print(p)
  if(save_plot){
    ggsave(
      plot = p,width = 6, height = 4,
      filename = paste0(
        plot_prefix,"mouse_embryo.pseudobulk_",
        C_motif,"_genebody.",
        strsplit(as.character(Sys.time()), " ")[[1]][1],".pdf"))    
  }else {
    print(p)
  }

}

plot_pseudobulk_profile("WCG")
plot_pseudobulk_profile("GCH")

## plot, QCed, celltype, meth level
stat_df <- read.table(
  "DNA/QC/basic_stat/nanoCool.mouse_embryo.merged_0708.basic_stat.tsv",
  header = T, stringsAsFactors = F
)

stat_df %>% 
  filter(QC=="Pass" & CellType!="M2") %>%
  dplyr::select(CellType, WCG_meth, GCH_meth, Sample) %>% 
  reshape2::melt(id.var = c("Sample", "CellType")) %>%
  ggplot(aes(x=fct_relevel(CellType, stage_order), y=value*100, color=CellType)) + 
  ggbeeswarm::geom_quasirandom() +
  scale_color_manual(values = celltype_col) + 
  facet_wrap(.~variable, scales = "free_y") + 
  labs(y= "Methylation level(%)") + 
  egg::theme_presentation() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        legend.position = "none")

plot_time <- strsplit(as.character(Sys.time()), " ")[[1]][1]
ggsave(
  paste0(plot_prefix,"mouse_embryo.meth_celltype.",plot_time,".pdf"),
  width = 9, height = 5
)

## stats, QCed, celltype, meth
stat_df %>% 
  # filter(QC=="Pass") %>%
  group_by(CellType) %>%
  # group_by(batch) %>% 
  summarise(
    pass=length(which(QC=="Pass")), 
    qc_ratio=pass/n(),
    avg_meth = mean(WCG_meth[QC=="Pass"]),
    avg_acc = mean(GCH_meth[QC=="Pass"])) 
