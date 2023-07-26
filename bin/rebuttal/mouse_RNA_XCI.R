source("/mnt/e/Project/nanoCool/Report/Manuscript/CellResearch/configure.R")

library(dplyr)
library(forcats)
library(ggplot2)

setwd(wkdir)

## read metadata
stat_df <- read.table(
  "DNA/QC/basic_stat/nanoCool.mouse_embryo.merged_0710.basic_stat.tsv",
  header = T, stringsAsFactors = F
)

## read raw data
read_umi_hap <- function(hap){
  umi_list <- lapply(
    list.files(path = "RNA/xci", full.names = T, pattern = hap),
    function(x){
      df <- read.table(x, header = T,stringsAsFactors = F)
    })
  
  merge2 <- function(x,y)base::merge(x,y,by=c("Gene"))
  umi_mtx <- umi_list %>% purrr::reduce(merge2) %>% tibble::column_to_rownames("Gene")
  return(umi_mtx)
}

umi_g1 <- read_umi_hap("genome1")
umi_g2 <- read_umi_hap("genome2")
all.equal(rownames(umi_g1), rownames(umi_g2))

## X or autosome-linked genes
x_genes <- read.table("RNA/mm10_chrX_genes.txt")
escapee <- read.table("DNA/anno/mouse_escapee.list")

plot_mat_frac <- function(chr, save_plot=F){
  
  switch(
    chr,
    "X"={genes <- rownames(umi_g1) %in% x_genes$V5 & !rownames(umi_g1) %in% escapee$V1},
    "auto"={genes <- !rownames(umi_g1) %in% x_genes$V5}
  )
  
  x_g1 <- apply(umi_g1, 2, function(x){return(sum(x[genes]))})
  x_g2 <- apply(umi_g2, 2, function(x){return(sum(x[genes]))})
  
  hap_df <- inner_join(
    data.frame(x_g1) %>% 
      tibble::rownames_to_column("Sample") %>% 
      mutate(sample=gsub("genome1[.]", "", Sample)) %>% 
      dplyr::select(-Sample),
    data.frame(x_g2) %>% 
      tibble::rownames_to_column("Sample") %>% 
      mutate(sample=gsub("genome2[.]", "", Sample)) %>% 
      dplyr::select(-Sample),
    by = "sample") %>% 
    mutate(x_ratio=x_g1/(x_g1+x_g2)) %>% 
    inner_join(
      stat_df %>% filter(QC=="Pass") %>% 
        mutate(sex2=ifelse(CellType!="M2", sex, Stage)) %>% 
        dplyr::select(RNA_ID, CellType, Stage, sex2),
      by=c("sample"="RNA_ID")) %>% 
    mutate(celltype=ifelse(grepl("^E", Stage), as.character(CellType), Stage))
  
  
  p <- hap_df %>% 
    filter(x_g1+x_g2>20) %>% 
    ggplot(aes(x=fct_relevel(celltype, c(stage_order)), y=x_ratio)) +
    geom_violin(scale = "width", linetype="longdash") + 
    ggbeeswarm::geom_quasirandom(aes(color=celltype),size=2) +
    facet_grid(.~sex2, scales = "free_x", space = "free") + 
    scale_color_manual(values=celltype_col) +
    labs(y=paste0("Maternal ",chr,"-fraction")) + 
    scale_y_continuous(
      limits = c(0,1),
      breaks = c(0, 0.5, 1),
      labels = c(0, 0.5, 1)) + 
    theme_bw(base_size = 15) +
    # egg::theme_presentation() +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1),
          strip.background = element_blank(), 
          panel.grid = element_blank())
  
  if(save_plot){
    plot_time <- strsplit(as.character(Sys.time()), " ")[[1]][1]
    ggsave(
      plot = p, width =10, height = 3, 
      filename = paste0(plot_prefix,"mouse_embryo.mat_",chr,"_frac.",plot_time,".pdf")
    )
  }else{
    print(p)
  }
}

plot_mat_frac(chr = "auto", save_plot=T)
plot_mat_frac(chr = "X", save_plot=T)

