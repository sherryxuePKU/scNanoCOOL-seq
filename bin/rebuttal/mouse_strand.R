source("/mnt/e/Project/nanoCool/Report/Manuscript/CellResearch/configure.R")

library(dplyr)
library(ggplot2)
library(stringr)
library(forcats)
library(ComplexHeatmap)
library(Seurat)
library(ggpubr)
library(ggsignif)

setwd(wkdir)

## load raw input
stat_df <- read.table(
  "DNA/QC/basic_stat/nanoCool.mouse_embryo.merged_0710.basic_stat.tsv",
  header = T, stringsAsFactors = F
)
rownames(stat_df) <- stat_df$Sample
mouse_em <- readRDS(seurat_obj_path)

## plot, heatmap, mat_hap
read_file_sb <- function(x){
  df <- read.table(x, header = F, stringsAsFactors = F)
  # df$V4 <- df$V4 - 0.5
  # df$V4[df$V4 > 0.5] <- 1
  # df$V4[df$V4 < 0.5] <- -1
  df <- df %>% filter(V1 %in% paste0("chr", 1:22))
  df <- df %>% group_by(V1) %>% 
    summarise(mean_f=mean(V4))
  x <- rev(strsplit(x, "/")[[1]])[1]
  x <- gsub(file_pattern, "", x)
  # colnames(df) <- c("chr", "start", "end", gsub(file_pattern, "", x))
  colnames(df) <- c("chr", x)
  
  return(df)
}

file_pattern <- ".g1.strand_bias.WCG_10win.bed"
fn <- list.files(pattern = file_pattern, path = "DNA/strand", full.names = T)
sb_df <- lapply(fn, read_file_sb)
sb_df <- sb_df %>% purrr::reduce(function(x,y)base::merge(x,y,by=c("chr")))
sb_df <- sb_df[,colnames(sb_df) %in% c("chr", stat_df$Sample[stat_df$CellType %in% c("Epi", "Late_TE")])]

f1 <- circlize::colorRamp2(
  seq(0, 1, length = 3),
  rev(RColorBrewer::brewer.pal(n = 9, name = "RdBu")[c(1,5,9)])
)

chr_col <- rep(c("black", "white"), round(length(unique(sb_df$chr))/2))
chr_col <- chr_col[1:length(unique(sb_df$chr))]
names(chr_col) <- stringr::str_sort(unique(sb_df$chr), numeric = T)

row_anno <- ComplexHeatmap::rowAnnotation(
  chr= sb_df$chr, col = list(chr=chr_col),
  show_annotation_name=F,
  border = TRUE,show_legend = FALSE
)

row_split <- factor(
  sb_df$chr, levels = stringr::str_sort(
    unique(sb_df$chr), numeric = T
  ))

col_anno <- HeatmapAnnotation(
  celltype= stat_df[colnames(sb_df)[2:ncol(sb_df)],"CellType"], 
  col = list(celltype=celltype_col),
  # gp = gpar(col = "black"), 
  show_annotation_name=F,
  border = TRUE,show_legend = FALSE)

col_split <- factor(
  stat_df[colnames(sb_df)[2:ncol(sb_df)],"CellType"], 
  levels = c("Epi", "Late_TE"))

col_idx <- stringr::str_sort(
  colnames(sb_df)[2:ncol(sb_df)], 
  numeric = TRUE)

col_sd <- apply(sb_df[,2:ncol(sb_df)], 2, sd)
col_sd <- col_sd[order(col_sd)]

ht <- Heatmap(
  as.matrix(sb_df[, col_idx]),
  cluster_rows = F,
  column_order = names(col_sd),
  col =  f1,
  row_gap = grid::unit(0, "mm"),
  border_gp = grid::gpar(col = "black"),
  name = "strand_bias",
  row_title_rot = 0,
  left_annotation = row_anno,
  show_column_names = F,
  row_split = row_split,
  top_annotation = col_anno,
  column_split = col_split
)

plot_time <- strsplit(as.character(Sys.time()), " ")[[1]][1]
pdf(paste0(plot_prefix,"mouse_emrbyo.strand_Dnmt1.heatmap.",plot_time,".pdf"),
    width = 12, height = 5)
draw(ht)
dev.off()

## cor, Dnmt1 & strand_bias_sd, FeatureScatter
sb_sd <- data.frame(strand_bias_sd=col_sd) %>% 
  tibble::rownames_to_column("DNA_ID") %>% 
  inner_join(
    mouse_em@meta.data %>% 
      dplyr::select(DNA_ID, RNA_ID),
    by = c("DNA_ID")) %>% 
  tibble::column_to_rownames("RNA_ID")

mouse_em <- AddMetaData(mouse_em, metadata = sb_sd)

p1 <- FeatureScatter(
  subset(
    mouse_em, 
    DNA_ID %in% names(col_sd) & 
      CellType %in% c("Epi", "Late_TE")),
  feature1 = "Dnmt1", feature2 = "strand_bias_sd",
  cols = celltype_col,
  pt.size = 2
) + stat_cor(
  aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")), 
  size=5, label.y = 0.20, label.x = 0.5) 
ratio <- diff(range(p1$data$Dnmt1))/diff(range(p1$data$strand_bias_sd))

p1 + coord_equal(ratio = ratio)

plot_time <- strsplit(as.character(Sys.time()), " ")[[1]][1]
ggsave( 
  paste0(plot_prefix,"mouse_emrbyo.strand_Dnmt1.",plot_time,".pdf"),
  width = 5, height = 5
)

## cor, Dnmt1 & strand_bias_sd, boxplot, test
p <- VlnPlot(
  subset(
    mouse_em, 
    DNA_ID %in% names(col_sd) & 
      CellType %in% c("Epi", "Late_TE")),
  features = c("Dnmt1","strand_bias_sd"),
  cols = celltype_col,
  pt.size = 0,adjust = .5
)

inner_join(
  p[[1]]$data %>% tibble::rownames_to_column("RNA_ID"), 
  p[[2]]$data %>% tibble::rownames_to_column("RNA_ID"),
  by = c("RNA_ID", "ident")) %>% 
  reshape2::melt(id.var=c("RNA_ID", "ident")) %>% 
  ggplot(aes(x=ident, y=value))+
  # ggbeeswarm::geom_beeswarm(aes(color=ident),cex = 2) +
  geom_boxplot(aes(fill=ident), outlier.shape = NA, width=.65) +
  geom_signif(
    comparisons = list(c("Epi", "Late_TE")),
    size = 0, map_signif_level = TRUE, vjust = 1) + 
  scale_fill_manual(values=celltype_col) + 
  facet_wrap(.~variable, scales = "free") +
  egg::theme_presentation() +
  theme(axis.title.x = element_blank())

plot_time <- strsplit(as.character(Sys.time()), " ")[[1]][1]
ggsave( 
  paste0(plot_prefix,"mouse_emrbyo.strand_Dnmt1_boxplot.",plot_time,".pdf"),
  width = 10, height = 4
)
