## library packages
library(dplyr)
library(data.table)
library(GenomicRanges)
library(ComplexHeatmap)

## meta info
tile_width <- 1e7

high_signal <- read.table("hg38-blacklist.v2.high_signal.bed",header = F, stringsAsFactors = F)
high_signal <- high_signal[,1:3]
colnames(high_signal) <- c("chr", "start", "end")
high_signal.gr <- makeGRangesFromDataFrame(high_signal)

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


cnv_mtx <- cnv_res %>% 
  select(Chromosome, Start, Cell, Ratio) %>%
  tidyr::spread(key = Cell, value=Ratio) 

rownames(cnv_mtx) <- paste0(cnv_mtx$Chromosome, ":", cnv_mtx$Start)
cnv_mtx <- cnv_mtx[,-c(1:2)]

cnv_mtx[cnv_mtx==-1] <- NA
cnv_mtx <- na.omit(cnv_mtx)

## filter high_signal region
cnv_mtx_region <- data.frame(
  chr=sapply(strsplit(rownames(cnv_mtx), ":"), "[[",1),
  start=as.integer(sapply(strsplit(rownames(cnv_mtx), ":"), "[[",2))
)
cnv_mtx_region$end <- cnv_mtx_region$start+tile_width
cnv_mtx.gr <- makeGRangesFromDataFrame(cnv_mtx_region)
cnv_mtx.gr$pos <- rownames(cnv_mtx)
cnv_mtx.gr <- cnv_mtx.gr[!overlapsAny(cnv_mtx.gr, high_signal.gr, minoverlap = 1e+5),]
cnv_mtx <- cnv_mtx[cnv_mtx.gr$pos,]

## prepare annotation for ComplexHeatmap

anno_col <- data.frame(chr=sapply(strsplit(rownames(cnv_mtx), ":"), "[[",1))
rownames(anno_col) <- rownames(cnv_mtx)

col_anno <- HeatmapAnnotation(
  chr = anno_col$chr,
  col = list(chr = chr_cols),
  border = T,
  show_annotation_name = F,
  show_legend = F
)

## colors
chr_cols <- rep(c("black", "white"),12)
names(chr_cols) <- paste0("chr", c(1:22, "X", "Y"))
col_func <- colorRamp2(c(0, 1,2),c("blue", "white", "red"))

ht <- Heatmap(
  t(cnv_mtx), 
  name = "CNV", 
  cluster_columns = F,
  cluster_rows = T,
  # clustering_distance_rows = "eudilicean",
  clustering_method_rows = "complete",
  show_row_names = F,
  row_names_gp = gpar(fontsize=5),
  show_column_names = F,
  show_row_dend = F,
  column_split = factor(anno_col$chr, levels = paste0("chr", c(1:22,"X", "Y"))),
  column_gap = unit(0, "mm"),
  top_annotation = col_anno,
  row_title_rot = 0,
  column_title_rot = 90,
  row_title_gp = gpar(fontsize=10),
  column_title_gp = gpar(fontsize=20),
  row_gap = unit(0, "mm"),
  cluster_row_slices = T,
  # height = unit(10*nrow(plot_mat)/nrow(remove_na_mar)+2,"inch"),
  # width = unit(8,"inch"),
  col = col_func,
  border = TRUE
)
draw(ht)

