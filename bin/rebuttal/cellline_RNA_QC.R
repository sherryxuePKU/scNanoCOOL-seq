source("/mnt/e/Project/nanoCool/Report/Manuscript/CellResearch/rebuttal_R/cellline_config.R")

library(dplyr)
library(ggplot2)
library(stringr)
library(Seurat)

setwd(wkdir)

#########################################
## preprocess RNA data for cell lines  ##
#########################################
stat_df <- read.table(stat_dir,header = T, stringsAsFactors = F)

res_rna <- read.table(
  "RNA/nanoCool_cellline.umi_counts.merge.tsv",
  header = T, row.names = 1, stringsAsFactors = F
)
rna_info <- read.table(
  "RNA/nanoCool_cellline.RNA.stat_info.txt", 
  header = T, stringsAsFactors = F
)
rownames(rna_info) <- rna_info$nanoCool_cellline.RNA_id

cellline <- CreateSeuratObject(
  counts = res_rna,
  min.cells = 3, 
  min.features = 1000,
  meta.data = rna_info,
  project = "nanoCool_cellline"
)
cellline <- subset(cellline, nanoCool_cellline.DNA_id %in% cellsPass)
cellline[["percent.mt"]]  <- PercentageFeatureSet(cellline, pattern = "^MT-")
cellline$Group <- "TGS_HFF1"
cellline$Group[grep("P13|P14", cellline$nanoCool_cellline.DNA_id)] <- "TGS_K562"
cellline <- NormalizeData(cellline, normalization.method = "LogNormalize", scale.factor = 10000)
cellline <- FindVariableFeatures(cellline, selection.method = "vst", nfeatures = 2000)
cellline <- ScaleData(cellline, features = rownames(cellline), verbose = F)

saveRDS(cellline, file = seurat_obj)

###############################
## stat quality of RNA data  ##
###############################

cellline <- readRDS(seurat_obj)
p <- cellline@meta.data %>% 
  dplyr::select("Group", "nFeature_RNA", "nCount_RNA") %>% 
  reshape2::melt(id.var="Group") %>% 
  ggplot(aes(x=Group, y=value, fill=Group))+
  geom_boxplot(outlier.shape = 21) + 
  scale_fill_manual(values = group_pal) + 
  facet_wrap(.~variable, scale="free")+
  theme_classic(base_size = 15)+
  theme(strip.background = element_blank())

plot_time <- strsplit(as.character(Sys.time()), " ")[[1]][1]
ggsave(
  plot = p, width = 8, height = 4, 
  filename = paste0(plot_prefix, "/cool_P7-P10.RNA_QC.",plot_time,".pdf")
)

table(cellline$Group)

