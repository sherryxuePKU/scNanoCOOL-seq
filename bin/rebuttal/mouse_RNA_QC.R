source("/mnt/e/Project/nanoCool/Report/Manuscript/CellResearch/configure.R")

library(Seurat)
library(stringr)
library(scCustomize)
library(ggplot2)
library(dplyr)

setwd(wkdir)

## merge count matrix
umi_list <- lapply(
  list.files(path = "RNA/umi_count/", full.names = T),
  function(x){
    df <- read.table(x, header = T,stringsAsFactors = F)
  })

umi_count <- umi_list %>% 
  purrr::reduce(function(x,y)base::merge(x,y,by=c("Gene"))) %>% 
  tibble::column_to_rownames("Gene")

## integrate with metadata
meta_info <- read.table("RNA/Mouse_embryo.RNA_meta.txt", header = T, stringsAsFactors = F)
rownames(meta_info) <- meta_info$RNA_ID
meta_info <- meta_info %>% 
  mutate(PLAT=str_split(RNA_ID, "[.]", simplify = T)[,3])
meta_info <- meta_info %>% filter(!grepl("cool_P35_", DNA_ID))

## sup stat
# meta_info <- meta_info %>% mutate(batch = str_split(DNA_ID, "_", simplify = T)[,2])
# 
# usable_batch <- paste0("P", c(23:25,28, 32,33, 34, 36, 37))
# mouse_em@meta.data %>% 
#   filter(batch %in% usable_batch) %>% 
#   summarise(n())
# 
# mouse_em@meta.data %>%
#   mutate(batch = str_split(DNA_ID, "_", simplify = T)[,2]) %>%
#   filter(batch %in% usable_batch) %>%
#   summarise(n(), ng=mean(nFeature_RNA))



common_cells <- intersect(meta_info$RNA_ID, colnames(umi_count))

## create seurat object
mouse_em <- CreateSeuratObject(
  counts = umi_count[,common_cells], 
  meta.data = meta_info[meta_info$RNA_ID %in% common_cells,], 
  min.cells = 3, min.features = 1000
)
mouse_em <- NormalizeData(mouse_em, normalization.method = "LogNormalize", scale.factor = 10000)
mouse_em <- FindVariableFeatures(mouse_em, selection.method = "vst", nfeatures = 2000)
mouse_em <- ScaleData(mouse_em, features = rownames(mouse_em), verbose = F)
mouse_em <- RunPCA(mouse_em, features = VariableFeatures(object = mouse_em))
mouse_em <- RunUMAP(mouse_em, dims = 1:50)
mouse_em <- RunTSNE(mouse_em, dims = 1:50)
mouse_em <- FindNeighbors(mouse_em, dims = 1:50, features = VariableFeatures(mouse_em))
mouse_em <- FindClusters(mouse_em, resolution = 0.8)

## preview
VlnPlot(mouse_em, features = c("nFeature_RNA", "nCount_RNA"), group.by = "CellType")
DimPlot(mouse_em, group.by = c("Stage", "PLAT", "seurat_clusters"), reduction = "tsne")
FeaturePlot(mouse_em, features = markers, reduction = "tsne")

p <- DimPlot(mouse_em, group.by = c("CellType"), reduction = "tsne")
cells.located <- CellSelector(plot = p)
cells.located <- mouse_em@meta.data %>% 
  filter(RNA_ID %in% cells.located & !CellType %in% c("ICM", "Early_TE"))
mouse_em <- subset(mouse_em, cells = setdiff(colnames(mouse_em), cells.located$RNA_ID))

## celltype assignment
Idents(mouse_em) <- mouse_em$seurat_clusters
new.cluster.ids <- c("ICM", "Late_TE", "Early_TE", "Epi", "PrE", "Epi", "M2")
names(new.cluster.ids) <- levels(mouse_em)
mouse_em <- RenameIdents(mouse_em, new.cluster.ids)

mouse_em$CellType <- Idents(mouse_em)

umap_df <- mouse_em@reductions$umap@cell.embeddings %>% 
  as.data.frame %>% 
  tibble::rownames_to_column("RNA_ID") %>% 
  inner_join(mouse_em@meta.data, by="RNA_ID")
outlier <- umap_df %>% 
  filter(UMAP_1  < -5 & CellType=="Epi")
mouse_em@meta.data[outlier$RNA_ID,"CellType"] <- "Late_TE"

## vlnplot, QC
VlnPlot(
  subset(mouse_em, CellType!="M2"), 
  features = c("nFeature_RNA", "nCount_RNA"), 
  group.by = "CellType",
  pt.size = 0,
  cols = celltype_col
)
ggsave(
  paste0(plot_prefix, "mouse_embryo.RNA_QC_celltype.0710.pdf"),
  width = 8, height = 4
)

mouse_em@meta.data %>% 
  filter(PLAT!="PLAT15") %>% 
  group_by(CellType) %>% 
  summarise(
    ncell = n(),
    avg_gene = mean(nFeature_RNA),
    avg_umi = mean(nCount_RNA) 
  )

# tsne, plat, cluster
DimPlot(
  subset(mouse_em, PLAT!="PLAT15" & CellType!="M2"),
  reduction = "tsne", group.by = c("PLAT", "seurat_clusters")
)

ggsave(
  paste0(plot_prefix, "mouse_embryo.tsne_QC.0710.pdf"),
  width = 10, height = 4
)

## tsne, celltype
# p <- DimPlot(
#   subset(mouse_em, PLAT!="PLAT15" & CellType!="M2"), 
#   reduction = "tsne", group.by = "CellType",
#   cols = celltype_col
# )
# cr <- diff(range(p$data$tSNE_1))/diff(range(p$data$tSNE_2))
# p + coord_equal(ratio = cr)

DimPlot_scCustom(
  subset(mouse_em, PLAT!="PLAT15" & CellType!="M2"), 
  reduction = "tsne", group.by = "CellType",
  pt.size = 1, figure_plot = TRUE,colors_use = celltype_col,
  split_seurat = T
) 

ggsave(
  paste0(plot_prefix, "mouse_embryo.tsne_celltype.0710.pdf"),
  width = 5, height = 4
)

## tsne, stage
# p <- DimPlot(
#   mouse_em, reduction = "tsne", group.by = "Stage",
#   cols = stage_col
# )
# 
# cr <- diff(range(p$data$tSNE_1))/diff(range(p$data$tSNE_2))
# p + coord_equal(ratio = cr)
DimPlot_scCustom(
  subset(mouse_em, PLAT!="PLAT15" & CellType!="M2"), 
  reduction = "tsne", group.by = "Stage",
  pt.size = 1, figure_plot = TRUE,colors_use = stage_col,
  split_seurat = T
) 

ggsave(
  paste0(plot_prefix, "mouse_embryo.tsne_stage.0710.pdf"),
  width = 5, height = 4
)

## tsne, featureplot
pal <- wesanderson::wes_palette(
  "Zissou1", length(viridis_light_high), 
  type = "continuous"
)

FeaturePlot_scCustom(
  subset(mouse_em, PLAT!="PLAT15" & CellType!="M2"), 
  colors_use = pal,
  features = markers, reduction = "tsne",
  na_cutoff = NULL
)

ggsave(
  paste0(plot_prefix,"mouse_embryo.markers_tsne.0710.pdf"),
  width = 10, height = 8
)

## stats, celltype
mouse_em@meta.data %>% 
  filter(PLAT!="PLAT15" & CellType!="M2") %>% 
  group_by(CellType) %>% 
  summarise(
    avg_gene = mean(nFeature_RNA),
    avg_umi = mean(nCount_RNA)
  )

saveRDS(mouse_em, file = "RNA/nanoCool_mouse.seu_obj.rds")

