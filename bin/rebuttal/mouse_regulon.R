source("/mnt/e/Project/nanoCool/Report/Manuscript/CellResearch/configure.R")

library(Seurat)
library(SCopeLoomR)
library(SCENIC)

setwd(wkdir)

## prepare input for pyscenic
mouse_em <- readRDS(seurat_obj_path)
exprMat <- as.matrix(mouse_em@assays$RNA@counts)
loom <- build_loom("RNA/pyscenic/Mouse_embryo.loom", dgem=exprMat)
close_loom()

## extract result from output of pyscenic
loom <- open_loom("RNA/pyscenic/Mouse_embryo.out_pyscenic.loom") 
regulons_incidMat <- get_regulons(loom, column.attr.name="Regulons")
# regulons_incidMat[1:4,1:4] 
regulons <- regulonsToGeneLists(regulons_incidMat)
regulonAUC <- get_regulons_AUC(loom,column.attr.name='RegulonsAUC')
regulonAucThresholds <- get_regulon_thresholds(loom)
tail(regulonAucThresholds[order(as.numeric(names(regulonAucThresholds)))])

embeddings <- get_embeddings(loom)  
close_loom(loom)

# rownames(regulonAUC)
# names(regulons)

sub_regulonAUC <- regulonAUC[,intersect(colnames(mouse_em),colnames(regulonAUC))]
sub_regulonAUC <- sub_regulonAUC[onlyNonDuplicatedExtended(rownames(sub_regulonAUC)),]
dim(sub_regulonAUC)

## RSS
rss <- calcRSS(
  AUC=getAUC(sub_regulonAUC), 
  cellAnnotation=mouse_em@meta.data[colnames(sub_regulonAUC),"CellType"]
) 
rss <- na.omit(rss)
write.table(
  as.data.frame(rss) %>% tibble::rownames_to_column("regulon"),
  file = "RNA/pyscenic/Mouse_embryo.pyscenic_out_0622.RSS.tsv",
  col.names = T, row.names = F, sep = "\t", quote = F)

plot_df <- rss %>% 
  as.data.frame()  %>% 
  tibble::rownames_to_column("regulon") %>% 
  reshape2::melt(id.var = "regulon") %>%
  mutate(sig = -log10(value+1)) %>% 
  dplyr::group_by(variable) %>% 
  # arrange(desc(value)) %>%
  # slice(1:3)
  # top_n(n = 10, wt = value) %>%
  # mutate(
  #   ranks = base::rank(sig, ties.method = "average"),
  #   label = ifelse(ranks %in% 1:5, regulon, ""))
  mutate(
    ranks = base::rank(sig, ties.method = "average"),
    label = if_else(
      regulon %in% paste0(markers, "(+)") & ranks %in% 1:50, regulon, ""))

plot_df %>% 
  ggplot(aes(x=ranks, y=value)) +
  geom_point(color="dodgerblue") + 
  ggrepel::geom_label_repel(
    aes(label=label),
    nudge_x = 0.5,
    segment.size = 0.2,
    box.padding = 0.5, 
    # direction = "x", hjust = "right", 
    max.overlaps = Inf) +
  facet_wrap(.~variable, scales = "free") + 
  labs(y = "Regulon specific score") +
  egg::theme_presentation(base_size = 20)

plot_time <- strsplit(as.character(Sys.time()), " ")[[1]][1]
ggsave( 
  paste0(plot_prefix,"mouse_emrbyo.pyscenic.RSS_celltype.",plot_time,".pdf"),
  width = 8, height = 6
)

## add assay to seurat object
regulon.as <- CreateAssayObject(getAUC(sub_regulonAUC))
mouse_em[["Regulon"]] <- CreateAssayObject(getAUC(sub_regulonAUC))

mouse_em <- subset(mouse_em, sex!="NA")
DefaultAssay(mouse_em) <- "Regulon"
mouse_em <- FindVariableFeatures(mouse_em)
mouse_em <- ScaleData(mouse_em, assay = "Regulon")
mouse_em <- RunPCA(mouse_em, assay = "Regulon", features = VariableFeatures(mouse_em))
mouse_em <- RunUMAP(
  mouse_em, dims = 1:50, assay = "Regulon", 
  reduction.name = "regulon.map", reduction.key = "regUMAP_"
)
DimPlot(
  mouse_em, group.by = "CellType", 
  reduction = "regulon.map", cols = celltype_col
)

## regulon & target
adj <- read.table('RNA/pyscenic/Mouse_embryo.adjacencies.tsv', header = T)
names(regulons) <- gsub("[(+)]", "", names(regulons))

saveRDS(regulons, file = "RNA/pyscenic/Mouse_embryo.regulons.rds")
