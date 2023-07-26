## configures for rebuttal.R

## paths
wkdir <- "/mnt/e/Project/nanoCool/Data/mouse_embryo"
plot_prefix <- "/mnt/e/Project/nanoCool/Report/Manuscript/CellResearch/rebuttal_plot/"
seurat_obj_path <- "RNA/nanoCool_mouse.seu_obj.rds"

## colors
hap_color <- c("C57"="#ee312b", "DBA"="#2e3192")

celltype_col <- c("#ea5050", "#329f5a", "#7893d3","#ecad0c", "#4aaae1", "grey", "#ee312b", "#2e3192")
names(celltype_col) <- c("Epi", "ICM", "PrE", "Late_TE", "Early_TE", "M2", "C57_M2", "DBA_M2")

stage_col <- c("#e29fc3", "#c8cdf6", "#ecbb7b","#f56264")
names(stage_col) <- c("C57_M2", "DBA_M2", "E3.5", "E4.5")

sex_col <- c("XX"="#ed6aab", "XY"="#70caee", "NA"="grey")

## markers
markers <- c(
  "Oct4", "Nanog", "Sox2","Sox17", "Pdgfra","Gata4", "Gata6","Cdx2", "Gata3","Tead4"
)
# markers <- c(
#   "Pou5f1", "Nanog", "Sox2",
#   "Sox17", "Pdgfra","Gata4", "Gata6",
#   "Cdx2", "Gata3", "Eomes", "Tead4",
#   "Klf4", "Ctcf", "Esrrb"
# )

stage_order <- c("ICM", "Early_TE", "Epi", "PrE", "Late_TE")

