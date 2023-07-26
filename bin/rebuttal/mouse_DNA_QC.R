source("/mnt/e/Project/nanoCool/Report/Manuscript/CellResearch/configure.R")

library(dplyr)
library(ggplot2)
library(stringr)
library(RColorBrewer)
library(forcats)

setwd(wkdir)

## load raw input
stat_df <- read.table(
  "DNA/QC/basic_stat/nanoCool.mouse_embryo.merged_0708.basic_stat.tsv",
  header = T, stringsAsFactors = F
)
mouse_em <- readRDS("RNA/nanoCool_mouse.seu_obj.rds")

## filter samples
stat_df <- stat_df %>% 
  mutate(
    batch=str_split(Sample, "_", simplify = T)[,2],
    QC = ifelse(
      Conversion_ratio >= 0.97 & Coverage >= 0.02 & 
        (!batch %in% paste0("P", c(26,27,29,31))) &
        GCH_meth < 0.6 & GCH_meth > 0.2, 
      "Pass", "Fail")) %>% 
  inner_join(mouse_em@meta.data, by=c("Sample"="DNA_ID"))

## plot, raw, batch
stat_df %>% 
  dplyr::select(
    Sample, batch,
    Duplication, Conversion_ratio,Mapping_ratio,
    Median_length, Coverage, WCG_site, GCH_site) %>%
  reshape2::melt(id.var=c("Sample", "batch")) %>%
  ggplot(aes(x=batch, y=value)) +
  ggbeeswarm::geom_quasirandom(aes(color=batch)) +
  facet_grid(variable~., scales = "free") + 
  egg::theme_presentation(base_size = 15) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text.y = element_text(angle = 0, hjust = 0),
        axis.title.x = element_blank())

plot_time <- strsplit(as.character(Sys.time()), " ")[[1]][1]
ggsave(
  paste0(plot_prefix,"mouse_emrbyo.DNA_QC.", plot_time, ".pdf"),
  width = 10, height = 8
)

## stats, QCed, cell type
stat_df %>% 
  filter(QC=="Pass" & CellType!="M2") %>% 
  dplyr::select(
    Sample, CellType,
    Duplication, Conversion_ratio,Mapping_ratio,
    Median_length, Coverage, WCG_site, GCH_site) %>%
  reshape2::melt(id.var=c("Sample", "CellType")) %>% 
  group_by(CellType, variable) %>% 
  summarise(avg_value = mean(value)) %>% 
  tidyr::spread(key = CellType, value = avg_value)

usable_batch <- paste0("P", c(23:25,28, 32,33, 34, 36, 37))
## stats, QCed, batch
stat_df %>% 
  filter(batch %in% usable_batch & CellType!="M2") %>%
  group_by(CellType) %>%
  # group_by(batch) %>% 
  summarise(
    pass=length(which(QC=="Pass")), 
    ngene = mean(nFeature_RNA),
    qc_ratio=pass/n(),
    avg_meth = mean(WCG_meth[QC=="Pass"]),
    avg_acc = mean(GCH_meth[QC=="Pass"])) 

## plot, QCed, celltype
stat_df %>% 
  filter(QC=="Pass" & CellType!="M2") %>% 
  dplyr::select(
    Sample, CellType,
    Duplication, Conversion_ratio,Mapping_ratio,
    Median_length, Coverage, WCG_site, GCH_site) %>%
  reshape2::melt(id.var=c("Sample", "CellType")) %>%
  ggplot(aes(x=fct_relevel(CellType, stage_order), y=value)) +
  # geom_boxplot(aes(color=CellType)) + 
  ggbeeswarm::geom_quasirandom(aes(color=CellType)) +
  scale_color_manual(values = celltype_col) +
  # scale_color_viridis_c() + 
  facet_wrap(variable~., scales = "free_y") + 
  egg::theme_presentation(base_size = 15) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text.y = element_text(angle = 0, hjust = 0),
        axis.title.x = element_blank())

plot_time <- strsplit(as.character(Sys.time()), " ")[[1]][1]
ggsave(
  paste0(plot_prefix,"mouse_emrbyo.DNA_QC_CellType.", plot_time, ".pdf"),
  width = 12, height = 10
)

## plot, QCed, celltype, tagging ratio
stat_df %>% 
  filter(QC=="Pass" & CellType!="M2") %>%
  dplyr::select(
    Sample, CellType,auto_tagratio) %>%
  ggplot(aes(x=fct_relevel(CellType, stage_order), y=auto_tagratio)) +
  # geom_boxplot(aes(color=CellType)) + 
  ggbeeswarm::geom_quasirandom(aes(color=CellType)) +
  scale_color_manual(values = celltype_col) +
  egg::theme_presentation(base_size = 15) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text.y = element_text(angle = 0, hjust = 0),
        axis.title.x = element_blank())

ggsave(
  paste0(plot_prefix,"mouse_emrbyo.DNA_QC_haplotag.0710.pdf"),
  width = 5, height = 4
)


write.table(
  stat_df, file = "DNA/QC/basic_stat/nanoCool.mouse_embryo.merged_0708.basic_stat.tsv",
  col.names = T, row.names = F, sep = "\t", quote = F
)
