source("/mnt/e/Project/nanoCool/Report/Manuscript/CellResearch/rebuttal_R/human_cellline/cellline_config.R")

library(dplyr)
library(ggplot2)
library(data.table)
library(stringr)
library(DescTools)

setwd(wkdir)

## metadata
stat_df <- read.table(stat_dir,header = T, stringsAsFactors = F)
cellsPass <- stat_df$Sample[stat_df$QC=="Pass"]

############################################
## GCH level flanking proximal NDR, K562 ##
############################################
#1) laod raw data
fname <- list.files("absoluteDistance/GCH_proximal_NDR/", pattern = "P13", full.names = T)
# fname <- fname[fname %like% paste(cellsPass, collapse="|")]
fname <- fname[fname %like any% paste0("%",cellsPass, ".GCH_NDR%")]

for(i in 1:length(fname)){
  methy <- read.table(fname[i], header = T)
  methy <- methy[-1,]
  methy$data <- methy$data * 100
  Sample <- gsub(".GCH_NDR.absoluteDistance.txt", "", methy$Sample)
  methy$Coord <- 1:nrow(methy)
  
  ## smooth by spline
  if(length(which(is.na(methy$data)))!=0) next
  # methy <- as.data.frame(spline(methy$Coord, methy$data))
  methy$Sample <- Sample
  # colnames(methy) <- c("Coord", "data", "Sample")
  # print(fname[i])
  
  if(i == 1)
    com.ave.methy <- methy else
      com.ave.methy <- rbind(com.ave.methy, methy)
}

#2) plot, mean & sd level  
p <- com.ave.methy %>% 
  filter(Sample %in% cellsPass) %>%
  group_by(Coord) %>% 
  summarise(data_avg=mean(data), data_sd=sd(data)) %>% 
  ggplot(aes(x=Coord))+
  geom_ribbon(
    aes(ymax=data_avg+data_sd, ymin=data_avg-data_sd), 
    alpha=.2, fill="dodgerblue")+
  geom_line(aes(y=data_avg, x=Coord),size=1.5, color="dodgerblue")+
  scale_x_continuous(
    breaks = c(1,40.5, 80), 
    labels = c("-1kb", "TSS","+1kb"))+
  labs(y="GCH methylation level (%)") + 
  # #coord_fixed(ratio = 0.8)+
  geom_vline(aes(xintercept = 40.5), linetype = "dashed", size = 1) +
  # facet_wrap(.~Sample)+
  theme_classic(base_size = 15)+
  theme(panel.grid.minor.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "bottom",
        axis.text.x = element_text(angle = 0, hjust = .5,vjust = 1),
        plot.title = element_text(vjust = 1, hjust = .5))

plot_time <- strsplit(as.character(Sys.time()), " ")[[1]][1]
ggsave(
  plot = p, width = 7, height = 5, 
  filename = paste0(plot_prefix, "/cool_P13.GCH_proximal_NDR.",plot_time,".pdf")
)

############################################
## signals of published data in NDR, K562 ##
############################################
library(ChIPseeker)
library(genomation)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

#1) load NDR & published data
ndr <- readBed(
  "NDR/nanoCool_5aza.0.CHlt80per.NOMe.genome.NDR_MethGC.3Depth_5GCH.tsv.len_ge240bp.bed"
)

## centralize 
ndr$mid <- (start(ndr) + end(ndr))/2
start(ndr) <- ndr$mid; end(ndr) <- ndr$mid

region <- makeBioRegionFromGranges(
  ndr, type = "start_site", by = "gene", 
  upstream = 3000, downstream = 3000
)

## ENCFF544LXB: H3K27ac; ENCFF616DLO: H3K4me3; ENCFF755HCK: EP300; ENCFF396BZQ: CTCF
## ENCFF333TAT: ATAC-seq; ENCFF163PXS:DNase-seq; ENCFF657CTC: GATA1; ENCFF389FLV: MAZ
encode_files <- list.files("../encode/", full.names = T)
names(encode_files) <- gsub(".bed.gz", "", list.files("../encode/", full.names = F))
encode_files <- encode_files[
  c("ENCFF333TAT", "ENCFF163PXS", "ENCFF616DLO", 
    "ENCFF544LXB", "ENCFF657CTC", "ENCFF389FLV")]

## computeMatrix & plot heatmap
tagMatrixList <- lapply(encode_files, getTagMatrix, windows=region)

plot_time <- strsplit(as.character(Sys.time()), " ")[[1]][1]
pdf(
  paste0(plot_prefix, "/cool_P13.NDR_public.heatmap.",plot_time,".pdf"),
  width = 6, height = 8
)
tagHeatmap(
  tagMatrixList[c("ENCFF333TAT", "ENCFF163PXS")], 
  xlim=c(-3000, 3000), color="#a6cde2",
  xlab = c("ATAC-seq", "DNase-seq")
)

tagHeatmap(
  tagMatrixList[c("ENCFF616DLO", "ENCFF544LXB")], 
  xlim=c(-3000, 3000), color="#38a232",
  xlab = c("H3K4me3","H3K27ac")
)

# tagHeatmap(
#   tagMatrixList[c("ENCFF396BZQ", "ENCFF755HCK")], 
#   xlim=c(-3000, 3000), color="#f37c74",
#   xlab = c("CTCF", "EP300")
# )

tagHeatmap(
  tagMatrixList[c("ENCFF657CTC", "ENCFF389FLV")],
  xlim=c(-3000, 3000), color="#f37c74",
  xlab = c("GATA1", "MAZ")
)

dev.off()

#############################################
## Differential NDR analysis, K562 vs HFF1 ##
#############################################
library(LOLA)
library(GenomicRanges)
library(genomation)
library(forcats)

#1) load raw input
regionDB <- loadRegionDB("/mnt/f/packages/LOLACore/hg38/", collections = "encode_tfbs")

ndr_k562 <- readBed(
  "NDR/nanoCool_5aza.0.CHlt80per.NOMe.genome.NDR_MethGC.3Depth_5GCH.tsv.len_ge240bp.bed"
)
ndr_hff1 <- readBed(
  "NDR/cool_P11_P12.NOMe.genome.NDR_MethGC.3Depth_5GCH.tsv.len_ge240bp.bed"
)

#2) enrichment analysis
useSets <- GRangesList(list(K562=ndr_k562, HFF1=ndr_hff1))
universe <- c(ndr_k562, ndr_hff1)
locResults <- runLOLA(
  useSets, universe, regionDB, 
  minOverlap = 1, direction = "enrichment"
)

#3) plot result, consensus description 
plot_df <- locResults %>% 
  filter(qValue < 1e-3 & treatment=="None" & !grepl("eGFP", antibody)) %>% 
  mutate(antibody=stringr::str_split(antibody, "_", simplify = T)[,1]) %>% 
  mutate(antibody=gsub("-", "", antibody)) %>%
  dplyr::select(userSet, antibody, oddsRatio) %>%
  group_by(userSet, antibody) %>% top_n(1, wt=oddsRatio) %>% 
  # dplyr::slice(order_by=oddsRatio, n=1)
  # distinct(antibody, .keep_all = T) %>% 
  group_by(userSet) %>%
  top_n(5, wt = oddsRatio) 

p <- plot_df %>%
  ggplot(aes(
    x=fct_reorder(antibody, oddsRatio, .desc = F), 
    y=oddsRatio, fill=userSet))+
  geom_bar(stat = "identity", position = "dodge", width = .75) + 
  coord_flip() + facet_wrap(userSet~., scale="free", ncol = 1) +
  scale_fill_manual(values = celltype_pal)+
  scale_y_continuous(expand = c(0,0)) +
  theme_classic(base_size = 15)+
  theme(axis.title.y = element_blank(),
        strip.background = element_blank(),
        # legend.position = "None",
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5),
        # axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(hjust = 1))

plot_time <- strsplit(as.character(Sys.time()), " ")[[1]][1]
ggsave(
  plot = p, width = 4, height = 4, 
  filename = paste0(plot_prefix, "/K562_HFF1.DNDRs_enrich.",plot_time,".pdf")
)

#########################################
##     Correlation between acc & expr  ##
#########################################
library(Seurat)
library(forecast)

#1) load raw input
cellline <- readRDS(seurat_obj)

wcg_gene_meth <- genomation::readBed("multi-omics/nanoCool_5aza.0.CHlt80per.WCG.TSS_2k.bed")
gch_gene_meth <- genomation::readBed("multi-omics/nanoCool_5aza.0.CHlt80per.GCH.TSS_2k.bed")

wcg_gene_meth <- wcg_gene_meth[!duplicated(wcg_gene_meth$name),]
gch_gene_meth <- gch_gene_meth[!duplicated(gch_gene_meth$name),]

meth_df <- inner_join(
  as.data.frame(wcg_gene_meth) %>% dplyr::select(name, score), 
  as.data.frame(gch_gene_meth) %>% dplyr::select(name, score),
  by = "name", suffix = c("_wcg", "_gch")
)
rna_expr_mtx <- GetAssayData(cellline, slot = "data") %>% as.data.frame()

cells_idx <- colnames(cellline)[cellline$Group=="TGS_K562"]
genes_ids <- intersect(meth_df$name, rownames(cellline))
# genes_ids <- Reduce(intersect, list(wcg_gene_meth$name, gch_gene_meth$name, rownames(cellline)))

rna_expr_mtx <- rna_expr_mtx[genes_ids,cells_idx]
rna_expr_df <- data.frame(
  name = rownames(rna_expr_mtx), Mean_expr = rowMeans(rna_expr_mtx)
)

#2) integrate omics & smoothing
mo_df <- inner_join(rna_expr_df, meth_df, by = "name") %>% 
  mutate(gene_rank = order(Mean_expr)) %>% 
  mutate(Expr_level = case_when(
    Mean_expr == 0 ~ "0",
    Mean_expr < 0.05 & Mean_expr > 0 ~ "(0, 0.05)",
    Mean_expr >= 0.05 & Mean_expr < .2 ~ "(0.05, 0.2)",
    Mean_expr >= .2 ~ "(0.2, Inf)")) %>% 
  mutate(expr_bin = cut(Mean_expr, breaks = 5))

mo_df$Expr_level <- factor(
  mo_df$Expr_level, 
  levels = c("0", "(0, 0.05)", "(0.05, 0.2)", "(0.2, Inf)")
)

#2) violin plot
p <- mo_df %>% 
  ggplot(aes(x=Expr_level, y=score_gch))+
  # geom_violin(fill="#b0ddde", adjust=1.5)+
  geom_violin(aes(fill=Expr_level), adjust=1.5)+
  scale_fill_manual(values = wesanderson::wes_palette("Zissou1"))+
  scale_x_discrete(labels = c("No","Low", "Median", "High")) + 
  stat_summary(fun.data = "mean_cl_boot", geom = "pointrange",colour = "black")+
  labs(y="Chromatin accessibility(GCH)", x="Expression level") +
  theme_classic(base_size = 15) 

plot_time <- strsplit(as.character(Sys.time()), " ")[[1]][1]
ggsave(
  plot = p, width = 5, height = 4, 
  filename = paste0(plot_prefix, "/K562_WT.RNA_GCH_cor.level.",plot_time,".pdf")
)

#3) classical line plot
ma_df <- inner_join(rna_expr_df, meth_df, by = "name")  %>%
  arrange(Mean_expr) %>% 
  mutate(
    expr_sm = ma(Mean_expr, order=10, centre=F),
    score_gch_sm = ma(score_gch, order=2000, centre=F),
    gene_rank = order(expr_sm, decreasing = F, na.last = F)
  ) %>% filter(!is.na(score_gch_sm)) 


a.diff <- max(ma_df$expr_sm) - min(ma_df$expr_sm)
b.diff <- max(ma_df$score_gch_sm) - min(ma_df$score_gch_sm)
# b.diff <- 0.1
a.min <- min(ma_df$expr_sm)
b.min <- min(ma_df$score_gch_sm)
# b.min <- 0.3

ratio <- diff(range(ma_df$gene_rank, na.rm=T))/diff(range(ma_df$expr_sm, na.rm=T))

test_res <- cor.test(ma_df$Mean_expr, ma_df$score_gch, method = "spearman")
tt <- paste0("rho=", round(test_res$estimate, 2))

p <- ma_df %>% 
  ggplot(aes(x=gene_rank))+
  geom_point(
    aes(y=Mean_expr), size=1, color="#1c7129")+
  geom_point(
    aes(y=(score_gch_sm - b.min)/b.diff * a.diff + a.min),
    size =1 , color = "#9515e5")+
  scale_y_continuous(
    name = "Mean expression level",
    sec.axis = sec_axis(
      trans = ~((. -a.min) * b.diff / a.diff) + b.min,
      name = "Chromatin accessibility(GCH)"))+
  labs(title = tt)+
  coord_equal(ratio = ratio) + 
  theme_bw(base_size = 15) +
  theme(panel.grid = element_blank())

plot_time <- strsplit(as.character(Sys.time()), " ")[[1]][1]
ggsave(
  plot = p, width = 5, height = 4, 
  filename = paste0(plot_prefix, "/K562_WT.RNA_GCH_cor.line.",plot_time,".pdf")
)

