source("/mnt/e/Project/nanoCool/Report/Manuscript/CellResearch/configure.R")

library(dplyr)
library(forcats)
library(ggplot2)
library(ggh4x)

setwd(wkdir)

## haplotagging
report <- read.table("DNA/haplotag/mouse_embryo.sc.haplotag_report.tsv", header = T)
stat_df <- read.table(
  "DNA/QC/basic_stat/nanoCool.mouse_embryo.merged_0708.basic_stat.tsv",
  header = T, stringsAsFactors = F
)

stat_df <- report %>% 
  inner_join(stat_df, by=c("Sample")) %>% 
  mutate(
    auto_tagratio=(C57B6_auto_reads+DBA_auto_reads)/Total_auto_reads,
    y_ratio=Total_chrY_reads/Total_auto_reads,
    dba_x_ratio=DBA_chrX_reads/DBA_auto_reads,
    sex = case_when(
       dba_x_ratio >= 0.7* y_ratio~"XX", 
       # dba_x_ratio < 0.5* y_ratio & dba_x_ratio >= 0.4* y_ratio~"NA",
       dba_x_ratio < 0.7* y_ratio ~ "XY")
  )

## cor between length & tagratio
# stat_df %>%
#   ggplot(aes(x=Median_length, y=auto_tagratio))+
#   geom_point(aes(color=batch)) +
#   stat_cor(
#     aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),
#     size=5, label.y = 0.10, label.x = 1000) +
#   theme_classic(base_size = 20)


## XCI & XCR

p <- stat_df %>% 
  filter(CellType!="M2" & QC=="Pass") %>% 
  ggplot(aes(x=y_ratio, y=dba_x_ratio)) + geom_density_2d() +
  # geom_point(color="#68c968") + 
  geom_point(aes(color=sex)) +
  scale_color_manual(values = sex_col) +
  geom_abline(slope = 0.7, linetype="longdash") +
  egg::theme_presentation()
# geom_point(aes(color=DBA_chrX_reads/DBA_auto_reads))

cr <- diff(range(p$data$y_ratio))/diff(range(p$data$dba_x_ratio))
p + coord_equal(ratio = cr)

plot_time <- strsplit(as.character(Sys.time()), " ")[[1]][1]
ggsave(
  paste0(plot_prefix,"mouse_emrbyo.diff_sex.",plot_time,".pdf"),
  width = 7, height = 5.5
)

write.table(
  stat_df, file = "DNA/QC/basic_stat/nanoCool.mouse_embryo.merged_0710.basic_stat.tsv",
  col.names = T, row.names = F, sep = "\t", quote = F
)

## pat vs mat methylation level

plot_hap_meth_chr <- function(C_motif, save_plot=F){
  
  switch(
    C_motif,
    "WCG"={
      ymin = 0; ymax = 0.6; ylabel = 0.55; 
      yaxis_t="DNA Methylation level"},
    "GCH"={
      ymin = 0.15; ymax = 0.7; ylabel = 0.65; 
      yaxis_t="Chromatin accessibility"}
  )
  
  hap_meth <- read.table(
    paste0("DNA/haplotag/mouse_embryo.sc.haplo_",C_motif,"_report.merged.0707.tsv"),
    header = T, stringsAsFactors = F, fill = T
  )
  
  hap_meth <- hap_meth %>% 
    inner_join(
      stat_df %>% filter(sex=="XX" & CellType!="M2" & QC=="Pass") %>% 
        select(Sample, CellType, sex), by="Sample")
  hap_meth$CellType <- factor(
    hap_meth$CellType, levels = stage_order
  )
  
  auto_test_res <- hap_meth %>% 
    # dplyr::select(Sample, C57_auto_meth, DBA_auto_meth, CellType) %>% 
    dplyr::group_by(CellType) %>% 
    dplyr::summarise(
      auto_p.val=wilcox.test(
        C57_auto_meth, DBA_auto_meth, paired = T)$p.value,
      x_p.val=wilcox.test(
        C57_chrX_meth, DBA_chrX_meth, paired = T)$p.value) %>% 
    # test_res <- test_res %>% 
    arrange(factor(CellType, levels=stage_order)) %>% 
    mutate(
      auto_padj = p.adjust(auto_p.val, method = "fdr"),
      x_padj = p.adjust(x_p.val, method = "fdr")) %>%
    mutate(
      auto_sig = gtools::stars.pval(auto_padj),
      x_sig = gtools::stars.pval(x_padj)
    )
  
  
  p <- hap_meth %>% 
    select(-sex) %>%
    reshape2::melt(id.var=c("Sample", "CellType")) %>% 
    tidyr::separate(
      col = variable, sep = "_", 
      into = c("Hap", "Chr", "Meth")) %>% 
    ggplot(aes(
      x=interaction(Hap,Chr,CellType),y=value))+
    annotate(
      "rect", xmin = c(4.5,12.5), xmax = c(8.5, 16.5),
      # ymin = 0.15, ymax = 0.7,
      # ymin = 0, ymax = 0.6,
      ymin = ymin, ymax = ymax, 
      alpha = .25,fill = "grey") +
    annotate(
      "text", y = ylabel, size=5,
      label = c(auto_test_res$auto_sig, auto_test_res$x_sig),
      x = c(seq(1.5, 17.5, 4), seq(3.5, 19.5, 4))) +
    geom_boxplot(outlier.shape = NA, fill=NA) + 
    geom_point(aes(color=Hap), size=1.5, shape=2)+
    scale_color_manual(values = hap_color) +
    scale_x_discrete(guide = "axis_nested") +
    scale_y_continuous(expand = c(0,0)) + 
    labs(y=yaxis_t) +
    # labs(y="Chromatin accessibility") +
    theme_classic(base_size = 15) +
    theme(axis.title.x = element_blank()) 
  
  if(save_plot){
    plot_time <- strsplit(as.character(Sys.time()), " ")[[1]][1]
    ggsave(plot = p, width = 12, height = 4,  
      filename = paste0(plot_prefix,"mouse_emrbyo.hap_meth_", C_motif,".",plot_time,".pdf")
    )}else{print(p)}
  
} 

plot_hap_meth_chr("WCG")
plot_hap_meth_chr("GCH")

## regional hap-reserved meth & acc
g1 <- read.table("DNA/haplotag/mouse_embryo.WCG_g1.auto_gene_meth.txt", header = T)
g2 <- read.table("DNA/haplotag/mouse_embryo.WCG_g2.auto_gene_meth.txt", header = T)

region_df <- inner_join(
  g1, g2, by=c("Sample"),
  suffix = c(".g1", ".g2")) %>% 
  inner_join(by = "Sample", 
  stat_df %>% 
    filter(sex=="XX" & QC=="Pass" & CellType!="M2") %>% 
    select(Sample, CellType))%>% 
  mutate(
    intra=Intragenic_meth.g1-Intragenic_meth.g2,
    intre=Intergenic_meth.g1-Intergenic_meth.g2) %>% 
  dplyr::select(Sample, intra, intre, CellType) %>%
  reshape2::melt(id.var=c("Sample", "CellType")) %>% 
  mutate(Hap=ifelse(value>0, "C57B6", "DBA")) 

region_df$CellType <- factor(region_df$CellType, levels = stage_order)
region_df <- region_df %>% 
  group_by(variable, CellType) %>%
  arrange(desc(value))

p2 <- region_df %>% 
  filter(variable=="intre") %>%
  ggplot(aes(x=forcats::fct_reorder(Sample,value,.desc = T), y=value, fill=Hap))+
  geom_col() + facet_grid2(~CellType, scales = "free_x") +
  labs(y="DNA Methylation(mat-pat, %)", x="Samples", title = "Intergenic") +
  scale_fill_manual(values = hap_color) +
  egg::theme_article() +
  theme(axis.text.x = element_blank(),
        plot.title = element_text(hjust=.5),
        axis.ticks.x = element_blank())

p1 <- region_df %>% 
  filter(variable=="intra") %>%
  ggplot(aes(x=forcats::fct_reorder(Sample,value,.desc = T), y=value, fill=Hap))+
  geom_col() + facet_grid2(~CellType, scales = "free_x") +
  labs(y="DNA Methylation(mat-pat, %)", title = "Intragenic") +
  scale_fill_manual(values = hap_color) +
  egg::theme_article() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        plot.title = element_text(hjust=.5),
        axis.ticks.x = element_blank())

cowplot::plot_grid(p1, p2, align = "hv", ncol = 1)

plot_time <- strsplit(as.character(Sys.time()), " ")[[1]][1]
ggsave( 
  paste0(plot_prefix,"mouse_emrbyo.hap_meth_region.",plot_time,".pdf"),
  width = 12, height = 5.5
)

## accumulative coverage curve
## WGBS_GSR: 0.05653135 (C57) & 0.065860875 (DBA)
plot_hap_cov <- function(x){
  switch(x,
         "g1"={ref_cov=0.05653135;color=hap_color["C57"];tt="C57B6/2N"},
         "g2"={ref_cov=0.065860875;color=hap_color["DBA"];tt="DBA/2N"}
  )
  # setwd("/mnt/e/Project/nanoCool/Data/mouse_embryo")
  # setwd(wkdir)
  fn <- paste0("DNA/haplotag/mouse_ICM.haplo_cov.sc_curve.", x, ".merged.tsv")
  
  hap_cov <- read.table(fn, header = F)
  col_names <- c("Sample", "hap", "frac", "ncell") 
  colnames(hap_cov) <- c(col_names, paste0("k_", 1:(ncol(hap_cov)-length(col_names))))
  hap_cov[nrow(hap_cov)+1,] <- hap_cov[nrow(hap_cov),]
  hap_cov[nrow(hap_cov),3:ncol(hap_cov)] <- 0
  p <- hap_cov %>% 
    reshape2::melt(id.var=col_names) %>%
    ggplot(aes(x=ncell, y=value))+
    geom_smooth(color=color, fill=color, alpha=.3, se=F) + 
    geom_hline(yintercept = ref_cov, linetype="longdash") + 
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    labs(x="No. Cell", y="Genomic coverage", title = tt) + 
    theme_classic(base_size = 15)+
    theme(plot.title = element_text(hjust = .5))
  print(p)
}

p1 <- plot_hap_cov("g1")
p2 <- plot_hap_cov("g2")

y_max <- max(range(p1$data$value), range(p2$data$value))
p1 <- p1 + ylim(0,y_max)
p2 <- p2 + ylim(0,y_max)

cowplot::plot_grid(p1,p2, rel_widths = c(1,1))

plot_time <- strsplit(as.character(Sys.time()), " ")[[1]][1]
ggsave( 
  paste0(plot_prefix,"mouse_emrbyo.hap_coverage.",plot_time,".pdf"),
  width = 8, height = 4
)

## benchmark, oocyte
# hap_report <- read.table("DNA/haplotag/cool_P37.sc.haplo_report.merged.tsv", header = T)
hap_report <- stat_df %>% 
  filter(CellType=="M2" & QC=="Pass" & Total_auto_reads > 2e5) %>% 
  mutate(
    hap_auto_total=C57B6_auto_reads+DBA_auto_reads,
    C57=C57B6_auto_reads/hap_auto_total,
    DBA=DBA_auto_reads/hap_auto_total) %>% 
  select(Stage, C57, DBA) %>% 
  reshape2::melt(id.var="Stage") 

hap_report %>%
  ggplot(aes(x=Stage, y=value, color=variable))+
  ggbeeswarm::geom_quasirandom()+ ylim(c(0,1)) +
  scale_color_manual(values = hap_color) +
  labs(y="Haplotype-ratio") +
  theme_bw(base_size = 20) +
  theme(axis.title.x = element_blank())

plot_time <- strsplit(as.character(Sys.time()), " ")[[1]][1]
ggsave( 
  paste0(plot_prefix,"mouse_emrbyo.oocyte_haptag_ratio.",plot_time,".pdf"),
  width = 6, height = 4
)

stat_df %>%
  filter(QC=="Pass") %>% 
  ggplot(aes(x=Median_length, y=auto_tagratio))+ 
  geom_point(color="dodgerblue") +
  stat_cor(
    aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),
    size=5, label.y = 0.10, label.x = 1000) +
  theme_classic(base_size = 20)

plot_time <- strsplit(as.character(Sys.time()), " ")[[1]][1]
ggsave( 
  paste0(plot_prefix,"mouse_emrbyo.haptag_length.",plot_time,".pdf"),
  width = 5.5, height = 5
)

# stat
hap_report %>% 
  group_by(Stage, variable) %>% 
  summarise(precision=mean(value))

## SNPs dist 
snp_df <- read.table("DNA/anno/all_DBA_2J_SNPs_C57BL_6NJ_reference.based_on_GRCm38.sorted.bed")
snp_df <- (snp_df$V2+snp_df$V3)/2
snp_diff <- diff(sort(snp_df), lag=1)
snp_diff <- data.frame(dist=log10(snp_diff+1))

snp_diff %>% 
  ggplot(aes(x=dist))+stat_ecdf() +
  geom_vline(xintercept = log10(900+1), linetype="longdash")

hist(snp_diff$dist[snp_diff>3], breaks = 100)


