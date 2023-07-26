source("/mnt/e/Project/nanoCool/Report/Manuscript/CellResearch/rebuttal_R/human_cellline/cellline_config.R")

library(dplyr)
library(ggplot2)
library(data.table)
library(stringr)

setwd(wkdir)

stat_df <- read.table(stat_dir,header = T, stringsAsFactors = F)

#####################################
## QC, WCG & genebody, pass & fail ##
#####################################

#1) load raw data
fname <- list.files("absoluteDistance/WCG_genebody/cool_P11_P12_P13_P14/")

for(i in 1:length(fname)){
  methy <- fread(paste0(indir, fname[i]), header = T)
  methy <- data.frame(methy)
  methy$data <- methy$data * 100
  
  if(i == 1)
    com.ave.methy <- methy else
      com.ave.methy <- rbind(com.ave.methy, methy)
}

colnames(com.ave.methy) <- c("methylation", "bin", "sample")

#2) plot pass profile
p <- com.ave.methy %>%
  filter(sample %in% stat_df$Sample[stat_df$WCG_genebody_profile=="Pass"]) %>% 
  filter(grepl("P11", sample)) %>%
  ggplot(aes(x = bin, y = methylation)) +
  geom_line(position = position_dodge(0.2), size = 1.5, color = "dodgerblue") +
  xlab("Region") + ylab("DNA Methylation Level (%)") +
  geom_vline(aes(xintercept = 50), linetype = "dashed", size = 1) +
  geom_vline(aes(xintercept = 250), linetype = "dashed", size = 1) +
  theme_classic() +
  facet_wrap(. ~ sample)
# dev.off()

plot_time <- strsplit(as.character(Sys.time()), " ")[[1]][1]
ggsave(plot = p, width = 20, height = 10, 
  filename = paste0(plot_prefix, "/cool_P11-P14.DNA_QC.WCG_genebody_profile.Pass.",plot_time,".pdf")
)

#3) plot fail profile
p <- com.ave.methy %>%
  filter(sample %in% stat_df$Sample[stat_df$WCG_genebody_profile=="Fail"]) %>% 
  filter(grepl("P11", sample)) %>%
  ggplot(aes(x = bin, y = methylation)) +
  geom_line(position = position_dodge(0.2), size = 1.5, color = "dodgerblue") +
  xlab("Region") + ylab("DNA Methylation Level (%)") +
  geom_vline(aes(xintercept = 50), linetype = "dashed", size = 1) +
  geom_vline(aes(xintercept = 250), linetype = "dashed", size = 1) +
  theme_classic() +
  facet_wrap(. ~ sample)
# dev.off()

plot_time <- strsplit(as.character(Sys.time()), " ")[[1]][1]
ggsave(plot = p, width = 10, height = 8, 
       filename = paste0(plot_prefix, "/cool_P11-P14.DNA_QC.WCG_genebody_profile.Fail.",plot_time,".pdf")
)

#4) automatic qc, WCG & genebody 
p <- com.ave.methy %>% 
  filter(bin>=100 & bin<=200) %>% 
  group_by(sample) %>% 
  summarise(cv=sd(methylation)/mean(methylation)) %>% 
  inner_join(
    stat_df %>% select(Sample, WCG_genebody_profile),
    by = c("sample"="Sample")
  ) %>% ggplot(aes(x=WCG_genebody_profile, y=cv, fill=WCG_genebody_profile)) +
  geom_violin(scale = "width", adjust = 2)+
  scale_fill_manual(values = qc_pal)+
  theme_classic(base_size = 15) +
  theme(legend.position = "none") + 
  labs(y="CV of the gene body", x="WCG_GENE_profile")

plot_time <- strsplit(as.character(Sys.time()), " ")[[1]][1]
ggsave(
  plot = p, width = 4, height = 4, 
  paste0(plot_prefix, "/cool_P11-P14.DNA_QC.WCG_GENE_CV.",plot_time,".pdf")
)

################################
## QC, GCH & TSS, pass & fail ##
################################

#1) load raw data
fname <- list.files(
  "absoluteDistance/GCH_TSS/cool_P11_P12_P13_P14/",
  pattern = ".GCH_TSS_avemethyrate.bed", full.names = T
)

for(i in 1:length(fname)){
  if(file.size(fname[i])==0) next
  methy <- fread(fname[i], header=F)
  methy <- data.frame(methy)
  colnames(methy) <- c("chr","bin.start","bin.end","gene","bin","methy.rate")
  
  # Calculate average methylation level in each bin #
  bin.methy <- aggregate(methy$methy.rate, list(methy$bin), mean)
  colnames(bin.methy) <- c("bin","methylation")
  # bin.methy$sample <- sample_name[i]
  fn <- rev(str_split(fname[i], "/")[[1]])[1]
  bin.methy$sample <- gsub(".GCH_TSS_avemethyrate.bed", "", fn)
  bin.methy$methylation <- bin.methy$methylation * 100
  
  if(i == 1)
    com.ave.methy <- bin.methy else
      com.ave.methy <- rbind(com.ave.methy, bin.methy)
}

#2) plot pass profile
p <- com.ave.methy %>% 
  filter(sample %in% stat_df$Sample[stat_df$GCH_TSS_profile=="Pass"]) %>% 
  filter(grepl("P11", sample)) %>% 
  ggplot(aes(x = bin, y = methylation)) +
  geom_line(position = position_dodge(0.2), size = 1, color = "dodgerblue") +
  xlab("Region") + ylab("DNA Methylation Level (%)") +
  geom_vline(aes(xintercept = 10.5), linetype = "dashed", size = 1) +
  # scale_x_continuous(breaks = c(0,10,20), labels = c("-20k", "TSS", "+20k"))+
  theme_classic() + facet_wrap(. ~ sample)+
  theme(axis.title.x = element_blank())

plot_time <- strsplit(as.character(Sys.time()), " ")[[1]][1]
ggsave(plot = p, width = 20, height = 10, 
       filename = paste0(plot_prefix, "/cool_P11-P14.DNA_QC.GCH_TSS_profile.Pass.",plot_time,".pdf")
)

#3) plot pass profile
p <- com.ave.methy %>% 
  filter(sample %in% stat_df$Sample[stat_df$GCH_TSS_profile=="Fail"]) %>% 
  filter(grepl("P11", sample)) %>% 
  ggplot(aes(x = bin, y = methylation)) +
  geom_line(position = position_dodge(0.2), size = 1, color = "dodgerblue") +
  xlab("Region") + ylab("DNA Methylation Level (%)") +
  geom_vline(aes(xintercept = 10.5), linetype = "dashed", size = 1) +
  # scale_x_continuous(breaks = c(0,10,20), labels = c("-20k", "TSS", "+20k"))+
  theme_classic() + facet_wrap(. ~ sample)+
  theme(axis.title.x = element_blank())

plot_time <- strsplit(as.character(Sys.time()), " ")[[1]][1]
ggsave(plot = p, width = 10, height = 8, 
       filename = paste0(plot_prefix, "/cool_P11-P14.DNA_QC.GCH_TSS_profile.Fail.",plot_time,".pdf")
)

#4) automatic qc, GCH & TSS 
p <- com.ave.methy %>% 
  filter(bin %in% c(1,10,20)) %>% 
  group_by(sample) %>% 
  filter(n()==3) %>% 
  summarise(
    tssenrich=methylation[bin==10]/mean(methylation[bin==1], methylation[bin==20])
  ) %>% inner_join(
    stat_df %>% dplyr::select(Sample, GCH_TSS_profile),
    by = c("sample"="Sample")
  ) %>% ggplot(aes(x=GCH_TSS_profile, y=tssenrich, fill=GCH_TSS_profile)) +
  geom_violin(scale = "width", adjust = 2) + 
  scale_fill_manual(values = qc_pal)+
  theme_classic(base_size = 15) +
  theme(legend.position = "none") + 
  labs(y="TSS enrichment", x="GCH_TSS_profile")

plot_time <- strsplit(as.character(Sys.time()), " ")[[1]][1]
ggsave(plot = p, width = 4, height = 4, 
       filename = paste0(plot_prefix, "/cool_P11-P14.DNA_QC.GCH_TSS_enrich.",plot_time,".pdf")
)

#####################################
## QC, conversion ratio & coverage ##
#####################################

#100cell/run
p <- stat_df %>% 
  filter(grepl("P11|P12|P13|P14", Sample)) %>% 
  mutate(QC=ifelse(
    WCG_genebody_profile=="Pass" & 
      GCH_TSS_profile=="Pass" &
      Conversion_ratio >=0.98,
    "Pass", "Fail"
  )) %>% ggplot(aes(x=Coverage, y=Conversion_ratio)) +
  geom_point(aes(shape=QC, color=Cell), size=1) + 
  scale_color_manual(values = celltype_pal) +
  scale_shape_manual(values = c("Pass"=0, "Fail"=8)) +
  geom_hline(yintercept = 98, linetype="longdash") +
  geom_vline(xintercept = 2, linetype="longdash") + 
  scale_x_continuous(breaks = c(0,2,10,20), labels = c(0,2,10,20))+
  theme_bw(base_size = 15) +
  theme(panel.grid = element_blank())

cr <- diff(range(p$data$Coverage))/diff(range(p$data$Conversion_ratio))
p <- p + coord_equal(ratio = cr)

plot_time <- strsplit(as.character(Sys.time()), " ")[[1]][1]
ggsave(
  plot = p, filename = paste0(plot_prefix,"/cool_P11-P14.WCG_genebody_profile.CV.",plot_time,".pdf"),
  width = 5, height = 5, units = "in"
)


###########################
## QC, depth & coverage  ##
###########################

## 10cell/run

stat_10_df <- read.table(
  "basic_stat/pre_res/cool_P7-P10.merged.basic_stat.txt",
  header = T, stringsAsFactors = F
)

stat_10_df <- stat_10_df %>% mutate(
  QC = ifelse(
    WCG_genebody_profile=="Pass" & 
      GCH_TSS_profile=="Pass" & 
      Conversion_ratio >= 0.98,
    yes = "Pass", no = "Fail"),
  CellType = case_when(
    grepl("P[78]", Sample) ~ "HFF1",
    grepl("P9|P10", Sample) ~ "K562"
  ))

stat_10_df %>% 
  group_by(QC, CellType) %>% 
  summarise(
    gch_site_n = mean(GCH_site),
    gch_site_p = gch_site_n/235107734
  )

## GCH density: 
2934860425/235107734


