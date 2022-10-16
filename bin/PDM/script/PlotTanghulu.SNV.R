library(dplyr)
library(ggplot2)

cmd <- commandArgs(T)

input <- cmd[1]
output <- cmd[2]
region <- cmd[3]
snv <- cmd[4]

snv <- as.character(snv)
# setwd("E:/OneDrive/VSCODE/pipeline/nanoCool_xxh/demo")
# input <- "cool_P12.CHlt80per.RB1_2.WCG_addCB.txt"
# cgi_df <- read.table("hg38.CGI.bed", stringsAsFactors = F)
# cgi <- strsplit(input, "[.]")[[1]][3]
# region <- paste0(cgi_df[cgi_df$V4==cgi,1:3], collapse = ":")

df <- read.table(input, header = T, stringsAsFactors = F)
# df <- read.table("F:/Project/nanoCool/Data/5aza/TRA/test/ABL1_Mut_GCH.v3.txt", header = T, stringsAsFactors = F)
#df$read <- paste0(df$cell_barcode, "_", rownames(df))
#df_long <- reshape2::melt(df, id.var=c("cell_barcode", "read"))
df$read <- rownames(df)

snv_col_idx <- which(colnames(df)==paste0("X", snv))
snv_df <- df[,c(snv_col_idx,ncol(df))]
# snv_df <- snv_df[snv_df[1,]!="NULL",]

df <- df[,-snv_col_idx]
df_long <- reshape2::melt(df, id.var=c("read")) %>% filter(value!="NULL")
df_long$value <- as.integer(df_long$value)

region <- strsplit(region, ":|-")[[1]]

df_long$variable <- gsub("X","",df_long$variable)

read_order <- df_long %>% 
  group_by(read) %>%
  # filter(value!="NULL") %>%
  dplyr::summarise(
    nsite=n(), meth=mean(value),
    start=min(as.integer(variable)),
    end=max(as.integer(variable))) %>%
  dplyr::mutate(meth_group=case_when(
    meth <=0.5 ~ "low",
    meth > 0.5 ~"high")) %>%
  dplyr::group_by(meth_group) %>%
  dplyr::arrange(nsite, meth) %>%
  dplyr::filter(
    nsite >=3 & 
      start <= as.integer(snv) & 
      end >= as.integer(snv)
  )

stopifnot(nrow(read_order)>=3)

df_long <- df_long %>% filter(value!="NULL") %>%
  filter(read %in% read_order$read)
snv_df <- snv_df %>% filter(read %in% read_order$read)
snv_df <- snv_df[snv_df[,1]!="NULL",]
# snv_df[,1] <- gsub("NULL", "", snv_df[,1])

df_long <- df_long %>% filter(read %in% snv_df$read)
read_order <- read_order %>% filter(read %in% snv_df$read)

df_long$read <- factor(df_long$read, levels = read_order$read)

if(nrow(df_long)!=0){
  title <- paste0(region)
  p <- df_long %>% 
    ggplot(aes(x=variable, y=read))+
    # annotate("segment", x = region[2], xend = region[2], 
    #          y = Inf,yend = -Inf, linetype='longdash')+
    # annotate("segment", x = region[3], xend = region[3], 
    #          y = Inf,yend = -Inf, linetype='longdash')+
    geom_line(aes(group=read))+
    annotate("segment", x = snv, xend = snv,
             y = Inf,yend = -Inf, size=5, color="white")+
    annotate("text", x = snv, y = snv_df$read, label=snv_df[,1], size=5)+
    geom_point(aes(fill=as.character(value)), shape=21, size=5)+
    # facet_grid(vars(cell_barcode), scales = "free_y", space = "free_y")+
    theme_bw()+labs(title=title)+
    scale_fill_manual(values = c("0"="white", "1"="black"))+
    theme(panel.grid = element_blank(),
          axis.text.y = element_blank(),
          plot.title = element_text(hjust = 0.5, size=10),
          legend.title = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          axis.title = element_blank(),
          axis.ticks.y = element_blank()) 
  print(p)

  plot_width <- length(unique(df_long$variable))*0.5+0.5
  plot_height <- length(unique(df_long$read))*0.1+2.5
  
  ggsave(filename = output, width = plot_width, height = plot_height)
}

#plot_width <- length(unique(df_long$variable))*1+1
#plot_height <- length(unique(df_long$read))*0.5+2

#ggsave(filename = output, width = plot_width, height = plot_height)
