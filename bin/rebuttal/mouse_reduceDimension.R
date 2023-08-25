source("/mnt/e/Project/nanoCool/Report/Manuscript/CellResearch/rebuttal_R/mouse_embryo/mouse_config.R")

library(dplyr)
library(pcaMethods)
library(ggplot2)
library(tsne)
library(umap)
library(stringr)

setwd(wkdir)
stat_df <- read.table(
  "DNA/QC/basic_stat/nanoCool.mouse_embryo.merged_0708.basic_stat.tsv",
  header = T, stringsAsFactors = F
)

##--------------------------raw imputation ----------------------

reduce_dim <- function(
    fn_input, reduction = "tsne", 
    per_na = 0.3, npc = 5, impute_method = "ppca", 
    na_as_zero = F, highly_var_feature = F, 
    fn_output = "None", save_p = FALSE){
  # percent <- 0.3
  # CG_df <- read.csv("DNA/reduce_dimension/win100k_embryo.mtx_WCG_5.csv",header=T)
  CG_df <- read.csv(fn_input)
  pos_df <- CG_df[,1]
  CG_df <- CG_df[,2:ncol(CG_df)]
  
  if(isFALSE(na_as_zero)){
    use_Pos <- apply(CG_df,1,function(x){sum(is.na(x))<per_na*ncol(CG_df)})
    
    use_GC_Merge <- CG_df[use_Pos,]
    use_GC_pos <- pos_df[use_Pos]
    use_GC_Merge <- as.matrix(use_GC_Merge)
    rownames(use_GC_Merge) <- use_GC_pos
    use_GC_Merge[use_GC_Merge=="NaN"] <- NA
  } else if(na_as_zero){
    CG_df[is.na(CG_df) | CG_df < 0.5] <- 0
    CG_df[CG_df >= 0.5] <- 1
    
    use_GC_Merge <- CG_df[!duplicated(pos_df),]
    rownames(use_GC_Merge) <- pos_df[!duplicated(pos_df)]
  }
  
  ## PCA
  set.seed(0)
  
  pc <- pca(use_GC_Merge,nPcs=npc,method=impute_method)
  imputed <- completeObs(pc)
  
  if(highly_var_feature){
    # stopifnot(n>2000)
    
    hvf <- top_frac(data.frame(sd=apply(imputed, 1, sd)), wt = sd, n = .2)
    imputed <- imputed[rownames(hvf),]
  }
  
  pca <- prcomp(t(imputed))
  
  n <- nrow(imputed)
  fn_tmp <- rev(strsplit(fn_input, "/")[[1]])[1]
  reg <- strsplit(fn_tmp, "[.]")[[1]][1]
  omic <- strsplit(fn_tmp, "[.]")[[1]][2]
  omic <- strsplit(omic, "_")[[1]][2]
  
  tt <- paste0(
    "Cmotif:", omic,"\nRegion:", reg, 
    "\nN:", n, "\tpercent_na:", per_na, 
    "\nnPC:", npc, "\tImpute:", impute_method
  )
  
  if(reduction=="tsne"){
    ## tSNE
    tsne_input <- as.data.frame(pca$x[,1:npc])
    mat_tsne<-tsne(tsne_input, max_iter = 1000)
    df_tsne<-data.frame(mat_tsne)
    colnames(df_tsne) <- paste0("tSNE_", 1:2)
    rownames(df_tsne) <- rownames(pca$x)
    
    df_tsne <- df_tsne %>% 
      tibble::rownames_to_column("Sample") %>% 
      inner_join(stat_df %>% select(Sample, CellType, Stage), by="Sample")
    
    p2 <- df_tsne %>% 
      ggplot(aes(x = tSNE_1, y = tSNE_2, color=CellType)) +
      geom_point() + labs(title = tt)  + 
      scale_color_manual(values = celltype_col) + 
      theme_bw(base_size = 10)
    axis_ratio <- diff(range(p2$data$tSNE_1))/diff(range(p2$data$tSNE_2))
    p2 <- p2 + coord_equal(ratio = axis_ratio) +
      theme(panel.grid = element_blank()) +
      theme(axis.ticks = element_blank(),
            axis.text = element_blank())
    
    if(save_p){
      plot_time <- strsplit(as.character(Sys.time()), " ")[[1]][1]
      ggsave(
        plot = p2, width = 7, height = 5, 
        filename = paste0(plot_prefix, fn_output,plot_time,".pdf")
      )
    } else {
      print(p2); return(p2)
    }}
  if(reduction=="umap"){
    ## UMAP
    # embryo_dist <- as.matrix(dist(t(imputed), method = "euclidean"))
    # mat_umap <- umap::umap(t(imputed))
    # df_umap <- as.data.frame(mat_umap$layout)
    # colnames(df_umap) <- paste0("UMAP_", 1:2)
    umap_input <- as.data.frame(pca$x[,1:npc])
    df_umap <- uwot::umap(umap_input, init = "pca")
    colnames(df_umap) <- paste0("UMAP_", 1:2)
    
    # library(uwot)
    # df_map <- uwot::umap(t(imputed), pca = 30, pca_center=T)
    # colnames(df_umap) <- paste0("UMAP_", 1:2)
    
    df_map <- df_umap %>% as.data.frame() %>% 
      tibble::rownames_to_column("Sample") %>% 
      inner_join(stat_df %>% select(Sample, CellType), by="Sample")
    
    p3 <- df_map %>% ggplot(aes(x = UMAP_1, y = UMAP_2)) +
      geom_point(aes(color=CellType)) + labs(title = tt)  + 
      scale_color_manual(values =celltype_col) +
      theme_bw(base_size = 10)
    axis_ratio <- diff(range(p3$data$UMAP_1))/diff(range(p3$data$UMAP_2))
    p3 <- p3 + coord_equal(ratio = axis_ratio) +
      theme(panel.grid = element_blank())+
      theme(axis.ticks = element_blank(),
            axis.text = element_blank())
    
    if(save_p){
      plot_time <- strsplit(as.character(Sys.time()), " ")[[1]][1]
      ggsave(
        plot = p3, width = 7, height = 5, 
        filename = paste0(plot_prefix, fn_output,plot_time,".pdf")
      )
    } else {
      print(p3); return(p3)
  }} 
  else {
    print("Un-recognized reduction method!")
    print("available methods: tsne, umap!")
    break
  }

  # gc()
}

# reduce_dim(
#   "DNA/reduce_dimension/Genebody.mtx_WCG_3.csv",
#   reduction = "tsne", per_na = 0.3, npc = 5,
#   impute_method = "ppca", fn_output = NA, save_p = F
# )

fns <- list.files("DNA/reduce_dimension/WCG/", pattern = "win", full.names = T) 
# fns <- str_sort(fns, numeric = T)
pl <- parallel::mclapply(fns, reduce_dim, mc.cores = 4)

pl <- purrr::discard(.x = pl, inherits, 'try-error')
pl <- purrr::discard(.x = pl, inherits, 'NULL')
cowplot::plot_grid(plotlist = pl)

plot_time <- strsplit(as.character(Sys.time()), " ")[[1]][1]
ggsave(
  width = 8, height = 4, 
  filename = paste0(plot_prefix, "mouse_embryo.5PC_30Na_tSNE.NDR.NaN_as_Zero.HVF.",plot_time,".pdf")
)
