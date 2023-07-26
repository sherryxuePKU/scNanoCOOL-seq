## paths
wkdir <- "/mnt/e/Project/nanoCool/Data/cell_line/"
stat_dir <- "basic_stat/TGS_NGS.merged.basic_stat.0711.txt"
plot_prefix <- "/mnt/e/Project/nanoCool/Report/Manuscript/CellResearch/rebuttal_plot/"
seurat_obj <- "RNA/nanoCool_cellline.RNA.seurat.rds"

## colors
celltype_pal <- c("K562"="#A6CDE2", "HFF1"="#1E78B4")

qc_pal <-  c("Fail"="#f28300", "Pass"="#68bcd6")

group_pal <- c("#A6CDE2", "#1E78B4", "#F59899", "#E11E26")
names(group_pal) <- c("TGS_K562", "TGS_HFF1", "NGS_K562", "NGS_HFF1")

cnv_color <- c("Gain"="#E41A1C","Loss"="#377EB8", "Normal"="darkgrey")

## functions
average_in_window <- function(window, gr, v, method = "absolute", empty_v = NA) {
  
  if(missing(v)) v = rep(1, length(gr))
  if(is.null(v)) v = rep(1, length(gr))
  if(is.atomic(v) && is.vector(v)) v = cbind(v)
  
  v = as.matrix(v)
  if(is.character(v) && ncol(v) > 1) {
    stop("`v` can only be a character vector.")
  }
  
  if(length(empty_v) == 1) {
    empty_v = rep(empty_v, ncol(v))
  }
  
  u = matrix(rep(empty_v, each = length(window)), nrow = length(window), ncol = ncol(v))
  
  mtch = as.matrix(findOverlaps(window, gr))
  intersect = pintersect(window[mtch[,1]], gr[mtch[,2]])
  w = width(intersect)
  v = v[mtch[,2], , drop = FALSE]
  n = nrow(v)
  
  ind_list = split(seq_len(n), mtch[, 1])
  window_index = as.numeric(names(ind_list))
  window_w = width(window)
  
  if(is.character(v)) {
    for(i in seq_along(ind_list)) {
      ind = ind_list[[i]]
      if(is.function(method)) {
        u[window_index[i], ] = method(v[ind], w[ind], window_w[i])
      } else {
        tb = tapply(w[ind], v[ind], sum)
        u[window_index[i], ] = names(tb[which.max(tb)])
      }
    }
  } else {
    if(method == "w0") {
      gr2 = reduce(gr, min.gapwidth = 0)
      mtch2 = as.matrix(findOverlaps(window, gr2))
      intersect2 = pintersect(window[mtch2[, 1]], gr2[mtch2[, 2]])
      
      width_intersect = tapply(width(intersect2), mtch2[, 1], sum)
      ind = unique(mtch2[, 1])
      width_setdiff = width(window[ind]) - width_intersect
      
      w2 = width(window[ind])
      
      for(i in seq_along(ind_list)) {
        ind = ind_list[[i]]
        x = colSums(v[ind, , drop = FALSE]*w[ind])/sum(w[ind])
        u[window_index[i], ] = (x*width_intersect[i] + empty_v*width_setdiff[i])/w2[i]
      }
      
    } else if(method == "absolute") {
      for(i in seq_along(ind_list)) {
        u[window_index[i], ] = colMeans(v[ind_list[[i]], , drop = FALSE])
      }
      
    } else if(method == "weighted") {
      for(i in seq_along(ind_list)) {
        ind = ind_list[[i]]
        u[window_index[i], ] = colSums(v[ind, , drop = FALSE]*w[ind])/sum(w[ind])
      }
    } else {
      if(is.function(method)) {
        for(i in seq_along(ind_list)) {
          ind = ind_list[[i]]
          u[window_index[i], ] = method(v[ind], w[ind], window_w[i])
        }
      } else {
        stop("wrong method.")
      }
    }
  }
  
  return(u)
}

## ggplot theme
cnv_theme <- theme_bw()+theme(
  panel.grid.major = element_blank(), 
  panel.grid.minor = element_blank(),
  panel.background = element_blank(), 
  axis.line = element_line(colour = "black"),
  strip.background = element_blank(),
  strip.text.y = element_text(angle = 0, size = 10),
  panel.border = element_rect(size = 1),
  axis.text = element_text(size = 10),
  legend.position = "none",
  axis.title = element_blank(),
  axis.ticks.x = element_blank())
