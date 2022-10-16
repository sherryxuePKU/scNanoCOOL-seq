library(EnrichedHeatmap)
library(data.table)
library(GenomicRanges)

set.seed(123)

cmd <- commandArgs(T)
ref <- cmd[1]
query <- cmd[2]
out <- cmd[3]

# dnase_df <- fread("Data/5aza/NDR/ENCFF185XRG.bed.gz", header = F, stringsAsFactors = F)
## prepare ref input
ref_df <- fread(ref, header = F, stringsAsFactors = F)
ref_df <- ref_df[,c(1:3, 7)]
colnames(ref_df) <- c("chr", "start", "end", "signal")
ref.gr <- makeGRangesFromDataFrame(
  ref_df,
   keep.extra.columns=T,
   ignore.strand=FALSE,
   seqinfo=NULL
)

## prepare query input
query_df <- read.table(query, header = F, stringsAsFactors = F)
query_df <- query_df[,1:3]
colnames(query_df) <- c("chr", "start", "end")

query_df$mid <- (query_df$start+query_df$end)/2
query.gr <- makeGRangesFromDataFrame(
  query_df,
   start.field = "mid",
   end.field = "mid",
   keep.extra.columns=FALSE,
   ignore.strand=FALSE,
   seqinfo=NULL
)

## calc matrix
mat1 <- normalizeToMatrix(
  ref.gr, query.gr, value_column = "signal", 
  extend = 1000, mean_mode = "coverage", w = 50
)
# EnrichedHeatmap(mat1, col = c("white", "red"), name = "DNase-seq")
saveRDS(mat1, file = out)

