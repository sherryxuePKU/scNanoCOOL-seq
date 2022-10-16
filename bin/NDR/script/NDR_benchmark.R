## library packages
library(dplyr)
library(ggplot2)

Intervene_pairwise_frac_matrix <- read.delim(
  "/path/to/Intervene_pairwise_frac_matrix.txt", row.names=1)

rownames(Intervene_pairwise_frac_matrix) <- c(
  "ENCODE_cCRE.K562","ENCODE.ATAC-seq", "ENCODE.DNase-seq",
  "scNanoCool-seq", "scCool-seq"
)
colnames(Intervene_pairwise_frac_matrix) <- rownames(Intervene_pairwise_frac_matrix)
Intervene_pairwise_frac_matrix.cp <- Intervene_pairwise_frac_matrix

## precision
Intervene_pairwise_frac_matrix$ID <- rownames(Intervene_pairwise_frac_matrix)
Intervene_pairwise_frac_matrix <- reshape2::melt(Intervene_pairwise_frac_matrix, id.var="ID")
colnames(Intervene_pairwise_frac_matrix) <- c("query", "reference", "precision")

p1 <- Intervene_pairwise_frac_matrix %>% 
  filter(reference %in% c("ENCODE.ATAC-seq", "ENCODE.DNase-seq","ENCODE_cCRE.K562") &
           query %in% c("scNanoCool-seq", "scCool-seq")) %>%
  ggplot(aes(x=query, y=precision, fill=reference))+
  geom_bar(stat = "identity", position = "dodge", width = 0.8)+
  geom_text(aes(label=round(precision, 2)), 
            position = position_dodge(0.75),vjust=0)+
  scale_fill_manual(values = c(
    "ENCODE.ATAC-seq"="#4DAF4A",
    "ENCODE.DNase-seq"="#984EA3",
    "ENCODE_cCRE.K562"="#FF7F00"))+
  theme_classic(base_size = 10)

## recall
Intervene_pairwise_frac_matrix.cp <- as.data.frame(t(Intervene_pairwise_frac_matrix.cp))

Intervene_pairwise_frac_matrix.cp$ID <- rownames(Intervene_pairwise_frac_matrix.cp)
Intervene_pairwise_frac_matrix.cp <- reshape2::melt(Intervene_pairwise_frac_matrix.cp, id.var="ID")
colnames(Intervene_pairwise_frac_matrix.cp) <- c("query", "reference", "recall")

p2 <- Intervene_pairwise_frac_matrix.cp %>% 
  filter(reference %in% c("ENCODE.ATAC-seq", "ENCODE.DNase-seq","ENCODE_cCRE.K562") &
           query %in% c("scNanoCool-seq", "scCool-seq")) %>%
  ggplot(aes(x=query, y=recall, fill=reference))+
  geom_bar(stat = "identity", position = "dodge", width = 0.8)+
  geom_text(aes(label=round(recall, 2)), 
            position = position_dodge(0.75),vjust=0)+
  scale_fill_manual(values = c(
    "ENCODE.ATAC-seq"="#4DAF4A",
    "ENCODE.DNase-seq"="#984EA3",
    "ENCODE_cCRE.K562"="#FF7F00"))+
  theme_classic(base_size = 10)

pr_df <- inner_join(
  Intervene_pairwise_frac_matrix,
  Intervene_pairwise_frac_matrix.cp,
  by=c("query"="query",
       "reference"="reference")
)

pr_df$F1 <- 2*pr_df$precision*pr_df$recall/(pr_df$precision+pr_df$recall)

p3 <- pr_df %>% 
  filter(reference %in% c("ENCODE.ATAC-seq", "ENCODE.DNase-seq","ENCODE_cCRE.K562") &
           query %in% c("scNanoCool-seq", "scCool-seq")) %>%
  ggplot(aes(x=query, y=F1, fill=reference))+
  geom_bar(stat = "identity", position = "dodge", width = 0.8)+
  geom_text(aes(label=round(recall, 2)), 
            position = position_dodge(0.75),vjust=0)+
  scale_fill_manual(values = c(
    "ENCODE.ATAC-seq"="#4DAF4A",
    "ENCODE.DNase-seq"="#984EA3",
    "ENCODE_cCRE.K562"="#FF7F00"))+
  theme_classic(base_size = 10)

p1+p2+p3
ggsave("Plot/Cellline.NDR_ge240bp.PR.pdf")
