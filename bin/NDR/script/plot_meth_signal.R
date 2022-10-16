## library packages
library(dplyr)
library(ggplot2)
library(scales)

## input
setwd("/path/to/input")

wcg_input <- "nanoCool_5aza.0.CHlt80per.3Depth_5GCH.proximal_NDR.WCG_methlevel.tsv"
gch_input <- "nanoCool_5aza.0.CHlt80per.3Depth_5GCH.proximal_NDR.GCH_methlevel.tsv"

wcg_pos <- read.table(wcg_input, header = F, stringsAsFactors = F)
gch_pos <- read.table(gch_input, header = F, stringsAsFactors = F)
wcg_pos_meth <- colMeans(wcg_pos[,-c(1,2)], na.rm = T)
gch_pos_meth <- colMeans(gch_pos[,-c(1,2)], na.rm = T)

ref_1 <- "nanoCool_5aza.0.CHlt80per.3Depth_5GCH.proximal_NDR.CTCF_signal.rds"
ref_2 <- "nanoCool_5aza.0.CHlt80per.3Depth_5GCH.proximal_NDR.POLR2X_signal.rds"

ref_1 <- readRDS(ref_1)
ref_1 <- rescale(colMeans(ref_1), to=c(0.5,1))
ref_2 <- readRDS(ref_2)
ref_2 <- rescale(colMeans(ref_2), to=c(0.5,1))

ndr_meth <- data.frame(
  WCG_meth=wcg_pos_meth,
  GCH_meth=gch_pos_meth,
  CTC_signal=ref_1,
  POLR2A_signal=ref_2,
  Pos=1:length(gch_pos_meth)
)

## smooth the data with spline
spline.d_1 <- as.data.frame(spline(ndr_meth$Pos, ndr_meth$WCG_meth))
spline.d_2 <- as.data.frame(spline(ndr_meth$Pos, ndr_meth$GCH_meth))
spline.d_3 <- as.data.frame(spline(ndr_meth$Pos, ndr_meth$CTC_signal))
spline.d_4 <- as.data.frame(spline(ndr_meth$Pos, ndr_meth$POLR2A_signal))

ggplot()+
  geom_line(data = spline.d_1, aes(x = x, y = y), 
            color="#E41A1C", size=1)+
  geom_line(data = spline.d_2, aes(x = x, y = y),
            color="#377EB8", size=1)+
  geom_line(data = spline.d_3, aes(x = x, y = y),
            color="#4DAF4A", size=1)+
  geom_line(data = spline.d_4, aes(x = x, y = y),
            color="#984EA3", size=1)+
  scale_x_continuous(
    breaks = c(1,10,20,30,40),
    labels = c("-1000","-500", "0","500", "1000"))+
  coord_fixed(ratio = 35)+
  scale_y_continuous(sec.axis = sec_axis(~. /1, name = "Normalized signal intensity"))+
  theme_classic(base_size = 20)+
  theme(legend.position = "bottom")+
  labs(y="Methylation level",x="Distance from NDR (bp)")
  