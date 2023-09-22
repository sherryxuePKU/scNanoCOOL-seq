#!/datc/houyu/software/miniconda3/envs/R4.0/bin/Rscript
arg <- commandArgs(T)
library(data.table)

input <- arg[1]
output <- arg[2]

filename <- rev(strsplit(input, "/")[[1]])[1]
sample <- strsplit(filename,".WCG")[[1]][1]
region <- sub("WCG_", "",strsplit(filename,"[.]")[[1]][2])

df <- fread(file = input, header = F, stringsAsFactors = F)
df <- df[,-1]
df <- colMeans(na.rm = T, df)
df <- data.frame(data=df)
df$Coord <- 1:nrow(df)
df$Sample <- sample

write.table(df, file = output, quote = F, sep = "\t", col.names = T, row.names = F)


# pdf(paste0(outdir, ".meth_",region,".pdf"))
# #ggplot(data = merge_file, aes(x=Coord, y=data, color=Sample))+geom_line()+scale_x_continuous(breaks = c(25, 125), labels = c("TSS", "TES"))+scale_color_manual(values = cols)+theme_bw()+theme(panel.grid.minor.x = element_blank(), axis.title = element_blank(), legend.position = "none")
# if(region=="genebody"){
#   ggplot(data = merge_file, aes(x=Coord, y=data, color=Sample))+geom_line()+scale_x_continuous(breaks = c(0,25, 125,150), labels = c("-15kb","TSS", "TES", "+15kb"))+ labs(y="DNA methylation level (%)")+theme_bw()+theme(panel.grid.minor.x = element_blank(), axis.title = element_blank(), legend.position = "none")
# } else if(region=="CGI"){
#   ggplot(data = merge_file, aes(x=Coord, y=data, color=Sample))+geom_line()+scale_x_continuous(breaks = c(0,25, 125,150), labels = c("-15kb","Start", "End", "+15kb"))+ labs(y="DNA methylation level (%)")+theme_bw()+theme(panel.grid.minor.x = element_blank(), axis.title = element_blank(), legend.position = "none")
# }
# dev.off()

