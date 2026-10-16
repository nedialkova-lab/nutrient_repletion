# ==============================================================================
# 2-1_targeted_metabolomics_4973.R
# Targeted metabolomics of amino acids (light/heavy ratios against 13C15N internal
# standards) along the 30-min glucose, nitrogen and phosphate repletions.
#
# Manuscript: Extended Data Fig. 5G (heatmaps), 5H (serine and leucine)
# Inputs    : data/metabolomics_4973/Orbi4973_AA_Replicates_Custome_quant.csv
#             data/metabolomics_4973/Orbi4973_231218_metabolomics_81samples.xlsx
# Outputs   : results/metabolomics_targeted (no extension, tab-separated)
#             source_data/targeted_metab_avg_{Glc,Nit,Pho}.tsv
#             figure_panels/heatmap_targeted_metab_{Glc,Nit,Pho}.pdf
#             figure_panels/targeted_{Ser,Leu}.pdf
# Run from the repository root (all paths are relative to it).
# ==============================================================================

library(ggplot2)
library(ggrepel)
library(ggpubr)
library(gplots)
library(RColorBrewer)
library(reshape)
library(openxlsx)
library(factoextra)
library(limma)
library(splines)
library(eulerr)
library(MicrobiomeProfiler)
library(clusterProfiler)

# Load data
targ_F5 = read.csv("data/metabolomics_4973/Orbi4973_AA_Replicates_Custome_quant.csv", sep=",", na.strings = c("#N/A"))
# Remove 1 contaminated sample during boiling water bath
targ_F5 = targ_F5[!(targ_F5$Replicate=="45"),]
sampletable = read.xlsx("data/metabolomics_4973/Orbi4973_231218_metabolomics_81samples.xlsx", rowNames = F)

# TICs, check blanks
tics_F5 = unique(targ_F5[,c("Replicate","Total.Ion.Current.Area")])
tics_F5[1,1] = "01" 
rownames(tics_F5) = sapply(tics_F5$Replicate, function(x) sampletable[sampletable$X6==as.numeric(x),"Sample.name"])
tics_F5$sample = factor(rownames(tics_F5), levels=rownames(tics_F5)[order(tics_F5$Total.Ion.Current.Area)])
ggplot(tics_F5, aes(x=Total.Ion.Current.Area, y=sample)) + geom_col() + theme_classic()

### TARGETED F5 ###
cleandf = as.data.frame(cast(targ_F5, Molecule.Name ~ Replicate, value="RatioLightToHeavy", mean))
rownames(cleandf) = cleandf$Molecule.Name; cleandf = cleandf[,-1]
colnames(cleandf)[1] <- "01" 
samplenames = sapply(colnames(cleandf), function(x) sampletable[sampletable$X6==as.numeric(x),"Sample.name"])
colnames(cleandf) = samplenames

# These samples only contain labeled amino acids, remove species not contained in standard used
standards = c("Alanine", "Lysine","Histidine","Arginine","Tyrosine","Phenylalanine","Methionine","Glutamic Acid",
              "Aspartic Acid","Leucine","Isoleucine","Valine","Threonine","Cysteine","Proline","Serine","Glycine")
cleandf = cleandf[rownames(cleandf) %in% standards,]

# Plot few aa
plotdf = as.data.frame(t(cleandf))
plotdf$cond = substr(rownames(plotdf),1,1)
plotdf$time = as.numeric(substr(rownames(plotdf),2,3))
p = ggplot(plotdf, aes(time,log2(Serine))) + facet_wrap(~cond) + geom_point(alpha=0.6) + stat_summary(geom="line", fun = "mean") + theme_test()
ggsave("figure_panels/targeted_Ser.pdf", p, width = 9, height = 2.5)
p = ggplot(plotdf, aes(time,log2(Leucine))) + facet_wrap(~cond) + geom_point(alpha=0.6) + stat_summary(geom="line", fun = "mean") + theme_test()
ggsave("figure_panels/targeted_Leu.pdf", p, width = 9, height = 2.5)
write.table(plotdf, "results/metabolomics_targeted.tsv", sep = "\t")

# Heatmap of log10 abundances of glucose deprived samples
Colors=brewer.pal(9,"YlOrRd")
Colors=colorRampPalette(Colors)(100)
heatmap.2(log2(t(cleandf)), trace="none", col=Colors, margins = c(8,7), Rowv = F, Colv = F)

# Differential expression of ctrls vs conditions
time = sapply(colnames(cleandf), function(x) sampletable[sampletable$Sample.name==x,"Time.point"])
avgCtrl_G = rowMeans(cleandf[,grepl("G0", colnames(cleandf))], na.rm = T)
avgCtrl_N = rowMeans(cleandf[,grepl("N0", colnames(cleandf))], na.rm = T)
avgCtrl_P = rowMeans(cleandf[,grepl("P0", colnames(cleandf))], na.rm = T)
fc_vs_ctrl_G = log2(cleandf[,grepl("G", colnames(cleandf))]/ avgCtrl_G)[,order(time[grepl("G",colnames(cleandf))])]
fc_vs_ctrl_N = log2(cleandf[,grepl("N", colnames(cleandf))]/ avgCtrl_N)[,order(time[grepl("N",colnames(cleandf))])]
fc_vs_ctrl_P = log2(cleandf[,grepl("P", colnames(cleandf))]/ avgCtrl_P)[,order(time[grepl("P",colnames(cleandf))])]
# Heatmap
Colors=rev(brewer.pal(11,"RdBu"))
Colors=colorRampPalette(Colors)(100)
heatmap.2(t(fc_vs_ctrl_G), trace="none", Rowv=F, Colv = F, col=Colors, margins = c(8,4),breaks=seq(-4,4,0.08),symkey=F)
heatmap.2(t(fc_vs_ctrl_N), trace="none", Rowv=F, Colv = F, col=Colors, margins = c(8,4),breaks=seq(-4,4,0.08),symkey=F)
heatmap.2(t(fc_vs_ctrl_P), trace="none", Rowv=F, Colv = F, col=Colors, margins = c(8,4),breaks=seq(-4,4,0.08),symkey=F)

# Plot averages
longdf = melt(as.matrix(cleandf))
longdf$condition = sapply(as.character(longdf$X2),function(x) strsplit(x," ")[[1]][1])
meandf = cast(longdf, X1 ~ condition, mean)
rownames(meandf)= meandf$X1
meandf= meandf[, -c(1)]
time = as.numeric(substr(colnames(meandf),2,4))
fc_vs_ctrl_G = log2(meandf[,grepl("G", colnames(meandf))]/ meandf[,"G0"])
fc_vs_ctrl_G = fc_vs_ctrl_G[order(rowMeans(fc_vs_ctrl_G)),order(time[grepl("G",colnames(meandf))])[2:9]]
fc_vs_ctrl_N = log2(meandf[,grepl("N", colnames(meandf))]/ meandf[,"N0"])
fc_vs_ctrl_N = fc_vs_ctrl_N[order(rowMeans(fc_vs_ctrl_N)),order(time[grepl("N",colnames(meandf))])[2:9]]
fc_vs_ctrl_P = log2(meandf[,grepl("P", colnames(meandf))]/ meandf[,"P0"])
fc_vs_ctrl_P = fc_vs_ctrl_P[order(rowMeans(fc_vs_ctrl_P)),order(time[grepl("P",colnames(meandf))])[2:9]]

pdf("figure_panels/heatmap_targeted_metab_Glc.pdf", width = 6, height = 4)
heatmap.2(t(fc_vs_ctrl_G), trace="none", Rowv=F, Colv = F, col=Colors, margins = c(8,4),breaks=seq(-3,3,0.06),symkey=F)
dev.off()
write.table(t(fc_vs_ctrl_G), "source_data/targeted_metab_avg_Glc.tsv", sep = "\t")

pdf("figure_panels/heatmap_targeted_metab_Nit.pdf", width = 6, height = 4)
heatmap.2(t(fc_vs_ctrl_N), trace="none", Rowv=F, Colv = F, col=Colors, margins = c(8,4),breaks=seq(-3,3,0.06),symkey=F)
dev.off()
write.table(t(fc_vs_ctrl_N), "source_data/targeted_metab_avg_Nit.tsv", sep = "\t")

pdf("figure_panels/heatmap_targeted_metab_Pho.pdf", width = 6, height = 4)
heatmap.2(t(fc_vs_ctrl_P), trace="none", Rowv=F, Colv = F, col=Colors, margins = c(8,4),breaks=seq(-3,3,0.06),symkey=F)
dev.off()
write.table(t(fc_vs_ctrl_P), "source_data/targeted_metab_avg_Pho.tsv", sep = "\t")
