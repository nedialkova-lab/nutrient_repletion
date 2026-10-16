# ==============================================================================
# 8-1_explore_tRNAseq.R
# mim-tRNAseq: anticodon-level tRNA abundance and aminoacylation (% CCA) along the
# glucose (0, 1, 4, 10 min) and nitrogen (0, 1, 2, 10 min) repletions.
#
# Manuscript: Fig. 4C; Extended Data Fig. 5C
# Inputs    : data/20251107_tRNA_sgRNA_RNA/{glucose,nitrogen}/Isodecoder_counts_raw.txt
#             data/20251107_tRNA_sgRNA_RNA/{glucose,nitrogen}/CCAprops.csv
# Outputs   : results/tRNA_{glucose,nitrogen}_anticodon_{abundance,CCAabundance,CCAperc}.tsv
#             figure_panels/dotplot_tRNA_{G,N}.pdf
#             figure_panels/lineplot_tRNA_G_Ser.pdf
#             figure_panels/lineplot_tRNA_N_Leu.pdf
# Run from the repository root (all paths are relative to it).
# ==============================================================================

library(edgeR)
library(ggplot2)
library(splines)
library(factoextra)
library(RColorBrewer)
library(eulerr)
library(gplots)
library(patchwork)
library(reshape2)
library(ggtree)
library(dplyr)
uniColors=brewer.pal(9,"YlOrRd")
uniColors=colorRampPalette(uniColors)(100)
biColors=rev(brewer.pal(11,"RdBu"))
biColors=colorRampPalette(biColors)(100)

#### Load data GLUCOSE ####
countsdf = read.csv("data/20251107_tRNA_sgRNA_RNA/glucose/Isodecoder_counts_raw.txt", sep="\t", row.names = 1)
ccalong = read.csv("data/20251107_tRNA_sgRNA_RNA/glucose/CCAprops.csv", sep=",")
# Explore abundance
countsdf = countsdf[,2:(ncol(countsdf)-2)]
colnames(countsdf) = sapply(colnames(countsdf), function(x) paste0(strsplit(strsplit(x,"\\.")[[1]][2],"_")[[1]][4:5], collapse="_"))
# Aminoacylation
ccadf = dcast(ccalong[ccalong$end=="CA",], formula = gene ~ sample, value.var = "percent")
rownames(ccadf) = ccadf$gene; ccadf = as.data.frame(ccadf[,-c(1)])
colnames(ccadf) = sapply(colnames(ccadf), function(x) paste0(strsplit(strsplit(x,"\\.")[[1]][1],"_")[[1]][4:5], collapse="_"))
# Keep only cyto tRNAs
countsdf = countsdf[!grepl("Escherichia_coli_", rownames(countsdf))&!grepl("mito_", rownames(countsdf)),]
ccadf = ccadf[!grepl("eColi", rownames(ccadf))&!grepl("mito", rownames(ccadf)),]
# map names from counts to cca dataframes
rownames(countsdf) = sapply(rownames(countsdf), function(x) paste0(sapply(strsplit(x,"-")[[1]][2:4], function(s) if(grepl("/",s)){sprintf("%s-%s",substr(s,1,1),"multi")}else{s}), collapse="-"))
commontrna = rownames(ccadf)[rownames(ccadf) %in% rownames(countsdf)]
commoncond = colnames(ccadf)[colnames(ccadf) %in% colnames(countsdf)]
# Calculate tRNA abundance (all)
anticodon = sapply(rownames(countsdf), function(x) paste0(strsplit(x,"-")[[1]][1:2], collapse="-")) 
antidf = aggregate(countsdf, by=list(anticodon), sum)
rownames(antidf) = antidf$Group.1; antidf=antidf[,-c(1)]
percdf = as.data.frame(t(apply(antidf,1,function(x) x*1e2/colSums(antidf))))
# Calculate aminoacylated tRNA abundances
chargedcounts = t(sapply(commontrna,function(x) as.numeric(countsdf[x,commoncond]*ccadf[x,commoncond]*1e-2)))
colnames(chargedcounts) = commoncond
anticodon = sapply(rownames(chargedcounts), function(x) paste0(strsplit(x,"-")[[1]][1:2], collapse="-")) 
chargedantidf = aggregate(chargedcounts, by=list(anticodon), sum)
rownames(chargedantidf) = chargedantidf$Group.1; chargedantidf=chargedantidf[,-c(1)]
chargedpercdf = t(apply(chargedantidf[,commoncond],1,function(x) x*1e2/colSums(antidf[,commoncond])))
# Anticodon CCA
anticca = chargedantidf[,commoncond]*1e2/antidf[,commoncond]
# Save
write.table(percdf,"results/tRNA_glucose_anticodon_abundance.tsv", sep = "\t")
write.table(chargedpercdf,"results/tRNA_glucose_anticodon_CCAabundance.tsv", sep = "\t")
write.table(anticca,"results/tRNA_glucose_anticodon_CCAperc.tsv", sep = "\t")

#### LIMMA DIFF EXP ####
timepointsG = sapply(colnames(antidf), function(x) as.numeric(substr(strsplit(x,"_")[[1]][1],2,4)))
# Remove genes with less than 1 cpm
d = DGEList(antidf)
dim(d) # number of tRNAs
# Define model to be fitted
X <- ns(timepointsG, knots=c(1.5, 4)) # represent time course as cubic spline
design = model.matrix(~ X)
# Voom transformation
y <- voom(d, design, plot = T)
# Fit model
fit <- lmFit(y, design)
fit.eb <- eBayes(fit)
res.abund.G <- topTable(fit.eb, coef=2:4, number = 100)
# Heatmap
colord = order(as.numeric(sapply(colnames(percdf), function(x) substr(strsplit(x, "_")[[1]][1], 2, 4))))
heatmap.2(t(percdf[,colord]), Rowv = F, Colv = T, col=uniColors, margins = c(7,7), trace = "none")

# Heatmaps
# raw
colord = order(as.numeric(sapply(colnames(anticca), function(x) substr(strsplit(x, "_")[[1]][1], 2, 4))))
heatmap.2(t(anticca[,colord]), Rowv = F, col=uniColors, margins = c(7,7),
          distfun = function(x) as.dist(1-cor(t(x), use = "pairwise", method = "pearson")),
          hclustfun = function(x) hclust(x, method="ward.D2"),trace = "none")
# delta
ctrl = rowMeans(anticca[,grepl("g0",colnames(anticca))])
deltadf = apply(anticca[,!grepl("g0",colnames(anticca))],2,function(x) x - ctrl)
deltadf = deltadf[rowSums(is.na(deltadf))==0,]
colord = order(as.numeric(sapply(colnames(deltadf), function(x) substr(strsplit(x, "_")[[1]][1], 2, 4))))
heatmap.2(t(deltadf[,colord]), Rowv = F, col=biColors, margins = c(7,7),
          distfun = function(x) as.dist(1-cor(t(x), use = "pairwise", method = "pearson")),
          hclustfun = function(x) hclust(x, method="ward.D2"),trace = "none")

## PLOT DOTPLOT
# Compute clustering of tRNAs
h_clust <- hclust(as.dist(1-cor(t(anticca), use = "pairwise", method = "pearson")))
ddgram_row <- as.dendrogram(h_clust)
ggtree_plot_col <- ggtree(ddgram_row, hang = -1) + layout_dendrogram()
deltadflong = melt(deltadf, varnames = c("tRNA","sample"), value.name = "AA")
perclong = melt(as.matrix(percdf), varnames = c("tRNA","sample"), value.name = "percent")
dotdf = merge(deltadflong,perclong, by=c("tRNA","sample"), all.x=T)
dotdf$sample = factor(dotdf$sample, levels=c(colnames(deltadf)[rev(colord)]))
dotdf$tRNA = factor(dotdf$tRNA, levels = h_clust$labels[h_clust$order])
dotplot <- ggplot(dotdf, aes(x=tRNA, y = sample, color = AA, size = percent)) + 
  geom_point() + 
  cowplot::theme_cowplot(font_size=12) + 
  theme(axis.line  = element_blank(), axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab('') + xlab('') +
  scale_color_gradientn(colours = biColors, limits = c(-40,40), oob = scales::squish, name = 'd%AA') +
  scale_size(limits = c(0,9)) +
  scale_y_discrete(position = "right")
ggtree_plot_col / plot_spacer() /dotplot + plot_layout(heights = c(2,-1.4,6))
ggsave("figure_panels/dotplot_tRNA_G.pdf", width=9, height=3)

# Serines
plotccadf = anticca[grepl("^Ser",rownames(anticca)),]
plotccadf$tRNA = rownames(plotccadf)
plotccadf = melt(plotccadf, value.name="AA")
plotabunddf = percdf[grepl("^Ser",rownames(anticca)),]
plotabunddf$tRNA = rownames(plotabunddf)
plotabunddf = melt(plotabunddf, value.name="percent")
plotmerged = merge(plotccadf,plotabunddf,by=c("variable","tRNA"))
plotmerged$time = sapply(as.character(plotmerged$variable), function(x) as.numeric(substr(strsplit(x,"_")[[1]][1],2,3)))
plotmerged$rep = sapply(as.character(plotmerged$variable), function(x) strsplit(x,"_")[[1]][2])
meanabund = plotmerged %>% group_by(tRNA) %>% summarise(avgabund=mean(percent))
plotmerged$tRNA = factor(plotmerged$tRNA, levels=unlist(meanabund[order(meanabund$avgabund, decreasing=T),"tRNA"]))
coef = 1e-1
ylim.prim = min(plotmerged$AA)-5
ggplot(plotmerged,aes(x=time)) +
  facet_grid(~tRNA) +
  geom_point(aes(y=AA), color="#d6604dff") + stat_summary(aes(y=AA), geom="line", fun = "mean", color="#d6604dff") +
  geom_point(aes(y=(percent/coef)+ylim.prim), color="#888888") + stat_summary(aes(y=(percent/coef)+ylim.prim), geom="line", fun = "mean", color="#888888") +
  scale_y_continuous(limits=c(ylim.prim,100), sec.axis = sec_axis(~(.-ylim.prim)*coef, name="percent")) +
  theme_test() + theme(axis.title.y.right = element_text(color = "#888888"), axis.text.y.right = element_text(colour = "#888888"))
ggsave("figure_panels/lineplot_tRNA_G_Ser.pdf", width=8, height=2)

#### Load data NITROGEN ####
countsdf = read.csv("data/20251107_tRNA_sgRNA_RNA/nitrogen/Isodecoder_counts_raw.txt", sep="\t", row.names = 1)
ccalong = read.csv("data/20251107_tRNA_sgRNA_RNA/nitrogen/CCAprops.csv", sep=",")
# Explore abundance
countsdf = countsdf[,2:(ncol(countsdf)-2)]
colnames(countsdf) = sapply(colnames(countsdf), function(x) paste0(strsplit(strsplit(x,"\\.")[[1]][2],"_")[[1]][4:5], collapse="_"))
# Aminoacylation
ccadf = dcast(ccalong[ccalong$end=="CA",], formula = gene ~ sample, value.var = "percent")
rownames(ccadf) = ccadf$gene; ccadf = as.data.frame(ccadf[,-c(1)])
colnames(ccadf) = sapply(colnames(ccadf), function(x) paste0(strsplit(strsplit(x,"\\.")[[1]][1],"_")[[1]][4:5], collapse="_"))
# Keep only cyto tRNAs
countsdf = countsdf[!grepl("Escherichia_coli_", rownames(countsdf))&!grepl("mito_", rownames(countsdf)),]
ccadf = ccadf[!grepl("eColi", rownames(ccadf))&!grepl("mito", rownames(ccadf)),]
# map names from counts to cca dataframes
rownames(countsdf) = sapply(rownames(countsdf), function(x) paste0(sapply(strsplit(x,"-")[[1]][2:4], function(s) if(grepl("/",s)){sprintf("%s-%s",substr(s,1,1),"multi")}else{s}), collapse="-"))
commontrna = rownames(ccadf)[rownames(ccadf) %in% rownames(countsdf)]
commoncond = colnames(ccadf)[colnames(ccadf) %in% colnames(countsdf)]
# Calculate tRNA abundance (all)
anticodon = sapply(rownames(countsdf), function(x) paste0(strsplit(x,"-")[[1]][1:2], collapse="-")) 
antidf = aggregate(countsdf, by=list(anticodon), sum)
rownames(antidf) = antidf$Group.1; antidf=antidf[,-c(1)]
percdf = as.data.frame(t(apply(antidf,1,function(x) x*1e2/colSums(antidf))))
# Calculate aminoacylated tRNA abundances
chargedcounts = t(sapply(commontrna,function(x) as.numeric(countsdf[x,commoncond]*ccadf[x,commoncond]*1e-2)))
colnames(chargedcounts) = commoncond
anticodon = sapply(rownames(chargedcounts), function(x) paste0(strsplit(x,"-")[[1]][1:2], collapse="-")) 
chargedantidf = aggregate(chargedcounts, by=list(anticodon), sum)
rownames(chargedantidf) = chargedantidf$Group.1; chargedantidf=chargedantidf[,-c(1)]
chargedpercdf = t(apply(chargedantidf[,commoncond],1,function(x) x*1e2/colSums(antidf[,commoncond])))
# Anticodon CCA
anticca = chargedantidf[,commoncond]*1e2/antidf[,commoncond]
# Save
write.table(percdf,"results/tRNA_nitrogen_anticodon_abundance.tsv", sep = "\t")
write.table(chargedpercdf,"results/tRNA_nitrogen_anticodon_CCAabundance.tsv", sep = "\t")
write.table(anticca,"results/tRNA_nitrogen_anticodon_CCAperc.tsv", sep = "\t")

#### LIMMA DIFF EXP ####
timepointsN = sapply(colnames(antidf), function(x) as.numeric(substr(strsplit(x,"_")[[1]][1],2,4)))
# Remove genes with less than 1 cpm
d = DGEList(antidf)
dim(d) # number of tRNAs
# Define model to be fitted
X <- ns(timepointsN, knots=c(1, 3)) # represent time course as cubic spline
design = model.matrix(~ X)
# Voom transformation
y <- voom(d, design, plot = T)
# Fit model
fit <- lmFit(y, design)
fit.eb <- eBayes(fit)
res.abund.N <- topTable(fit.eb, coef=2:4, number = 100)
# RPMs
percdf = t(apply(antidf,1,function(x) x*1e2/colSums(antidf)))
colord = order(as.numeric(sapply(colnames(percdf), function(x) substr(strsplit(x, "_")[[1]][1], 2, 4))))
heatmap.2(t(percdf[,colord]), Rowv = F, Colv = T, col=uniColors, margins = c(7,7), trace = "none")

# Heatmaps
# raw
colord = order(as.numeric(sapply(colnames(anticca), function(x) substr(strsplit(x, "_")[[1]][1], 2, 4))))
heatmap.2(t(anticca[,colord]), Rowv = F, col=uniColors, margins = c(7,7),
          distfun = function(x) as.dist(1-cor(t(x), use = "pairwise", method = "pearson")),
          hclustfun = function(x) hclust(x, method="ward.D2"),trace = "none")
# delta
ctrl = rowMeans(anticca[,grepl("n0",colnames(anticca)), drop=F])
deltadf = apply(anticca[,!grepl("n0",colnames(anticca))],2,function(x) x - ctrl)
deltadf = deltadf[rowSums(is.na(deltadf))==0,]
colord = order(as.numeric(sapply(colnames(deltadf), function(x) substr(strsplit(x, "_")[[1]][1], 2, 4))))
heatmap.2(t(deltadf[,colord]), Rowv = F, col=biColors, margins = c(7,7),
          distfun = function(x) as.dist(1-cor(t(x), use = "pairwise", method = "pearson")),
          hclustfun = function(x) hclust(x, method="ward.D2"),trace = "none")

## PLOT DOTPLOT
# Compute clustering of tRNAs
h_clust <- hclust(as.dist(1-cor(t(anticca), use = "pairwise", method = "pearson")))
ddgram_row <- as.dendrogram(h_clust)
ggtree_plot_col <- ggtree(ddgram_row, hang = -1) + layout_dendrogram()
deltadflong = melt(deltadf, varnames = c("tRNA","sample"), value.name = "AA")
perclong = melt(as.matrix(percdf), varnames = c("tRNA","sample"), value.name = "percent")
dotdf = merge(deltadflong,perclong, by=c("tRNA","sample"), all.x=T)
dotdf$sample = factor(dotdf$sample, levels=c(colnames(deltadf)[rev(colord)]))
dotdf$tRNA = factor(dotdf$tRNA, levels = h_clust$labels[h_clust$order])
dotplot <- ggplot(dotdf, aes(x=tRNA, y = sample, color = AA, size = percent)) + 
  geom_point() + 
  cowplot::theme_cowplot(font_size=12) + 
  theme(axis.line  = element_blank(), axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab('') + xlab('') +
  scale_color_gradientn(colours = biColors, limits = c(-40,40), oob = scales::squish, name = 'd%AA') +
  scale_size(limits = c(0,9)) +
  scale_y_discrete(position = "right")
ggtree_plot_col / plot_spacer() /dotplot + plot_layout(heights = c(2,-1.4,6))
ggsave("figure_panels/dotplot_tRNA_N.pdf", width=9, height=3)

# Leucines
plotccadf = anticca[grepl("^Leu",rownames(anticca)),]
plotccadf$tRNA = rownames(plotccadf)
plotccadf = melt(plotccadf, value.name="AA")
plotabunddf = as.data.frame(percdf[grepl("^Leu",rownames(anticca)),])
plotabunddf$tRNA = rownames(plotabunddf)
plotabunddf = melt(plotabunddf, value.name="percent")
plotmerged = merge(plotccadf,plotabunddf,by=c("variable","tRNA"))
plotmerged$time = sapply(as.character(plotmerged$variable), function(x) as.numeric(substr(strsplit(x,"_")[[1]][1],2,3)))
plotmerged$rep = sapply(as.character(plotmerged$variable), function(x) strsplit(x,"_")[[1]][2])
meanabund = plotmerged %>% group_by(tRNA) %>% summarise(avgabund=mean(percent))
plotmerged$tRNA = factor(plotmerged$tRNA, levels=unlist(meanabund[order(meanabund$avgabund, decreasing=T),"tRNA"]))
coef = 7e-2
ylim.prim = min(plotmerged$AA)-5
ggplot(plotmerged,aes(x=time)) +
  facet_grid(~tRNA) +
  geom_point(aes(y=AA), color="#d6604dff") + stat_summary(aes(y=AA), geom="line", fun = "mean", color="#d6604dff") +
  geom_point(aes(y=(percent/coef)+ylim.prim), color="#888888") + stat_summary(aes(y=(percent/coef)+ylim.prim), geom="line", fun = "mean", color="#888888") +
  scale_y_continuous(limits=c(ylim.prim,100), sec.axis = sec_axis(~(.-ylim.prim)*coef, name="percent")) +
  theme_test() + theme(axis.title.y.right = element_text(color = "#888888"), axis.text.y.right = element_text(colour = "#888888"))
ggsave("figure_panels/lineplot_tRNA_N_Leu.pdf", width=8, height=2)

