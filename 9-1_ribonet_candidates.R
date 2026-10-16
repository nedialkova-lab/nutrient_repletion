# ==============================================================================
# 9-1_ribonet_candidates.R
# Exploratory (not shown in the manuscript): RxT, RPF and RNA time courses along the
# glucose repletion for selected candidate genes (e.g. YNR028W/CPR8, YPL207W/TYW1).
# Plots are only displayed, nothing is saved.
#
# Manuscript: exploratory analysis, no figure panel
# Requires  : outputs of 6-1, 7-1, 7-3
# Inputs    : results/TPM_RPF.tsv
#             results/TPM_RNA.tsv
#             results/RxT.tsv
# Outputs   : none (interactive plots only)
# Run from the repository root (all paths are relative to it).
# ==============================================================================

library(ggplot2)
library(RColorBrewer)
library(UniProt.ws)

# Data
rpfdf = read.csv("results/TPM_RPF.tsv",sep="\t")
rnadf = read.csv("results/TPM_RNA.tsv",sep="\t")
rxtdf = read.csv("results/RxT.tsv",sep="\t")

# Subset data to G
rpfdf = rpfdf[grepl("^g",colnames(rpfdf))]
rnadf = rnadf[grepl("^g",colnames(rnadf))]
rxtdf = rxtdf[grepl("^g",colnames(rxtdf))]

# Plot genes
gene = "YNR028W"
genedf = as.data.frame(t(cbind(rxtdf[gene,],rpfdf[gene,],rnadf[gene,])))
genedf$sample = c(colnames(rxtdf),colnames(rpfdf),colnames(rnadf))
genedf$metric = c(rep("RxT",ncol(rxtdf)),rep("RPF",ncol(rpfdf)),rep("RNA",ncol(rnadf)))
genedf$time = sapply(genedf$sample, function(x) as.numeric(substr(strsplit(x,"_")[[1]][1],2,4)))
ggplot(genedf[genedf$metric=="RxT",], aes(x=time, y=get(gene))) +
  geom_point() +
  stat_summary(geom="line") +
  ylab(gene) +
  theme_bw()
ggplot(genedf, aes(x=time, y=get(gene), color=metric)) +
  geom_point() +
  stat_summary(geom="line") +
  ylab(gene) +
  theme_bw()

gene = "YPL207W"
genedf = as.data.frame(t(cbind(rxtdf[gene,],rpfdf[gene,],rnadf[gene,])))
genedf$sample = c(colnames(rxtdf),colnames(rpfdf),colnames(rnadf))
genedf$metric = c(rep("RxT",ncol(rxtdf)),rep("RPF",ncol(rpfdf)),rep("RNA",ncol(rnadf)))
genedf$time = sapply(genedf$sample, function(x) as.numeric(substr(strsplit(x,"_")[[1]][1],2,4)))
ggplot(genedf[genedf$metric=="RxT",], aes(x=time, y=get(gene))) +
  geom_point() +
  stat_summary(geom="line") +
  ylab(gene) +
  theme_bw()
ggplot(genedf, aes(x=time, y=get(gene), color=metric)) +
  geom_point() +
  stat_summary(geom="line") +
  ylab(gene) +
  theme_bw()
