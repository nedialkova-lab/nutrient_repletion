# ==============================================================================
# 3-1_explore_ctrl_lipms.R
# Exploratory (not shown in the manuscript): comparison of the time-0 total proteome
# (Orbi4696) with the trypsin-only (TrP) controls of the LiP-MS experiment (Venn
# diagram of detected proteins, PCA). Plots are only displayed, nothing is saved.
#
# Manuscript: exploratory analysis, no figure panel
# Inputs    : data/proteomics_test/Orbi4696_ProteinGroups.tsv
#             data/20240507_lipms_upshifts_5199/Orbi5199_LiP+Ctrl_HybridLib_Report.tsv
#             data/20240507_lipms_upshifts_5199/240507_lipms.xlsx
# Outputs   : none (interactive plots only)
# Run from the repository root (all paths are relative to it).
# ==============================================================================

library(ggplot2)
library(gplots)
library(RColorBrewer)
library(factoextra)
library(reshape)
library(clusterProfiler)
library(msigdbr)
library(limma)
library(openxlsx)
library(eulerr)
#https://www.biostat.jhsph.edu/~kkammers/software/CVproteomics/R_guide.html

# Load proteomics data
protdata = read.csv("data/proteomics_test/Orbi4696_ProteinGroups.tsv", sep="\t", row.names = "PG.ProteinGroups")
cleandf = as.matrix(protdata[,grepl("Quantity",colnames(protdata))])
# Remove NAs and low confidence proteins
cleandf = cleandf[(rowSums(is.na(cleandf))==0)&(!apply(protdata[,grepl("IsSingleHit",colnames(protdata))],1,any)),]
samplenames = c("C1","C2","C3","G1","G2","G3","N1","N2","N3","P1","P2","P3")
colnames(cleandf) = samplenames

# Load lipms data
lipdata = read.csv("data/20240507_lipms_upshifts_5199/Orbi5199_LiP+Ctrl_HybridLib_Report.tsv", sep="\t")
protdf = cast(lipdata, formula = PG.ProteinGroups ~ R.FileName, value = "PG.Quantity", fun.aggregate = mean)
rownames(protdf) = protdf$PG.ProteinGroups; protdf = as.data.frame(protdf[,-c(1)])
sampletable = read.xlsx("data/20240507_lipms_upshifts_5199/240507_lipms.xlsx", rowNames = F)
# Keep only TrP
samplenums = as.numeric(gsub('.*Sample([0-9]*).*','\\1',colnames(protdf)))
islip = grepl("LiP",colnames(protdf))
protdf = protdf[,!islip]
samplenums = samplenums[!islip]
colnames(protdf) = sapply(samplenums, function(x) sampletable[sampletable$LiP_MSnumber==x,"Sample.name"])
# Subset to controls (time 0)
protdf = protdf[,substr(colnames(protdf),2,2)==0]
protdf = as.data.frame(protdf[(rowSums(is.na(protdf))==0),])

# Merge
alldf = cbind(cleandf[rownames(cleandf) %in% rownames(protdf),],protdf[rownames(cleandf)[rownames(cleandf) %in% rownames(protdf)],sort(colnames(protdf))])
cordf = (cor(alldf, use = "pairwise", method = "spearman")^2)[1:12,13:21]
Colors=rev(brewer.pal(9,"YlOrRd"))
Colors=colorRampPalette(Colors)(100)
heatmap.2(cordf, trace="none", col=Colors, margins = c(6,6), Rowv = F, Colv = F, breaks=seq(0.5,0.55,0.0005))

# Venn diagram
overlap = list(Proteomics = unique(unlist(sapply(rownames(cleandf), function(x) unlist(strsplit(x,";"))))),
               LiPMSctrl = unique(unlist(sapply(rownames(protdf), function(x) unlist(strsplit(x,";"))))))
fit <- euler(overlap)
plot(fit, fills = c("#f4d7d7","#d7e6f4"), edges=T, legend=F, labels=T, quantities=T,main="Proteins")

# Principal component analysis
res.pca = prcomp(t(cleandf), scale = T)
fviz_eig(res.pca)
fviz_pca_ind(res.pca, repel = TRUE, pointsize=3, invisible="quali",
             habillage=substr(samplenames,1,1)) +
  scale_color_brewer(palette="Dark2") +
  theme_minimal()

res.pca = prcomp(t(protdf), center = T, scale. = T)
fviz_eig(res.pca)
fviz_pca_ind(res.pca, repel = TRUE, pointsize=3, invisible="quali",
             habillage=substr(colnames(protdf),1,1)) +
  scale_color_brewer(palette="Dark2") +
  theme_minimal()

res.pca = prcomp(t(alldf), center = T, scale. = T)
fviz_eig(res.pca)
fviz_pca_ind(res.pca, repel = TRUE, pointsize=3, invisible="quali",
             habillage=substr(colnames(alldf),1,1)) +
  scale_color_brewer(palette="Dark2") +
  theme_minimal()
