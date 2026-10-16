# ==============================================================================
# 3-2_explore_LiPMS_5199.R
# Exploratory (not shown in the manuscript): quality control of the LiP-MS data
# (PCA per repletion, correlation between replicates, batch effects) before the
# batch correction of 3-3. Plots are only displayed, nothing is saved.
#
# Manuscript: exploratory analysis, no figure panel
# Inputs    : data/20240507_lipms_upshifts_5199/Orbi5199_LiP+Ctrl_HybridLib_Report.tsv
#             data/20240507_lipms_upshifts_5199/240507_lipms.xlsx
# Outputs   : none (interactive plots only)
# Run from the repository root (all paths are relative to it).
# ==============================================================================

library(ggplot2)
library(RColorBrewer)
library(factoextra)
library(reshape)
library(limma)
library(openxlsx)
library(splines)
#https://www.biostat.jhsph.edu/~kkammers/software/CVproteomics/R_guide.html

# Load data
lipdata = read.csv("data/20240507_lipms_upshifts_5199/Orbi5199_LiP+Ctrl_HybridLib_Report.tsv", sep="\t")
protdf = cast(lipdata, formula = PG.ProteinGroups ~ R.FileName, value = "PG.Quantity", fun.aggregate = mean)
rownames(protdf) = protdf$PG.ProteinGroups; protdf = as.data.frame(protdf[,-c(1)])
pepdf = cast(lipdata, formula = PEP.StrippedSequence ~ R.FileName, value = "PEP.Quantity", fun.aggregate = mean)
rownames(pepdf) = pepdf$PG.ProteinGroups; pepdf = as.data.frame(pepdf[,-c(1)])
sampletable = read.xlsx("data/20240507_lipms_upshifts_5199/240507_lipms.xlsx", rowNames = F)

# Remove NAs and low confidence proteins
protdf = as.data.frame(protdf[(rowSums(is.na(protdf))==0),])
pepdf = as.data.frame(pepdf[(rowSums(is.na(pepdf))==0),])
# Rename columns
samplenums = as.numeric(gsub('.*Sample([0-9]*).*','\\1',colnames(protdf)))
ctrllip = rep("TrP",144)
islip = grepl("LiP",colnames(protdf))
ctrllip[islip] = "LiP"
samplenames = sprintf("%s %s", sapply(samplenums, function(x) sampletable[sampletable$LiP_MSnumber==x,"Sample.name"]), ctrllip)
colnames(protdf) = samplenames
colnames(pepdf) = samplenames

# Principal component analysis
res.pca = prcomp(t(protdf), center = T, scale. = T)
fviz_eig(res.pca)
fviz_pca_ind(res.pca, repel = TRUE, pointsize=3, invisible="quali",
             habillage=sapply(samplenums, function(x) sampletable[sampletable$LiP_MSnumber==x,"Condition"])) +
  scale_color_brewer(palette="Dark2") +
  theme_minimal()

res.pca = prcomp(t(pepdf), center = T, scale. = T)
fviz_eig(res.pca)
fviz_pca_ind(res.pca, repel = TRUE, pointsize=3, invisible="quali",
             habillage=sapply(samplenums, function(x) sampletable[sampletable$LiP_MSnumber==x,"Condition"])) +
  scale_color_brewer(palette="Dark2") +
  theme_minimal()

# Glucose
res.pca = prcomp(t(protdf[,grepl("G",colnames(protdf))]), center = T, scale. = T)
fviz_eig(res.pca)
fviz_pca_ind(res.pca, repel = TRUE, pointsize=3, invisible="quali",
             habillage=sapply(samplenums[grepl("G",colnames(protdf))], function(x) sampletable[sampletable$LiP_MSnumber==x,"Time.point"])) +
  scale_color_manual(values = colorRampPalette(brewer.pal(9, "Blues"))(12)[4:12]) +
  theme_minimal()

res.pca = prcomp(t(pepdf[,grepl("G",colnames(pepdf))]), center = T, scale. = T)
fviz_eig(res.pca)
fviz_pca_ind(res.pca, repel = TRUE, pointsize=3, invisible="quali",
             habillage=sapply(samplenums[grepl("G",colnames(pepdf))], function(x) sampletable[sampletable$LiP_MSnumber==x,"Time.point"])) +
  scale_color_manual(values = colorRampPalette(brewer.pal(9, "Blues"))(12)[4:12]) +
  theme_minimal()

# Nitrogen
res.pca = prcomp(t(protdf[,grepl("N",colnames(protdf))]), center = T, scale. = T)
fviz_eig(res.pca)
fviz_pca_ind(res.pca, repel = TRUE, pointsize=3, invisible="quali",
             habillage=sapply(samplenums[grepl("N",colnames(protdf))], function(x) sampletable[sampletable$LiP_MSnumber==x,"Time.point"])) +
  scale_color_manual(values = colorRampPalette(brewer.pal(9, "Blues"))(12)[4:12]) +
  theme_minimal()

res.pca = prcomp(t(pepdf[,grepl("N",colnames(pepdf))]), center = T, scale. = T)
fviz_eig(res.pca)
fviz_pca_ind(res.pca, repel = TRUE, pointsize=3, invisible="quali",
             habillage=sapply(samplenums[grepl("N",colnames(pepdf))], function(x) sampletable[sampletable$LiP_MSnumber==x,"Time.point"])) +
  scale_color_manual(values = colorRampPalette(brewer.pal(9, "Blues"))(12)[4:12]) +
  theme_minimal()

# Phosphate
res.pca = prcomp(t(protdf[,grepl("P",substr(colnames(protdf),1,1))]), center = T, scale. = T)
fviz_eig(res.pca)
fviz_pca_ind(res.pca, repel = TRUE, pointsize=3, invisible="quali",
             habillage=sapply(samplenums[grepl("P",substr(colnames(protdf),1,1))], function(x) sampletable[sampletable$LiP_MSnumber==x,"Time.point"])) +
  scale_color_manual(values = colorRampPalette(brewer.pal(9, "Blues"))(12)[4:12]) +
  theme_minimal()

res.pca = prcomp(t(pepdf[,grepl("P",substr(colnames(protdf),1,1))]), center = T, scale. = T)
fviz_eig(res.pca)
fviz_pca_ind(res.pca, repel = TRUE, pointsize=3, invisible="quali",
             habillage=sapply(samplenums[grepl("P",substr(colnames(protdf),1,1))], function(x) sampletable[sampletable$LiP_MSnumber==x,"Time.point"])) +
  scale_color_manual(values = colorRampPalette(brewer.pal(9, "Blues"))(12)[4:12]) +
  theme_minimal()


# Corplot between replicates
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}
orderdf = protdf[,sort(colnames(protdf)[grepl("TrP", colnames(protdf))])]
upper_tri <- melt(get_upper_tri(round(cor(orderdf, method = "spearman"),2)), na.rm = TRUE)
ggplot(data = upper_tri, aes(X2, X1, fill = value))+
  geom_tile(color = "white") +
  scale_fill_gradientn(colors = brewer.pal(9,"YlOrRd")) +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 9, hjust = 1),
        axis.title.x=element_blank(), axis.title.y=element_blank())+
  coord_fixed()

orderdf = pepdf[,sort(colnames(pepdf)[grepl("TrP", colnames(pepdf))])]
upper_tri <- melt(get_upper_tri(round(cor(orderdf, method = "spearman"),2)), na.rm = TRUE)
ggplot(data = upper_tri, aes(X2, X1, fill = value))+
  geom_tile(color = "white") +
  scale_fill_gradientn(colors = brewer.pal(9,"YlOrRd")) +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 9, hjust = 1),
        axis.title.x=element_blank(), axis.title.y=element_blank())+
  coord_fixed()

orderdf = protdf[,sort(colnames(protdf)[grepl("LiP", colnames(protdf))])]
upper_tri <- melt(get_upper_tri(round(cor(orderdf, method = "spearman"),2)), na.rm = TRUE)
ggplot(data = upper_tri, aes(X2, X1, fill = value))+
  geom_tile(color = "white") +
  scale_fill_gradientn(colors = brewer.pal(9,"YlOrRd")) +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 9, hjust = 1),
        axis.title.x=element_blank(), axis.title.y=element_blank())+
  coord_fixed()

orderdf = pepdf[,sort(colnames(pepdf)[grepl("LiP", colnames(pepdf))])]
upper_tri <- melt(get_upper_tri(round(cor(orderdf, method = "spearman"),2)), na.rm = TRUE)
ggplot(data = upper_tri, aes(X2, X1, fill = value))+
  geom_tile(color = "white") +
  scale_fill_gradientn(colors = brewer.pal(9,"YlOrRd")) +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 9, hjust = 1),
        axis.title.x=element_blank(), axis.title.y=element_blank())+
  coord_fixed()

# Check Batch effects
# Glucose
idx = grepl("G",colnames(pepdf))

res.pca = prcomp(t(pepdf[,idx]), center = T, scale. = T)
fviz_pca_ind(res.pca, repel = TRUE, pointsize=3, invisible="quali",
             habillage=sapply(samplenums[idx], function(x) sampletable[sampletable$LiP_MSnumber==x,"Date.LiP"])) +
  scale_color_brewer(palette="Dark2") +
  theme_minimal()
fviz_pca_ind(res.pca, repel = TRUE, pointsize=3, invisible="quali",
             habillage=sapply(samplenums[idx], function(x) sampletable[sampletable$LiP_MSnumber==x,"Time.point"])) +
  scale_color_manual(values = colorRampPalette(brewer.pal(9, "Blues"))(12)[4:12]) +
  theme_minimal()

# Nitrogen
idx = grepl("N",colnames(pepdf))

res.pca = prcomp(t(pepdf[,idx]), center = T, scale. = T)
fviz_pca_ind(res.pca, repel = TRUE, pointsize=3, invisible="quali",
             habillage=sapply(samplenums[idx], function(x) sampletable[sampletable$LiP_MSnumber==x,"Date.LiP"])) +
  scale_color_brewer(palette="Dark2") +
  theme_minimal()
fviz_pca_ind(res.pca, repel = TRUE, pointsize=3, invisible="quali",
             habillage=sapply(samplenums[idx], function(x) sampletable[sampletable$LiP_MSnumber==x,"Time.point"])) +
  scale_color_manual(values = colorRampPalette(brewer.pal(9, "Blues"))(12)[4:12]) +
  theme_minimal()

# Phosphate
idx = grepl("P",substr(colnames(protdf),1,1))

res.pca = prcomp(t(pepdf[,idx]), center = T, scale. = T)
fviz_pca_ind(res.pca, repel = TRUE, pointsize=3, invisible="quali",
             habillage=sapply(samplenums[idx], function(x) sampletable[sampletable$LiP_MSnumber==x,"Date.LiP"])) +
  scale_color_brewer(palette="Dark2") +
  theme_minimal()
fviz_pca_ind(res.pca, repel = TRUE, pointsize=3, invisible="quali",
             habillage=sapply(samplenums[idx], function(x) sampletable[sampletable$LiP_MSnumber==x,"Time.point"])) +
  scale_color_manual(values = colorRampPalette(brewer.pal(9, "Blues"))(12)[4:12]) +
  theme_minimal()

