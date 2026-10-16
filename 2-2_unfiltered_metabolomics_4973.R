# ==============================================================================
# 2-2_unfiltered_metabolomics_4973.R
# Untargeted metabolomics (MASSter + tima identification): normalisation to internal
# standards + PQN, batch correction (limma::removeBatchEffect) and differential
# analysis along the 30-min time course (limma, cubic spline).
#
# Manuscript: Fig. 1C (AMP), 1D (Venn); Extended Data Fig. 2C (adenosine, guanosine); Asp panel not shown
# Inputs    : data/metabolomics_4973/Orbi4973_AA_Replicates_Custome_quant.csv
#             data/metabolomics_4973/study_identified_tima.xlsx
#             data/metabolomics_4973/study_identified_ms1.xlsx
#             data/metabolomics_4973/Orbi4973_231218_metabolomics_81samples.xlsx
# Outputs   : results/metabolomics_unfiltered_normalized.tsv
#             results/metabolomics_unfiltered_normalized_batchcorrected.tsv
#             results/diffMetabolomics_unfiltered_{G,N,P}.tsv
#             figure_panels/venn_metabolomics.pdf
#             figure_panels/metabolite_{AMP,Asp,guanosine,adenosine}.pdf
# Run from the repository root (all paths are relative to it).
# ==============================================================================

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
library(MetaboAnalystR)

# Load data
targ_F5 = read.csv("data/metabolomics_4973/Orbi4973_AA_Replicates_Custome_quant.csv", sep=",", na.strings = c("#N/A"))
untarg_features = read.xlsx("data/metabolomics_4973/study_identified_tima.xlsx", sheet = 2, rowNames = T)
untarg_id = read.xlsx("data/metabolomics_4973/study_identified_tima.xlsx", sheet = 3, rowNames = F)
untarg_ms1 = read.xlsx("data/metabolomics_4973/study_identified_ms1.xlsx", sheet = 2, rowNames = T)
untarg_ms1$id_source = sapply(untarg_ms1$id_top_name, function(x) if(is.na(x)){NA}else{"masster-ms1-ChEBI"})
untarg_F5 = read.xlsx("data/metabolomics_4973/study_identified_tima.xlsx", sheet = 5, rowNames = T)
untarg_F5[untarg_F5==0] = NA
untarg_F5["feature"] = sapply(rownames(untarg_F5), function(x) sprintf("%.5f__%.3f", as.numeric(untarg_features[x,"mz"]), as.numeric(untarg_features[x,"rt"])/60))
sampletable = read.xlsx("data/metabolomics_4973/Orbi4973_231218_metabolomics_81samples.xlsx", rowNames = F)

# Normalize to internal standards
is_df = targ_F5[, c("Replicate", "Molecule.Name", "heavy.Area")]
is_df = is_df[!is.na(is_df$heavy.Area),]
is_df[grepl("01re",is_df$Replicate),1] = "01"
# Map replicate → sample name
is_df$Sample.name = sapply(is_df$Replicate, function(x) {
  sampletable[sampletable$X6 == as.numeric(x), "Sample.name"]
})
is_wide = reshape2::dcast(is_df, Molecule.Name ~ Sample.name, value.var = "heavy.Area")
standards = c("Alanine", "Lysine","Histidine","Arginine","Tyrosine","Phenylalanine","Methionine","Glutamic Acid",
              "Aspartic Acid","Leucine","Isoleucine","Valine","Threonine","Cysteine","Proline","Serine","Glycine")
is_wide = is_wide[is_wide$Molecule.Name %in% standards,]
# compute normalization factor (median across IS per sample)
is_matrix = log2(as.matrix(is_wide[,-1]))
colnames(is_matrix) = colnames(is_wide)[-1]
is_normfactor = apply(is_matrix, 2, median, na.rm = TRUE)
is_normfactor = is_normfactor - mean(is_normfactor)

### UNTARGETED F5 ###
cleandf = untarg_F5[,grep("Sample",colnames(untarg_F5))]
# Sample N30 R2 got contaminated with water bath water during extraction, remove
cleandf = cleandf[,!grepl("Sample45", colnames(cleandf))]
rownames(cleandf) = untarg_F5$feature
samplenums = gsub('.*Sample([0-9]*).*','\\1',colnames(cleandf))
samplenames = sapply(samplenums, function(x) sampletable[sampletable$X6==as.numeric(x),"Sample.name"])
colnames(cleandf) = samplenames
cleandf = log2(cleandf)

# normalize to IS
cleandf = sweep(cleandf, 2, is_normfactor[colnames(cleandf)], "-")

# PQN normalization
# reference = median spectrum
ref = apply(cleandf, 1, median, na.rm = TRUE)
quotients = sweep(cleandf, 1, ref, "-")
pqn_factors = apply(quotients, 2, median, na.rm = TRUE)
pqn_factors = pqn_factors - mean(pqn_factors)
cleandf_norm = sweep(cleandf, 2, pqn_factors, "-")

# Save
annot = sapply(rownames(cleandf_norm), function(x) unlist(untarg_features[rownames(untarg_F5)[untarg_F5$feature==x],]))
annot = as.data.frame(t(sapply(colnames(annot), function(x) if(is.na(annot["id_top_name",x])){unlist(untarg_ms1[rownames(untarg_F5)[untarg_F5$feature==x],])}else{annot[,x]})))
annot$fullname = sapply(rownames(annot), function(x) if(is.na(annot[x,"id_top_name"])|annot[x,"id_source"]=="masster-ms1-ChEBI"){NA}else{untarg_id[(untarg_id$consensus_id %in% rownames(untarg_F5)[untarg_F5$feature==x])&(untarg_id$shortname %in% annot[x,"id_top_name"]),"name"]})
export_col = c("rt","mz","number_samples","number_ms2","inty_mean","quality",
               "adduct_top","adduct_charge_top","adduct_mass_neutral_top",
               "id_top_name","fullname","id_top_class","id_top_adduct","id_top_score","id_source")
mergeddf = cbind(2^cleandf_norm, annot[,export_col])
write.table(mergeddf, "results/metabolomics_unfiltered_normalized.tsv", sep="\t", quote = F)

# Remove replicate effects for visualization
reps = as.character(sapply(samplenums, function(x) paste0(sampletable[sampletable$X6==as.numeric(x),c("Condition","Replicate")],collapse="")))
condtime = as.character(sapply(samplenums, function(x) paste0(sampletable[sampletable$X6==as.numeric(x),c("Condition","Time.point")], collapse="")))
repcorrected = removeBatchEffect(cleandf_norm, batch = reps, design = model.matrix(~ condtime))
# Annotate and save
mergeddf = cbind(2^repcorrected, annot[,export_col])
write.table(mergeddf, "results/metabolomics_unfiltered_normalized_batchcorrected.tsv", sep="\t", quote = F)

# Principal component analysis
hasna_feat = rowSums(is.na(repcorrected)) > 0
res.pca = prcomp(t(repcorrected[!hasna_feat,]), center = T, scale. = T)
fviz_eig(res.pca)
fviz_pca_ind(res.pca, repel = TRUE, pointsize=3, invisible="quali",
             habillage=sapply(samplenums, function(x) sampletable[sampletable$X6==as.numeric(x),"Condition"])) +
  scale_color_brewer(palette="Dark2") +
  theme_minimal()
# Glucose
res.pca = prcomp(t(repcorrected[!hasna_feat,grepl("G",colnames(repcorrected))]), center = T, scale. = T)
fviz_eig(res.pca)
fviz_pca_ind(res.pca, repel = TRUE, pointsize=3, invisible="quali",
             habillage=sapply(samplenums[grepl("G",colnames(repcorrected))], function(x) sampletable[sampletable$X6==as.numeric(x),"Time.point"])) +
  scale_color_manual(values = colorRampPalette(brewer.pal(9, "Blues"))(12)[4:12]) +
  theme_minimal()
# Nitrogen
res.pca = prcomp(t(repcorrected[!hasna_feat,grepl("N",colnames(repcorrected))]), center = T, scale. = T)
fviz_eig(res.pca)
fviz_pca_ind(res.pca, repel = TRUE, pointsize=3, invisible="quali",
             habillage=sapply(samplenums[grepl("N",colnames(repcorrected))], function(x) sampletable[sampletable$X6==as.numeric(x),"Time.point"])) +
  scale_color_manual(values = colorRampPalette(brewer.pal(9, "Blues"))(12)[4:12]) +
  theme_minimal()
# Phosphate
res.pca = prcomp(t(repcorrected[!hasna_feat,grepl("P",colnames(repcorrected))]), center = T, scale. = T)
fviz_eig(res.pca)
fviz_pca_ind(res.pca, repel = TRUE, pointsize=3, invisible="quali",
             habillage=sapply(samplenums[grepl("P",colnames(repcorrected))], function(x) sampletable[sampletable$X6==as.numeric(x),"Time.point"])) +
  scale_color_manual(values = colorRampPalette(brewer.pal(9, "Blues"))(12)[4:12]) +
  theme_minimal()

# Heatmap of log2 abundances of samples
Colors=brewer.pal(9,"YlOrRd")
Colors=colorRampPalette(Colors)(100)
heatmap.2(t(repcorrected[!hasna_feat,]), trace="none", col=Colors, margins = c(6,5),
          distfun = function(x) dist(x, method="euclidean"),hclustfun = function(x) hclust(x, method="ward.D2"))

# Corplot between replicates
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}
orderdf = repcorrected[,order(colnames(repcorrected))]
upper_tri <- melt(get_upper_tri(round(cor(orderdf, method = "spearman", use = "pairwise"),2)), na.rm = TRUE)

ggplot(data = upper_tri, aes(X2, X1, fill = value))+
  geom_tile(color = "white") +
  scale_fill_gradientn(colors = brewer.pal(9,"YlOrRd")) +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 9, hjust = 1),
        axis.title.x=element_blank(), axis.title.y=element_blank())+
  coord_fixed()

# Average and corplot
longdf = melt(repcorrected)
longdf$condition = sapply(as.character(longdf$X2),function(x) strsplit(x," ")[[1]][1])
meandf = cast(longdf, X1 ~ condition, mean, na.rm=T)
rownames(meandf)= meandf$X1
meandf= 2^meandf[, -c(1)]
upper_tri <- melt(get_upper_tri(round(cor(as.data.frame(meandf), method = "spearman", use = "pairwise"),2)), na.rm = TRUE)

ggplot(data = upper_tri, aes(X2, X1, fill = value))+
  geom_tile(color = "white") +
  scale_fill_gradientn(colors = brewer.pal(9,"YlOrRd")) +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 9, hjust = 1),
        axis.title.x=element_blank(), axis.title.y=element_blank())+
  coord_fixed()

# ANOVA with limma
time = sapply(samplenums, function(x) sampletable[sampletable$X6==as.numeric(x),"Time.point"])
# Glucose
idx = grepl("G",colnames(repcorrected))
keepmet = (apply(!is.na(repcorrected[,idx]),1, function(x) any(table(substr(names(x)[x],1,3))>=2)) & rowSums(!is.na(repcorrected[,idx]))>=4)
subsetdf = repcorrected[keepmet,idx]
X <- ns(time[idx], knots=c(0.5, 2.5, 5)) # represent time course as spline
design = model.matrix(~ X)
fit <- lmFit(subsetdf, design)
fit.eb <- eBayes(fit)
res.an.G <- topTable(fit.eb, coef=2:5, number = 5000)
sum(res.an.G[!is.na(res.an.G$adj.P.Val),"adj.P.Val"]<0.05)
# Nitrogen
idx = grepl("N",colnames(repcorrected))
keepmet = (apply(!is.na(repcorrected[,idx]),1, function(x) any(table(substr(names(x)[x],1,3))>=2)) & rowSums(!is.na(repcorrected[,idx]))>=4)
subsetdf = repcorrected[keepmet,idx]
X <- ns(time[idx], knots=c(0.5, 2.5, 5)) # represent time course as spline
design = model.matrix(~ X)
fit <- lmFit(subsetdf, design)
fit.eb <- eBayes(fit)
res.an.N <- topTable(fit.eb, coef=2:5, number = 5000)
sum(res.an.N[!is.na(res.an.N$adj.P.Val),"adj.P.Val"]<0.05)
# Phosphate
idx = grepl("P",colnames(repcorrected))
keepmet = (apply(!is.na(repcorrected[,idx]),1, function(x) any(table(substr(names(x)[x],1,3))>=2)) & rowSums(!is.na(repcorrected[,idx]))>=4)
subsetdf = repcorrected[keepmet,idx]
X <- ns(time[idx], knots=c(0.5, 2.5, 5)) # represent time course as spline
design = model.matrix(~ X)
fit <- lmFit(subsetdf, design)
fit.eb <- eBayes(fit)
res.an.P <- topTable(fit.eb, coef=2:5, number = 5000)
sum(res.an.P[!is.na(res.an.P$adj.P.Val),"adj.P.Val"]<0.05)
# Annotate diff expression results
annotG = sapply(rownames(res.an.G), function(x) unlist(untarg_features[rownames(untarg_F5)[untarg_F5$feature==x], grepl("id_",colnames(untarg_features))]))
annotG = as.data.frame(t(sapply(colnames(annotG), function(x) if(is.na(annotG["id_top_name",x])){unlist(untarg_ms1[rownames(untarg_F5)[untarg_F5$feature==x], grepl("id_",colnames(untarg_ms1))])}else{annotG[,x]})))
annotG$name = sapply(rownames(annotG), function(x) if(is.na(annotG[x,"id_top_name"])|annotG[x,"id_source"]=="masster-ms1-ChEBI"){NA}else{untarg_id[(untarg_id$consensus_id %in% rownames(untarg_F5)[untarg_F5$feature==x])&(untarg_id$shortname %in% annotG[x,"id_top_name"]),"name"]})
annotN = sapply(rownames(res.an.N), function(x) unlist(untarg_features[rownames(untarg_F5)[untarg_F5$feature==x], grepl("id_",colnames(untarg_features))]))
annotN = as.data.frame(t(sapply(colnames(annotN), function(x) if(is.na(annotN["id_top_name",x])){unlist(untarg_ms1[rownames(untarg_F5)[untarg_F5$feature==x], grepl("id_",colnames(untarg_ms1))])}else{annotN[,x]})))
annotN$name = sapply(rownames(annotN), function(x) if(is.na(annotN[x,"id_top_name"])|annotN[x,"id_source"]=="masster-ms1-ChEBI"){NA}else{untarg_id[(untarg_id$consensus_id %in% rownames(untarg_F5)[untarg_F5$feature==x])&(untarg_id$shortname %in% annotN[x,"id_top_name"]),"name"]})
annotP = sapply(rownames(res.an.P), function(x) unlist(untarg_features[rownames(untarg_F5)[untarg_F5$feature==x], grepl("id_",colnames(untarg_features))]))
annotP = as.data.frame(t(sapply(colnames(annotP), function(x) if(is.na(annotP["id_top_name",x])){unlist(untarg_ms1[rownames(untarg_F5)[untarg_F5$feature==x], grepl("id_",colnames(untarg_ms1))])}else{annotP[,x]})))
annotP$name = sapply(rownames(annotP), function(x) if(is.na(annotP[x,"id_top_name"])|annotP[x,"id_source"]=="masster-ms1-ChEBI"){NA}else{untarg_id[(untarg_id$consensus_id %in% rownames(untarg_F5)[untarg_F5$feature==x])&(untarg_id$shortname %in% annotP[x,"id_top_name"]),"name"]})
outresG = cbind(res.an.G, annotG, meandf[rownames(res.an.G),grepl("G", colnames(meandf))])
outresN = cbind(res.an.N, annotN, meandf[rownames(res.an.N),grepl("N", colnames(meandf))])
outresP = cbind(res.an.P, annotP, meandf[rownames(res.an.P),grepl("P", colnames(meandf))])
# Save diff expression results
write.table(outresG, "results/diffMetabolomics_unfiltered_G.tsv", sep="\t", quote = F)
write.table(outresN, "results/diffMetabolomics_unfiltered_N.tsv", sep="\t", quote = F)
write.table(outresP, "results/diffMetabolomics_unfiltered_P.tsv", sep="\t", quote = F)

# Venn diagram
overlap = list(G_diff = rownames(res.an.G)[!is.na(res.an.G$adj.P.Val)&(res.an.G$adj.P.Val<0.05)],
               N_diff = rownames(res.an.N)[!is.na(res.an.N$adj.P.Val)&(res.an.N$adj.P.Val<0.05)],
               P_diff = rownames(res.an.P)[!is.na(res.an.P$adj.P.Val)&(res.an.P$adj.P.Val<0.05)])
fit <- euler(overlap)
pdf("figure_panels/venn_metabolomics.pdf", width = 4, height = 4)
plot(fit, fills = c("#f4d7d7","#f4f4d7","#d7e6f4"), edges=T, legend=F, labels=T, quantities=T,main="Metabolites")
dev.off()

# Plot few metabolites
plotdf = as.data.frame(t(2^repcorrected))
plotdf$cond = substr(rownames(plotdf),1,1)
plotdf$time = as.numeric(substr(rownames(plotdf),2,3))
metab = "348.07040__2.226" #AMP
p = ggplot(plotdf, aes(time,get(metab))) + facet_wrap(~cond) + geom_point(alpha=0.6) + stat_summary(geom="line", fun = "mean") + theme_test()
ggsave("figure_panels/metabolite_AMP.pdf", p, width = 9, height = 2.5)
metab = "134.04470__1.513" #Asp
p = ggplot(plotdf, aes(time,get(metab))) + facet_wrap(~cond) + geom_point(alpha=0.6) + stat_summary(geom="line", fun = "mean") + theme_test()
ggsave("figure_panels/metabolite_Asp.pdf", p, width = 9, height = 2.5)
metab = "284.09890__4.351" #guanosine
p = ggplot(plotdf, aes(time,get(metab))) + facet_wrap(~cond) + geom_point(alpha=0.6) + stat_summary(geom="line", fun = "mean") + theme_test()
ggsave("figure_panels/metabolite_guanosine.pdf", p, width = 9, height = 2.5)
metab = "268.10400__3.272" # adenosine
p = ggplot(plotdf, aes(time,get(metab))) + facet_wrap(~cond) + geom_point(alpha=0.6) + stat_summary(geom="line", fun = "mean") + theme_test()
ggsave("figure_panels/metabolite_adenosine.pdf", p, width = 9, height = 2.5)