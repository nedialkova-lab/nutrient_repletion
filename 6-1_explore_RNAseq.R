# ==============================================================================
# 6-1_explore_RNAseq.R
# Total RNA-seq: TPM normalisation, differential expression along the 10-min
# repletions (limma-voom, cubic spline).
#
# Manuscript: Fig. 2A (RNA part of the Venn diagrams)
# Inputs    : data/20250204_tRNA_RNA/featureCounts.rna.txt
# Outputs   : results/diffRNA_{G,N,P}.tsv
#             results/TPM_RNA.tsv
#             figure_panels/venn_RNAseq.pdf
# Run from the repository root (all paths are relative to it).
# ==============================================================================

library(edgeR)
library(splines)
library(factoextra)
library(RColorBrewer)
library(eulerr)
library(clusterProfiler)
library(org.Sc.sgd.db)
library(gplots)
library(UniProt.ws)
library(RUVSeq)
library(reshape)
timecolors = colorRampPalette(brewer.pal(9, "Blues"))(10)[5:10]
Colors=brewer.pal(11,"YlOrRd")
Colors=colorRampPalette(Colors)(100)

# Load data
countsdf = read.csv("data/20250204_tRNA_RNA/featureCounts.rna.txt", sep="\t", skip = 1)
cleandf = countsdf[,7:ncol(countsdf)]
rownames(cleandf) = sapply(countsdf$Geneid, function(x) strsplit(x,"_")[[1]][1])
colnames(cleandf) = sapply(colnames(cleandf), function(x) paste0(strsplit(strsplit(x,"Aligned")[[1]][1],"_")[[1]][5:6], collapse="_"))

#### LIMMA DIFF EXP ####
#### Glucose ####
idx = grepl("g",substr(colnames(cleandf),1,1))
subsetdf = cleandf[,idx]
heatmap.2(cor(subsetdf, method = "spearman", use = "pairwise"), Rowv = F, Colv = F, col=Colors, margins = c(7,7))
timepointsG = sapply(colnames(subsetdf), function(x) as.numeric(substr(strsplit(x,"_")[[1]][1],2,4)))
# Remove genes with less than 1 cpm
d0 = DGEList(subsetdf)
cutoff <- 1
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop,] 
dim(d) # number of genes left
# Normalize by spikeins for TPM
spike_in_factor = colSums(subsetdf[c(nrow(subsetdf),nrow(subsetdf)-1),]) / colSums(subsetdf)
counts_without_spike_in = colSums(subsetdf[-c(nrow(subsetdf),nrow(subsetdf)-1),])
norm.factors <- spike_in_factor / counts_without_spike_in
norm.factors <- norm.factors / prod(norm.factors)^(1/length(norm.factors))
# Define model to be fitted
X <- ns(timepointsG, knots=c(0.5, 2.5, 5)) # represent time course as cubic spline
design = model.matrix(~ X)
# Voom transformation
y <- voom(d, design, lib.size = counts_without_spike_in, plot = T)
# Fit model
fit <- lmFit(y, design)
fit.eb <- eBayes(fit)
res.an.G <- topTable(fit.eb, coef=2:5, number = 6000)
# TPMs
rpk = apply(subsetdf,2,function(x) x/(countsdf$Length/1000))
kbnorm = colSums(rpk)/1e6
tpmdf = t(apply(rpk,1,function(x) x/kbnorm))

#### Nitrogen ####
idx = grepl("n",substr(colnames(cleandf),1,1))
subsetdf = cleandf[,idx]
heatmap.2(cor(subsetdf, method = "spearman", use = "pairwise"), Rowv = F, Colv = F, col=Colors, margins = c(7,7))
timepointsN = sapply(colnames(subsetdf), function(x) as.numeric(substr(strsplit(x,"_")[[1]][1],2,4)))
# Remove genes with less than 1 cpm
d0 = DGEList(subsetdf)
cutoff <- 1
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop,] 
dim(d) # number of genes left
# Normalize by spikeins
spike_in_factor = colSums(subsetdf[c(nrow(subsetdf),nrow(subsetdf)-1),]) / colSums(subsetdf)
counts_without_spike_in = colSums(subsetdf[-c(nrow(subsetdf),nrow(subsetdf)-1),])
norm.factors <- spike_in_factor / counts_without_spike_in
norm.factors <- norm.factors / prod(norm.factors)^(1/length(norm.factors))
# Define model to be fitted
X <- ns(timepointsN, knots=c(0.5, 2.5, 5)) # represent time course as cubic spline
design = model.matrix(~ X)
# Voom transformation
y <- voom(d, design, lib.size = counts_without_spike_in, plot = T)
# Fit model
fit <- lmFit(y, design)
fit.eb <- eBayes(fit)
res.an.N <- topTable(fit.eb, coef=2:5, number = 6000)
# TPMs
rpk = apply(subsetdf,2,function(x) x/(countsdf$Length/1000))
kbnorm = colSums(rpk)/1e6
tpmdf = cbind(tpmdf,t(apply(rpk,1,function(x) x/kbnorm)))

#### Phosphate ####
idx = grepl("p",substr(colnames(cleandf),1,1))
subsetdf = cleandf[,idx]
heatmap.2(cor(subsetdf, method = "spearman", use = "pairwise"), Rowv = F, Colv = F, col=Colors, margins = c(7,7))
timepointsP = sapply(colnames(subsetdf), function(x) as.numeric(substr(strsplit(x,"_")[[1]][1],2,4)))
# Remove genes with less than 1 cpm
d0 = DGEList(subsetdf)
cutoff <- 1
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop,] 
dim(d) # number of genes left
# Normalize by spikeins
spike_in_factor = colSums(subsetdf[c(nrow(subsetdf),nrow(subsetdf)-1),]) / colSums(subsetdf)
counts_without_spike_in = colSums(subsetdf[-c(nrow(subsetdf),nrow(subsetdf)-1),])
norm.factors <- spike_in_factor / counts_without_spike_in
norm.factors <- norm.factors / prod(norm.factors)^(1/length(norm.factors))
# Define model to be fitted
X <- ns(timepointsP, knots=c(0.5, 2.5, 5)) # represent time course as cubic spline
design = model.matrix(~ X)
# Voom transformation
y <- voom(d, design, lib.size = counts_without_spike_in, plot = T)
# Fit model
fit <- lmFit(y, design)
fit.eb <- eBayes(fit)
res.an.P <- topTable(fit.eb, coef=2:5, number = 6000)
# TPMs
rpk = apply(subsetdf,2,function(x) x/(countsdf$Length/1000))
kbnorm = colSums(rpk)/1e6
tpmdf = cbind(tpmdf,t(apply(rpk,1,function(x) x/kbnorm)))

# GO ontology
gotab = select(org.Sc.sgd.db, rownames(cleandf), columns=c("GOALL","ORF"), keytype="ORF")
mapgo = go2term(unique(gotab$GO))
gotab$goterm = as.character(sapply(gotab$GO, function(x) mapgo[mapgo$go_id==x,"Term"]))
gotab = unique(gotab[gotab$ONTOLOGY=="BP",c("goterm","ORF")])

# Venn diagram
gotrans = gotab[gotab$goterm=="translation","ORF"]
overlap = list(G_diff = rownames(res.an.G[(res.an.G$adj.P.Val<0.05)&!is.na(res.an.G$adj.P.Val),]),
               N_diff = rownames(res.an.N[(res.an.N$adj.P.Val<0.05)&!is.na(res.an.N$adj.P.Val),]),
               P_diff = rownames(res.an.P[(res.an.P$adj.P.Val<0.05)&!is.na(res.an.P$adj.P.Val),]))
               #GO_Translation = gotrans)
fit <- euler(overlap)
pdf("figure_panels/venn_RNAseq.pdf", width = 4, height = 4)
plot(fit, fills = c("#f4d7d7","#f4f4d7","#d7e6f4"), edges=T, legend=F, labels=T, quantities=T,main="RNA")
dev.off()

# Write results
res.an.G = cbind(res.an.G, tpmdf[rownames(res.an.G),grepl("g", colnames(tpmdf))],data.frame(gene=mapIds(org.Sc.sgd.db, rownames(res.an.G), column="GENENAME", keytype="ORF")))
res.an.N = cbind(res.an.N, tpmdf[rownames(res.an.N),grepl("n", colnames(tpmdf))],data.frame(gene=mapIds(org.Sc.sgd.db, rownames(res.an.N), column="GENENAME", keytype="ORF")))
res.an.P = cbind(res.an.P, tpmdf[rownames(res.an.P),grepl("^p", colnames(tpmdf))],data.frame(gene=mapIds(org.Sc.sgd.db, rownames(res.an.P), column="GENENAME", keytype="ORF")))
write.table(res.an.G, "results/diffRNA_G.tsv", sep="\t", quote = F)
write.table(res.an.N, "results/diffRNA_N.tsv", sep="\t", quote = F)
write.table(res.an.P, "results/diffRNA_P.tsv", sep="\t", quote = F)
write.table(tpmdf, "results/TPM_RNA.tsv", sep="\t", quote = F)

#### Principal component analysis ####
tpm_filt = tpmdf
tpm_filt[tpm_filt<1] <- NA
tpm_filt = tpm_filt[rowSums(is.na(tpm_filt))==0,]
# Glucose
res.pca = prcomp(t(tpm_filt[,grepl("g",colnames(tpm_filt))]), center = T, scale. = T)
fviz_eig(res.pca)
fviz_pca_ind(res.pca, repel = TRUE, pointsize=3, invisible="quali",
             habillage=sapply(colnames(tpm_filt[,grepl("g",colnames(tpm_filt))]), function(x) as.numeric(substr(strsplit(x,"_")[[1]][1],2,4)))) +
  scale_color_manual(values = timecolors) +
  theme_minimal()
# Nitrogen
res.pca = prcomp(t(tpm_filt[,grepl("n",colnames(tpm_filt))]), center = T, scale. = T)
fviz_eig(res.pca)
fviz_pca_ind(res.pca, repel = TRUE, pointsize=3, invisible="quali",
             habillage=sapply(colnames(tpm_filt[,grepl("n",colnames(tpm_filt))]), function(x) as.numeric(substr(strsplit(x,"_")[[1]][1],2,4)))) +
  scale_color_manual(values = timecolors) +
  theme_minimal()
# Phosphate
res.pca = prcomp(t(tpm_filt[,grepl("p",substr(colnames(cleandf),1,1))]), center = T, scale. = T)
fviz_eig(res.pca)
fviz_pca_ind(res.pca, repel = TRUE, pointsize=3, invisible="quali",
             habillage=sapply(colnames(tpm_filt[,grepl("p",substr(colnames(cleandf),1,1))]), function(x) as.numeric(substr(strsplit(x,"_")[[1]][1],2,4)))) +
  scale_color_manual(values = timecolors) +
  theme_minimal()
