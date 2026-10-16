# ==============================================================================
# 7-3_spikein_normalization.R
# Spike-in (Fluc trisomes, Rluc disomes) calibration of absolute ribosomes per
# transcript (RxT, Tomuro et al.), replicate correlation and differential RxT
# analysis (limma, cubic spline).
#
# Manuscript: Fig. 2B; Extended Data Fig. 3A
# Inputs    : data/20240529_RNA_RPF/counts_psite.txt
#             data/20241210_RNA_RPF/counts_psite.txt
#             data/20250204_tRNA_RNA/featureCounts.rna.txt
# Outputs   : results/RxT.tsv
#             results/diffRxT_{G,N,P}.tsv
#             figure_panels/spikein_calibration.pdf
#             figure_panels/cormatrix_RxT.pdf
#             figure_panels/venn_RxT.pdf
# Run from the repository root (all paths are relative to it).
# ==============================================================================

library(edgeR)
library(splines)
library(factoextra)
library(RColorBrewer)
library(eulerr)
library(clusterProfiler)
library(msigdbr)
library(gplots)
library(UniProt.ws)
library(RUVSeq)
library(reshape)
library(ggrepel)
library(org.Sc.sgd.db)

# Load data
countsdf1 = read.csv("data/20240529_RNA_RPF/counts_psite.txt", sep="\t", row.names = 1)
countsdf2 = read.csv("data/20241210_RNA_RPF/counts_psite.txt", sep="\t", row.names = 1)
rpfdf = cbind(countsdf1[,2:ncol(countsdf1)],countsdf2[rownames(countsdf1),2:ncol(countsdf2)])
colnames(rpfdf) = sapply(colnames(rpfdf), function(x) paste0(strsplit(x,"_")[[1]][4:5], collapse="_"))

rna.countsdf = read.csv("data/20250204_tRNA_RNA/featureCounts.rna.txt", sep="\t", skip = 1)
rnadf = rna.countsdf[,7:ncol(rna.countsdf)]
rownames(rnadf) = sapply(rna.countsdf$Geneid, function(x) strsplit(x,"_")[[1]][1])
colnames(rnadf) = sapply(colnames(rnadf), function(x) paste0(strsplit(strsplit(x,"Aligned")[[1]][1],"_")[[1]][5:6], collapse="_"))

# Check percentage spikeins
summary(rowSums(apply(rpfdf[c("Rluc","Fluc"),],1,function(x) x*100/colSums(rpfdf))))
summary(rowSums(apply(rnadf[c("Rluc","Fluc"),],1,function(x) x*100/colSums(rnadf))))

#### TPM ####
rpk = apply(rnadf,2,function(x) x/(rna.countsdf$Length))
kbnorm = colSums(rpk)/1e6
rna.tpmdf = t(apply(rpk,1,function(x) x/kbnorm))
rna.tpmdf[rna.tpmdf<1] <- NA

rpf.rpmdf = t(apply(rpfdf[rownames(rna.tpmdf),],1,function(x) x*1e6/colSums(rpfdf)))
rpf.rpmdf[rpf.rpmdf<1] <- NA

rnasamples = sapply(colnames(rpf.rpmdf), function(x) grep(substr(x,1,nchar(x)-3),colnames(rna.tpmdf)))
ggplot(data.frame(RPF=as.numeric(rpf.rpmdf["Rluc",]), RNA=as.numeric(rna.tpmdf["Rluc",rnasamples])), aes(x=RNA,y=RPF,label=colnames(rpf.rpmdf))) +
  geom_point() +
  geom_text_repel()

#### RIBOSOMES PER TRANSCRIPT ####
# Get normalization factors from spikeins
nf = as.numeric(5/((rpf.rpmdf["Rluc",]/rna.tpmdf["Rluc",rnasamples]) + (rpf.rpmdf["Fluc",]/rna.tpmdf["Fluc",rnasamples])))
# Normalize
rpfxtrans = t(apply(rpf.rpmdf/rna.tpmdf[,rnasamples],1,function(x) x*nf))
# Plot normalized samples
plotspikes = melt(as.matrix(rpfxtrans[c("Rluc","Fluc"),]))
plotspikes$cond = substr(plotspikes$X2,0,1)
plotspikes$expected = sapply(plotspikes$X1,function(x) if(x=="Rluc"){2}else{3})
ggplot(plotspikes, aes(x=expected,y=value,label=X2, color=X1)) +
  geom_point() +
  geom_text_repel() +
  expand_limits(x = 0, y = 0) +
  theme_bw()

p = ggplot(plotspikes, aes(x=X1,y=value,label=X2)) +
  facet_wrap(~cond) +
  geom_jitter(alpha=0.5) +
  stat_summary(fun = "mean", geom = "crossbar", width = 0.4, linewidth=0.4, color="#D6604D") +
  expand_limits(y = 0) +
  theme_bw()
ggsave("figure_panels/spikein_calibration.pdf", p, width = 3.5, height = 3)

# Check overall trend over time
longdf = melt(rpfxtrans)
longdf["sample"] = sapply(as.character(longdf$X2), function(x) strsplit(x, "_")[[1]][1])
meandf = cast(longdf, "X1 ~ sample", fun.aggregate = "mean", na.rm=T)
longmeandf = melt(meandf)
longmeandf["cond"] = factor(sapply(as.character(longmeandf$sample), function(x) substr(x, 1,1)), levels=c("g","n","p"))
longmeandf["time"] = as.factor(as.numeric(sapply(as.character(longmeandf$sample), function(x) substr(x, 2,3))))
ggplot(longmeandf, aes(x=time,y=value)) +
  facet_wrap(~cond) +
  geom_jitter(alpha=0.1, size=0.5, aes(color=cond)) +
  geom_boxplot(outliers = F, aes(fill=cond)) +
  scale_color_manual(values=c("#ca403c", "#909136", "#377dbf")) +
  scale_fill_manual(values=c("#f4d7d7", "#f4f4d8", "#d7e6f3")) +
  scale_y_log10() +
  theme_bw()

#### Correlation analysis ####
Colors=brewer.pal(11,"YlOrRd")
Colors=colorRampPalette(Colors)(100)
pdf("figure_panels/cormatrix_RxT.pdf", width = 5, height = 5)
# Glucose
subsetsamples = colnames(rpfxtrans)[grepl("g",colnames(rpfxtrans))]
subsetsamples = subsetsamples[order(as.numeric(sapply(subsetsamples, function(x) substr(strsplit(x,"_")[[1]][1],2,3))))]
heatmap.2(cor(rpfxtrans[,subsetsamples], method = "spearman", use = "pairwise"), 
          Rowv = F, Colv = F, col=Colors, trace = "none", breaks = seq(0.85, 1, 0.0015))
# Nitrogen
subsetsamples = colnames(rpfxtrans)[grepl("n",colnames(rpfxtrans))]
subsetsamples = subsetsamples[order(as.numeric(sapply(subsetsamples, function(x) substr(strsplit(x,"_")[[1]][1],2,3))))]
heatmap.2(cor(rpfxtrans[,subsetsamples], method = "spearman", use = "pairwise"), 
          Rowv = F, Colv = F, col=Colors, trace = "none", breaks = seq(0.85, 1, 0.0015))
# Phosphate
subsetsamples = colnames(rpfxtrans)[grepl("^p",colnames(rpfxtrans))]
subsetsamples = subsetsamples[order(as.numeric(sapply(subsetsamples, function(x) substr(strsplit(x,"_")[[1]][1],2,3))))]
heatmap.2(cor(rpfxtrans[,subsetsamples], method = "spearman", use = "pairwise"), 
          Rowv = F, Colv = F, col=Colors, trace = "none", breaks = seq(0.85, 1, 0.0015))
dev.off()

# ANOVA with limma
# Glucose
idx = grepl("^g",colnames(rpfxtrans))
subsetdf = rpfxtrans[,idx]
timepointsG = sapply(colnames(subsetdf), function(x) as.numeric(substr(strsplit(x,"_")[[1]][1],2,4)))
X <- ns(timepointsG, knots=c(0.5, 2.5, 5)) # represent time course as spline
design = model.matrix(~ X)
fit <- lmFit(log10(subsetdf), design)
fit.eb <- eBayes(fit)
res.an.G <- topTable(fit.eb, number = 6000)
sum(res.an.G$adj.P.Val<0.05, na.rm = T)
# Nitrogen
idx = grepl("^n",colnames(rpfxtrans))
subsetdf = rpfxtrans[,idx]
timepointsN = sapply(colnames(subsetdf), function(x) as.numeric(substr(strsplit(x,"_")[[1]][1],2,4)))
X <- ns(timepointsN, knots=c(0.5, 2.5, 5)) # represent time course as spline
design = model.matrix(~ X)
fit <- lmFit(log10(subsetdf), design)
fit.eb <- eBayes(fit)
res.an.N <- topTable(fit.eb, number = 6000)
sum(res.an.N$adj.P.Val<0.05, na.rm = T)
# Phosphate
idx = grepl("^p",colnames(rpfxtrans))
subsetdf = rpfxtrans[,idx]
timepointsP = sapply(colnames(subsetdf), function(x) as.numeric(substr(strsplit(x,"_")[[1]][1],2,4)))
X <- ns(timepointsP, knots=c(0.5, 2.5, 5)) # represent time course as spline
design = model.matrix(~ X)
fit <- lmFit(log10(subsetdf), design)
fit.eb <- eBayes(fit)
res.an.P <- topTable(fit.eb, number = 6000)
sum(res.an.P$adj.P.Val<0.05, na.rm = T)

# GO ontology
gotab = select(org.Sc.sgd.db, rownames(rpfxtrans), columns=c("GOALL","ORF"), keytype="ORF")
mapgo = go2term(unique(gotab$GO))
gotab$goterm = as.character(sapply(gotab$GO, function(x) mapgo[mapgo$go_id==x,"Term"]))
gotab = unique(gotab[gotab$ONTOLOGY=="BP",c("goterm","ORF")])

# Venn diagram
gotrans = gotab[gotab$goterm=="translation","ORF"]
overlap = list(G_diff = rownames(res.an.G[(res.an.G$adj.P.Val<0.05)&!is.na(res.an.G$adj.P.Val),]),
               N_diff = rownames(res.an.N[(res.an.N$adj.P.Val<0.05)&!is.na(res.an.N$adj.P.Val),]),
               P_diff = rownames(res.an.P[(res.an.P$adj.P.Val<0.05)&!is.na(res.an.P$adj.P.Val),]))
fit <- euler(overlap)
pdf("figure_panels/venn_RxT.pdf", width = 4, height = 4)
plot(fit, fills = c("#f4d7d7","#f4f4d7","#d7e6f4"), edges=T, legend=F, labels=T, quantities=T,main="RxT")
dev.off()

# Differential expression of ctrls vs conditions
longdf = melt(as.matrix(rpfxtrans))
longdf$condition = sapply(as.character(longdf$X2),function(x) strsplit(x,"_")[[1]][1])
meandf = cast(longdf, X1 ~ condition, mean, na.rm=T)
rownames(meandf)= meandf$X1
meandf= meandf[, -c(1)]

# Write results
res.an.G = cbind(res.an.G, meandf[rownames(res.an.G),grepl("g", colnames(meandf))],data.frame(gene=mapIds(org.Sc.sgd.db, rownames(res.an.G), column="GENENAME", keytype="ORF")))
res.an.N = cbind(res.an.N, meandf[rownames(res.an.N),grepl("n", colnames(meandf))],data.frame(gene=mapIds(org.Sc.sgd.db, rownames(res.an.N), column="GENENAME", keytype="ORF")))
res.an.P = cbind(res.an.P, meandf[rownames(res.an.P),grepl("p", colnames(meandf))],data.frame(gene=mapIds(org.Sc.sgd.db, rownames(res.an.P), column="GENENAME", keytype="ORF")))
write.table(res.an.G, "results/diffRxT_G.tsv", sep="\t", quote = F)
write.table(res.an.N, "results/diffRxT_N.tsv", sep="\t", quote = F)
write.table(res.an.P, "results/diffRxT_P.tsv", sep="\t", quote = F)
write.table(rpfxtrans, "results/RxT.tsv", sep="\t", quote = F)
