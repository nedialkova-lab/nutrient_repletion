# ==============================================================================
# 3-4_abundances_timecourse_5199.R
# Differential protein abundance (TrP samples) along the 20-min repletions
# (limma, cubic spline).
#
# Manuscript: Fig. 5B (abundance counts, see also 4-1)
# Requires  : outputs of 3-3_LiPMS_batchcorrection.R
# Inputs    : results/proteomics_batchcorrected.tsv
#             data/20240507_lipms_upshifts_5199/240507_lipms.xlsx
# Outputs   : results/diffProtAbund_{G,N,P}.tsv
#             figure_panels/venn_prot_abundance.pdf
# Run from the repository root (all paths are relative to it).
# ==============================================================================

library(openxlsx)
library(reshape)
library(factoextra)
library(limma)
library(splines)
library(eulerr)
library(UniProt.ws)

# Load data
normdf = read.delim("results/proteomics_batchcorrected.tsv", sep="\t", row.names = 1)
colnames(normdf) = gsub("\\."," ",colnames(normdf))
sampletable = read.xlsx("data/20240507_lipms_upshifts_5199/240507_lipms.xlsx", rowNames = F)

# ANOVA with limma
# Glucose
idx = grepl("G",colnames(normdf))
keepprot = (apply(!is.na(normdf[,idx]),1, function(x) any(table(substr(names(x)[x],1,3))>=2)) & rowSums(!is.na(normdf[,idx]))>=4)
subsetdf = normdf[keepprot,idx]
timepointsG = sapply(colnames(subsetdf), function(x) sampletable[sampletable$Sample.name==x,"Time.point"])
X <- ns(timepointsG, knots=c(0.5, 2.5, 5)) # represent time course as spline
design = model.matrix(~ X)
fit <- lmFit(log2(subsetdf), design)
fit.eb <- eBayes(fit)
res.an.G <- topTable(fit.eb, number = 5000)
sum(res.an.G$adj.P.Val<0.05, na.rm = T)
# Nitrogen
idx = grepl("N",colnames(normdf))
keepprot = (apply(!is.na(normdf[,idx]),1, function(x) any(table(substr(names(x)[x],1,3))>=2)) & rowSums(!is.na(normdf[,idx]))>=4)
subsetdf = normdf[keepprot,idx]
timepointsN = sapply(colnames(subsetdf), function(x) sampletable[sampletable$Sample.name==x,"Time.point"])
X <- ns(timepointsN, knots=c(0.5, 2.5, 5)) # represent time course as spline
design = model.matrix(~ X)
fit <- lmFit(log2(subsetdf), design)
fit.eb <- eBayes(fit)
res.an.N <- topTable(fit.eb, number = 5000)
sum(res.an.N$adj.P.Val<0.05, na.rm = T)
# Phosphate
idx = grepl("P",colnames(normdf))
keepprot = (apply(!is.na(normdf[,idx]),1, function(x) any(table(substr(names(x)[x],1,3))>=2)) & rowSums(!is.na(normdf[,idx]))>=4)
subsetdf = normdf[keepprot,idx]
timepointsP = sapply(colnames(subsetdf), function(x) sampletable[sampletable$Sample.name==x,"Time.point"])
X <- ns(timepointsP, knots=c(0.5, 2.5, 5)) # represent time course as spline
design = model.matrix(~ X)
fit <- lmFit(log2(subsetdf), design)
fit.eb <- eBayes(fit)
res.an.P <- topTable(fit.eb, number = 5000)
sum(res.an.P$adj.P.Val<0.05, na.rm = T)
# Save diff expression results
flatprots = unlist(sapply(unique(c(rownames(res.an.G),rownames(res.an.N),rownames(res.an.P))),function(x) unlist(strsplit(x, ";"))))
mapping = mapUniProt(query=flatprots, from = "UniProtKB_AC-ID", to = "Gene_Name")
res.an.G["gene"] = sapply(rownames(res.an.G), function(x) paste0(mapping[mapping$From %in% unlist(strsplit(x,";")),"To"], collapse = ";"))
res.an.N["gene"] = sapply(rownames(res.an.N), function(x) paste0(mapping[mapping$From %in% unlist(strsplit(x,";")),"To"], collapse = ";"))
res.an.P["gene"] = sapply(rownames(res.an.P), function(x) paste0(mapping[mapping$From %in% unlist(strsplit(x,";")),"To"], collapse = ";"))

# Venn diagram
overlap = list(G_diff = unique(unlist(sapply(rownames(res.an.G)[res.an.G$adj.P.Val<0.05&!is.na(res.an.G$adj.P.Val)], function(x) strsplit(x,";")))),
               N_diff = unique(unlist(sapply(rownames(res.an.N)[res.an.N$adj.P.Val<0.05&!is.na(res.an.N$adj.P.Val)], function(x) strsplit(x,";")))),
               P_diff = unique(unlist(sapply(rownames(res.an.P)[res.an.P$adj.P.Val<0.05&!is.na(res.an.P$adj.P.Val)], function(x) strsplit(x,";")))))
fit <- euler(overlap)
pdf("figure_panels/venn_prot_abundance.pdf", width = 4, height = 4)
plot(fit, fills = c("#f4d7d7","#f4f4d7","#d7e6f4"), edges=T, legend=F, labels=T, quantities=T,main="Protein Abundance")
dev.off()

# Average per cond
longdf = melt(as.matrix(log2(normdf)))
longdf$condition = sapply(as.character(longdf$X2),function(x) strsplit(x," ")[[1]][1])
meandf = cast(longdf, X1 ~ condition, mean, na.rm=T)
rownames(meandf)= meandf$X1
meandf= 2^meandf[, -c(1)]

# Save diff results
res.an.G = cbind(res.an.G, meandf[rownames(res.an.G),grepl("G", colnames(meandf))])
res.an.N = cbind(res.an.N, meandf[rownames(res.an.N),grepl("N", colnames(meandf))])
res.an.P = cbind(res.an.P, meandf[rownames(res.an.P),grepl("P", colnames(meandf))])
write.table(res.an.G, "results/diffProtAbund_G.tsv", sep="\t", quote = F)
write.table(res.an.N, "results/diffProtAbund_N.tsv", sep="\t", quote = F)
write.table(res.an.P, "results/diffProtAbund_P.tsv", sep="\t", quote = F)
