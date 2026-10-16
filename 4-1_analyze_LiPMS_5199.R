# ==============================================================================
# 4-1_analyze_LiPMS_5199.R
# Structural changes from LiP-MS: LiPAnalyzeR 'LiPPep ~ TrPProt' model (FTHTjoin),
# filtering and differential analysis of the residuals (limma, cubic spline).
#
# Manuscript: Fig. 5B; Extended Data Fig. 6D
# Requires  : outputs of 3-3_LiPMS_batchcorrection.R and 3-4_abundances_timecourse_5199.R
# Inputs    : data/20240507_lipms_upshifts_5199/Orbi5199_LiP+Ctrl_HybridLib_Report.tsv
#             data/20240507_lipms_upshifts_5199/240507_lipms.xlsx
#             results/{TrPProt,LiPPep}.txt
#             results/diffProtAbund_{G,N,P}.tsv
#             data/interaction_promis_S11.txt (exploratory comparison at the end of the script)
# Outputs   : results/resLiP.tsv
#             results/diffStructurome_{G,N,P}.tsv
#             figure_panels/venn_resLiP_proteins.pdf
#             figure_panels/prot_abundance_vs_resLiP.pdf
#             figure_panels/GOtranslation_prot_abundance_vs_resLiP.pdf
# Run from the repository root (all paths are relative to it).
# ==============================================================================

library(LiPAnalyzeR)
library(openxlsx)
library(UniProt.ws)
library(clusterProfiler)
library(org.Sc.sgd.db)
library(factoextra)
library(ggrepel)
library(gplots)
library(reshape)
library(RColorBrewer)
library(limma)
library(splines)
library(eulerr)

# Load data
hybriddata = read.csv("data/20240507_lipms_upshifts_5199/Orbi5199_LiP+Ctrl_HybridLib_Report.tsv", sep="\t")
sampletable = read.xlsx("data/20240507_lipms_upshifts_5199/240507_lipms.xlsx", rowNames = F)

# Separate lip and trp
trpdata = hybriddata[hybriddata$R.Replicate %in% 73:144,]
lipdata = hybriddata[hybriddata$R.Replicate %in% 1:72,]

# Annotate sample names
lipdata["Sample"] = apply(lipdata["R.Replicate"],1,function(x) sampletable[sampletable$LiP_MSnumber==x,"Sample.name"])
trpdata["Sample"] = apply(trpdata["R.Replicate"],1,function(x) sampletable[sampletable$TrP_MSnumber==x,"Sample.name"])
lipdata["Condition"] = apply(lipdata["R.Replicate"],1,function(x) sampletable[sampletable$LiP_MSnumber==x,"Condition"])
trpdata["Condition"] = apply(trpdata["R.Replicate"],1,function(x) sampletable[sampletable$TrP_MSnumber==x,"Condition"])

#### LiPAnalyzeR RUV ####

# Create annotation files
annotPepProt <- getPepProtAnnot(trpdata, lipdata, "PEP.StrippedSequence", 
                                "PEP.StrippedSequence", "PG.ProteinGroups", "PEP.DigestType....Trypsin.P..LysC.", 
                                "PEP.PeptidePosition", "PEP.AllOccurringProteinAccessions", 
                                "PEP.NrOfMissedCleavages", "PEP.IsProteotypic")
annotSample <- getSampleAnnot(reportOut = lipdata, sampleName = "Sample", sampleCondition = "Condition", contrastCoding=NULL)
annotSample["timepoint"] = as.numeric(sapply(rownames(annotSample),function(x) sampletable[sampletable$Sample.name==x,"Time.point"]))
annotSample["strip"] = as.factor(sapply(rownames(annotSample),function(x) sampletable[sampletable$Sample.name==x,"Strip"]))
annotSample["rep"] = as.factor(sapply(rownames(annotSample),function(x) paste0(sampletable[sampletable$Sample.name==x,c("Condition","Replicate")], collapse="")))

# Fetch gene names
flatprots = unlist(apply(annotPepProt["Protein"],1,function(x) unlist(strsplit(x, ";"))))
mapping = mapUniProt(query=flatprots, from = "UniProtKB_AC-ID", to = "Gene_Name")
annotPepProt["gene_name"] = apply(annotPepProt["Protein"], 1, function(x) paste0(mapping[mapping$From %in% unlist(strsplit(x,";")),"To"], collapse = ";"))

# Filter and log transform
TrPProt = read.csv("results/TrPProt.txt", sep = "\t")
LiPPep = read.csv("results/LiPPep.txt", sep = "\t")
colnames(TrPProt) = gsub("\\."," ", colnames(TrPProt)); colnames(LiPPep) = gsub("\\."," ", colnames(LiPPep))
QuantityList = list(LiPPep=as.matrix(LiPPep), TrPProt=as.matrix(TrPProt))

# Correct for protein abundance
modelRUV <- runModel(quantityList = QuantityList,
                     annotS = annotSample,
                     formulaRUV = "Y~XProt",
                     formulaContrast = NULL,
                     lowRUV=c(-Inf, 0), upRUV=c(Inf, Inf))
normdf = modelRUV$modelResid
write.table(normdf,"results/resLiP.tsv", sep="\t", quote = F)

# Corplot between replicates
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}
orderdf = normdf[,order(colnames(normdf))]
upper_tri <- melt(get_upper_tri(round(cor(orderdf, method = "spearman", use = "pairwise"),2)), na.rm = TRUE)
ggplot(data = upper_tri, aes(X2, X1, fill = value))+
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "#053061", high = "#67001F", mid = "white", midpoint = 0, limit = c(-1,1), space = "Lab", name="Spearman\nR", na.value = "lightgrey") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 9, hjust = 1),
        axis.title.x=element_blank(), axis.title.y=element_blank())+
  coord_fixed()

# PCA
# Remove NAs
idx = grepl("G",colnames(normdf))
subsetdf = normdf[,idx]
pcadata = subsetdf[(rowSums(is.na(subsetdf))==0),]
# Plot
res.pca = prcomp(t(pcadata))
#fviz_eig(res.pca)
fviz_pca_ind(res.pca, repel = TRUE, pointsize=3, invisible="quali",
             habillage=annotSample[colnames(normdf)[idx],"timepoint"]) +
  scale_color_manual(values = colorRampPalette(brewer.pal(9, "Blues"))(12)[4:12]) +
  theme_minimal()

# ANOVA with limma
# Glucose
idx = grepl("G",colnames(normdf))
keeppept = (apply(!is.na(normdf[,idx]),1, function(x) any(table(substr(names(x)[x],1,3))>=2)) & rowSums(!is.na(normdf[,idx]))>=4)
subsetdf = normdf[keeppept,idx]
timepointsG = annotSample[colnames(normdf)[idx],"timepoint"]
X <- ns(timepointsG, knots=c(0.5, 2.5, 5)) # represent time course as spline
design = model.matrix(~ X)
fit <- lmFit(subsetdf, design)
fit.eb <- eBayes(fit)
res.an.G <- topTable(fit.eb, number = 60000)
sum(res.an.G$adj.P.Val<0.05, na.rm = T)
# Nitrogen
idx = grepl("N",colnames(normdf))
keeppept = (apply(!is.na(normdf[,idx]),1, function(x) any(table(substr(names(x)[x],1,3))>=2)) & rowSums(!is.na(normdf[,idx]))>=4)
subsetdf = normdf[keeppept,idx]
timepointsN = annotSample[colnames(normdf)[idx],"timepoint"]
X <- ns(timepointsN, knots=c(0.5, 2.5, 5)) # represent time course as spline
design = model.matrix(~ X)
fit <- lmFit(subsetdf, design)
fit.eb <- eBayes(fit)
res.an.N <- topTable(fit.eb, number = 60000)
sum(res.an.N$adj.P.Val<0.05, na.rm = T)
# Phosphate
idx = grepl("P",colnames(normdf))
keeppept = (apply(!is.na(normdf[,idx]),1, function(x) any(table(substr(names(x)[x],1,3))>=2)) & rowSums(!is.na(normdf[,idx]))>=4)
subsetdf = normdf[keeppept,idx]
timepointsP = annotSample[colnames(normdf)[idx],"timepoint"]
X <- ns(timepointsP, knots=c(0.5, 2.5, 5)) # represent time course as spline
design = model.matrix(~ X)
fit <- lmFit(subsetdf, design)
fit.eb <- eBayes(fit)
res.an.P <- topTable(fit.eb, number = 60000)
sum(res.an.P$adj.P.Val<0.05, na.rm = T)
# Annotate diff expression results
res.an.G = cbind(res.an.G, annotPepProt[rownames(res.an.G),])
res.an.N = cbind(res.an.N, annotPepProt[rownames(res.an.N),])
res.an.P = cbind(res.an.P, annotPepProt[rownames(res.an.P),])

# Count significant
sort(table(res.an.G[res.an.G$adj.P.Val<0.05,"gene_name"]), decreasing = T)[1:50]
sort(table(res.an.N[res.an.N$adj.P.Val<0.05,"gene_name"]), decreasing = T)[1:50]
sort(table(res.an.P[res.an.P$adj.P.Val<0.05,"gene_name"]), decreasing = T)[1:50]

# Venn diagram
overlap = list(G_diff = rownames(res.an.G)[res.an.G$adj.P.Val<0.05&!is.na(res.an.G$adj.P.Val)],
               N_diff = rownames(res.an.N)[res.an.N$adj.P.Val<0.05&!is.na(res.an.N$adj.P.Val)],
               P_diff = rownames(res.an.P)[res.an.P$adj.P.Val<0.05&!is.na(res.an.P$adj.P.Val)])
fit <- euler(overlap)
plot(fit, fills = c("#f4d7d7","#f4f4d7","#d7e6f4"), edges=T, legend=F, labels=T, quantities=T,main="LiPRes Peptides")

overlap = list(G_diff = unique(unlist(sapply(res.an.G[res.an.G$adj.P.Val<0.05&!is.na(res.an.G$adj.P.Val),"Protein"], function(x) strsplit(x,";")))),
               N_diff = unique(unlist(sapply(res.an.N[res.an.N$adj.P.Val<0.05&!is.na(res.an.N$adj.P.Val),"Protein"], function(x) strsplit(x,";")))),
               P_diff = unique(unlist(sapply(res.an.P[res.an.P$adj.P.Val<0.05&!is.na(res.an.P$adj.P.Val),"Protein"], function(x) strsplit(x,";")))))
fit <- euler(overlap)
pdf("figure_panels/venn_resLiP_proteins.pdf", width = 4, height = 4)
plot(fit, fills = c("#f4d7d7","#f4f4d7","#d7e6f4"), edges=T, legend=F, labels=T, quantities=T,main="LiPRes Proteins")
dev.off()

# Differential expression of ctrls vs conditions
longdf = melt(as.matrix(normdf))
longdf$condition = sapply(as.character(longdf$X2),function(x) strsplit(x," ")[[1]][1])
meandf = cast(longdf, X1 ~ condition, mean, na.rm=T)
rownames(meandf)= meandf$X1
meandf= meandf[, -c(1)]

# Save results
outresG = cbind(res.an.G, meandf[rownames(res.an.G),grepl("G", colnames(meandf))])
outresN = cbind(res.an.N, meandf[rownames(res.an.N),grepl("N", colnames(meandf))])
outresP = cbind(res.an.P, meandf[rownames(res.an.P),grepl("P", colnames(meandf))])
write.table(outresG,"results/diffStructurome_G.tsv", sep="\t", row.names = F, quote = F)
write.table(outresN,"results/diffStructurome_N.tsv", sep="\t", row.names = F, quote = F)
write.table(outresP,"results/diffStructurome_P.tsv", sep="\t", row.names = F, quote = F)

# Compare with abundances
diffG = read.csv("results/diffProtAbund_G.tsv", sep = "\t")
diffN = read.csv("results/diffProtAbund_N.tsv", sep = "\t")
diffP = read.csv("results/diffProtAbund_P.tsv", sep = "\t")
overlap = list(abundance = unique(unlist(sapply(diffG[diffG$adj.P.Val<0.05&!is.na(diffG$adj.P.Val),"gene"], function(x) unlist(strsplit(x,";"))))),
               LiPMS = unique(unlist(sapply(res.an.G[res.an.G$adj.P.Val<0.05&!is.na(res.an.G$adj.P.Val),"gene_name"], function(x) unlist(strsplit(x,";"))))))
fit <- euler(overlap)
plot(fit, fills = c("#f4d7d7","#f4f4d7"), edges=T, legend=F, labels=T, quantities=T,main="Glucose")

overlap = list(abundance = unique(unlist(sapply(diffN[diffN$adj.P.Val<0.05&!is.na(diffN$adj.P.Val),"gene"], function(x) unlist(strsplit(x,";"))))),
               LiPMS = unique(unlist(sapply(res.an.N[res.an.N$adj.P.Val<0.05&!is.na(res.an.N$adj.P.Val),"gene_name"], function(x) unlist(strsplit(x,";"))))))
fit <- euler(overlap)
plot(fit, fills = c("#f4d7d7","#f4f4d7"), edges=T, legend=F, labels=T, quantities=T,main="Nitrogen")

overlap = list(abundance = unique(unlist(sapply(diffP[diffP$adj.P.Val<0.05&!is.na(diffP$adj.P.Val),"gene"], function(x) unlist(strsplit(x,";"))))),
               LiPMS = unique(unlist(sapply(res.an.P[res.an.P$adj.P.Val<0.05&!is.na(res.an.P$adj.P.Val),"gene_name"], function(x) unlist(strsplit(x,";"))))))
fit <- euler(overlap)
plot(fit, fills = c("#f4d7d7","#f4f4d7"), edges=T, legend=F, labels=T, quantities=T,main="Phosphate")

#barplot
bardf = data.frame(cond=c("Glucose","Glucose","Nitrogen","Nitrogen","Phosphate","Phosphate"),
                   metric=c("abundance","LiPMS","abundance","LiPMS","abundance","LiPMS"),
                   prots=c(length(unique(unlist(sapply(rownames(diffG)[diffG$adj.P.Val<0.05&!is.na(diffG$adj.P.Val)], function(x) unlist(strsplit(x,";")))))),
                           length(unique(unlist(sapply(res.an.G[res.an.G$adj.P.Val<0.05&!is.na(res.an.G$adj.P.Val),"Protein"], function(x) unlist(strsplit(x,";")))))),
                           length(unique(unlist(sapply(rownames(diffN)[diffN$adj.P.Val<0.05&!is.na(diffN$adj.P.Val)], function(x) unlist(strsplit(x,";")))))),
                           length(unique(unlist(sapply(res.an.N[res.an.N$adj.P.Val<0.05&!is.na(res.an.N$adj.P.Val),"Protein"], function(x) unlist(strsplit(x,";")))))),
                           length(unique(unlist(sapply(rownames(diffP)[diffP$adj.P.Val<0.05&!is.na(diffP$adj.P.Val)], function(x) unlist(strsplit(x,";")))))),
                           length(unique(unlist(sapply(res.an.P[res.an.P$adj.P.Val<0.05&!is.na(res.an.P$adj.P.Val),"Protein"], function(x) unlist(strsplit(x,";"))))))))
p = ggplot(bardf, aes(x=cond, y=prots, fill=metric)) +
  geom_bar(position="dodge", stat="identity") +
  scale_fill_manual(values=c("#9fccba","#548c78")) +
  theme_bw() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggsave("figure_panels/prot_abundance_vs_resLiP.pdf", p, width = 3, height = 2.5)

# GO translation
includedgenes = unique(unlist(sapply(c(diffG[,"gene"],diffN[,"gene"],diffP[,"gene"],res.an.G[,"gene_name"],res.an.N[,"gene_name"],res.an.P[,"gene_name"]), function(x) unlist(strsplit(x,";")))))
gotab = select(org.Sc.sgd.db, includedgenes, columns=c("GOALL","GENENAME"), keytype="GENENAME")
mapgo = go2term(unique(gotab$GO))
gotab$goterm = as.character(sapply(gotab$GO, function(x) mapgo[mapgo$go_id==x,"Term"]))
gotab = unique(gotab[gotab$ONTOLOGY=="BP",c("goterm","GENENAME")])

compute_recovery <- function(result_df, gene_col, gotab, goterm = "translation") {
  
  # All detected genes in this dataset (universe)
  univ_genes <- unique(unlist(strsplit(result_df[, gene_col], ";")))
  
  # Significant hits
  sig_genes <- unique(unlist(strsplit(
    result_df[!is.na(result_df$adj.P.Val) & result_df$adj.P.Val < 0.05, gene_col], ";"
  )))
  
  go_genes <- gotab[gotab$goterm == goterm, "GENENAME"]
  
  detected_in_go <- sum(univ_genes %in% go_genes)   # denominator: detected translation genes
  hits_in_go     <- sum(sig_genes  %in% go_genes)   # numerator: significant translation genes
  
  list(
    recovery         = hits_in_go / detected_in_go,
    hits_translation = hits_in_go,
    det_translation  = detected_in_go
  )
}

datasets <- list(
  list(res = diffG,    gene_col = "gene",      label = "Glucose",   metric = "abundance"),
  list(res = res.an.G, gene_col = "gene_name", label = "Glucose",   metric = "LiPMS"),
  list(res = diffN,    gene_col = "gene",      label = "Nitrogen",  metric = "abundance"),
  list(res = res.an.N, gene_col = "gene_name", label = "Nitrogen",  metric = "LiPMS"),
  list(res = diffP,    gene_col = "gene",      label = "Phosphate", metric = "abundance"),
  list(res = res.an.P, gene_col = "gene_name", label = "Phosphate", metric = "LiPMS")
)

gotransdf <- do.call(rbind, lapply(datasets, function(d) {
  r <- compute_recovery(d$res, d$gene_col, gotab)
  data.frame(
    cond             = d$label,
    metric           = d$metric,
    recovery         = r$recovery,
    hits_translation = r$hits_translation,
    det_translation  = r$det_translation
  )
}))

p = ggplot(gotransdf, aes(x=cond, y=recovery*100, fill=metric)) +
  geom_bar(position="dodge", stat="identity") +
  scale_fill_manual(values=c("#9fccba","#548c78")) +
  theme_bw() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggsave("figure_panels/GOtranslation_prot_abundance_vs_resLiP.pdf", p, width = 3, height = 2.5)


# Wood plots
res.an.G = cbind(res.an.G, meandf[rownames(res.an.G),grepl("G", colnames(meandf))] - meandf[rownames(res.an.G),"G0"])
res.an.N = cbind(res.an.N, meandf[rownames(res.an.N),grepl("N", colnames(meandf))] - meandf[rownames(res.an.N),"N0"])
res.an.P = cbind(res.an.P, meandf[rownames(res.an.P),grepl("P", colnames(meandf))] - meandf[rownames(res.an.P),"P0"])
makeWoodsPlotSingleProtein(res.an.G, annotPepProt, deltaColorIsTryptic=TRUE,
                           protName = "SES1", nameProtQuant="gene_name",
                           coefCol="G10", pvalCol="adj.P.Val")
makeWoodsPlotSingleProtein(res.an.G, annotPepProt, deltaColorIsTryptic=TRUE,
                           protName = "TEF1", nameProtQuant="gene_name",
                           coefCol="G1", pvalCol="adj.P.Val")
makeWoodsPlotSingleProtein(res.an.G, annotPepProt, deltaColorIsTryptic=TRUE,
                           protName = "STM1", nameProtQuant="gene_name",
                           coefCol="G20", pvalCol="adj.P.Val")

# Plot few peptides
plotdf = as.data.frame(t(normdf))
plotdf$cond = substr(rownames(plotdf),1,1)
plotdf$time = as.numeric(substr(rownames(plotdf),2,3))
pept = "ELVSCSNCTDYQSR" # SES1
p = ggplot(plotdf, aes(time,get(pept))) + facet_wrap(~cond, scales = "free_y") + geom_point(alpha=0.6) + stat_summary(geom="line", fun = "mean") + theme_test()
pept = "IGTVPVGR" #TEF1/2
p = ggplot(plotdf, aes(time,get(pept))) + facet_wrap(~cond, scales = "free_y") + geom_point(alpha=0.6) + stat_summary(geom="line", fun = "mean") + theme_test()
pept = "EAYVPATK" #SMT1
p = ggplot(plotdf, aes(time,get(pept))) + facet_wrap(~cond, scales = "free_y") + geom_point(alpha=0.6) + stat_summary(geom="line", fun = "mean") + theme_test()

