# ==============================================================================
# 1-1_explore_proteomics.R
# Total proteomics of the nutrient-depleted starting cultures (glucose, nitrogen,
# phosphate depletion) compared to cells grown in SD medium. limma differential
# abundance and GO enrichment (clusterProfiler) of up/down-regulated proteins.
#
# Manuscript: Extended Data Fig. 1A (Venn), 1B (heatmap), 1C (GO dotplot)
# Inputs    : data/proteomics_test/Orbi4696_ProteinGroups.tsv
# Outputs   : results/diffTotalProteome.tsv
#             results/GOanalysis_{Glc,Nit,Pho}_totalproteome.tsv
#             source_data/heatmap_totalproteome.tsv
#             figure_panels/{venn,heatmap,GOanalysis}_totalproteome.pdf
# Run from the repository root (all paths are relative to it).
# ==============================================================================

library(ggplot2)
library(gplots)
library(RColorBrewer)
library(factoextra)
library(reshape)
library(clusterProfiler)
library(eulerr)
library(limma)
library(org.Sc.sgd.db)
library(dplyr)
#https://www.biostat.jhsph.edu/~kkammers/software/CVproteomics/R_guide.html

# Load data
protdata = read.csv("data/proteomics_test/Orbi4696_ProteinGroups.tsv", sep="\t", row.names = "PG.ProteinGroups")
cleandf = as.matrix(protdata[,grepl("Quantity",colnames(protdata))])
# Remove NAs and low confidence proteins
cleandf = cleandf[(rowSums(is.na(cleandf))==0)&(!apply(protdata[,grepl("IsSingleHit",colnames(protdata))],1,any)),]
samplenames = c("C1","C2","C3","G1","G2","G3","N1","N2","N3","P1","P2","P3")
colnames(cleandf) = samplenames

# Heatmap of log10 abundances
Colors=brewer.pal(9,"YlOrRd")
Colors=colorRampPalette(Colors)(100)
heatmap.2(log10(t(cleandf)), trace="none", col=Colors,
          distfun = function(x) dist(x, method="euclidean"),hclustfun = function(x) hclust(x, method="ward.D2"))

# PCA analysis
res.pca = prcomp(t(cleandf), scale = T)
fviz_eig(res.pca)
fviz_pca_ind(res.pca, repel = TRUE, pointsize=3, invisible="quali",
             habillage=substr(samplenames,1,1)) +
  scale_color_brewer(palette="Dark2") +
  theme_minimal()

# Differential expression
avgCtrl = rowMeans(cleandf[,c("C1","C2","C3")])
fc_vs_C = log2(cleandf[,c("G1","G2","G3","N1","N2","N3","P1","P2","P3")] / avgCtrl)

# ANOVA with limma
cond = as.factor(substr(colnames(cleandf),1,1))
design = model.matrix(~ cond)
fit <- lmFit(log10(cleandf), design)
fit.eb <- eBayes(fit)
res.anova <- topTable(fit.eb, number=5000)
# Save diff expression results
res.anova_G <- topTable(fit.eb, coef="condG", number=5000)
res.anova_N <- topTable(fit.eb, coef="condN", number=5000)
res.anova_P <- topTable(fit.eb, coef="condP", number=5000)
colnames(res.anova_G) = sprintf("%s_%s","G",colnames(res.anova_G))
colnames(res.anova_N) = sprintf("%s_%s","N",colnames(res.anova_N))
colnames(res.anova_P) = sprintf("%s_%s","P",colnames(res.anova_P))
outdf = cbind(res.anova_G[,c(1,4,5)],
              res.anova_N[rownames(res.anova_G),c(1,4,5)],
              res.anova_P[rownames(res.anova_G),c(1,4,5)])
# Add gene names
outdf["gene"] = protdata[rownames(res.anova_G),"PG.Genes"]
write.table(outdf,"results/diffTotalProteome.tsv", sep="\t", quote = F)

# Heatmap
subsetgenes = unique(c(rownames(outdf)[outdf$G_adj.P.Val<0.05],rownames(outdf)[outdf$N_adj.P.Val<0.05],rownames(outdf)[outdf$P_adj.P.Val<0.05]))
Colors=rev(brewer.pal(9,"RdBu"))
Colors=colorRampPalette(Colors)(100)
pdf("figure_panels/heatmap_totalproteome.pdf", width = 10, height = 4)
heatmap.2(t(fc_vs_C[subsetgenes,]), trace="none", col=Colors,labCol = protdata[subsetgenes,"PG.Genes"],
          distfun = function(x) dist(x, method="euclidean"),hclustfun = function(x) hclust(x, method="ward.D2"),
          breaks=seq(-2,2,0.04))
dev.off()
write.table(t(fc_vs_C[subsetgenes,]), "source_data/heatmap_totalproteome.tsv", sep = "\t")

# Venn diagram
overlap = list(G_diff = unique(unlist(sapply(rownames(outdf)[outdf$G_adj.P.Val<0.05&!is.na(outdf$G_adj.P.Val)], function(x) strsplit(x,";")))),
               N_diff = unique(unlist(sapply(rownames(outdf)[outdf$N_adj.P.Val<0.05&!is.na(outdf$N_adj.P.Val)], function(x) strsplit(x,";")))),
               P_diff = unique(unlist(sapply(rownames(outdf)[outdf$P_adj.P.Val<0.05&!is.na(outdf$P_adj.P.Val)], function(x) strsplit(x,";")))))
fit <- euler(overlap)
pdf("figure_panels/venn_totalproteome.pdf", width = 4, height = 4)
plot(fit, fills = c("#f4d7d7","#f4f4d7","#d7e6f4"), edges=T, legend=F, labels=T, quantities=T,main="Protein Groups")
dev.off()

# GO analysis
includedgenes = unique(unlist(sapply(unique(outdf$gene), function(x) strsplit(x,";")[[1]])))

gotab = select(org.Sc.sgd.db, includedgenes, columns=c("GOALL","GENENAME"), keytype="GENENAME")
mapgo = go2term(unique(gotab$GO))
gotab$goterm = as.character(sapply(gotab$GO, function(x) mapgo[mapgo$go_id==x,"Term"]))
gotab = unique(gotab[gotab$ONTOLOGY=="BP",c("goterm","GENENAME")])

keggtab <- getGeneKEGGLinks(species="sce")
keggtab$gene <- mapIds(org.Sc.sgd.db, keggtab$GeneID, column="GENENAME", keytype="ORF")
mappaths = getKEGGPathwayNames(species="sce")
keggtab$pathname = sapply(keggtab$PathwayID, function(x) substr(mappaths[mappaths$PathwayID==x,"Description"],0, nchar(mappaths[mappaths$PathwayID==x,"Description"])-43))
keggtab = unique(keggtab[,c("pathname","gene")])

genesG = data.frame(genes=unique(unlist(sapply(outdf[outdf$G_adj.P.Val<0.05&!is.na(outdf$G_adj.P.Val),"gene"], function(x) unlist(strsplit(x,";"))))))
genesG$FC = factor(sapply(genesG$genes, function(x) if(mean(outdf[grepl(x,outdf$gene),"G_logFC"])>0){"UP"}else{"DOWN"}), levels = c("UP","DOWN"))
res.goG <- compareCluster(genes~FC, data=genesG, fun="enricher", TERM2GENE = gotab, pvalueCutoff=0.05, universe=includedgenes)
dotplot(res.goG, showCategory=10, title = "GO - Glc")
write.table(res.goG@compareClusterResult, "results/GOanalysis_Glc_totalproteome.tsv", sep = "\t")
res.keggG <- compareCluster(genes~FC, data=genesG, fun="enricher", TERM2GENE = keggtab, pvalueCutoff=0.05, universe=includedgenes)
dotplot(res.keggG, showCategory=10, title = "KEGG - Glc")

genesN = data.frame(genes=unique(unlist(sapply(outdf[outdf$N_adj.P.Val<0.05&!is.na(outdf$N_adj.P.Val),"gene"], function(x) unlist(strsplit(x,";"))))))
genesN$FC = factor(sapply(genesN$genes, function(x) if(mean(outdf[grepl(x,outdf$gene),"N_logFC"])>0){"UP"}else{"DOWN"}), levels = c("UP","DOWN"))
res.goN <- compareCluster(genes~FC, data=genesN, fun="enricher", TERM2GENE = gotab, pvalueCutoff=0.05, universe=includedgenes)
dotplot(res.goN, showCategory=10, title = "GO - Nit")
write.table(res.goN@compareClusterResult, "results/GOanalysis_Nit_totalproteome.tsv", sep = "\t")
res.keggN <- compareCluster(genes~FC, data=genesN, fun="enricher", TERM2GENE = keggtab, pvalueCutoff=0.05, universe=includedgenes)
dotplot(res.keggN, showCategory=10, title = "KEGG - Nit")

genesP = data.frame(genes=unique(unlist(sapply(outdf[outdf$P_adj.P.Val<0.05&!is.na(outdf$P_adj.P.Val),"gene"], function(x) unlist(strsplit(x,";"))))))
genesP$FC = factor(sapply(genesP$genes, function(x) if(mean(outdf[grepl(x,outdf$gene),"P_logFC"])>0){"UP"}else{"DOWN"}), levels = c("UP","DOWN"))
res.goP <- compareCluster(genes~FC, data=genesP, fun="enricher", TERM2GENE = gotab, pvalueCutoff=0.05, universe=includedgenes)
dotplot(res.goP, showCategory=10, title = "GO - Pho")
write.table(res.goP@compareClusterResult, "results/GOanalysis_Pho_totalproteome.tsv", sep = "\t")
res.keggP <- compareCluster(genes~FC, data=genesP, fun="enricher", TERM2GENE = keggtab, pvalueCutoff=0.05, universe=includedgenes)
dotplot(res.keggP, showCategory=10, title = "KEGG - Pho")

# Merged dotplot
genesG$FC = sapply(genesG$FC, function(x) sprintf("Glc-%s", x))
genesN$FC = sapply(genesN$FC, function(x) sprintf("Nit-%s", x))
genesP$FC = sapply(genesP$FC, function(x) sprintf("Pho-%s", x))
genesall = rbind(genesG, genesN, genesP)
res.go <- compareCluster(genes~FC, data=genesall, fun="enricher", TERM2GENE = gotab, pvalueCutoff=0.05, universe=includedgenes)
p = dotplot(res.go, showCategory=3, title = "GO", by="p.adjust", size="GeneRatio")
ggsave("figure_panels/GOanalysis_totalproteome.pdf", p, width = 6.5, height = 6)
