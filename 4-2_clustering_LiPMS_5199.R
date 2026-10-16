# ==============================================================================
# 4-2_clustering_LiPMS_5199.R
# Clustering (1 - |Pearson R|, Ward.D2) and GO enrichment of significantly changing
# LiP-MS peptides, and median absolute change of translation-annotated peptides.
#
# Manuscript: Fig. 5C (top, LiP-MS); Extended Data Fig. 6A-C
# Requires  : outputs of 4-1_analyze_LiPMS_5199.R
# Inputs    : results/diffStructurome_{G,N,P}.tsv
#             results/resLiP.tsv
# Outputs   : results/GOanalysis_{Glc,Pho}_resLiP.tsv
#             source_data/heatmap_resLiP_{Glc,Nit,Pho}.tsv
#             figure_panels/heatmap_resLiP_{Glc,Nit,Pho}{,_dendrogram,_trends}.pdf
#             figure_panels/lineplot_resLiP_GOtrans.pdf
# Run from the repository root (all paths are relative to it).
# ==============================================================================

library(gplots)
library(cluster)
library(ggplot2)
library(RColorBrewer)
library(clusterProfiler)
library(msigdbr)
library(dendextend)
library(limma)
library(org.Sc.sgd.db)

# Record info about translation
gotransdf = data.frame(cond=c(rep("G",8),rep("N",8),rep("P",8)),
                       time=rep(c(0,1,2,4,7,10,14,20),3),
                       fcmedian=NA, fcquantmin=NA, fcquantmax=NA)

#### GLUCOSE ####
# Load data
difflip = read.csv("results/diffStructurome_G.tsv",sep="\t", row.names = "quantID")
meandf = difflip[,(ncol(difflip)-7):ncol(difflip)]
resLiP = read.csv("results/resLiP.tsv",sep="\t")

# Clustering by correlation
X = resLiP[rownames(difflip)[(difflip$adj.P.Val<0.05)&(!is.na(difflip$adj.P.Val))],grepl("G",colnames(resLiP))]
# Make clusters
cordist = 1-abs(cor(t(X), use = "pairwise", method = "pearson"))
notna = rowSums(is.na(cordist))==0
clust_aggl = hclust(as.dist(cordist[notna,notna]), method = "ward.D2")
cutgroups = cutree(clust_aggl, k=3)
dend = color_branches(as.dendrogram(clust_aggl), k=3)
pdf("figure_panels/heatmap_resLiP_Glc_dendrogram.pdf", width = 6, height = 2.5); plot(dend); dev.off()
colmap = setNames(get_leaves_branches_col(dend), names(cutgroups)[order.dendrogram(dend)])
# Check clusters
time = order(as.numeric(substr(colnames(meandf),2,4)))
heatdf = meandf[rownames(X)[notna],]
signcordf = heatdf  - meandf[rownames(X)[notna],"G0"]
for (x in unique(cutgroups)){
  tempclust = cutree(hclust(as.dist(1-cor(t(heatdf[cutgroups==x,]), use = "pairwise", method = "pearson")), method = "ward.D2"), k=2)
  signcordf[names(tempclust)[tempclust==2],] = -signcordf[names(tempclust)[tempclust==2],]
}
Colors=rev(brewer.pal(11,"RdBu"))
Colors=colorRampPalette(Colors)(100)
heatdf = heatdf  - meandf[rownames(X)[notna],"G0"]
heatdf = heatdf[rowSums(!is.na(heatdf))>0,]; signcordf = signcordf[rowSums(!is.na(signcordf))>0,]
orderedpept = unlist(sapply(unique(cutgroups), function(x) rownames(heatdf[rownames(heatdf) %in% names(cutgroups)[cutgroups==x],])[hclust(dist(heatdf[rownames(heatdf) %in% names(cutgroups)[cutgroups==x],], method="euclidean"), method = "ward.D2")$order]))
pdf("figure_panels/heatmap_resLiP_Glc.pdf", width = 7, height = 5)
heatmap.2(t(heatdf[orderedpept,time[2:8]]), trace="none", Rowv=F, Colv=F, col=Colors, margins = c(4,4),breaks=seq(-2,2,0.04),symkey=F,
          labCol = difflip[orderedpept,"gene_name"], ColSideColors = colmap[orderedpept])
dev.off()
write.table(t(heatdf[orderedpept,time[2:8]]),"source_data/heatmap_resLiP_Glc.tsv", sep = "\t")

# Plot trends
trendsdf = data.frame(row.names = 1:(length(time)*length(unique(cutgroups))))
trendsdf$timegroup = c(sapply(unique(cutgroups), function(x) colnames(signcordf[rownames(signcordf) %in% names(cutgroups)[cutgroups==x],])))
trendsdf$time = apply(trendsdf,1, function(x) as.numeric(substr(x["timegroup"],2,4)))
trendsdf$cluster = as.character(c(sapply(unique(cutgroups),function(x) rep(x,length(time)))))
trendsdf$fcmedian = apply(trendsdf,1, function(x) median(signcordf[rownames(signcordf) %in% names(cutgroups)[cutgroups==x["cluster"]],x["timegroup"]], na.rm=T))
trendsdf$fcquantmin = apply(trendsdf,1, function(x) quantile(signcordf[rownames(signcordf) %in% names(cutgroups)[cutgroups==x["cluster"]],x["timegroup"]], probs=0.25, na.rm=T))
trendsdf$fcquantmax = apply(trendsdf,1, function(x) quantile(signcordf[rownames(signcordf) %in% names(cutgroups)[cutgroups==x["cluster"]],x["timegroup"]], probs=0.75, na.rm=T))
p = ggplot(data=trendsdf, aes(x=time, y=fcmedian, ymin=fcquantmin, ymax=fcquantmax)) +
  facet_wrap(~cluster, ncol = 1) + geom_hline(yintercept=0, color="grey") +
  geom_line(aes(color=cluster)) + 
  geom_ribbon(alpha=0.5, aes(fill=cluster)) +
  scale_color_manual(values=unique(colmap[orderedpept])) + scale_fill_manual(values=unique(colmap[orderedpept])) +
  theme_bw()
ggsave("figure_panels/heatmap_resLiP_Glc_trends.pdf", p, width = 2.5, height = 3)

# GO enrichment
includedgenes = unique(unlist(sapply(unique(difflip$gene_name), function(x) strsplit(x,";")[[1]])))
gotab = select(org.Sc.sgd.db, includedgenes, columns=c("GOALL","GENENAME"), keytype="GENENAME")
mapgo = go2term(unique(gotab$GO))
gotab$goterm = as.character(sapply(gotab$GO, function(x) mapgo[mapgo$go_id==x,"Term"]))
gotab = unique(gotab[gotab$ONTOLOGY=="BP",c("goterm","GENENAME")])
countgs = table(gotab$goterm) # Keep only GS with at least 10 genes
gotab = gotab[gotab$goterm %in% names(countgs)[countgs>=10],]
gotab_pept = unlist(sapply(rownames(difflip), function(x) gotab[gotab$GENENAME %in% strsplit(difflip[x,"gene_name"],";")[[1]],"goterm"]))
peptnames = unlist(sapply(rownames(difflip), function(x) rep(x,sum(gotab$GENENAME %in% strsplit(difflip[x,"gene_name"],";")[[1]]))))
gotab_pept = unique(data.frame(gs_name=gotab_pept, gene_symbol=peptnames, row.names = NULL))

keggtab <- getGeneKEGGLinks(species="sce")
keggtab$gene <- mapIds(org.Sc.sgd.db, keggtab$GeneID, column="GENENAME", keytype="ORF")
mappaths = getKEGGPathwayNames(species="sce")
keggtab$pathname = sapply(keggtab$PathwayID, function(x) substr(mappaths[mappaths$PathwayID==x,"Description"],0, nchar(mappaths[mappaths$PathwayID==x,"Description"])-43))
keggtab = unique(keggtab[,c("pathname","gene")])
countgs = table(keggtab$pathname) # Keep only GS with at least 10 genes
keggtab = keggtab[keggtab$pathname %in% names(countgs)[countgs>=10],]
keggtab_pept = unlist(sapply(rownames(difflip), function(x) keggtab[keggtab$gene %in% strsplit(difflip[x,"gene_name"],";")[[1]],"pathname"]))
peptnames = unlist(sapply(rownames(difflip), function(x) rep(x,sum(keggtab$gene %in% strsplit(difflip[x,"gene_name"],";")[[1]]))))
keggtab_pept = unique(data.frame(gs_name=keggtab_pept, gene_symbol=peptnames, row.names = NULL))

gene2pept = data.frame(gs_name=difflip[,"gene_name"],gene_symbol=rownames(difflip))

clusterdf = data.frame(gene_name=difflip[names(cutgroups),"gene_name"],cluster=cutgroups, peptide = names(cutgroups))

gotest <- compareCluster(peptide~cluster, data=clusterdf, fun="enricher",
                         TERM2GENE = gene2pept, pvalueCutoff=0.05)
dotplot(gotest, showCategory=10)

gotest <- compareCluster(peptide~cluster, data=clusterdf, fun="enricher",
                         TERM2GENE = gotab_pept, pvalueCutoff=0.05,
                         minGSSize=100, maxGSSize=10000)
write.table(gotest@compareClusterResult, "results/GOanalysis_Glc_resLiP.tsv", sep = "\t")
dotplot(gotest, showCategory=10)

gotest <- compareCluster(peptide~cluster, data=clusterdf, fun="enricher",
                         TERM2GENE = keggtab_pept, pvalueCutoff=0.05,
                         minGSSize=100, maxGSSize=10000)
dotplot(gotest, showCategory=10)

# Translation
transdf = abs(meandf - meandf[,"G0"])
gotransdf[gotransdf$cond=="G","fcmedian"] = sapply(c(0,1,2,4,7,10,14,20), function(x) median(transdf[rownames(transdf) %in% gotab_pept[gotab_pept$gs_name=="translation","gene_symbol"],sprintf("G%i",x)], na.rm=T))
gotransdf[gotransdf$cond=="G","fcquantmin"] = sapply(c(0,1,2,4,7,10,14,20), function(x) quantile(transdf[rownames(transdf) %in% gotab_pept[gotab_pept$gs_name=="translation","gene_symbol"],sprintf("G%i",x)], probs=0.25, na.rm=T))
gotransdf[gotransdf$cond=="G","fcquantmax"] = sapply(c(0,1,2,4,7,10,14,20), function(x) quantile(transdf[rownames(transdf) %in% gotab_pept[gotab_pept$gs_name=="translation","gene_symbol"],sprintf("G%i",x)], probs=0.75, na.rm=T))

# At prot level
clusterdf = unique(data.frame(cluster = unlist(sapply(names(cutgroups), function(x) rep(cutgroups[x],length(strsplit(difflip[x,"gene_name"],";")[[1]])))),
                              gene = unlist(sapply(names(cutgroups), function(x) strsplit(difflip[x,"gene_name"],";")[[1]]))))

gotest <- compareCluster(gene~cluster, data=clusterdf, fun="enricher",
                         TERM2GENE = gotab, pvalueCutoff=0.05, universe=includedgenes)
dotplot(gotest, showCategory=10)

gotest <- compareCluster(gene~cluster, data=clusterdf, fun="enricher",
                         TERM2GENE = keggtab, pvalueCutoff=0.05, universe=includedgenes)
dotplot(gotest, showCategory=10)

#### NITROGEN ####
# Load data
difflip = read.csv("results/diffStructurome_N.tsv",sep="\t", row.names = "quantID")
meandf = difflip[,(ncol(difflip)-7):ncol(difflip)]
resLiP = read.csv("results/resLiP.tsv",sep="\t")

# Clustering by correlation
X = resLiP[rownames(difflip)[(difflip$adj.P.Val<0.05)&(!is.na(difflip$adj.P.Val))],grepl("N",colnames(resLiP))]
# Make clusters
cordist = 1-abs(cor(t(X), use = "pairwise", method = "pearson"))
notna = rowSums(is.na(cordist))==0
clust_aggl = hclust(as.dist(cordist[notna,notna]), method = "ward.D2")
cutgroups = cutree(clust_aggl, k=2)
dend = color_branches(as.dendrogram(clust_aggl), k=2)
pdf("figure_panels/heatmap_resLiP_Nit_dendrogram.pdf", width = 6, height = 2.5); plot(dend); dev.off()
colmap = setNames(get_leaves_branches_col(dend), names(cutgroups)[order.dendrogram(dend)])
# Check clusters
time = order(as.numeric(substr(colnames(meandf),2,4)))
heatdf = meandf[rownames(X)[notna],]
signcordf = heatdf  - meandf[rownames(X)[notna],"N0"]
for (x in unique(cutgroups)){
  tempclust = cutree(hclust(as.dist(1-cor(t(heatdf[cutgroups==x,]), use = "pairwise", method = "pearson")), method = "ward.D2"), k=2)
  signcordf[names(tempclust)[tempclust==2],] = -signcordf[names(tempclust)[tempclust==2],]
}
Colors=rev(brewer.pal(11,"RdBu"))
Colors=colorRampPalette(Colors)(100)
heatdf = heatdf  - meandf[rownames(X)[notna],"N0"]
heatdf = heatdf[rowSums(!is.na(heatdf))>0,]; signcordf = signcordf[rowSums(!is.na(signcordf))>0,]
orderedpept = unlist(sapply(unique(cutgroups), function(x) rownames(heatdf[rownames(heatdf) %in% names(cutgroups)[cutgroups==x],])[hclust(dist(heatdf[rownames(heatdf) %in% names(cutgroups)[cutgroups==x],], method="euclidean"), method = "ward.D2")$order]))
pdf("figure_panels/heatmap_resLiP_Nit.pdf", width = 7, height = 5)
heatmap.2(t(heatdf[orderedpept,time[2:8]]), trace="none", Rowv=F, Colv=F, col=Colors, margins = c(4,4),breaks=seq(-2,2,0.04),symkey=F,
          labCol = difflip[orderedpept,"gene_name"], ColSideColors = colmap[orderedpept])
dev.off()
write.table(t(heatdf[orderedpept,time[2:8]]),"source_data/heatmap_resLiP_Nit.tsv", sep = "\t")

# Plot trends
trendsdf = data.frame(row.names = 1:(length(time)*length(unique(cutgroups))))
trendsdf$timegroup = c(sapply(unique(cutgroups), function(x) colnames(signcordf[rownames(signcordf) %in% names(cutgroups)[cutgroups==x],])))
trendsdf$time = apply(trendsdf,1, function(x) as.numeric(substr(x["timegroup"],2,4)))
trendsdf$cluster = as.character(c(sapply(unique(cutgroups),function(x) rep(x,length(time)))))
trendsdf$fcmedian = apply(trendsdf,1, function(x) median(signcordf[rownames(signcordf) %in% names(cutgroups)[cutgroups==x["cluster"]],x["timegroup"]], na.rm=T))
trendsdf$fcquantmin = apply(trendsdf,1, function(x) quantile(signcordf[rownames(signcordf) %in% names(cutgroups)[cutgroups==x["cluster"]],x["timegroup"]], probs=0.25, na.rm=T))
trendsdf$fcquantmax = apply(trendsdf,1, function(x) quantile(signcordf[rownames(signcordf) %in% names(cutgroups)[cutgroups==x["cluster"]],x["timegroup"]], probs=0.75, na.rm=T))
p = ggplot(data=trendsdf, aes(x=time, y=fcmedian, ymin=fcquantmin, ymax=fcquantmax)) +
  facet_wrap(~cluster, ncol = 1) + geom_hline(yintercept=0, color="grey") +
  geom_line(aes(color=cluster)) + 
  geom_ribbon(alpha=0.5, aes(fill=cluster)) +
  scale_color_manual(values=unique(colmap[orderedpept])) + scale_fill_manual(values=unique(colmap[orderedpept])) +
  theme_bw()
ggsave("figure_panels/heatmap_resLiP_Nit_trends.pdf", p, width = 2.5, height = 2.5)

# GO enrichment
includedgenes = unique(unlist(sapply(unique(difflip$gene_name), function(x) strsplit(x,";")[[1]])))
gotab = select(org.Sc.sgd.db, includedgenes, columns=c("GOALL","GENENAME"), keytype="GENENAME")
mapgo = go2term(unique(gotab$GO))
gotab$goterm = as.character(sapply(gotab$GO, function(x) mapgo[mapgo$go_id==x,"Term"]))
gotab = unique(gotab[gotab$ONTOLOGY=="BP",c("goterm","GENENAME")])
countgs = table(gotab$goterm) # Keep only GS with at least 10 genes
gotab = gotab[gotab$goterm %in% names(countgs)[countgs>=10],]
gotab_pept = unlist(sapply(rownames(difflip), function(x) gotab[gotab$GENENAME %in% strsplit(difflip[x,"gene_name"],";")[[1]],"goterm"]))
peptnames = unlist(sapply(rownames(difflip), function(x) rep(x,sum(gotab$GENENAME %in% strsplit(difflip[x,"gene_name"],";")[[1]]))))
gotab_pept = unique(data.frame(gs_name=gotab_pept, gene_symbol=peptnames, row.names = NULL))

keggtab <- getGeneKEGGLinks(species="sce")
keggtab$gene <- mapIds(org.Sc.sgd.db, keggtab$GeneID, column="GENENAME", keytype="ORF")
mappaths = getKEGGPathwayNames(species="sce")
keggtab$pathname = sapply(keggtab$PathwayID, function(x) substr(mappaths[mappaths$PathwayID==x,"Description"],0, nchar(mappaths[mappaths$PathwayID==x,"Description"])-43))
keggtab = unique(keggtab[,c("pathname","gene")])
countgs = table(keggtab$pathname) # Keep only GS with at least 10 genes
keggtab = keggtab[keggtab$pathname %in% names(countgs)[countgs>=10],]
keggtab_pept = unlist(sapply(rownames(difflip), function(x) keggtab[keggtab$gene %in% strsplit(difflip[x,"gene_name"],";")[[1]],"pathname"]))
peptnames = unlist(sapply(rownames(difflip), function(x) rep(x,sum(keggtab$gene %in% strsplit(difflip[x,"gene_name"],";")[[1]]))))
keggtab_pept = unique(data.frame(gs_name=keggtab_pept, gene_symbol=peptnames, row.names = NULL))

gene2pept = data.frame(gs_name=difflip[,"gene_name"],gene_symbol=rownames(difflip))

clusterdf = data.frame(gene_name=difflip[names(cutgroups),"gene_name"],cluster=cutgroups, peptide = names(cutgroups))

gotest <- compareCluster(peptide~cluster, data=clusterdf, fun="enricher",
                         TERM2GENE = gene2pept, pvalueCutoff=0.05)
dotplot(gotest, showCategory=10)

gotest <- compareCluster(peptide~cluster, data=clusterdf, fun="enricher",
                         TERM2GENE = gotab_pept, pvalueCutoff=0.05)
dotplot(gotest, showCategory=10)

gotest <- compareCluster(peptide~cluster, data=clusterdf, fun="enricher",
                         TERM2GENE = keggtab_pept, pvalueCutoff=0.05)
dotplot(gotest, showCategory=10)

# Translation
transdf = abs(meandf - meandf[,"N0"])
gotransdf[gotransdf$cond=="N","fcmedian"] = sapply(c(0,1,2,4,7,10,14,20), function(x) median(transdf[rownames(transdf) %in% gotab_pept[gotab_pept$gs_name=="translation","gene_symbol"],sprintf("N%i",x)], na.rm=T))
gotransdf[gotransdf$cond=="N","fcquantmin"] = sapply(c(0,1,2,4,7,10,14,20), function(x) quantile(transdf[rownames(transdf) %in% gotab_pept[gotab_pept$gs_name=="translation","gene_symbol"],sprintf("N%i",x)], probs=0.25, na.rm=T))
gotransdf[gotransdf$cond=="N","fcquantmax"] = sapply(c(0,1,2,4,7,10,14,20), function(x) quantile(transdf[rownames(transdf) %in% gotab_pept[gotab_pept$gs_name=="translation","gene_symbol"],sprintf("N%i",x)], probs=0.75, na.rm=T))

# At prot level
clusterdf = unique(data.frame(cluster = unlist(sapply(names(cutgroups), function(x) rep(cutgroups[x],length(strsplit(difflip[x,"gene_name"],";")[[1]])))),
                              gene = unlist(sapply(names(cutgroups), function(x) strsplit(difflip[x,"gene_name"],";")[[1]]))))

gotest <- compareCluster(gene~cluster, data=clusterdf, fun="enricher",
                         TERM2GENE = gotab, pvalueCutoff=0.05, universe=includedgenes)
dotplot(gotest, showCategory=10)

gotest <- compareCluster(gene~cluster, data=clusterdf, fun="enricher",
                         TERM2GENE = keggtab, pvalueCutoff=0.05, universe=includedgenes)
dotplot(gotest, showCategory=10)

#### PHOSPHATE ####
# Load data
difflip = read.csv("results/diffStructurome_P.tsv",sep="\t", row.names = "quantID")
meandf = difflip[,(ncol(difflip)-7):ncol(difflip)]
resLiP = read.csv("results/resLiP.tsv",sep="\t")

# Clustering by correlation
X = resLiP[rownames(difflip)[(difflip$adj.P.Val<0.05)&(!is.na(difflip$adj.P.Val))],grepl("P",colnames(resLiP))]
# Make clusters
cordist = 1-abs(cor(t(X), use = "pairwise", method = "pearson"))
notna = rowSums(is.na(cordist))==0
clust_aggl = hclust(as.dist(cordist[notna,notna]), method = "ward.D2")
cutgroups = cutree(clust_aggl, k=3)
dend = color_branches(as.dendrogram(clust_aggl), k=3)
pdf("figure_panels/heatmap_resLiP_Pho_dendrogram.pdf", width = 6, height = 2.5); plot(dend); dev.off()
colmap = setNames(get_leaves_branches_col(dend), names(cutgroups)[order.dendrogram(dend)])
# Check clusters
time = order(as.numeric(substr(colnames(meandf),2,4)))
heatdf = meandf[rownames(X)[notna],]
signcordf = heatdf  - meandf[rownames(X)[notna],"P0"]
for (x in unique(cutgroups)){
  tempclust = cutree(hclust(as.dist(1-cor(t(heatdf[cutgroups==x,]), use = "pairwise", method = "pearson")), method = "ward.D2"), k=2)
  signcordf[names(tempclust)[tempclust==2],] = -signcordf[names(tempclust)[tempclust==2],]
}
Colors=rev(brewer.pal(11,"RdBu"))
Colors=colorRampPalette(Colors)(100)
heatdf = heatdf  - meandf[rownames(X)[notna],"P0"]
heatdf = heatdf[rowSums(!is.na(heatdf))>0,]; signcordf = signcordf[rowSums(!is.na(signcordf))>0,]
orderedpept = unlist(sapply(unique(cutgroups), function(x) rownames(heatdf[rownames(heatdf) %in% names(cutgroups)[cutgroups==x],])[hclust(dist(heatdf[rownames(heatdf) %in% names(cutgroups)[cutgroups==x],], method="euclidean"), method = "ward.D2")$order]))
pdf("figure_panels/heatmap_resLiP_Pho.pdf", width = 7, height = 5)
heatmap.2(t(heatdf[orderedpept,time[2:8]]), trace="none", Rowv=F, Colv=F, col=Colors, margins = c(4,4),breaks=seq(-2,2,0.04),symkey=F,
          labCol = difflip[orderedpept,"gene_name"], ColSideColors = colmap[orderedpept])
dev.off()
write.table(t(heatdf[orderedpept,time[2:8]]),"source_data/heatmap_resLiP_Pho.tsv", sep = "\t")

# Plot trends
trendsdf = data.frame(row.names = 1:(length(time)*length(unique(cutgroups))))
trendsdf$timegroup = c(sapply(unique(cutgroups), function(x) colnames(signcordf[rownames(signcordf) %in% names(cutgroups)[cutgroups==x],])))
trendsdf$time = apply(trendsdf,1, function(x) as.numeric(substr(x["timegroup"],2,4)))
trendsdf$cluster = as.character(c(sapply(unique(cutgroups),function(x) rep(x,length(time)))))
trendsdf$fcmedian = apply(trendsdf,1, function(x) median(signcordf[rownames(signcordf) %in% names(cutgroups)[cutgroups==x["cluster"]],x["timegroup"]], na.rm=T))
trendsdf$fcquantmin = apply(trendsdf,1, function(x) quantile(signcordf[rownames(signcordf) %in% names(cutgroups)[cutgroups==x["cluster"]],x["timegroup"]], probs=0.25, na.rm=T))
trendsdf$fcquantmax = apply(trendsdf,1, function(x) quantile(signcordf[rownames(signcordf) %in% names(cutgroups)[cutgroups==x["cluster"]],x["timegroup"]], probs=0.75, na.rm=T))
p = ggplot(data=trendsdf, aes(x=time, y=fcmedian, ymin=fcquantmin, ymax=fcquantmax)) +
  facet_wrap(~cluster, ncol = 1) + geom_hline(yintercept=0, color="grey") +
  geom_line(aes(color=cluster)) + 
  geom_ribbon(alpha=0.5, aes(fill=cluster)) +
  scale_color_manual(values=unique(colmap[orderedpept])) + scale_fill_manual(values=unique(colmap[orderedpept])) +
  theme_bw()
ggsave("figure_panels/heatmap_resLiP_Pho_trends.pdf", p, width = 2.5, height = 3)

# GO enrichment
includedgenes = unique(unlist(sapply(unique(difflip$gene_name), function(x) strsplit(x,";")[[1]])))
gotab = select(org.Sc.sgd.db, includedgenes, columns=c("GOALL","GENENAME"), keytype="GENENAME")
mapgo = go2term(unique(gotab$GO))
gotab$goterm = as.character(sapply(gotab$GO, function(x) mapgo[mapgo$go_id==x,"Term"]))
gotab = unique(gotab[gotab$ONTOLOGY=="BP",c("goterm","GENENAME")])
countgs = table(gotab$goterm) # Keep only GS with at least 10 genes
gotab = gotab[gotab$goterm %in% names(countgs)[countgs>=10],]
gotab_pept = unlist(sapply(rownames(difflip), function(x) gotab[gotab$GENENAME %in% strsplit(difflip[x,"gene_name"],";")[[1]],"goterm"]))
peptnames = unlist(sapply(rownames(difflip), function(x) rep(x,sum(gotab$GENENAME %in% strsplit(difflip[x,"gene_name"],";")[[1]]))))
gotab_pept = unique(data.frame(gs_name=gotab_pept, gene_symbol=peptnames, row.names = NULL))

keggtab <- getGeneKEGGLinks(species="sce")
keggtab$gene <- mapIds(org.Sc.sgd.db, keggtab$GeneID, column="GENENAME", keytype="ORF")
mappaths = getKEGGPathwayNames(species="sce")
keggtab$pathname = sapply(keggtab$PathwayID, function(x) substr(mappaths[mappaths$PathwayID==x,"Description"],0, nchar(mappaths[mappaths$PathwayID==x,"Description"])-43))
keggtab = unique(keggtab[,c("pathname","gene")])
countgs = table(keggtab$pathname) # Keep only GS with at least 10 genes
keggtab = keggtab[keggtab$pathname %in% names(countgs)[countgs>=10],]
keggtab_pept = unlist(sapply(rownames(difflip), function(x) keggtab[keggtab$gene %in% strsplit(difflip[x,"gene_name"],";")[[1]],"pathname"]))
peptnames = unlist(sapply(rownames(difflip), function(x) rep(x,sum(keggtab$gene %in% strsplit(difflip[x,"gene_name"],";")[[1]]))))
keggtab_pept = unique(data.frame(gs_name=keggtab_pept, gene_symbol=peptnames, row.names = NULL))

gene2pept = data.frame(gs_name=difflip[,"gene_name"],gene_symbol=rownames(difflip))

clusterdf = data.frame(gene_name=difflip[names(cutgroups),"gene_name"],cluster=cutgroups, peptide = names(cutgroups))

gotest <- compareCluster(peptide~cluster, data=clusterdf, fun="enricher",
                         TERM2GENE = gene2pept, pvalueCutoff=0.05)
dotplot(gotest, showCategory=10)

gotest <- compareCluster(peptide~cluster, data=clusterdf, fun="enricher",
                         TERM2GENE = gotab_pept, pvalueCutoff=0.05)
write.table(gotest@compareClusterResult, "results/GOanalysis_Pho_resLiP.tsv", sep = "\t")
dotplot(gotest, showCategory=10)

gotest <- compareCluster(peptide~cluster, data=clusterdf, fun="enricher",
                         TERM2GENE = keggtab_pept, pvalueCutoff=0.05)
dotplot(gotest, showCategory=10)

# Translation
transdf = abs(meandf - meandf[,"P0"])
gotransdf[gotransdf$cond=="P","fcmedian"] = sapply(c(0,1,2,4,7,10,14,20), function(x) median(transdf[rownames(transdf) %in% gotab_pept[gotab_pept$gs_name=="translation","gene_symbol"],sprintf("P%i",x)], na.rm=T))
gotransdf[gotransdf$cond=="P","fcquantmin"] = sapply(c(0,1,2,4,7,10,14,20), function(x) quantile(transdf[rownames(transdf) %in% gotab_pept[gotab_pept$gs_name=="translation","gene_symbol"],sprintf("P%i",x)], probs=0.25, na.rm=T))
gotransdf[gotransdf$cond=="P","fcquantmax"] = sapply(c(0,1,2,4,7,10,14,20), function(x) quantile(transdf[rownames(transdf) %in% gotab_pept[gotab_pept$gs_name=="translation","gene_symbol"],sprintf("P%i",x)], probs=0.75, na.rm=T))

# At prot level
clusterdf = unique(data.frame(cluster = unlist(sapply(names(cutgroups), function(x) rep(cutgroups[x],length(strsplit(difflip[x,"gene_name"],";")[[1]])))),
                              gene = unlist(sapply(names(cutgroups), function(x) strsplit(difflip[x,"gene_name"],";")[[1]]))))

gotest <- compareCluster(gene~cluster, data=clusterdf, fun="enricher",
                         TERM2GENE = gotab, pvalueCutoff=0.05, universe=includedgenes)
dotplot(gotest, showCategory=10)

gotest <- compareCluster(gene~cluster, data=clusterdf, fun="enricher",
                         TERM2GENE = keggtab, pvalueCutoff=0.05, universe=includedgenes)
dotplot(gotest, showCategory=10)

##### TRANSLATION #####
ggplot(data=gotransdf, aes(x=time, y=fcmedian, ymin=fcquantmin, ymax=fcquantmax)) +
  facet_wrap(~cond, ncol = 1) + geom_hline(yintercept=0, color="grey") +
  geom_line(aes(color=cond)) + 
  geom_ribbon(alpha=0.5, aes(fill=cond)) +
  scale_color_manual(values=c("#ca3d3d","#919128","#347cbf")) + scale_fill_manual(values=c("#f4d7d7","#f4f4d7","#d7e6f4")) +
  theme_bw()
p = ggplot(data=gotransdf, aes(x=time, y=fcmedian, color=cond)) +
  geom_hline(yintercept=0, color="grey") +
  geom_line(aes(color=cond), linewidth=1) + 
  scale_color_manual(values=c("#ca3d3d","#919128","#347cbf")) +
  coord_cartesian(ylim = c(0, 0.3)) +
  theme_bw()
ggsave("figure_panels/lineplot_resLiP_GOtrans.pdf", p, width = 4, height = 2.5)
