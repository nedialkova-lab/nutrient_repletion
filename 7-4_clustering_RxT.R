# ==============================================================================
# 7-4_clustering_RxT.R
# Clustering and GO enrichment of transcripts with significantly changing ribosomes
# per transcript.
#
# Manuscript: Fig. 2C; Extended Data Fig. 3B-C
# Requires  : outputs of 7-3_spikein_normalization.R
# Inputs    : results/diffRxT_{G,N,P}.tsv
#             results/RxT.tsv
# Outputs   : results/GOanalysis_{Glc,Nit,Pho}_RxT.tsv
#             source_data/heatmap_RxT_{Glc,Nit,Pho}.tsv
#             figure_panels/heatmap_RxT_{Glc,Nit,Pho}{,_dendrogram,_trends}.pdf
# Run from the repository root (all paths are relative to it).
# ==============================================================================

library(gplots)
library(cluster)
library(ggplot2)
library(RColorBrewer)
library(clusterProfiler)
library(org.Sc.sgd.db)
library(stringr)
library(UniProt.ws)
library(dendextend)
library(limma)

#### GLUCOSE ####
# Load data
diffrxt = read.csv("results/diffRxT_G.tsv",sep="\t", row.names = 1)
meandf = diffrxt[,(ncol(diffrxt)-6):(ncol(diffrxt)-1)]
normdf = read.csv("results/RxT.tsv",sep="\t")

# Clustering by correlation
X = normdf[rownames(diffrxt)[(diffrxt$adj.P.Val<0.05)&(!is.na(diffrxt$adj.P.Val))],grepl("g",colnames(normdf))]
# Make clusters
cordist = 1-cor(t(X), use = "pairwise", method = "pearson")
notna = rowSums(is.na(cordist))==0
clust_aggl = hclust(as.dist(cordist[notna,notna]), method = "ward.D2")
cutgroups = cutree(clust_aggl, k=4)
dend = color_branches(as.dendrogram(clust_aggl), k=4)
pdf("figure_panels/heatmap_RxT_Glc_dendrogram.pdf", width = 6, height = 2.5); plot(dend); dev.off()
colmap = setNames(get_leaves_branches_col(dend), names(cutgroups)[order.dendrogram(dend)])
# Check clusters
time = order(as.numeric(substr(colnames(meandf),2,4)))
heatdf = log2(meandf[rownames(X)[notna],] / meandf[rownames(X)[notna],"g0"])
heatdf = heatdf[(rowSums(!is.na(heatdf))>0)&(rowSums(is.infinite(as.matrix(heatdf)))==0),]
Colors=rev(brewer.pal(11,"RdBu"))
Colors=colorRampPalette(Colors)(100)
orderedgene = unlist(sapply(unique(cutgroups), function(x) rownames(heatdf[rownames(heatdf) %in% names(cutgroups)[cutgroups==x],])[hclust(dist(heatdf[rownames(heatdf) %in% names(cutgroups)[cutgroups==x],], method="euclidean"), method = "ward.D2")$order]))
pdf("figure_panels/heatmap_RxT_Glc.pdf", width = 7, height = 5)
heatmap.2(t(heatdf[orderedgene,time[2:6]]), trace="none", Rowv=F, Colv=F, col=Colors, margins = c(6,4),breaks=seq(-3,3,0.06),symkey=F,
          labCol = diffrxt[orderedgene,"gene"], ColSideColors = colmap[orderedgene])
dev.off()
write.table(t(heatdf[orderedgene,time[2:6]]),"source_data/heatmap_RxT_Glc.tsv", sep = "\t")

# Plot trends
trendsdf = data.frame(row.names = 1:(length(time)*length(unique(cutgroups))))
trendsdf$timegroup = c(sapply(unique(cutgroups), function(x) colnames(heatdf[rownames(heatdf) %in% names(cutgroups)[cutgroups==x],])))
trendsdf$time = apply(trendsdf,1, function(x) as.numeric(substr(x["timegroup"],2,4)))
trendsdf$cluster = as.character(c(sapply(unique(cutgroups),function(x) rep(x,length(time)))))
trendsdf$fcmedian = apply(trendsdf,1, function(x) median(heatdf[rownames(heatdf) %in% names(cutgroups)[cutgroups==x["cluster"]],x["timegroup"]], na.rm=T))
trendsdf$fcquantmin = apply(trendsdf,1, function(x) quantile(heatdf[rownames(heatdf) %in% names(cutgroups)[cutgroups==x["cluster"]],x["timegroup"]], probs=0.25, na.rm=T))
trendsdf$fcquantmax = apply(trendsdf,1, function(x) quantile(heatdf[rownames(heatdf) %in% names(cutgroups)[cutgroups==x["cluster"]],x["timegroup"]], probs=0.75, na.rm=T))
p = ggplot(data=trendsdf, aes(x=time, y=fcmedian, ymin=fcquantmin, ymax=fcquantmax)) +
  facet_wrap(~cluster, ncol = 1) + geom_hline(yintercept=0, color="grey") +
  geom_line(aes(color=cluster)) + 
  geom_ribbon(alpha=0.5, aes(fill=cluster)) +
  scale_color_manual(values=unique(colmap[orderedgene])) + scale_fill_manual(values=unique(colmap[orderedgene])) +
  theme_bw()
ggsave("figure_panels/heatmap_RxT_Glc_trends.pdf", p, width = 2.5, height = 4)


# GO enrichment
includedgenes = rownames(diffrxt)
gotab = select(org.Sc.sgd.db, includedgenes, columns=c("GOALL","ORF"), keytype="ORF")
mapgo = go2term(unique(gotab$GO))
gotab$goterm = as.character(sapply(gotab$GO, function(x) mapgo[mapgo$go_id==x,"Term"]))
gotab = unique(gotab[gotab$ONTOLOGY=="BP",c("goterm","ORF")])

keggtab <- getGeneKEGGLinks(species="sce")
mappaths = getKEGGPathwayNames(species="sce")
keggtab$pathname = sapply(keggtab$PathwayID, function(x) substr(mappaths[mappaths$PathwayID==x,"Description"],0, nchar(mappaths[mappaths$PathwayID==x,"Description"])-43))
keggtab = unique(keggtab[,c("pathname","GeneID")])

clusterdf = data.frame(cluster = cutgroups,
                       gene = names(cutgroups))

gotest <- compareCluster(gene~cluster, data=clusterdf, fun="enricher",
                         TERM2GENE = gotab, pvalueCutoff=0.05, universe = includedgenes)
dotplot(gotest, showCategory=10)
write.table(gotest@compareClusterResult, "results/GOanalysis_Glc_RxT.tsv", sep = "\t")
gotest <- compareCluster(gene~cluster, data=clusterdf, fun="enricher",
                         TERM2GENE = keggtab, pvalueCutoff=0.05, universe = includedgenes)
dotplot(gotest, showCategory=10)

#### NITROGEN ####
# Load data
diffrxt = read.csv("results/diffRxT_N.tsv",sep="\t", row.names = 1)
meandf = diffrxt[,(ncol(diffrxt)-6):(ncol(diffrxt)-1)]
normdf = read.csv("results/RxT.tsv",sep="\t")

# Clustering by correlation
X = normdf[rownames(diffrxt)[(diffrxt$adj.P.Val<0.05)&(!is.na(diffrxt$adj.P.Val))],grepl("n",colnames(normdf))]
# Make clusters
cordist = 1-cor(t(X), use = "pairwise", method = "pearson")
notna = rowSums(is.na(cordist))==0
clust_aggl = hclust(as.dist(cordist[notna,notna]), method = "ward.D2")
cutgroups = cutree(clust_aggl, k=1)
dend = color_branches(as.dendrogram(clust_aggl), k=1)
pdf("figure_panels/heatmap_RxT_Nit_dendrogram.pdf", width = 6, height = 2.5); plot(dend); dev.off()
colmap = setNames(get_leaves_branches_col(dend), names(cutgroups)[order.dendrogram(dend)])
# Check clusters
time = order(as.numeric(substr(colnames(meandf),2,4)))
heatdf = log2(meandf[rownames(X)[notna],] / meandf[rownames(X)[notna],"n0"])
heatdf = heatdf[(rowSums(!is.na(heatdf))>0)&(rowSums(is.infinite(as.matrix(heatdf)))==0),]
Colors=rev(brewer.pal(11,"RdBu"))
Colors=colorRampPalette(Colors)(100)
orderedgene = unlist(sapply(unique(cutgroups), function(x) rownames(heatdf[rownames(heatdf) %in% names(cutgroups)[cutgroups==x],])[hclust(dist(heatdf[rownames(heatdf) %in% names(cutgroups)[cutgroups==x],], method="euclidean"), method = "ward.D2")$order]))
pdf("figure_panels/heatmap_RxT_Nit.pdf", width = 7, height = 5)
heatmap.2(t(heatdf[orderedgene,time[2:6]]), trace="none", Rowv=F, Colv=F, col=Colors, margins = c(6,4),breaks=seq(-3,3,0.06),symkey=F,
          labCol = diffrxt[orderedgene,"gene"], ColSideColors = colmap[orderedgene])
dev.off()
write.table(t(heatdf[orderedgene,time[2:6]]),"source_data/heatmap_RxT_Nit.tsv", sep = "\t")

# Plot trends
trendsdf = data.frame(row.names = 1:(length(time)*length(unique(cutgroups))))
trendsdf$timegroup = c(sapply(unique(cutgroups), function(x) colnames(heatdf[rownames(heatdf) %in% names(cutgroups)[cutgroups==x],])))
trendsdf$time = apply(trendsdf,1, function(x) as.numeric(substr(x["timegroup"],2,4)))
trendsdf$cluster = as.character(c(sapply(unique(cutgroups),function(x) rep(x,length(time)))))
trendsdf$fcmedian = apply(trendsdf,1, function(x) median(heatdf[rownames(heatdf) %in% names(cutgroups)[cutgroups==x["cluster"]],x["timegroup"]], na.rm=T))
trendsdf$fcquantmin = apply(trendsdf,1, function(x) quantile(heatdf[rownames(heatdf) %in% names(cutgroups)[cutgroups==x["cluster"]],x["timegroup"]], probs=0.25, na.rm=T))
trendsdf$fcquantmax = apply(trendsdf,1, function(x) quantile(heatdf[rownames(heatdf) %in% names(cutgroups)[cutgroups==x["cluster"]],x["timegroup"]], probs=0.75, na.rm=T))
p = ggplot(data=trendsdf, aes(x=time, y=fcmedian, ymin=fcquantmin, ymax=fcquantmax)) +
  facet_wrap(~cluster, ncol = 1) + geom_hline(yintercept=0, color="grey") +
  geom_line(aes(color=cluster)) + 
  geom_ribbon(alpha=0.5, aes(fill=cluster)) +
  scale_color_manual(values=unique(colmap[orderedgene])) + scale_fill_manual(values=unique(colmap[orderedgene])) +
  theme_bw()
ggsave("figure_panels/heatmap_RxT_Nit_trends.pdf", p, width = 2.5, height = 1.5)

# GO enrichment
includedgenes = rownames(diffrxt)
gotab = select(org.Sc.sgd.db, includedgenes, columns=c("GOALL","ORF"), keytype="ORF")
mapgo = go2term(unique(gotab$GO))
gotab$goterm = as.character(sapply(gotab$GO, function(x) mapgo[mapgo$go_id==x,"Term"]))
gotab = unique(gotab[gotab$ONTOLOGY=="BP",c("goterm","ORF")])

keggtab <- getGeneKEGGLinks(species="sce")
mappaths = getKEGGPathwayNames(species="sce")
keggtab$pathname = sapply(keggtab$PathwayID, function(x) substr(mappaths[mappaths$PathwayID==x,"Description"],0, nchar(mappaths[mappaths$PathwayID==x,"Description"])-43))
keggtab = unique(keggtab[,c("pathname","GeneID")])

clusterdf = data.frame(cluster = cutgroups,
                       gene = names(cutgroups))

gotest <- compareCluster(gene~cluster, data=clusterdf, fun="enricher",
                         TERM2GENE = gotab, pvalueCutoff=0.05, universe = includedgenes)
dotplot(gotest, showCategory=10)
write.table(gotest@compareClusterResult, "results/GOanalysis_Nit_RxT.tsv", sep = "\t")
gotest <- compareCluster(gene~cluster, data=clusterdf, fun="enricher",
                         TERM2GENE = keggtab, pvalueCutoff=0.05, universe = includedgenes)
dotplot(gotest, showCategory=10)

#### PHOSPHATE ####
# Load data
diffrxt = read.csv("results/diffRxT_P.tsv",sep="\t", row.names = 1)
meandf = diffrxt[,(ncol(diffrxt)-6):(ncol(diffrxt)-1)]
normdf = read.csv("results/RxT.tsv",sep="\t")

# Clustering by correlation
X = normdf[rownames(diffrxt)[(diffrxt$adj.P.Val<0.05)&(!is.na(diffrxt$adj.P.Val))],grepl("^p",colnames(normdf))]
# Make clusters
cordist = 1-cor(t(X), use = "pairwise", method = "pearson")
notna = rowSums(is.na(cordist))==0
clust_aggl = hclust(as.dist(cordist[notna,notna]), method = "ward.D2")
cutgroups = cutree(clust_aggl, k=2)
dend = color_branches(as.dendrogram(clust_aggl), k=2)
pdf("figure_panels/heatmap_RxT_Pho_dendrogram.pdf", width = 6, height = 2.5); plot(dend); dev.off()
colmap = setNames(get_leaves_branches_col(dend), names(cutgroups)[order.dendrogram(dend)])
# Check clusters
time = order(as.numeric(substr(colnames(meandf),2,4)))
heatdf = log2(meandf[rownames(X)[notna],] / meandf[rownames(X)[notna],"p0"])
heatdf = heatdf[(rowSums(!is.na(heatdf))>0)&(rowSums(is.infinite(as.matrix(heatdf)))==0),]
Colors=rev(brewer.pal(11,"RdBu"))
Colors=colorRampPalette(Colors)(100)
orderedgene = unlist(sapply(unique(cutgroups), function(x) rownames(heatdf[rownames(heatdf) %in% names(cutgroups)[cutgroups==x],])[hclust(dist(heatdf[rownames(heatdf) %in% names(cutgroups)[cutgroups==x],], method="euclidean"), method = "ward.D2")$order]))
pdf("figure_panels/heatmap_RxT_Pho.pdf", width = 7, height = 5)
heatmap.2(t(heatdf[orderedgene,time[2:6]]), trace="none", Rowv=F, Colv=F, col=Colors, margins = c(6,4),breaks=seq(-3,3,0.06),symkey=F,
          labCol = diffrxt[orderedgene,"gene"], ColSideColors = colmap[orderedgene])
dev.off()
write.table(t(heatdf[orderedgene,time[2:6]]),"source_data/heatmap_RxT_Pho.tsv", sep = "\t")

# Plot trends
trendsdf = data.frame(row.names = 1:(length(time)*length(unique(cutgroups))))
trendsdf$timegroup = c(sapply(unique(cutgroups), function(x) colnames(heatdf[rownames(heatdf) %in% names(cutgroups)[cutgroups==x],])))
trendsdf$time = apply(trendsdf,1, function(x) as.numeric(substr(x["timegroup"],2,4)))
trendsdf$cluster = as.character(c(sapply(unique(cutgroups),function(x) rep(x,length(time)))))
trendsdf$fcmedian = apply(trendsdf,1, function(x) median(heatdf[rownames(heatdf) %in% names(cutgroups)[cutgroups==x["cluster"]],x["timegroup"]], na.rm=T))
trendsdf$fcquantmin = apply(trendsdf,1, function(x) quantile(heatdf[rownames(heatdf) %in% names(cutgroups)[cutgroups==x["cluster"]],x["timegroup"]], probs=0.25, na.rm=T))
trendsdf$fcquantmax = apply(trendsdf,1, function(x) quantile(heatdf[rownames(heatdf) %in% names(cutgroups)[cutgroups==x["cluster"]],x["timegroup"]], probs=0.75, na.rm=T))
p = ggplot(data=trendsdf, aes(x=time, y=fcmedian, ymin=fcquantmin, ymax=fcquantmax)) +
  facet_wrap(~cluster, ncol = 1) + geom_hline(yintercept=0, color="grey") +
  geom_line(aes(color=cluster)) + 
  geom_ribbon(alpha=0.5, aes(fill=cluster)) +
  scale_color_manual(values=unique(colmap[orderedgene])) + scale_fill_manual(values=unique(colmap[orderedgene])) +
  theme_bw()
ggsave("figure_panels/heatmap_RxT_Pho_trends.pdf", p, width = 2.5, height = 2.5)

# GO enrichment
includedgenes = rownames(diffrxt)
gotab = select(org.Sc.sgd.db, includedgenes, columns=c("GOALL","ORF"), keytype="ORF")
mapgo = go2term(unique(gotab$GO))
gotab$goterm = as.character(sapply(gotab$GO, function(x) mapgo[mapgo$go_id==x,"Term"]))
gotab = unique(gotab[gotab$ONTOLOGY=="BP",c("goterm","ORF")])

keggtab <- getGeneKEGGLinks(species="sce")
mappaths = getKEGGPathwayNames(species="sce")
keggtab$pathname = sapply(keggtab$PathwayID, function(x) substr(mappaths[mappaths$PathwayID==x,"Description"],0, nchar(mappaths[mappaths$PathwayID==x,"Description"])-43))
keggtab = unique(keggtab[,c("pathname","GeneID")])

clusterdf = data.frame(cluster = cutgroups,
                       gene = names(cutgroups))

gotest <- compareCluster(gene~cluster, data=clusterdf, fun="enricher",
                         TERM2GENE = gotab, pvalueCutoff=0.05, universe = includedgenes)
dotplot(gotest, showCategory=10)
write.table(gotest@compareClusterResult, "results/GOanalysis_Pho_RxT.tsv", sep = "\t")
gotest <- compareCluster(gene~cluster, data=clusterdf, fun="enricher",
                         TERM2GENE = keggtab, pvalueCutoff=0.05, universe = includedgenes)
dotplot(gotest, showCategory=10)
