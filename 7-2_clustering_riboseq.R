# ==============================================================================
# 7-2_clustering_riboseq.R
# Clustering and GO enrichment of differentially changing ribosome footprints.
#
# Manuscript: (Supplementary Table only; heatmaps not shown in the manuscript)
# Requires  : outputs of 7-1_explore_riboseq.R
# Inputs    : results/diffRPF_{G,N,P}.tsv
#             results/TPM_RPF.tsv
# Outputs   : results/GOanalysis_{Glc,Nit,Pho}_riboseq.tsv
#             source_data/heatmap_riboseq_{Glc,Nit,Pho}.tsv
#             figure_panels/heatmap_riboseq_{Glc,Nit,Pho}{,_dendrogram,_trends}.pdf
# Run from the repository root (all paths are relative to it).
# ==============================================================================

library(edgeR)
library(gplots)
library(cluster)
library(ggplot2)
library(RColorBrewer)
library(clusterProfiler)
library(org.Sc.sgd.db)
library(stringr)
library(UniProt.ws)
library(dendextend)

#### GLUCOSE ####
# Load data
diffrpf = read.csv("results/diffRPF_G.tsv",sep="\t", row.names = 1)
meandf = diffrpf[,(ncol(diffrpf)-6):(ncol(diffrpf)-1)]
normdf = read.csv("results/TPM_RPF.tsv",sep="\t")

# Clustering by correlation
X = normdf[rownames(diffrpf)[(diffrpf$adj.P.Val<0.05)&(!is.na(diffrpf$adj.P.Val))],grepl("g",colnames(normdf))]
# Make clusters
cordist = 1-cor(t(X), use = "pairwise", method = "pearson")
notna = rowSums(is.na(cordist))==0
clust_aggl = hclust(as.dist(cordist[notna,notna]), method = "ward.D2")
cutgroups = cutree(clust_aggl, k=3)
dend = color_branches(as.dendrogram(clust_aggl), k=3)
pdf("figure_panels/heatmap_riboseq_Glc_dendrogram.pdf", width = 6, height = 2.5); plot(dend); dev.off()
colmap = setNames(get_leaves_branches_col(dend), names(cutgroups)[order.dendrogram(dend)])
# Check clusters
time = order(as.numeric(substr(colnames(meandf),2,4)))
heatdf = log2(meandf[rownames(X)[notna],] / meandf[rownames(X)[notna],"g0"])
heatdf = heatdf[(rowSums(!is.na(heatdf))>0)&(rowSums(is.infinite(as.matrix(heatdf)))==0),]
Colors=rev(brewer.pal(11,"RdBu"))
Colors=colorRampPalette(Colors)(100)
orderedgene = unlist(sapply(unique(cutgroups), function(x) rownames(heatdf[rownames(heatdf) %in% names(cutgroups)[cutgroups==x],])[hclust(dist(heatdf[rownames(heatdf) %in% names(cutgroups)[cutgroups==x],], method="euclidean"), method = "ward.D2")$order]))
pdf("figure_panels/heatmap_riboseq_Glc.pdf", width = 7, height = 5)
heatmap.2(t(heatdf[orderedgene,time[2:6]]), trace="none", Rowv=F, Colv=F, col=Colors, margins = c(6,4),breaks=seq(-3,3,0.06),symkey=F,
          labCol = diffrpf[orderedgene,"gene"], ColSideColors = colmap[orderedgene])
dev.off()
write.table(t(heatdf[orderedgene,time[2:6]]),"source_data/heatmap_riboseq_Glc.tsv", sep = "\t")

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
ggsave("figure_panels/heatmap_riboseq_Glc_trends.pdf", p, width = 2.5, height = 3)

# GO enrichment
includedgenes = rownames(diffrpf)
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
write.table(gotest@compareClusterResult, "results/GOanalysis_Glc_riboseq.tsv", sep = "\t")
dotplot(gotest, showCategory=10)
gotest <- compareCluster(gene~cluster, data=clusterdf, fun="enricher",
                         TERM2GENE = keggtab, pvalueCutoff=0.05, universe = includedgenes)
dotplot(gotest, showCategory=10)

#### NITROGEN ####
# Load data
diffrpf = read.csv("results/diffRPF_N.tsv",sep="\t", row.names = 1)
meandf = diffrpf[,(ncol(diffrpf)-6):(ncol(diffrpf)-1)]
normdf = read.csv("results/TPM_RPF.tsv",sep="\t")

# Clustering by correlation
X = normdf[rownames(diffrpf)[(diffrpf$adj.P.Val<0.05)&(!is.na(diffrpf$adj.P.Val))],grepl("n",colnames(normdf))]
# Make clusters
cordist = 1-cor(t(X), use = "pairwise", method = "pearson")
notna = rowSums(is.na(cordist))==0
clust_aggl = hclust(as.dist(cordist[notna,notna]), method = "ward.D2")
cutgroups = cutree(clust_aggl, k=3)
dend = color_branches(as.dendrogram(clust_aggl), k=3)
pdf("figure_panels/heatmap_riboseq_Nit_dendrogram.pdf", width = 6, height = 2.5); plot(dend); dev.off()
colmap = setNames(get_leaves_branches_col(dend), names(cutgroups)[order.dendrogram(dend)])
# Check clusters
time = order(as.numeric(substr(colnames(meandf),2,4)))
heatdf = log2(meandf[rownames(X)[notna],] / meandf[rownames(X)[notna],"n0"])
heatdf = heatdf[(rowSums(!is.na(heatdf))>0)&(rowSums(is.infinite(as.matrix(heatdf)))==0),]
Colors=rev(brewer.pal(11,"RdBu"))
Colors=colorRampPalette(Colors)(100)
orderedgene = unlist(sapply(unique(cutgroups), function(x) rownames(heatdf[rownames(heatdf) %in% names(cutgroups)[cutgroups==x],])[hclust(dist(heatdf[rownames(heatdf) %in% names(cutgroups)[cutgroups==x],], method="euclidean"), method = "ward.D2")$order]))
pdf("figure_panels/heatmap_riboseq_Nit.pdf", width = 7, height = 5)
heatmap.2(t(heatdf[orderedgene,time[2:6]]), trace="none", Rowv=F, Colv=F, col=Colors, margins = c(6,4),breaks=seq(-3,3,0.06),symkey=F,
          labCol = diffrpf[orderedgene,"gene"], ColSideColors = colmap[orderedgene])
dev.off()
write.table(t(heatdf[orderedgene,time[2:6]]),"source_data/heatmap_riboseq_Nit.tsv", sep = "\t")

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
ggsave("figure_panels/heatmap_riboseq_Nit_trends.pdf", p, width = 2.5, height = 3)

# GO enrichment
includedgenes = rownames(diffrpf)
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
write.table(gotest@compareClusterResult, "results/GOanalysis_Nit_riboseq.tsv", sep = "\t")
dotplot(gotest, showCategory=10)
gotest <- compareCluster(gene~cluster, data=clusterdf, fun="enricher",
                         TERM2GENE = keggtab, pvalueCutoff=0.05, universe = includedgenes)
dotplot(gotest, showCategory=10)

#### PHOSPHATE ####
# Load data
diffrpf = read.csv("results/diffRPF_P.tsv",sep="\t", row.names = 1)
meandf = diffrpf[,(ncol(diffrpf)-6):(ncol(diffrpf)-1)]
normdf = read.csv("results/TPM_RPF.tsv",sep="\t")

# Clustering by correlation
X = normdf[rownames(diffrpf)[(diffrpf$adj.P.Val<0.05)&(!is.na(diffrpf$adj.P.Val))],grepl("^p",colnames(normdf))]
# Make clusters
cordist = 1-cor(t(X), use = "pairwise", method = "pearson")
notna = rowSums(is.na(cordist))==0
clust_aggl = hclust(as.dist(cordist[notna,notna]), method = "ward.D2")
cutgroups = cutree(clust_aggl, k=2)
dend = color_branches(as.dendrogram(clust_aggl), k=2)
pdf("figure_panels/heatmap_riboseq_Pho_dendrogram.pdf", width = 6, height = 2.5); plot(dend); dev.off()
colmap = setNames(get_leaves_branches_col(dend), names(cutgroups)[order.dendrogram(dend)])
# Check clusters
time = order(as.numeric(substr(colnames(meandf),2,4)))
heatdf = log2(meandf[rownames(X)[notna],] / meandf[rownames(X)[notna],"p0"])
heatdf = heatdf[(rowSums(!is.na(heatdf))>0)&(rowSums(is.infinite(as.matrix(heatdf)))==0),]
Colors=rev(brewer.pal(11,"RdBu"))
Colors=colorRampPalette(Colors)(100)
orderedgene = unlist(sapply(unique(cutgroups), function(x) rownames(heatdf[rownames(heatdf) %in% names(cutgroups)[cutgroups==x],])[hclust(dist(heatdf[rownames(heatdf) %in% names(cutgroups)[cutgroups==x],], method="euclidean"), method = "ward.D2")$order]))
pdf("figure_panels/heatmap_riboseq_Pho.pdf", width = 7, height = 5)
heatmap.2(t(heatdf[orderedgene,time[2:6]]), trace="none", Rowv=F, Colv=F, col=Colors, margins = c(6,4),breaks=seq(-3,3,0.06),symkey=F,
          labCol = diffrpf[orderedgene,"gene"], ColSideColors = colmap[orderedgene])
dev.off()
write.table(t(heatdf[orderedgene,time[2:6]]),"source_data/heatmap_riboseq_Pho.tsv", sep = "\t")

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
ggsave("figure_panels/heatmap_riboseq_Pho_trends.pdf", p, width = 2.5, height = 2.5)

# GO enrichment
includedgenes = rownames(diffrpf)
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
write.table(gotest@compareClusterResult, "results/GOanalysis_Pho_riboseq.tsv", sep = "\t")
dotplot(gotest, showCategory=10)
gotest <- compareCluster(gene~cluster, data=clusterdf, fun="enricher",
                         TERM2GENE = keggtab, pvalueCutoff=0.05, universe = includedgenes)
dotplot(gotest, showCategory=10)
