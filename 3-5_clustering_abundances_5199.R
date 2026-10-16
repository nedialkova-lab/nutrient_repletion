# ==============================================================================
# 3-5_clustering_abundances_5199.R
# Clustering and GO enrichment of significantly changing protein abundances, and
# median absolute change of translation-annotated proteins.
#
# Manuscript: Fig. 5C (bottom, abundance)
# Requires  : outputs of 3-4_abundances_timecourse_5199.R
# Inputs    : results/diffProtAbund_{G,N,P}.tsv
#             results/proteomics_batchcorrected.tsv
# Outputs   : results/GOanalysis_{Glc,Nit}_prot_abundance.tsv
#             source_data/heatmap_prot_abundance_{Glc,Nit,Pho}.tsv
#             figure_panels/heatmap_prot_abundance_{Glc,Nit,Pho}{,_dendrogram,_trends}.pdf
#             figure_panels/lineplot_prot_abundance_GOtrans.pdf
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

# Record info about translation
gotransdf = data.frame(cond=c(rep("G",8),rep("N",8),rep("P",8)),
                       time=rep(c(0,1,2,4,7,10,14,20),3),
                       fcmedian=NA, fcquantmin=NA, fcquantmax=NA)

#### GLUCOSE ####
# Load data
diffprot = read.csv("results/diffProtAbund_G.tsv",sep="\t", row.names = 1)
meandf = diffprot[,(ncol(diffprot)-7):ncol(diffprot)]
normdf = read.csv("results/proteomics_batchcorrected.tsv", sep="\t", row.names = 1)

# Clustering by correlation
X = normdf[rownames(diffprot)[(diffprot$adj.P.Val<0.05)&(!is.na(diffprot$adj.P.Val))],grepl("G",colnames(normdf))]
# Make clusters
cordist = 1-cor(t(X), use = "pairwise", method = "pearson")
notna = rowSums(is.na(cordist))==0
clust_aggl = hclust(as.dist(cordist[notna,notna]), method = "ward.D2")
cutgroups = cutree(clust_aggl, k=3)
dend = color_branches(as.dendrogram(clust_aggl), k=3)
pdf("figure_panels/heatmap_prot_abundance_Glc_dendrogram.pdf", width = 6, height = 2.5); plot(dend); dev.off()
colmap = setNames(get_leaves_branches_col(dend), names(cutgroups)[order.dendrogram(dend)])
# Check clusters
time = order(as.numeric(substr(colnames(meandf),2,4)))
heatdf = log2(meandf[rownames(X)[notna],] / meandf[rownames(X)[notna],"G0"])
heatdf = heatdf[rowSums(!is.na(heatdf))>0,]
Colors=rev(brewer.pal(11,"RdBu"))
Colors=colorRampPalette(Colors)(100)
orderedprot = unlist(sapply(unique(cutgroups), function(x) rownames(heatdf[rownames(heatdf) %in% names(cutgroups)[cutgroups==x],])[hclust(dist(heatdf[rownames(heatdf) %in% names(cutgroups)[cutgroups==x],], method="euclidean"), method = "ward.D2")$order]))
pdf("figure_panels/heatmap_prot_abundance_Glc.pdf", width = 7, height = 5)
heatmap.2(t(heatdf[orderedprot,time[2:8]]), trace="none", Rowv=F, Colv=F, col=Colors, margins = c(6,4),breaks=seq(-2,2,0.04),symkey=F,
          labCol = diffprot[orderedprot,"gene"], ColSideColors = colmap[orderedprot])
dev.off()
write.table(t(heatdf[orderedprot,time[2:8]]),"source_data/heatmap_prot_abundance_Glc.tsv", sep = "\t")

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
  scale_color_manual(values=unique(colmap[orderedprot])) + scale_fill_manual(values=unique(colmap[orderedprot])) +
  theme_bw()
ggsave("figure_panels/heatmap_prot_abundance_Glc_trends.pdf", p, width = 2.5, height = 3)

# GO enrichment
includedgenes = unique(unlist(sapply(unique(diffprot$gene), function(x) strsplit(x,";")[[1]])))
gotab = select(org.Sc.sgd.db, includedgenes, columns=c("GOALL","GENENAME"), keytype="GENENAME")
mapgo = go2term(unique(gotab$GO))
gotab$goterm = as.character(sapply(gotab$GO, function(x) mapgo[mapgo$go_id==x,"Term"]))
gotab = unique(gotab[gotab$ONTOLOGY=="BP",c("goterm","GENENAME")])

keggtab <- getGeneKEGGLinks(species="sce")
keggtab$gene <- mapIds(org.Sc.sgd.db, keggtab$GeneID, column="GENENAME", keytype="ORF")
mappaths = getKEGGPathwayNames(species="sce")
keggtab$pathname = sapply(keggtab$PathwayID, function(x) substr(mappaths[mappaths$PathwayID==x,"Description"],0, nchar(mappaths[mappaths$PathwayID==x,"Description"])-43))
keggtab = unique(keggtab[,c("pathname","gene")])

flatprots = unique(unlist(sapply(rownames(normdf),function(x) unlist(strsplit(x, ";")))))
mapping = mapUniProt(query=flatprots, from = "UniProtKB_AC-ID", to = "Gene_Name")
clusterdf = data.frame(cluster = unlist(sapply(names(cutgroups), function(x) rep(cutgroups[x],str_count(x,";")+1))),
                       protein = unlist(sapply(names(cutgroups), function(x) strsplit(x,";")[[1]])))
clusterdf$gene_name = as.character(sapply(clusterdf$protein, function(x) mapping[mapping$From %in% x, "To"]))

gotest <- compareCluster(gene_name~cluster, data=clusterdf, fun="enricher",
                         TERM2GENE = gotab, pvalueCutoff=0.05, universe = includedgenes)
write.table(gotest@compareClusterResult, "results/GOanalysis_Glc_prot_abundance.tsv", sep = "\t")
dotplot(gotest, showCategory=10)
gotest <- compareCluster(gene_name~cluster, data=clusterdf, fun="enricher",
                         TERM2GENE = keggtab, pvalueCutoff=0.05, universe = includedgenes)
dotplot(gotest, showCategory=10)

# Translation
transdf = abs(log2(meandf / meandf[,"G0"]))
gotransdf[gotransdf$cond=="G","fcmedian"] = sapply(c(0,1,2,4,7,10,14,20), function(x) median(transdf[sapply(diffprot[rownames(transdf),"gene"], function(y) any(strsplit(y,",")[[1]] %in% gotab[gotab$goterm=="translation","GENENAME"])),sprintf("G%i",x)], na.rm=T))
gotransdf[gotransdf$cond=="G","fcquantmin"] = sapply(c(0,1,2,4,7,10,14,20), function(x) quantile(transdf[sapply(diffprot[rownames(transdf),"gene"], function(y) any(strsplit(y,",")[[1]] %in% gotab[gotab$goterm=="translation","GENENAME"])),sprintf("G%i",x)], probs=0.25, na.rm=T))
gotransdf[gotransdf$cond=="G","fcquantmax"] = sapply(c(0,1,2,4,7,10,14,20), function(x) quantile(transdf[sapply(diffprot[rownames(transdf),"gene"], function(y) any(strsplit(y,",")[[1]] %in% gotab[gotab$goterm=="translation","GENENAME"])),sprintf("G%i",x)], probs=0.75, na.rm=T))

#### NITROGEN ####
# Load data
diffprot = read.csv("results/diffProtAbund_N.tsv",sep="\t", row.names = 1)
meandf = diffprot[,(ncol(diffprot)-7):ncol(diffprot)]
normdf = read.csv("results/proteomics_batchcorrected.tsv", sep="\t", row.names = 1)

# Clustering by correlation
X = normdf[rownames(diffprot)[(diffprot$adj.P.Val<0.05)&(!is.na(diffprot$adj.P.Val))],grepl("N",colnames(normdf))]
# Make clusters
cordist = 1-cor(t(X), use = "pairwise", method = "pearson")
notna = rowSums(is.na(cordist))==0
clust_aggl = hclust(as.dist(cordist[notna,notna]), method = "ward.D2")
cutgroups = cutree(clust_aggl, k=2)
dend = color_branches(as.dendrogram(clust_aggl), k=2)
pdf("figure_panels/heatmap_prot_abundance_Nit_dendrogram.pdf", width = 6, height = 2.5); plot(dend); dev.off()
colmap = setNames(get_leaves_branches_col(dend), names(cutgroups)[order.dendrogram(dend)])
# Check clusters
time = order(as.numeric(substr(colnames(meandf),2,4)))
heatdf = log2(meandf[rownames(X)[notna],] / meandf[rownames(X)[notna],"N0"])
heatdf = heatdf[rowSums(!is.na(heatdf))>0,]
Colors=rev(brewer.pal(11,"RdBu"))
Colors=colorRampPalette(Colors)(100)
orderedprot = unlist(sapply(unique(cutgroups), function(x) rownames(heatdf[rownames(heatdf) %in% names(cutgroups)[cutgroups==x],])[hclust(dist(heatdf[rownames(heatdf) %in% names(cutgroups)[cutgroups==x],], method="euclidean"), method = "ward.D2")$order]))
pdf("figure_panels/heatmap_prot_abundance_Nit.pdf", width = 7, height = 5)
heatmap.2(t(heatdf[orderedprot,time[2:8]]), trace="none", Rowv=F, Colv=F, col=Colors, margins = c(6,4),breaks=seq(-2,2,0.04),symkey=F,
          labCol = diffprot[orderedprot,"gene"], ColSideColors = colmap[orderedprot])
dev.off()
write.table(t(heatdf[orderedprot,time[2:8]]),"source_data/heatmap_prot_abundance_Nit.tsv", sep = "\t")

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
  scale_color_manual(values=unique(colmap[orderedprot])) + scale_fill_manual(values=unique(colmap[orderedprot])) +
  theme_bw()
ggsave("figure_panels/heatmap_prot_abundance_Nit_trends.pdf", p, width = 2.5, height = 2.5)

# GO enrichment
includedgenes = unique(unlist(sapply(unique(diffprot$gene), function(x) strsplit(x,";")[[1]])))
gotab = select(org.Sc.sgd.db, includedgenes, columns=c("GOALL","GENENAME"), keytype="GENENAME")
mapgo = go2term(unique(gotab$GO))
gotab$goterm = as.character(sapply(gotab$GO, function(x) mapgo[mapgo$go_id==x,"Term"]))
gotab = unique(gotab[gotab$ONTOLOGY=="BP",c("goterm","GENENAME")])

keggtab <- getGeneKEGGLinks(species="sce")
keggtab$gene <- mapIds(org.Sc.sgd.db, keggtab$GeneID, column="GENENAME", keytype="ORF")
mappaths = getKEGGPathwayNames(species="sce")
keggtab$pathname = sapply(keggtab$PathwayID, function(x) substr(mappaths[mappaths$PathwayID==x,"Description"],0, nchar(mappaths[mappaths$PathwayID==x,"Description"])-43))
keggtab = unique(keggtab[,c("pathname","gene")])

clusterdf = data.frame(cluster = unlist(sapply(names(cutgroups), function(x) rep(cutgroups[x],str_count(x,";")+1))),
                       protein = unlist(sapply(names(cutgroups), function(x) strsplit(x,";")[[1]])))
clusterdf$gene_name = as.character(sapply(clusterdf$protein, function(x) mapping[mapping$From %in% x, "To"]))

gotest <- compareCluster(gene_name~cluster, data=clusterdf, fun="enricher",
                         TERM2GENE = gotab, pvalueCutoff=0.05, universe = includedgenes)
write.table(gotest@compareClusterResult, "results/GOanalysis_Nit_prot_abundance.tsv", sep = "\t")
dotplot(gotest, showCategory=10)
gotest <- compareCluster(gene_name~cluster, data=clusterdf, fun="enricher",
                         TERM2GENE = keggtab, pvalueCutoff=0.05, universe = includedgenes)
dotplot(gotest, showCategory=10)

# Translation
transdf = abs(log2(meandf / meandf[,"N0"]))
gotransdf[gotransdf$cond=="N","fcmedian"] = sapply(c(0,1,2,4,7,10,14,20), function(x) median(transdf[sapply(diffprot[rownames(transdf),"gene"], function(y) any(strsplit(y,",")[[1]] %in% gotab[gotab$goterm=="translation","GENENAME"])),sprintf("N%i",x)], na.rm=T))
gotransdf[gotransdf$cond=="N","fcquantmin"] = sapply(c(0,1,2,4,7,10,14,20), function(x) quantile(transdf[sapply(diffprot[rownames(transdf),"gene"], function(y) any(strsplit(y,",")[[1]] %in% gotab[gotab$goterm=="translation","GENENAME"])),sprintf("N%i",x)], probs=0.25, na.rm=T))
gotransdf[gotransdf$cond=="N","fcquantmax"] = sapply(c(0,1,2,4,7,10,14,20), function(x) quantile(transdf[sapply(diffprot[rownames(transdf),"gene"], function(y) any(strsplit(y,",")[[1]] %in% gotab[gotab$goterm=="translation","GENENAME"])),sprintf("N%i",x)], probs=0.75, na.rm=T))

#### PHOSPHATE ####
# Load data
diffprot = read.csv("results/diffProtAbund_P.tsv",sep="\t", row.names = 1)
meandf = diffprot[,(ncol(diffprot)-7):ncol(diffprot)]
normdf = read.csv("results/proteomics_batchcorrected.tsv", sep="\t", row.names = 1)

# Clustering by correlation
X = normdf[rownames(diffprot)[(diffprot$adj.P.Val<0.05)&(!is.na(diffprot$adj.P.Val))],grepl("P",colnames(normdf))]
# Make clusters
cordist = 1-cor(t(X), use = "pairwise", method = "pearson")
notna = rowSums(is.na(cordist))==0
clust_aggl = hclust(as.dist(cordist[notna,notna]), method = "ward.D2")
cutgroups = cutree(clust_aggl, k=2)
dend = color_branches(as.dendrogram(clust_aggl), k=2)
pdf("figure_panels/heatmap_prot_abundance_Pho_dendrogram.pdf", width = 6, height = 2.5); plot(dend); dev.off()
colmap = setNames(get_leaves_branches_col(dend), names(cutgroups)[order.dendrogram(dend)])
# Check clusters
time = order(as.numeric(substr(colnames(meandf),2,4)))
heatdf = log2(meandf[rownames(X)[notna],] / meandf[rownames(X)[notna],"P0"])
heatdf = heatdf[rowSums(!is.na(heatdf))>0,]
Colors=rev(brewer.pal(11,"RdBu"))
Colors=colorRampPalette(Colors)(100)
orderedprot = unlist(sapply(unique(cutgroups), function(x) rownames(heatdf[rownames(heatdf) %in% names(cutgroups)[cutgroups==x],])[hclust(dist(heatdf[rownames(heatdf) %in% names(cutgroups)[cutgroups==x],], method="euclidean"), method = "ward.D2")$order]))
pdf("figure_panels/heatmap_prot_abundance_Pho.pdf", width = 7, height = 5)
heatmap.2(t(heatdf[orderedprot,time[2:8]]), trace="none", Rowv=F, Colv=F, col=Colors, margins = c(6,4),breaks=seq(-2,2,0.04),symkey=F,
          labCol = diffprot[orderedprot,"gene"], ColSideColors = colmap[orderedprot])
dev.off()
write.table(t(heatdf[orderedprot,time[2:8]]),"source_data/heatmap_prot_abundance_Pho.tsv", sep = "\t")

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
  scale_color_manual(values=unique(colmap[orderedprot])) + scale_fill_manual(values=unique(colmap[orderedprot])) +
  theme_bw()
ggsave("figure_panels/heatmap_prot_abundance_Pho_trends.pdf", p, width = 2.5, height = 2.5)

# Translation
includedgenes = unique(unlist(sapply(unique(diffprot$gene), function(x) strsplit(x,";")[[1]])))
gotab = select(org.Sc.sgd.db, includedgenes, columns=c("GOALL","GENENAME"), keytype="GENENAME")
mapgo = go2term(unique(gotab$GO))
gotab$goterm = as.character(sapply(gotab$GO, function(x) mapgo[mapgo$go_id==x,"Term"]))
gotab = unique(gotab[gotab$ONTOLOGY=="BP",c("goterm","GENENAME")])

transdf = abs(log2(meandf / meandf[,"P0"]))
gotransdf[gotransdf$cond=="P","fcmedian"] = sapply(c(0,1,2,4,7,10,14,20), function(x) median(transdf[sapply(diffprot[rownames(transdf),"gene"], function(y) any(strsplit(y,",")[[1]] %in% gotab[gotab$goterm=="translation","GENENAME"])),sprintf("P%i",x)], na.rm=T))
gotransdf[gotransdf$cond=="P","fcquantmin"] = sapply(c(0,1,2,4,7,10,14,20), function(x) quantile(transdf[sapply(diffprot[rownames(transdf),"gene"], function(y) any(strsplit(y,",")[[1]] %in% gotab[gotab$goterm=="translation","GENENAME"])),sprintf("P%i",x)], probs=0.25, na.rm=T))
gotransdf[gotransdf$cond=="P","fcquantmax"] = sapply(c(0,1,2,4,7,10,14,20), function(x) quantile(transdf[sapply(diffprot[rownames(transdf),"gene"], function(y) any(strsplit(y,",")[[1]] %in% gotab[gotab$goterm=="translation","GENENAME"])),sprintf("P%i",x)], probs=0.75, na.rm=T))

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
ggsave("figure_panels/lineplot_prot_abundance_GOtrans.pdf", p, width = 4, height = 2.5)
