# ==============================================================================
# 2-3_clustering_metabolomics_4973.R
# Clustering of significantly changing metabolic features (1 - Pearson R, Ward.D2)
# and peak set enrichment analysis (mummichog, MetaboAnalystR) per cluster.
# NOTE: MetaboAnalystR downloads KEGG libraries, so internet access is required.
#
# Manuscript: Fig. 1B; Extended Data Fig. 2A-B
# Requires  : outputs of 2-2_unfiltered_metabolomics_4973.R
# Inputs    : results/diffMetabolomics_unfiltered_{G,N,P}.tsv
#             results/metabolomics_unfiltered_normalized_batchcorrected.tsv
# Outputs   : source_data/heatmap_untargeted_metab_{Glc,Nit,Pho}.tsv
#             source_data/PSEAclusters_untargeted_metab_{Glc,Nit,Pho}.tsv
#             figure_panels/heatmap_untargeted_metab_{Glc,Nit,Pho}{,_dendrogram,_trends}.pdf
#             figure_panels/PSEAclusters_untargeted_metab_{Glc,Nit,Pho}.pdf
# Run from the repository root (all paths are relative to it).
# ==============================================================================

library(gplots)
library(cluster)
library(ggplot2)
library(RColorBrewer)
library(dendextend)
library(MetaboAnalystR)

#### GLUCOSE ####
# Load data
diffmet = read.csv("results/diffMetabolomics_unfiltered_G.tsv",sep="\t", row.names = 1)
meandf = diffmet[,(ncol(diffmet)-8):ncol(diffmet)]
cleandf = read.csv("results/metabolomics_unfiltered_normalized_batchcorrected.tsv",sep="\t")
cleandf =  cleandf[,1:(ncol(cleandf)-15)]

# Clustering by correlation
X = cleandf[rownames(diffmet)[(diffmet$adj.P.Val<0.05)&(!is.na(diffmet$adj.P.Val))],grepl("G",colnames(cleandf))]
# Make clusters
cordist = 1-cor(t(X), use = "pairwise", method = "pearson")
notna = rowSums(is.na(cordist))==0
clust_aggl = hclust(as.dist(cordist[notna,notna]), method = "ward.D2")
cutgroups = cutree(clust_aggl, k=3)
dend = color_branches(as.dendrogram(clust_aggl), k=3)
pdf("figure_panels/heatmap_untargeted_metab_Glc_dendrogram.pdf", width = 6, height = 2.5); plot(dend); dev.off()
colmap = setNames(get_leaves_branches_col(dend), names(cutgroups)[order.dendrogram(dend)])
# Check clusters
time = order(as.numeric(substr(colnames(meandf),2,4)))
heatdf = log2(meandf[rownames(X)[notna],] / meandf[rownames(X)[notna],"G0"])
heatdf = heatdf[rowSums(!is.na(heatdf))>0,]
Colors=rev(brewer.pal(11,"RdBu"))
Colors=colorRampPalette(Colors)(100)
orderedmet = unlist(sapply(unique(cutgroups), function(x) rownames(heatdf[rownames(heatdf) %in% names(cutgroups)[cutgroups==x],])[hclust(dist(heatdf[rownames(heatdf) %in% names(cutgroups)[cutgroups==x],], method="euclidean"), method = "ward.D2")$order]))
pdf("figure_panels/heatmap_untargeted_metab_Glc.pdf", width = 7, height = 5)
heatmap.2(t(heatdf[orderedmet,time[2:9]]), trace="none", Rowv=F, Colv=F, col=Colors, margins = c(8,4),breaks=seq(-10,10,0.2),symkey=F,
          labCol = orderedmet, ColSideColors = colmap[orderedmet])
dev.off()
write.table(t(heatdf[orderedmet,time[2:9]]),"source_data/heatmap_untargeted_metab_Glc.tsv", sep = "\t")

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
  scale_color_manual(values=unique(colmap[orderedmet])) + scale_fill_manual(values=unique(colmap[orderedmet])) +
  theme_bw() 
ggsave("figure_panels/heatmap_untargeted_metab_Glc_trends.pdf", p, width = 2.5, height = 4)

# mummichog
PSEA <- function(res.an){
  # Guard: check sufficient significant features before proceeding
  n_sig <- sum(res.an[,"adj.P.Val"] <= 0.05, na.rm = TRUE)
  if (n_sig < 20) {
    warning("Cluster skipped: only ", n_sig, " significant features (< 20 threshold)")
    return(NULL)
  }
  # Prepare input tables
  write.table(data.frame(mz=sapply(rownames(res.an), function(x) as.numeric(strsplit(x,"__")[[1]][1])),
                         rt=sapply(rownames(res.an), function(x) as.numeric(strsplit(x,"__")[[1]][2])),
                         p.value=res.an[,"adj.P.Val"]), "mummichog_input.txt", sep="\t", quote = F, row.names = F)
  # Remove adducts for positive mode and currency metabolites
  curr.vec <- c("Water (C00001)", "Proton (C00080)", "Oxygen (C00007)", "NADPH (C00005)",
                "NADP (C00006)", "NADH (C00004)", "NAD (C00003)", "Adenosine triphosphate (C00002)",
                "Pyrophosphate (C00013)","Phosphate (C00009)","Carbon dioxide (C00011)",
                "Hydrogen (C00282)","Hydrogen peroxide (C00027)","Sodium (C01330)")
  add.vec <- c("M [1+]","M+H [1+]","M+2H [2+]","M+3H [3+]","M+Na [1+]",
               "M+H+Na [2+]","M+K [1+]","M+H2O+H [1+]","M-H2O+H [1+]",
               "M-H4O2+H [1+]","M(C13)+H [1+]","M(C13)+2H [2+]","M(C13)+3H [3+]",
               "M(Cl37)+H [1+]","M-NH3+H [1+]","M-CO+H [1+]","M-CO2+H [1+]",
               "M-HCOOH+H [1+]","M+HCOONa [1+]","M-HCOONa+H [1+]","M+NaCl [1+]",
               "M-C3H4O2+H [1+]","M+HCOOK [1+]","M-HCOOK+H [1+]")
  # Enrichment of peaks
  mSet <- InitDataObjects("mass_all", "mummichog", FALSE, default.dpi = 72)
  mSet <- SetPeakFormat(mSet, "mrp") # mz, rt, pval
  mSet <- UpdateInstrumentParameters(mSet, 3.0, "positive", "no", 0.02) # mass accuracy 3 ppm, positive mode, rt tolerance is 2% of max RT
  mSet <- Read.PeakListData(mSet, "mummichog_input.txt")
  mSet <- SetRTincluded(mSet, "minutes")
  mSet <- SanityCheckMummichogData(mSet)
  mSet <- Setup.MapData(mSet, curr.vec) # Map selected currency to data
  mSet <- PerformCurrencyMapping(mSet)
  mSet <- Setup.AdductData(mSet, add.vec) # Set up the selected adducts
  mSet <- PerformAdductMapping(mSet, "positive")
  mSet <- SetPeakEnrichMethod(mSet, "mum", "v2")
  mSet <- SetMummichogPval(mSet, 0.05) # p-value cutoff: 0.05
  mSet <- PerformPSEA(mSet, "sce_kegg", "current", 3 , 100)
  resPSEA = as.data.frame(mSet$mummi.resmat)
  resPSEA["EnrichFact"] = resPSEA["Sig.EC"]/resPSEA["Expected"]
  resPSEA["Pathway"] = rownames(resPSEA)
  # Clean temporary files
  file.remove("mummichog_input.txt","currency_cmpd.qs","initial_ecs.qs","mum_raw.qs","mum_res.qs","mummichog_matched_compound_all.csv",
              "mummichog_pathway_enrichment_mummichog.csv","mummichog_query_mummichog.json","pos_adduct.qs","raw_dataview.csv","sce_kegg.qs")
  return(resPSEA)
}
clustersenrich = data.frame(cluster=NULL, pathway=NULL, FET=NULL, EnrichFact= NULL)
for (c in unique(cutgroups)){
  tempdiff = diffmet[!is.na(diffmet$adj.P.Val),]; tempdiff[names(cutgroups)[cutgroups!=c],"adj.P.Val"] = 1
  PSEAres = PSEA(tempdiff)
  if (is.null(PSEAres)) next  # skip clusters with too few features
  clustersenrich = rbind(clustersenrich, cbind(data.frame(cluster=rep(as.character(c), nrow(PSEAres))), PSEAres[c("Pathway","FET","EnrichFact")]), make.row.names=F)
}
p = ggplot(clustersenrich[clustersenrich$FET<0.05,], aes(x=cluster, y=reorder(reorder(Pathway,EnrichFact),-as.numeric(cluster)), size=EnrichFact, fill=FET)) +
  geom_point(shape = 21, color = "black", stroke = 0.5) + scale_fill_gradient(low = "#e06663", high = "#327eba") + theme_bw()
ggsave("figure_panels/PSEAclusters_untargeted_metab_Glc.pdf", p, width = 5.5, height = 7)
write.table(clustersenrich[clustersenrich$FET<0.05,],"source_data/PSEAclusters_untargeted_metab_Glc.tsv", sep = "\t")

#### NITROGEN ####
# Load data
diffmet = read.csv("results/diffMetabolomics_unfiltered_N.tsv",sep="\t", row.names = 1)
meandf = diffmet[,(ncol(diffmet)-8):ncol(diffmet)]
cleandf = read.csv("results/metabolomics_unfiltered_normalized_batchcorrected.tsv",sep="\t")
cleandf =  cleandf[,1:(ncol(cleandf)-15)]

# Clustering by correlation
X = cleandf[rownames(diffmet)[(diffmet$adj.P.Val<0.05)&(!is.na(diffmet$adj.P.Val))],grepl("N",colnames(cleandf))]
# Make clusters
cordist = 1-cor(t(X), use = "pairwise", method = "pearson")
notna = rowSums(is.na(cordist))==0
clust_aggl = hclust(as.dist(cordist[notna,notna]), method = "ward.D2")
cutgroups = cutree(clust_aggl, k=3)
dend = color_branches(as.dendrogram(clust_aggl), k=3)
pdf("figure_panels/heatmap_untargeted_metab_Nit_dendrogram.pdf", width = 6, height = 2.5); plot(dend); dev.off()
colmap = setNames(get_leaves_branches_col(dend), names(cutgroups)[order.dendrogram(dend)])
# Check clusters
time = order(as.numeric(substr(colnames(meandf),2,4)))
heatdf = log2(meandf[rownames(X)[notna],] / meandf[rownames(X)[notna],"N0"])
heatdf = heatdf[rowSums(!is.na(heatdf))>0,]
Colors=rev(brewer.pal(11,"RdBu"))
Colors=colorRampPalette(Colors)(100)
orderedmet = unlist(sapply(unique(cutgroups), function(x) rownames(heatdf[rownames(heatdf) %in% names(cutgroups)[cutgroups==x],])[hclust(dist(heatdf[rownames(heatdf) %in% names(cutgroups)[cutgroups==x],], method="euclidean"), method = "ward.D2")$order]))
pdf("figure_panels/heatmap_untargeted_metab_Nit.pdf", width = 7, height = 5)
heatmap.2(t(heatdf[orderedmet,time[2:9]]), trace="none", Rowv=F, Colv=F, col=Colors, margins = c(8,4),breaks=seq(-10,10,0.2),symkey=F,
          labCol = orderedmet, ColSideColors = colmap[orderedmet])
dev.off()
write.table(t(heatdf[orderedmet,time[2:9]]),"source_data/heatmap_untargeted_metab_Nit.tsv", sep = "\t")

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
  scale_color_manual(values=unique(colmap[orderedmet])) + scale_fill_manual(values=unique(colmap[orderedmet])) +
  theme_bw() 
ggsave("figure_panels/heatmap_untargeted_metab_Nit_trends.pdf", p, width = 2.5, height = 4)

# mummichog
clustersenrich = data.frame(cluster=NULL, pathway=NULL, FET=NULL, EnrichFact= NULL)
for (c in unique(cutgroups)){
  tempdiff = diffmet[!is.na(diffmet$adj.P.Val),]; tempdiff[names(cutgroups)[cutgroups!=c],"adj.P.Val"] = 1
  PSEAres = PSEA(tempdiff)
  if (is.null(PSEAres)) next  # skip clusters with too few features
  clustersenrich = rbind(clustersenrich, cbind(data.frame(cluster=rep(as.character(c), nrow(PSEAres))), PSEAres[c("Pathway","FET","EnrichFact")]), make.row.names=F)
}
p = ggplot(clustersenrich[clustersenrich$FET<0.05,], aes(x=cluster, y=reorder(reorder(Pathway,EnrichFact),-as.numeric(cluster)), size=EnrichFact, fill=FET)) +
  geom_point(shape = 21, color = "black", stroke = 0.5) + scale_fill_gradient(low = "#e06663", high = "#327eba") + theme_bw()
ggsave("figure_panels/PSEAclusters_untargeted_metab_Nit.pdf", p, width = 5.5, height = 5)
write.table(clustersenrich[clustersenrich$FET<0.05,],"source_data/PSEAclusters_untargeted_metab_Nit.tsv", sep = "\t")

#### PHOSPHATE ####
# Load data
diffmet = read.csv("results/diffMetabolomics_unfiltered_P.tsv",sep="\t", row.names = 1)
meandf = diffmet[,(ncol(diffmet)-8):ncol(diffmet)]
cleandf = read.csv("results/metabolomics_unfiltered_normalized_batchcorrected.tsv",sep="\t")
cleandf =  cleandf[,1:(ncol(cleandf)-15)]

# Clustering by correlation
X = cleandf[rownames(diffmet)[(diffmet$adj.P.Val<0.05)&(!is.na(diffmet$adj.P.Val))],grepl("P",colnames(cleandf))]
# Make clusters
cordist = 1-cor(t(X), use = "pairwise", method = "pearson")
notna = rowSums(is.na(cordist))==0
clust_aggl = hclust(as.dist(cordist[notna,notna]), method = "ward.D2")
cutgroups = cutree(clust_aggl, k=3)
dend = color_branches(as.dendrogram(clust_aggl), k=3)
pdf("figure_panels/heatmap_untargeted_metab_Pho_dendrogram.pdf", width = 6, height = 2.5); plot(dend); dev.off()
colmap = setNames(get_leaves_branches_col(dend), names(cutgroups)[order.dendrogram(dend)])
# Check clusters
time = order(as.numeric(substr(colnames(meandf),2,4)))
heatdf = log2(meandf[rownames(X)[notna],] / meandf[rownames(X)[notna],"P0"])
heatdf = heatdf[rowSums(!is.na(heatdf))>0,]
Colors=rev(brewer.pal(11,"RdBu"))
Colors=colorRampPalette(Colors)(100)
orderedmet = unlist(sapply(unique(cutgroups), function(x) rownames(heatdf[rownames(heatdf) %in% names(cutgroups)[cutgroups==x],])[hclust(dist(heatdf[rownames(heatdf) %in% names(cutgroups)[cutgroups==x],], method="euclidean"), method = "ward.D2")$order]))
pdf("figure_panels/heatmap_untargeted_metab_Pho.pdf", width = 7, height = 5)
heatmap.2(t(heatdf[orderedmet,time[2:9]]), trace="none", Rowv=F, Colv=F, col=Colors, margins = c(8,4),breaks=seq(-10,10,0.2),symkey=F,
          labCol = orderedmet, ColSideColors = colmap[orderedmet])
dev.off()
write.table(t(heatdf[orderedmet,time[2:9]]),"source_data/heatmap_untargeted_metab_Pho.tsv", sep = "\t")


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
  scale_color_manual(values=unique(colmap[orderedmet])) + scale_fill_manual(values=unique(colmap[orderedmet])) +
  theme_bw() 
ggsave("figure_panels/heatmap_untargeted_metab_Pho_trends.pdf", p, width = 2.5, height = 4)

# mummichog
clustersenrich = data.frame(cluster=NULL, pathway=NULL, FET=NULL, EnrichFact= NULL)
for (c in unique(cutgroups)){
  tempdiff = diffmet[!is.na(diffmet$adj.P.Val),]; tempdiff[names(cutgroups)[cutgroups!=c],"adj.P.Val"] = 1
  PSEAres = PSEA(tempdiff)
  if (is.null(PSEAres)) next  # skip clusters with too few features
  clustersenrich = rbind(clustersenrich, cbind(data.frame(cluster=rep(as.character(c), nrow(PSEAres))), PSEAres[c("Pathway","FET","EnrichFact")]), make.row.names=F)
}
p = ggplot(clustersenrich[clustersenrich$FET<0.05,], aes(x=cluster, y=reorder(reorder(Pathway,EnrichFact),-as.numeric(cluster)), size=EnrichFact, fill=FET)) +
  geom_point(shape = 21, color = "black", stroke = 0.5) + scale_fill_gradient(low = "#e06663", high = "#327eba") + theme_bw()
ggsave("figure_panels/PSEAclusters_untargeted_metab_Pho.pdf", p, width = 5.5, height = 5)
write.table(clustersenrich[clustersenrich$FET<0.05,],"source_data/PSEAclusters_untargeted_metab_Pho.tsv", sep = "\t")
