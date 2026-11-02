# ==============================================================================
# 5-1_LiPMS_vs_Metab.R
# Correlation of significantly changing LiP-MS peptides with significantly changing
# identified metabolites (Pearson, >= 8 paired observations, fdrtool local FDR),
# clustering of the correlation matrix and selected peptide-metabolite examples.
# NOTE: the export of source_data/LiPMS_vs_Metab.tsv (read by the 5-2 notebooks) is
# commented out because the annotation step is slow; the pre-computed file is
# provided on Zenodo. Uncomment the '#corsig[...]' and '#write.table' lines to
# regenerate it.
#
# Manuscript: Fig. 5D, 5F, 5H-L; Extended Data Fig. 6E-I
# Requires  : outputs of 2-2, 4-1
# Inputs    : results/diffMetabolomics_unfiltered_{G,N,P}.tsv
#             results/metabolomics_unfiltered_normalized_batchcorrected.tsv
#             results/diffStructurome_{G,N,P}.tsv
#             results/resLiP.tsv
# Outputs   : figure_panels/cormatrix_heatmap{,_translation,_translation_snapshots}.pdf
#             figure_panels/cormatrix_proteins.pdf
#             figure_panels/GOtrans_top15_cor{LiP,Metab}.pdf
#             figure_panels/{MES1vsAMP,SES1vsAMP,TEF1-2vsAsp,STM1vsAsp}_{woodplot,correlation}.pdf
#             figure_panels/CDC60_Gua_woodplot.pdf
#             figure_panels/CDC60vsGuan_{correlation,N_lineplot}.pdf
#             figure_panels/{MES1vsAMP,SES1vsAMP,TEF1-2vsAsp}_G_lineplot.pdf
#             plots/cormatrix_heatmap{,_translation}.html (interactive)
#             source_data/LiPMS_vs_Metab_filteredsignificant.tsv
# Run from the repository root (all paths are relative to it).
# ==============================================================================

library(GeneNet)
library(graph)
library(Rgraphviz)
library(reshape)
library(gplots)
library(RColorBrewer)
library(ggplot2)
library(patchwork)
library(ggpubr)
library(ggrepel)
library(eulerr)
library(LiPAnalyzeR)
library(clusterProfiler)
library(org.Sc.sgd.db)
library(heatmaply)
library(htmlwidgets)

plotcor <- function(df,m,p){
  if(grepl("_",p)){p = annotprot[apply(annotprot[,c("gene_name","startPosition","endPosition")],1, paste0, collapse="_") %in% p,"Peptide"]}
  plotdf = as.data.frame(t(as.data.frame(df[c(m,p),])))
  plotdf$labels = rownames(plotdf)
  plotdf$cond = substr(plotdf$labels,1,1)
  p1 <- ggplot(plotdf, aes(x=get(m),y=get(p))) +
    geom_smooth(color="black", method="lm", se=T) +
    geom_point(aes(color=cond), alpha=0.6) +
    stat_cor() +
    geom_text_repel(aes(label=labels)) +
    xlab(m) + ylab(p) +
    theme_classic()
  longdf = melt(plotdf[,1:(ncol(plotdf)-1)],id.vars = "labels")
  longdf$cond = substr(longdf$labels,1,1)
  longdf$time = as.numeric(substr(longdf$labels,2,3))
  p2 <- ggplot(longdf[longdf$cond=="G",], aes(x=time, y=value)) +
    facet_wrap(~ variable, nrow=2, scales = "free_y")+
    geom_line() +
    geom_point() +
    ggtitle("Glucose") +
    theme_bw()
  p3 <- ggplot(longdf[longdf$cond=="N",], aes(x=time, y=value)) +
    facet_wrap(~ variable, nrow=2, scales = "free_y")+
    geom_line() +
    geom_point() +
    ggtitle("Nitrogen") +
    theme_bw()
  p4 <- ggplot(longdf[longdf$cond=="P",], aes(x=time, y=value)) +
    facet_wrap(~ variable, nrow=2, scales = "free_y")+
    geom_line() +
    geom_point() +
    ggtitle("Phosphate") +
    theme_bw()
  (p1 | p2) / (p3 | p4)
}
Colors=brewer.pal(9,"YlOrRd")
Colors=colorRampPalette(Colors)(100)

#### LOAD DATA ####
# Metabolomics
diffmetG = read.csv("results/diffMetabolomics_unfiltered_G.tsv",sep="\t", row.names = 1)
diffmetN = read.csv("results/diffMetabolomics_unfiltered_N.tsv",sep="\t", row.names = 1)
diffmetP = read.csv("results/diffMetabolomics_unfiltered_P.tsv",sep="\t", row.names = 1)
annotmet = rbind(diffmetG[,c("id_top_name","id_top_class","id_top_adduct","id_top_score","id_source","name")],
                 diffmetN[,c("id_top_name","id_top_class","id_top_adduct","id_top_score","id_source","name")],
                 diffmetP[,c("id_top_name","id_top_class","id_top_adduct","id_top_score","id_source","name")])
annotmet = annotmet[rownames(annotmet) %in% c(rownames(diffmetG),rownames(diffmetN),rownames(diffmetP)),]

cleandf = read.csv("results/metabolomics_unfiltered_normalized_batchcorrected.tsv",sep="\t")
cleandf =  cleandf[,1:(ncol(cleandf)-15)]
longdf = melt(as.matrix(cleandf))
longdf$condition = sapply(as.character(longdf$X2),function(x) strsplit(x,"\\.")[[1]][1])
metdf = cast(longdf, X1 ~ condition, median, na.rm=T)
rownames(metdf)= metdf$X1
metdf= metdf[, -c(1)]
# Center per condition
metdf = t(apply(as.data.frame(metdf),1,function(x) x - c(rep(mean(x[grepl("G",colnames(metdf))],na.rm=T),sum(grepl("G",colnames(metdf)))),
                                                         rep(mean(x[grepl("N",colnames(metdf))],na.rm=T),sum(grepl("N",colnames(metdf)))),
                                                         rep(mean(x[grepl("P",colnames(metdf))],na.rm=T),sum(grepl("P",colnames(metdf)))))))

# LiP-MS
difflipG = read.csv("results/diffStructurome_G.tsv",sep="\t", row.names = "quantID")
difflipN = read.csv("results/diffStructurome_N.tsv",sep="\t", row.names = "quantID")
difflipP = read.csv("results/diffStructurome_P.tsv",sep="\t", row.names = "quantID")
annotprot = rbind(difflipG[,c("Peptide","Protein","gene_name","startPosition","endPosition","isTryptic")],
                  difflipN[,c("Peptide","Protein","gene_name","startPosition","endPosition","isTryptic")],
                  difflipP[,c("Peptide","Protein","gene_name","startPosition","endPosition","isTryptic")])
annotprot = annotprot[rownames(annotprot) %in% c(rownames(difflipG),rownames(difflipN),rownames(difflipP)),]

resLiP = read.csv("results/resLiP.tsv",sep="\t")
longdf = melt(as.matrix(resLiP))
longdf$condition = sapply(as.character(longdf$X2),function(x) strsplit(x,"\\.")[[1]][1])
lipdf = cast(longdf, X1 ~ condition, median, na.rm=T)
rownames(lipdf)= lipdf$X1
lipdf= lipdf[, -c(1)]
lipdf[is.na(lipdf)] <- NA
# Center per condition
lipdf = t(apply(as.data.frame(lipdf),1,function(x) x - c(rep(mean(x[grepl("G",colnames(lipdf))],na.rm=T),sum(grepl("G",colnames(lipdf)))),
                                                         rep(mean(x[grepl("N",colnames(lipdf))],na.rm=T),sum(grepl("N",colnames(lipdf)))),
                                                         rep(mean(x[grepl("P",colnames(lipdf))],na.rm=T),sum(grepl("P",colnames(lipdf)))))))

# Significant hits
# keep only identified metabolites with score >0.2
sigmet = unique(c(rownames(diffmetG)[(diffmetG$adj.P.Val<0.05)&!is.na(diffmetG$adj.P.Val)&!is.na(diffmetG$id_top_name)],
                  rownames(diffmetN)[(diffmetN$adj.P.Val<0.05)&!is.na(diffmetN$adj.P.Val)&!is.na(diffmetN$id_top_name)],
                  rownames(diffmetP)[(diffmetP$adj.P.Val<0.05)&!is.na(diffmetP$adj.P.Val)&!is.na(diffmetP$id_top_name)]))
siglip = unique(c(rownames(difflipG)[(difflipG$adj.P.Val<0.05)&!is.na(difflipG$adj.P.Val)],
                  rownames(difflipN)[(difflipN$adj.P.Val<0.05)&!is.na(difflipN$adj.P.Val)],
                  rownames(difflipP)[(difflipP$adj.P.Val<0.05)&!is.na(difflipP$adj.P.Val)]))

#### CORRELATION TABLE ####

jointdf = rbind(lipdf[siglip,],metdf[sigmet,colnames(lipdf)])
cordf= cor(t(jointdf), use = "pairwise", method = "pearson")[1:length(siglip),(length(siglip)+1):nrow(jointdf)]
rownames(cordf) = siglip; colnames(cordf) = sigmet
corsigdf = sapply(sigmet, function(x) sapply(siglip, function(y)
  if(sum(!is.na(as.numeric(metdf[x, colnames(lipdf)]))&!is.na(as.numeric(lipdf[y, ])))>=8){
    cor.test(as.numeric(metdf[x, colnames(lipdf)]), 
             as.numeric(lipdf[y, ]), 
             use = "pairwise", 
             method = "pearson")$p.value
    }else{NA}))

# Compute significance
corsiglong = as.data.frame(melt(corsigdf))
fdr.out = fdrtool(corsiglong[!is.na(corsiglong$value), "value"], statistic = "pvalue", plot = T)
pval = fdr.out$pval
qval = fdr.out$qval
lfdr = fdr.out$lfdr
corsig = corsiglong[,1:2]
corsig[!is.na(corsiglong$value),"pval"] = pval
corsig[!is.na(corsiglong$value),"qval"] = qval
corsig[!is.na(corsiglong$value),"lfdr"] = lfdr
corsigdf = cast(corsig, formula = "X1 ~ X2", value = "lfdr")
rownames(corsigdf) = corsigdf$X1; corsigdf = corsigdf[,-c(1)]
# Show only significant correlations
cordffilt <- cordf
cordffilt[(corsigdf>0.05)|is.na(corsigdf)] <- NA
sum(corsigdf<0.05, na.rm = T)
# Annotate and save
# NOTE: commented out because it is slow; uncomment these lines (and the two below
# the GO 'translation' block) to regenerate source_data/LiPMS_vs_Metab.tsv, which is
# otherwise provided pre-computed on Zenodo and read by the 5-2 notebooks.
#corsig["R2"] = apply(corsig,1,function(x) if(is.na(x["pval"])){NA}else{cordf[x["X1"],x["X2"]]})
#corsig["annot_pept"] = apply(corsig,1,function(x) sprintf("%s_%s-%s", annotprot[x["X1"],"gene_name"],annotprot[x["X1"],"startPosition"],annotprot[x["X1"],"endPosition"]))
#corsig["annot_metab"] = apply(corsig,1,function(x) annotmet[x["X2"],"id_top_name"])

includedgenes = unique(unlist(sapply(annotprot[,"gene_name"], function(x) unlist(strsplit(x,";")))))
gotab = AnnotationDbi::select(org.Sc.sgd.db, includedgenes, columns=c("GOALL","GENENAME"), keytype="GENENAME")
mapgo = go2term(unique(gotab$GO))
gotab$goterm = as.character(sapply(gotab$GO, function(x) mapgo[mapgo$go_id==x,"Term"]))
gotab = unique(gotab[gotab$ONTOLOGY=="BP",c("goterm","GENENAME")])

goterm = "translation"
go_genes <- gotab[gotab$goterm == goterm, "GENENAME"]
#corsig["GO_translation"] = apply(corsig,1,function(x) any(strsplit(annotprot[x["X1"],"gene_name"],";")[[1]] %in% go_genes))
#write.table(corsig[!is.na(corsig$pval),],"source_data/LiPMS_vs_Metab.tsv", sep = "\t")
# Keep only significant correlations (local FDR < 0.05) for the Supplementary Table
corsigfull = read.csv("source_data/LiPMS_vs_Metab.tsv", sep = "\t")
write.table(corsigfull[(corsigfull$lfdr<0.05)&!is.na(corsigfull$lfdr),], "source_data/LiPMS_vs_Metab_filteredsignificant.tsv", sep = "\t")

# Heatmap
# Clean all nan rows
keeppept = rowSums(!is.na(cordffilt))>0 
keepmet = colSums(!is.na(cordffilt))>0
cordffilt_clean = cordffilt[keeppept,keepmet]
cordf_clean = cordf[keeppept,keepmet]
# Plot
row_dist <- dist(t(cordf_clean^2), method = "euclidean")
col_dist <- dist(cordf_clean^2, method = "euclidean")
row_hc <- hclust(row_dist, method = "ward.D2")
col_hc <- hclust(col_dist, method = "ward.D2")
pdf("figure_panels/cormatrix_heatmap.pdf", width = 10, height = 4)
heatmap.2(t(cordffilt_clean^2), 
          trace = "none", col = Colors, margins = c(4, 4),
          breaks = seq(0, 1, 0.01), symkey = F, na.color = "white",
          Rowv = as.dendrogram(row_hc),
          Colv = as.dendrogram(col_hc), dendrogram = "both")
dev.off()


p = heatmaply(t(cordffilt_clean^2),
          Rowv = as.dendrogram(row_hc),
          Colv = as.dendrogram(col_hc),
          colors = Colors,
          na.value = "white",
          limits = c(0, 1),
          margins = c(40, 40),
          #showticklabels= F,
          labCol = sapply(rownames(cordffilt_clean), function(x) sprintf("%s_%s-%s", annotprot[x,"gene_name"], annotprot[x,"startPosition"], annotprot[x,"endPosition"])),
          labRow = sapply(colnames(cordffilt_clean), function(x) sprintf("%s|%s", x, annotmet[x,"id_top_name"])))
saveWidget(p, "plots/cormatrix_heatmap.html", selfcontained = TRUE)

# Find most regulated prots
subdf = cbind(cordffilt,data.frame(pept=rownames(cordffilt)))
melteddf = melt(subdf, id.vars = "pept")
melteddf = melteddf[!is.na(melteddf$value),]
top30 = names(sort(table(melteddf[!is.na(melteddf$value),"pept"]), decreasing=T))[1:30]
melteddf = melteddf[melteddf$pept %in% top30, ]
melteddf$annot = factor(apply(annotprot[melteddf$pept,c("gene_name","startPosition","endPosition")],1, paste0, collapse="_"),
                        levels = apply(annotprot[top30,c("gene_name","startPosition","endPosition")],1, paste0, collapse="_"))
ggplot(melteddf, aes(x=annot)) +
  geom_bar(stat = "count") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

# Find most regulatory metab
subdf = cbind(t(cordffilt),data.frame(metab=sapply(colnames(cordffilt), function(x) sprintf("%s\n%s",x,annotmet[x,"id_top_name"]))))
melteddf = melt(subdf, id.vars = "metab")
melteddf = melteddf[!is.na(melteddf$value),]
top30 = names(sort(table(melteddf[!is.na(melteddf$value),"metab"]), decreasing=T))[1:30]
melteddf = melteddf[melteddf$metab %in% top30, ]
melteddf$metab = factor(melteddf$metab, levels = top30)

ggplot(melteddf, aes(x=metab)) +
  geom_bar(stat = "count") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

# Significant associations in translation
goterm = "translation"
go_genes <- gotab[gotab$goterm == goterm, "GENENAME"]
sublip = siglip[sapply(annotprot[siglip,"gene_name"], function(x) any(strsplit(x,";")[[1]] %in% go_genes))]
# lip
subdf = cbind(cordffilt,data.frame(pept=rownames(cordffilt)))
melteddf = melt(subdf[sublip,], id.vars = "pept")
melteddf = melteddf[!is.na(melteddf$value),]
top15 = names(sort(table(melteddf[!is.na(melteddf$value),"pept"]), decreasing=T))[1:15]
melteddf = melteddf[melteddf$pept %in% top15, ]
melteddf$annot = factor(apply(annotprot[melteddf$pept,c("gene_name","startPosition","endPosition")],1, paste0, collapse="_"),
                        levels = apply(annotprot[top15,c("gene_name","startPosition","endPosition")],1, paste0, collapse="_"))
p = ggplot(melteddf, aes(x=annot)) +
  geom_bar(stat = "count") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggsave("figure_panels/GOtrans_top15_corLiP.pdf", p, height = 3, width = 4)
#metab
subdf = cbind(t(cordffilt),data.frame(metab=sapply(colnames(cordffilt), function(x) sprintf("%s\n%s",x,annotmet[x,"id_top_name"]))))
melteddf = melt(subdf[,c(sublip,"metab")], id.vars = "metab")
melteddf = melteddf[!is.na(melteddf$value),]
top15 = names(sort(table(melteddf[!is.na(melteddf$value),"metab"]), decreasing=T))[1:15]
melteddf = melteddf[melteddf$metab %in% top15, ]
melteddf$metab = factor(melteddf$metab, levels = top15)

p = ggplot(melteddf, aes(x=metab)) +
  geom_bar(stat = "count") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggsave("figure_panels/GOtrans_top15_corMetab.pdf", p, height = 5, width = 6)

# Heatmap
# Clean all nan rows
keeppept = rowSums(!is.na(cordffilt[sublip,]))>0
keepmet = colSums(!is.na(cordffilt[sublip,]))>0
cordffilt_clean = cordffilt[sublip,][keeppept,keepmet]
cordf_clean = cordf[sublip,][keeppept,keepmet]
# Plot
row_dist <- dist(t(cordf_clean^2), method = "euclidean")
col_dist <- dist(cordf_clean^2, method = "euclidean")
row_hc <- hclust(row_dist, method = "ward.D2")
col_hc <- hclust(col_dist, method = "ward.D2")
pdf("figure_panels/cormatrix_heatmap_translation.pdf", width = 6, height = 4)
heatmap.2(t(cordffilt_clean^2), 
          trace = "none", col = Colors, margins = c(4, 4),
          breaks = seq(0, 1, 0.01), symkey = F, na.color = "white",
          Rowv = as.dendrogram(row_hc),
          Colv = as.dendrogram(col_hc), dendrogram = "both")
heatmap.2(t(cordf_clean^2), 
          trace = "none", col = Colors, margins = c(4, 4),
          breaks = seq(0, 1, 0.01), symkey = F, na.color = "white",
          Rowv = as.dendrogram(row_hc),
          Colv = as.dendrogram(col_hc), dendrogram = "both")
dev.off()

pdf("figure_panels/cormatrix_heatmap_translation_snapshots.pdf", width = 10, height = 10)
heatmap.2(t(cordffilt_clean^2)[row_hc$order, col_hc$order][37:1,124:128],
          trace = "none", col = Colors, margins = c(4, 4),
          breaks = seq(0, 1, 0.01), symkey = F, na.color = "white",
          Rowv = F, Colv = F, dendrogram = "none",
          labCol = sapply(colnames(t(cordffilt_clean^2)[row_hc$order, col_hc$order][37:1,124:128]), function(x) sprintf("%s_%s-%s", annotprot[x,"gene_name"], annotprot[x,"startPosition"], annotprot[x,"endPosition"])),
          labRow = sapply(rownames(t(cordffilt_clean^2)[row_hc$order, col_hc$order][37:1,124:128]), function(x) sprintf("%s|%s", x, annotmet[x,"id_top_name"])))
heatmap.2(t(cordffilt_clean^2)[row_hc$order, col_hc$order][37:1,179:210],
          trace = "none", col = Colors, margins = c(4, 4),
          breaks = seq(0, 1, 0.01), symkey = F, na.color = "white",
          Rowv = F, Colv = F, dendrogram = "none",
          labCol = sapply(colnames(t(cordffilt_clean^2)[row_hc$order, col_hc$order][37:1,179:210]), function(x) sprintf("%s_%s-%s", annotprot[x,"gene_name"], annotprot[x,"startPosition"], annotprot[x,"endPosition"])),
          labRow = sapply(rownames(t(cordffilt_clean^2)[row_hc$order, col_hc$order][37:1,179:210]), function(x) sprintf("%s|%s", x, annotmet[x,"id_top_name"])))
heatmap.2(t(cordffilt_clean^2)[row_hc$order, col_hc$order][93:87,14:60],
          trace = "none", col = Colors, margins = c(4, 4),
          breaks = seq(0, 1, 0.01), symkey = F, na.color = "white",
          Rowv = F, Colv = F, dendrogram = "none",
          labCol = sapply(colnames(t(cordffilt_clean^2)[row_hc$order, col_hc$order][93:87,14:60]), function(x) sprintf("%s_%s-%s", annotprot[x,"gene_name"], annotprot[x,"startPosition"], annotprot[x,"endPosition"])),
          labRow = sapply(rownames(t(cordffilt_clean^2)[row_hc$order, col_hc$order][93:87,14:60]), function(x) sprintf("%s|%s", x, annotmet[x,"id_top_name"])))
heatmap.2(t(cordffilt_clean^2)[row_hc$order, col_hc$order][170:139,1:13],
          trace = "none", col = Colors, margins = c(4, 4),
          breaks = seq(0, 1, 0.01), symkey = F, na.color = "white",
          Rowv = F, Colv = F, dendrogram = "none",
          labCol = sapply(colnames(t(cordffilt_clean^2)[row_hc$order, col_hc$order][170:139,1:13]), function(x) sprintf("%s_%s-%s", annotprot[x,"gene_name"], annotprot[x,"startPosition"], annotprot[x,"endPosition"])),
          labRow = sapply(rownames(t(cordffilt_clean^2)[row_hc$order, col_hc$order][170:139,1:13]), function(x) sprintf("%s|%s", x, annotmet[x,"id_top_name"])))
dev.off()

p = heatmaply(t(cordffilt_clean^2),
              Rowv = as.dendrogram(row_hc),
              Colv = as.dendrogram(col_hc),
              colors = Colors,
              na.value = "white",
              limits = c(0, 1),
              margins = c(40, 40),
              #showticklabels= F,
              labCol = sapply(rownames(cordffilt_clean), function(x) sprintf("%s_%s-%s", annotprot[x,"gene_name"], annotprot[x,"startPosition"], annotprot[x,"endPosition"])),
              labRow = sapply(colnames(cordffilt_clean), function(x) sprintf("%s|%s", x, annotmet[x,"id_top_name"])))
saveWidget(p, "plots/cormatrix_heatmap_translation.html", selfcontained = TRUE)

# Check AA synthetases
cyto_tRNAsyn = c("SES1","MES1","DPS1","GLN4","GRS1","HTS1","GUS1","CDC60","THS1","ILS1","KRS1",
                 "YHR020W","WRS1","ALA1","YDR341C","DED81","YNL247W","FRS1","FRS2","TYS1","VAS1")
sprintf("All gene ids (20 aaRS): %i. Detected: %s. Significant: %s. Correlating: %s", length(cyto_tRNAsyn),
        paste0(unique(annotprot[sapply(annotprot$gene_name, function(x) any(strsplit(x,";")[[1]] %in% cyto_tRNAsyn)),"gene_name"]),collapse=","),
        paste0(unique(annotprot[(sapply(annotprot$gene_name, function(x) any(strsplit(x,";")[[1]] %in% cyto_tRNAsyn)))&(rownames(annotprot) %in% siglip),"gene_name"]),collapse=","),
        paste0(unique(annotprot[(sapply(annotprot$gene_name, function(x) any(strsplit(x,";")[[1]] %in% cyto_tRNAsyn)))&(rownames(annotprot) %in% corsigfull[corsigfull$lfdr<0.05,"X1"]),"gene_name"]),collapse=","))

# Find protein regulated by metab
metab = "348.07040__2.226" #AMP
metab = "134.04470__1.513" #Asp
corpept = names(cordffilt[!is.na(cordffilt[,metab]),metab])
sort(table(annotprot[corpept,"gene_name"]), decreasing=T)
# Subset for GO term
sort(table(annotprot[corpept[corpept %in% sublip],"gene_name"]), decreasing=T)

# Check 1 protein
pdf("figure_panels/cormatrix_proteins.pdf", width = 8, height = 7)
gene="MES1"
subdf = cordffilt[rownames(cordffilt) %in% rownames(annotprot)[annotprot["gene_name"]==gene],, drop=F]^2
subdf = subdf[rowSums(!is.na(subdf))>0, colSums(!is.na(subdf))>0, drop=F]
subdf = subdf[order(sapply(rownames(subdf), function(x) as.numeric(strsplit(annotprot[x,"startPosition"],";")[[1]][1]))), # order by gene name
              order(sapply(colnames(subdf), function(x) as.numeric(strsplit(x,"__")[[1]][1]))), drop=F] # order by mass
heatmap.2(subdf, trace="none", col=Colors, margins = c(20,12),breaks=seq(0,1,0.01),
          labCol = sapply(colnames(subdf), function(x) sprintf("%s|%s",x,annotmet[x,"id_top_name"])),
          labRow = apply(annotprot[rownames(subdf),c("gene_name","startPosition","endPosition")],1, paste0, collapse="_"), Rowv = F, Colv = F)
gene="SES1"
subdf = cordffilt[rownames(cordffilt) %in% rownames(annotprot)[annotprot["gene_name"]==gene],, drop=F]^2
subdf = subdf[rowSums(!is.na(subdf))>0, colSums(!is.na(subdf))>0, drop=F]
subdf = subdf[order(sapply(rownames(subdf), function(x) as.numeric(strsplit(annotprot[x,"startPosition"],";")[[1]][1]))), # order by gene name
              order(sapply(colnames(subdf), function(x) as.numeric(strsplit(x,"__")[[1]][1]))), drop=F] # order by mass
heatmap.2(subdf, trace="none", col=Colors, margins = c(20,12),breaks=seq(0,1,0.01),
          labCol = sapply(colnames(subdf), function(x) sprintf("%s|%s",x,annotmet[x,"id_top_name"])),
          labRow = apply(annotprot[rownames(subdf),c("gene_name","startPosition","endPosition")],1, paste0, collapse="_"), Rowv = F, Colv = F)
gene="TEF1;TEF2"
subdf = cordffilt[rownames(cordffilt) %in% rownames(annotprot)[annotprot["gene_name"]==gene],]^2
subdf = subdf[rowSums(!is.na(subdf))>0, colSums(!is.na(subdf))>0, drop=F]
subdf = subdf[order(sapply(rownames(subdf), function(x) as.numeric(strsplit(annotprot[x,"startPosition"],";")[[1]][1]))), # order by gene name
              order(sapply(colnames(subdf), function(x) as.numeric(strsplit(x,"__")[[1]][1]))), drop=F] # order by mass
heatmap.2(subdf, trace="none", col=Colors, margins = c(20,12),breaks=seq(0,1,0.01),
          labCol = sapply(colnames(subdf), function(x) sprintf("%s|%s",x,annotmet[x,"id_top_name"])),
          labRow = apply(annotprot[rownames(subdf),c("gene_name","startPosition","endPosition")],1, paste0, collapse="_"), Rowv = F, Colv = F)
gene = "STM1"
subdf = cordffilt[rownames(cordffilt) %in% rownames(annotprot)[annotprot["gene_name"]==gene],]^2
subdf = subdf[rowSums(!is.na(subdf))>0, colSums(!is.na(subdf))>0, drop=F]
subdf = subdf[order(sapply(rownames(subdf), function(x) as.numeric(strsplit(annotprot[x,"startPosition"],";")[[1]][1]))), # order by gene name
              order(sapply(colnames(subdf), function(x) as.numeric(strsplit(x,"__")[[1]][1]))), drop=F] # order by mass
heatmap.2(subdf, trace="none", col=Colors, margins = c(20,12),breaks=seq(0,1,0.01),
          labCol = sapply(colnames(subdf), function(x) sprintf("%s|%s",x,annotmet[x,"id_top_name"])),
          labRow = apply(annotprot[rownames(subdf),c("gene_name","startPosition","endPosition")],1, paste0, collapse="_"), Rowv = F, Colv = F)
gene = "CDC60"
subdf = cordffilt[rownames(cordffilt) %in% rownames(annotprot)[annotprot["gene_name"]==gene],]^2
subdf = subdf[rowSums(!is.na(subdf))>0, colSums(!is.na(subdf))>0, drop=F]
subdf = subdf[order(sapply(rownames(subdf), function(x) as.numeric(strsplit(annotprot[x,"startPosition"],";")[[1]][1]))), # order by gene name
              order(sapply(colnames(subdf), function(x) as.numeric(strsplit(x,"__")[[1]][1]))), drop=F] # order by mass
heatmap.2(subdf, trace="none", col=Colors, margins = c(20,12),breaks=seq(0,1,0.01),
          labCol = sapply(colnames(subdf), function(x) sprintf("%s|%s",x,annotmet[x,"id_top_name"])),
          labRow = apply(annotprot[rownames(subdf),c("gene_name","startPosition","endPosition")],1, paste0, collapse="_"), Rowv = F, Colv = F)
dev.off()

# Wood's plot
jointdfall = rbind(lipdf,metdf[,colnames(lipdf)])
metab = "348.07040__2.226" #AMP
prot="MES1"
protall = rownames(annotprot)[grepl(prot,annotprot$gene_name)]
wooddf = as.data.frame(t(cor(t(jointdfall[c(metab,protall),]), use = "pairwise", method = "pearson")[1,2:(length(protall)+1), drop=F]^2))
wooddf["pv"] = sapply(rownames(wooddf), function(x) if(x %in% rownames(corsigdf)){corsigdf[x,metab]}else{1})
makeWoodsPlotSingleProtein(wooddf, annotprot, deltaColorIsTryptic=TRUE,
                           protName = prot, nameProtQuant="gene_name",
                           coefCol=metab, pvalCol="pv")
ggsave("figure_panels/MES1vsAMP_woodplot.pdf", width=8, height = 2)

prot="SES1"
protall = rownames(annotprot)[grepl(prot,annotprot$gene_name)]
wooddf = as.data.frame(t(cor(t(jointdfall[c(metab,protall),]), use = "pairwise", method = "pearson")[1,2:(length(protall)+1), drop=F]^2))
wooddf["pv"] = sapply(rownames(wooddf), function(x) if(x %in% rownames(corsigdf)){corsigdf[x,metab]}else{1})
makeWoodsPlotSingleProtein(wooddf, annotprot, deltaColorIsTryptic=TRUE,
                           protName = prot, nameProtQuant="gene_name",
                           coefCol=metab, pvalCol="pv")
ggsave("figure_panels/SES1vsAMP_woodplot.pdf", width=8, height = 2)

metab = "134.04470__1.513" #Asp
prot="TEF1"
protall = rownames(annotprot)[grepl(prot,annotprot$gene_name)]
wooddf = as.data.frame(t(cor(t(jointdfall[c(metab,protall),]), use = "pairwise", method = "pearson")[1,2:(length(protall)+1), drop=F]^2))
wooddf["pv"] = sapply(rownames(wooddf), function(x) if(x %in% rownames(corsigdf)){corsigdf[x,metab]}else{1})
makeWoodsPlotSingleProtein(wooddf, annotprot, deltaColorIsTryptic=TRUE,
                           protName = prot, nameProtQuant="gene_name",
                           coefCol=metab, pvalCol="pv")
ggsave("figure_panels/TEF1-2vsAsp_woodplot.pdf", width=8, height = 2)

metab = "134.04470__1.513" #Asp
prot="STM1"
protall = rownames(annotprot)[grepl(prot,annotprot$gene_name)]
wooddf = as.data.frame(t(cor(t(jointdfall[c(metab,protall),]), use = "pairwise", method = "pearson")[1,2:(length(protall)+1), drop=F]^2))
wooddf["pv"] = sapply(rownames(wooddf), function(x) if(x %in% rownames(corsigdf)){corsigdf[x,metab]}else{1})
makeWoodsPlotSingleProtein(wooddf, annotprot, deltaColorIsTryptic=TRUE,
                           protName = prot, nameProtQuant="gene_name",
                           coefCol=metab, pvalCol="pv")
ggsave("figure_panels/STM1vsAsp_woodplot.pdf", width=8, height = 2)

metab = "284.09890__4.351" #Guanosine
prot="CDC60"
protall = rownames(annotprot)[grepl(prot,annotprot$gene_name)]
wooddf = as.data.frame(t(cor(t(jointdfall[c(metab,protall),]), use = "pairwise", method = "pearson")[1,2:(length(protall)+1), drop=F]^2))
wooddf["pv"] = sapply(rownames(wooddf), function(x) if(x %in% rownames(corsigdf)){corsigdf[x,metab]}else{1})
makeWoodsPlotSingleProtein(wooddf, annotprot, deltaColorIsTryptic=TRUE,
                           protName = prot, nameProtQuant="gene_name",
                           coefCol=metab, pvalCol="pv")
ggsave("figure_panels/CDC60_Gua_woodplot.pdf", width=8, height = 2)

# Check specific peptide
#most correlating
metab = "348.07040__2.226" #AMP
prot="MES1"
pept = names(which.max(abs(cordffilt[rownames(cordffilt) %in% annotprot[annotprot$gene_name==prot,"Peptide"],metab])))
plotcor(jointdf,metab,pept)
ggsave("figure_panels/MES1vsAMP_correlation.pdf", width=8, height = 7)

# lineplot
plotdf = as.data.frame(t(as.data.frame(jointdf[c(metab,pept),])))
plotdf$labels = rownames(plotdf)
plotdf$cond = substr(plotdf$labels,1,1)
plotdf$time = as.numeric(substr(plotdf$labels,2,3))
plotdf = plotdf[plotdf$cond=="G",]
coef = 3e-8
ylim.prim = min(plotdf[,metab])-100
ylim.sec = min(plotdf[,pept], na.rm = T)-0.2
ggplot(plotdf, aes(x=time)) +
  geom_point(aes(y=get(metab)), color="#d6604dff") + stat_summary(aes(y=get(metab)), geom="line", fun = "mean", color="#d6604dff") +
  geom_point(aes(y=((get(pept)-ylim.sec)/coef)+ylim.prim), color="#558d78") + stat_summary(aes(y=((get(pept)-ylim.sec)/coef)+ylim.prim), geom="line", fun = "mean", color="#558d78") +
  scale_y_continuous(limits=c(ylim.prim,NA), sec.axis = sec_axis(~(.-ylim.prim)*coef+ylim.sec, name="pept")) +
  theme_test() + theme(axis.title.y.right = element_text(color = "#558d78"), axis.text.y.right = element_text(colour = "#558d78"))
ggsave("figure_panels/MES1vsAMP_G_lineplot.pdf", width=4, height = 2.5)

metab = "348.07040__2.226" #AMP
prot="SES1"
pept = names(which.max(abs(cordffilt[rownames(cordffilt) %in% annotprot[annotprot$gene_name==prot,"Peptide"],metab])))
plotcor(jointdf,metab,pept)
ggsave("figure_panels/SES1vsAMP_correlation.pdf", width=8, height = 7)

# lineplot
plotdf = as.data.frame(t(as.data.frame(jointdf[c(metab,pept),])))
plotdf$labels = rownames(plotdf)
plotdf$cond = substr(plotdf$labels,1,1)
plotdf$time = as.numeric(substr(plotdf$labels,2,3))
plotdf = plotdf[plotdf$cond=="G",]
coef = 7e-7
ylim.prim = min(plotdf[,metab])-100
ylim.sec = min(plotdf[,pept], na.rm = T)-0.2
ggplot(plotdf, aes(x=time)) +
  geom_point(aes(y=get(metab)), color="#d6604dff") + stat_summary(aes(y=get(metab)), geom="line", fun = "mean", color="#d6604dff") +
  geom_point(aes(y=((get(pept)-ylim.sec)/coef)+ylim.prim), color="#558d78") + stat_summary(aes(y=((get(pept)-ylim.sec)/coef)+ylim.prim), geom="line", fun = "mean", color="#558d78") +
  scale_y_continuous(limits=c(ylim.prim,NA), sec.axis = sec_axis(~(.-ylim.prim)*coef+ylim.sec, name="pept")) +
  theme_test() + theme(axis.title.y.right = element_text(color = "#558d78"), axis.text.y.right = element_text(colour = "#558d78"))
ggsave("figure_panels/SES1vsAMP_G_lineplot.pdf", width=4, height = 2.5)

metab = "134.04470__1.513" #Asp
prot="TEF1;TEF2"
pept = names(which.max(abs(cordffilt[rownames(cordffilt) %in% annotprot[annotprot$gene_name==prot,"Peptide"],metab])))
plotcor(jointdf,metab,pept)
ggsave("figure_panels/TEF1-2vsAsp_correlation.pdf", width=8, height = 7)

# lineplot
plotdf = as.data.frame(t(as.data.frame(jointdf[c(metab,pept),])))
plotdf$labels = rownames(plotdf)
plotdf$cond = substr(plotdf$labels,1,1)
plotdf$time = as.numeric(substr(plotdf$labels,2,3))
plotdf = plotdf[plotdf$cond=="G",]
coef = 1.2e-8
ylim.prim = min(plotdf[,metab])-100
ylim.sec = min(plotdf[,pept], na.rm = T)-0.1
ggplot(plotdf, aes(x=time)) +
  geom_point(aes(y=get(metab)), color="#d6604dff") + stat_summary(aes(y=get(metab)), geom="line", fun = "mean", color="#d6604dff") +
  geom_point(aes(y=((get(pept)-ylim.sec)/coef)+ylim.prim), color="#558d78") + stat_summary(aes(y=((get(pept)-ylim.sec)/coef)+ylim.prim), geom="line", fun = "mean", color="#558d78") +
  scale_y_continuous(limits=c(ylim.prim,NA), sec.axis = sec_axis(~(.-ylim.prim)*coef+ylim.sec, name="pept")) +
  theme_test() + theme(axis.title.y.right = element_text(color = "#558d78"), axis.text.y.right = element_text(colour = "#558d78"))
ggsave("figure_panels/TEF1-2vsAsp_G_lineplot.pdf", width=4, height = 2.5)


metab = "134.04470__1.513" #Asp
prot = "STM1"
pept = names(which.max(abs(cordffilt[rownames(cordffilt) %in% annotprot[annotprot$gene_name==prot,"Peptide"],metab])))
plotcor(jointdf,metab,pept)
ggsave("figure_panels/STM1vsAsp_correlation.pdf", width=8, height = 7)


metab = "284.09890__4.351" #Guanosine
prot="CDC60"
pept = names(which.max(abs(cordffilt[rownames(cordffilt) %in% annotprot[annotprot$gene_name==prot,"Peptide"],metab])))
plotcor(jointdf,metab,pept)
ggsave("figure_panels/CDC60vsGuan_correlation.pdf", width=8, height = 7)

# lineplot
plotdf = as.data.frame(t(as.data.frame(jointdf[c(metab,pept),])))
plotdf$labels = rownames(plotdf)
plotdf$cond = substr(plotdf$labels,1,1)
plotdf$time = as.numeric(substr(plotdf$labels,2,3))
plotdf = plotdf[plotdf$cond=="N",]
coef = 1e-6
ylim.prim = min(plotdf[,metab])-100
ylim.sec = min(plotdf[,pept], na.rm = T)-0.2
ggplot(plotdf, aes(x=time)) +
  geom_point(aes(y=get(metab)), color="#d6604dff") + stat_summary(aes(y=get(metab)), geom="line", fun = "mean", color="#d6604dff") +
  geom_point(aes(y=((get(pept)-ylim.sec)/coef)+ylim.prim), color="#558d78") + stat_summary(aes(y=((get(pept)-ylim.sec)/coef)+ylim.prim), geom="line", fun = "mean", color="#558d78") +
  scale_y_continuous(limits=c(ylim.prim,NA), sec.axis = sec_axis(~(.-ylim.prim)*coef+ylim.sec, name="pept")) +
  theme_test() + theme(axis.title.y.right = element_text(color = "#558d78"), axis.text.y.right = element_text(colour = "#558d78"))
ggsave("figure_panels/CDC60vsGuan_N_lineplot.pdf", width=4, height = 2.5)
