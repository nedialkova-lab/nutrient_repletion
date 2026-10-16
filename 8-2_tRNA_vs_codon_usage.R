# ==============================================================================
# 8-2_tRNA_vs_codon_usage.R
# Correlation of codon occupancies (short and long RPFs) with cognate tRNA abundance
# and aminoacylation (fdrtool q-values), and GSEA of genes ranked by codon frequency.
#
# Manuscript: Fig. 4D-E; Extended Data Fig. 5D-F
# Requires  : outputs of 7-3 (RxT) and 8-1 (tRNA tables)
# Inputs    : data/sacCer3.no_dubious_ORFs_coding.fa
#             data/pairing_rules_Johansson2008.txt
#             results/RxT.tsv
#             data/20241210_RNA_RPF/{glucose,nitrogen}_asite_occupancy.tsv
#             results/tRNA_{glucose,nitrogen}_anticodon_{abundance,CCAabundance,CCAperc}.tsv
# Outputs   : figure_panels/tRNA_vs_codons.pdf
#             figure_panels/<CODON>_scatters.pdf
#             figure_panels/GSEA_codons.pdf
#             results/codon_vs_tRNA_correlations.tsv
#             results/GSEA_codons_SerLeu.tsv
# Run from the repository root (all paths are relative to it).
# ==============================================================================

library(gplots)
library(ggplot2)
library(seqinr)
library(RColorBrewer)
library(reshape2)
library(ggpubr)
library(ggrepel)
library(patchwork)
library(clusterProfiler)
library(org.Sc.sgd.db)
library(limma)
uniColors=brewer.pal(9,"YlOrRd")
uniColors=colorRampPalette(uniColors)(100)
biColors=rev(brewer.pal(11,"RdBu"))
biColors=colorRampPalette(biColors)(100)

# Load codon usage
cds = read.fasta("data/sacCer3.no_dubious_ORFs_coding.fa")
codons = lapply(cds, function(x) sapply(seq(1,length(x),3), function(p) toupper(paste0(x[p:(p+2)], collapse = ""))))
codonlabels = c("ATG","GTC","AAA","TTA","ACT","TCA","ATC","GCC","GCT","GGT","TCT","GCA","ACC","CTA","CAA","GAC","GAA","AGA",
                "AAC","TTG","GTG","TAC","GAT","CAC","TTC","CCA","GTT","ATT","TGG","TCC","AGC","AAG","AGT","ATA","GGG","CCG",
                "AGG","TGC","CTT","CCC","TTT","CCT","ACA","TCG","AAT","CGA","CGT","TAT","TGT","CTG","CGG","CAT","CTC","GGA",
                "GTA","GAG","CAG","ACG","GGC","GCG","CGC")
codus = data.frame(row.names = codonlabels)
glen = c()
for (g in names(codons)){
  genename = strsplit(g,"_")[[1]][1]
  codus[genename] = sapply(codonlabels, function(x) sum(codons[[g]]==x))
  glen[genename] = length(codons[[g]])*3/1000
}
codus = t(codus)
codfreq = codus/rowSums(codus)

# Paríring rules
pairing = read.csv("data/pairing_rules_Johansson2008.txt", sep = "\t")

# Load RxT
rxts = read.csv("results/RxT.tsv", sep = "\t")

# Loop through each condition in rxts
coefdf = data.frame(row.names = colnames(codfreq))
conds = unique(sapply(colnames(rxts), function(x) strsplit(x,"_")[[1]][1]))
for (cond in conds) {
  # Repeat traincodus twice and convert to data frame
  traintemp <- rbind(codfreq,codfreq)
  # Subtract log of genelength from rxts and handle missing values
  rxttemp <- log(rxts[rownames(codfreq), grepl(sprintf("%s_",cond), colnames(rxts))]) - log(glen[rownames(codfreq)])
  rxttemp <- unlist(rxttemp)
  isnan <- is.na(rxttemp)
  rxttemp <- rxttemp[!isnan]
  traintemp <- traintemp[!isnan, ]
  # Fit linear regression model
  for (cod in colnames(traintemp)){
    reg <- lm(rxttemp ~ ., data = as.data.frame(traintemp[,cod]))
    # Store coefficients
    coefdf[cod,cond] <- coef(reg)[-1]
  }
}

# Genomic codon usage
codorder = names(sort(colSums(codus)))
ggplot(data.frame(codon=factor(codorder,levels=codorder), codfreq=sort(colSums(codus)),row.names = codorder), aes(x=codon,y=codfreq)) +
  geom_bar(stat="identity") + theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# Codon pause scores
AsiteG = read.csv("data/20241210_RNA_RPF/glucose_asite_occupancy.tsv", sep="\t", row.names = 1)
samplecond = sapply(rownames(AsiteG), function(x) substr(x,1,(nchar(x)-1)))
AsiteG = aggregate(log2(AsiteG), list(samplecond), mean)
rownames(AsiteG) = AsiteG$Group.1; AsiteG = AsiteG[,-c(1)]
AsiteG = AsiteG[c("G0_short","G1_short","G2_short","G4_short","G7_short","G10_short",
                  "G0_long","G1_long","G2_long","G4_long","G7_long","G10_long" ),]


AsiteN = read.csv("data/20241210_RNA_RPF/nitrogen_asite_occupancy.tsv", sep="\t", row.names = 1)
samplecond = sapply(rownames(AsiteN), function(x) substr(x,1,(nchar(x)-1)))
AsiteN = aggregate(log2(AsiteN), list(samplecond), mean)
rownames(AsiteN) = AsiteN$Group.1; AsiteN= AsiteN[,-c(1)]
AsiteN = AsiteN[c("N0_short","N1_short","N2_short","N4_short","N7_short","N10_short",
                  "N0_long","N1_long","N2_long","N4_long","N7_long","N10_long" ),]

#AsiteG = log2(read.csv("data/20241210_RNA_RPF/glucose_pause_Asite_all.tsv", sep="\t", row.names = 1))
#AsiteG = AsiteG[c(grep("short",rownames(AsiteG))[c(1,3,4,5,6,2)],grep("long",rownames(AsiteG))[c(1,3,4,5,6,2)]),]
#AsiteN = log2(read.csv("data/20241210_RNA_RPF/nitrogen_pause_Asite_all.tsv", sep="\t", row.names = 1))
#AsiteN = AsiteN[c(grep("short",rownames(AsiteN))[c(1,3,4,5,6,2)],grep("long",rownames(AsiteN))[c(1,3,4,5,6,2)]),]

# tRNAs
percdfG = read.csv("results/tRNA_glucose_anticodon_abundance.tsv", sep = "\t")
chargedpercdfG = read.csv("results/tRNA_glucose_anticodon_CCAabundance.tsv", sep = "\t")
anticcaG = read.csv("results/tRNA_glucose_anticodon_CCAperc.tsv", sep = "\t")
percdfN = read.csv("results/tRNA_nitrogen_anticodon_abundance.tsv", sep = "\t")
chargedpercdfN = read.csv("results/tRNA_nitrogen_anticodon_CCAabundance.tsv", sep = "\t")
anticcaN = read.csv("results/tRNA_nitrogen_anticodon_CCAperc.tsv", sep = "\t")
# Remove iMet
percdfG = percdfG[rownames(percdfG)!="iMet-CAT",]
chargedpercdfG = chargedpercdfG[rownames(chargedpercdfG)!="iMet-CAT",]
anticcaG = anticcaG[rownames(anticcaG)!="iMet-CAT",]
percdfN = percdfN[rownames(percdfN)!="iMet-CAT",]
chargedpercdfN = chargedpercdfN[rownames(chargedpercdfN)!="iMet-CAT",]
anticcaN = anticcaN[rownames(anticcaN)!="iMet-CAT",]

# ALL CONDITIONS TOGETHER
percdf = cbind(percdfG,percdfN[rownames(percdfG),])
chargedpercdf = cbind(chargedpercdfG,chargedpercdfN[rownames(chargedpercdfG),])
anticca = cbind(anticcaG,anticcaN[rownames(anticcaG),])
Asite = rbind(AsiteG,AsiteN[,colnames(AsiteG)])
# Create long df
condtRNA = sapply(colnames(percdf), function(x) strsplit(x,"_")[[1]][1])
allcond = conds[grepl("g",conds)|grepl("n",conds)]
commoncond = allcond[allcond %in% condtRNA]
pairs = c("short","long","coef","codfreq")
ncors = length(allcond)*length(pairs)
longdf = data.frame()
for (cod in codorder){
  cognates = strsplit(pairing[pairing$CODON==cod,"ANTICODON"],",")[[1]]
  tempdf = data.frame(codon = rep(cod,ncors*length(cognates)),
                      anticodon = c(sapply(cognates, function(x) rep(x,ncors))),
                      cond = c(sapply(cognates, function(y) sapply(allcond, function(x) rep(x,length(pairs))))),
                      pair = c(sapply(cognates, function(y) rep(pairs,length(allcond)))))
  for (row in rownames(tempdf)){
    tempdf[row,"anti_abund"] = mean(as.numeric(percdf[grepl(tempdf[row,"anticodon"],rownames(percdf)),
                                                       grepl(sprintf("%s_",tempdf[row,"cond"]),colnames(percdf))]))
    tempdf[row,"anti_corabund"] = mean(as.numeric(chargedpercdf[grepl(tempdf[row,"anticodon"],rownames(chargedpercdf)),
                                                                 grepl(sprintf("%s_",tempdf[row,"cond"]),colnames(chargedpercdf))]))
    tempdf[row,"anti_cca"] = mean(as.numeric(anticca[grepl(tempdf[row,"anticodon"],rownames(anticca)),
                                                      grepl(sprintf("%s_",tempdf[row,"cond"]),colnames(anticca))]))
    if (tempdf[row,"pair"]=="short"){
      tempdf[row,"codval"] = Asite[grepl(sprintf("%s_",tempdf[row,"cond"]),tolower(rownames(Asite)))&grepl("short",tolower(rownames(Asite))),
                                    tempdf[row,"codon"]]
    } else if (tempdf[row,"pair"]=="long"){
      tempdf[row,"codval"] = Asite[grepl(sprintf("%s_",tempdf[row,"cond"]),tolower(rownames(Asite)))&grepl("long",tolower(rownames(Asite))),
                                    tempdf[row,"codon"]]
    } else if (tempdf[row,"pair"]=="coef"){
      tempdf[row,"codval"] = coefdf[tempdf[row,"codon"], tempdf[row,"cond"]]
    } else if (tempdf[row,"pair"]=="codfreq"){
      tempdf[row,"codval"] = colSums(codus)[tempdf[row,"codon"]]
    }
  }
  # Center per cond
  #tempdf$anti_abund = apply(tempdf,1,function(x) as.numeric(x["anti_abund"]) - mean(tempdf[grepl(substr(x["cond"],1,1),tempdf$cond)&grepl(x["anticodon"],tempdf$anticodon),"anti_abund"], na.rm=T))
  #tempdf$anti_corabund = apply(tempdf,1,function(x) as.numeric(x["anti_corabund"]) - mean(tempdf[grepl(substr(x["cond"],1,1),tempdf$cond)&grepl(x["anticodon"],tempdf$anticodon),"anti_corabund"], na.rm=T))
  #tempdf$anti_cca = apply(tempdf,1,function(x) as.numeric(x["anti_cca"]) - mean(tempdf[grepl(substr(x["cond"],1,1),tempdf$cond)&grepl(x["anticodon"],tempdf$anticodon),"anti_cca"], na.rm=T))
  #tempdf$codval = apply(tempdf,1,function(x) as.numeric(x["codval"]) - mean(tempdf[grepl(substr(x["cond"],1,1),tempdf$cond)&grepl(x["codon"],tempdf$codon)&grepl(x["pair"],tempdf$pair),"codval"], na.rm=T))
  # Record
  longdf = rbind(longdf,tempdf)
}
longdf$codanti = apply(longdf,1,function(x) sprintf("%s-%s",x["codon"],x["anticodon"]))
longdf$time = as.numeric(substr(longdf$cond,2,4))
nonandf = longdf[rowSums(is.na(longdf))==0,]

# Correlate per codon
cordf1 = as.data.frame(t(sapply(unique(nonandf$codanti), function(x) sapply(c("short","long","coef"), function(y) cor(nonandf[(nonandf$codanti==x)&(nonandf$pair==y),"anti_abund"],
                                                                                                                    nonandf[(nonandf$codanti==x)&(nonandf$pair==y),"codval"])))))
cordf1test = as.data.frame(t(sapply(unique(nonandf$codanti), function(x) sapply(c("short","long","coef"), function(y) cor.test(nonandf[(nonandf$codanti==x)&(nonandf$pair==y),"anti_abund"],
                                                                                                                             nonandf[(nonandf$codanti==x)&(nonandf$pair==y),"codval"])$p.value))))
cordf2 = as.data.frame(t(sapply(unique(nonandf$codanti), function(x) sapply(c("short","long","coef"), function(y) cor(nonandf[(nonandf$codanti==x)&(nonandf$pair==y),"anti_cca"],
                                                                                                                    nonandf[(nonandf$codanti==x)&(nonandf$pair==y),"codval"])))))
cordf2test = as.data.frame(t(sapply(unique(nonandf$codanti), function(x) sapply(c("short","long","coef"), function(y) cor.test(nonandf[(nonandf$codanti==x)&(nonandf$pair==y),"anti_cca"],
                                                                                                                             nonandf[(nonandf$codanti==x)&(nonandf$pair==y),"codval"])$p.value))))
cordf1$tRNA = "abund"; cordf1$codanti = rownames(cordf1)
cordf1test$tRNA = "abund"; cordf1test$codanti = rownames(cordf1test)
cordf2$tRNA = "CCA"; cordf2$codanti = rownames(cordf2)
cordf2test$tRNA = "CCA"; cordf2test$codanti = rownames(cordf2test)
cordf = rbind(cordf1,cordf2)
cordftest = rbind(cordf1test,cordf2test)

cordftest[,c("short","long")] = fdrtool::fdrtool(c(cordftest$short,cordftest$long), statistic = "pvalue", plot = T)$qval
cordflong = dcast(melt(cordf), "codanti+variable ~ tRNA")
cordftestlong = dcast(melt(cordftest), "codanti+variable ~ tRNA")
# Export codon-anticodon correlations (Pearson R and fdrtool q-values) for the Supplementary Table
write.table(merge(cordflong, cordftestlong, by = c("codanti","variable"), suffixes = c("_R","_qval")), "results/codon_vs_tRNA_correlations.tsv", sep="\t", quote = F, row.names = F)
cordflong$sig = apply(cordftestlong[,c("abund","CCA")],1,function(x) any(x<0.1))
cordflong$aacod = sapply(cordflong$codanti, function(x) sprintf("%s%s",pairing[pairing$CODON==substr(x,1,3),"AMINOACID"],x))

ggplot(cordflong,aes(x=abund,y=CCA,label=aacod,color=sig)) +
  facet_wrap("variable") +
  geom_point() + 
  scale_color_manual(values=c("lightgrey","#558d78")) +
  geom_text_repel(color="black", data = cordflong[cordflong$sig,], size=2) +
  xlim(-1,1) + ylim(-1,1) +
  theme_bw()
ggsave("figure_panels/tRNA_vs_codons.pdf", width=15, height=5)

# Plot 1 codon
codon = "TCC"
plotdf = nonandf[(nonandf$codon==codon)&(nonandf$pair=="short"),]
plotdf$cset = substr(plotdf$cond,1,1)
p1 = ggplot(plotdf, aes(x=codval,y=anti_abund, shape = anticodon)) +
  geom_smooth(color="black", method="lm", se=T) +
  geom_point(aes(color=cset), alpha=0.6) +
  stat_cor() +
  geom_text_repel(aes(label=cond)) +
  theme_classic()
p2 = ggplot(plotdf, aes(x=codval,y=anti_cca, shape = anticodon)) +
  geom_smooth(color="black", method="lm", se=T) +
  geom_point(aes(color=cset), alpha=0.6) +
  stat_cor() +
  geom_text_repel(aes(label=cond)) +
  theme_classic()
plotdf = nonandf[(nonandf$codon==codon)&(nonandf$pair=="long"),]
plotdf$cset = substr(plotdf$cond,1,1)
p3 = ggplot(plotdf, aes(x=codval,y=anti_abund, shape = anticodon)) +
  geom_smooth(color="black", method="lm", se=T) +
  geom_point(aes(color=cset), alpha=0.6) +
  stat_cor() +
  geom_text_repel(aes(label=cond)) +
  theme_classic()
p4 = ggplot(plotdf, aes(x=codval,y=anti_cca, shape = anticodon)) +
  geom_smooth(color="black", method="lm", se=T) +
  geom_point(aes(color=cset), alpha=0.6) +
  stat_cor() +
  geom_text_repel(aes(label=cond)) +
  theme_classic()
(p1 | p2)/(p3 | p4)
ggsave("figure_panels/TCC_scatters.pdf", width = 8, height = 5.5)

codon = "TTA"
plotdf = nonandf[(nonandf$codon==codon)&(nonandf$pair=="short"),]
plotdf$cset = substr(plotdf$cond,1,1)
p1 = ggplot(plotdf, aes(x=codval,y=anti_abund, shape = anticodon)) +
  geom_smooth(color="black", method="lm", se=T) +
  geom_point(aes(color=cset), alpha=0.6) +
  stat_cor() +
  geom_text_repel(aes(label=cond)) +
  theme_classic()
p2 = ggplot(plotdf, aes(x=codval,y=anti_cca, shape = anticodon)) +
  geom_smooth(color="black", method="lm", se=T) +
  geom_point(aes(color=cset), alpha=0.6) +
  stat_cor() +
  geom_text_repel(aes(label=cond)) +
  theme_classic()
plotdf = nonandf[(nonandf$codon==codon)&(nonandf$pair=="long"),]
plotdf$cset = substr(plotdf$cond,1,1)
p3 = ggplot(plotdf, aes(x=codval,y=anti_abund, shape = anticodon)) +
  geom_smooth(color="black", method="lm", se=T) +
  geom_point(aes(color=cset), alpha=0.6) +
  stat_cor() +
  geom_text_repel(aes(label=cond)) +
  theme_classic()
p4 = ggplot(plotdf, aes(x=codval,y=anti_cca, shape = anticodon)) +
  geom_smooth(color="black", method="lm", se=T) +
  geom_point(aes(color=cset), alpha=0.6) +
  stat_cor() +
  geom_text_repel(aes(label=cond)) +
  theme_classic()
(p1 | p2)/(p3 | p4)
ggsave("figure_panels/TTA_scatters.pdf", width = 8, height = 5.5)

codon = "CGG"
plotdf = nonandf[(nonandf$codon==codon)&(nonandf$pair=="short"),]
plotdf$cset = substr(plotdf$cond,1,1)
p1 = ggplot(plotdf, aes(x=codval,y=anti_abund, shape = anticodon)) +
  geom_smooth(color="black", method="lm", se=T) +
  geom_point(aes(color=cset), alpha=0.6) +
  stat_cor() +
  geom_text_repel(aes(label=cond)) +
  theme_classic()
p2 = ggplot(plotdf, aes(x=codval,y=anti_cca, shape = anticodon)) +
  geom_smooth(color="black", method="lm", se=T) +
  geom_point(aes(color=cset), alpha=0.6) +
  stat_cor() +
  geom_text_repel(aes(label=cond)) +
  theme_classic()
plotdf = nonandf[(nonandf$codon==codon)&(nonandf$pair=="long"),]
plotdf$cset = substr(plotdf$cond,1,1)
p3 = ggplot(plotdf, aes(x=codval,y=anti_abund, shape = anticodon)) +
  geom_smooth(color="black", method="lm", se=T) +
  geom_point(aes(color=cset), alpha=0.6) +
  stat_cor() +
  geom_text_repel(aes(label=cond)) +
  theme_classic()
p4 = ggplot(plotdf, aes(x=codval,y=anti_cca, shape = anticodon)) +
  geom_smooth(color="black", method="lm", se=T) +
  geom_point(aes(color=cset), alpha=0.6) +
  stat_cor() +
  geom_text_repel(aes(label=cond)) +
  theme_classic()
(p1 | p2)/(p3 | p4)
ggsave("figure_panels/CGG_scatters.pdf", width = 8, height = 5.5)

codon = "CAA"
plotdf = nonandf[(nonandf$codon==codon)&(nonandf$pair=="short"),]
plotdf$cset = substr(plotdf$cond,1,1)
p1 = ggplot(plotdf, aes(x=codval,y=anti_abund, shape = anticodon)) +
  geom_smooth(color="black", method="lm", se=T) +
  geom_point(aes(color=cset), alpha=0.6) +
  stat_cor() +
  geom_text_repel(aes(label=cond)) +
  theme_classic()
p2 = ggplot(plotdf, aes(x=codval,y=anti_cca, shape = anticodon)) +
  geom_smooth(color="black", method="lm", se=T) +
  geom_point(aes(color=cset), alpha=0.6) +
  stat_cor() +
  geom_text_repel(aes(label=cond)) +
  theme_classic()
plotdf = nonandf[(nonandf$codon==codon)&(nonandf$pair=="long"),]
plotdf$cset = substr(plotdf$cond,1,1)
p3 = ggplot(plotdf, aes(x=codval,y=anti_abund, shape = anticodon)) +
  geom_smooth(color="black", method="lm", se=T) +
  geom_point(aes(color=cset), alpha=0.6) +
  stat_cor() +
  geom_text_repel(aes(label=cond)) +
  theme_classic()
p4 = ggplot(plotdf, aes(x=codval,y=anti_cca, shape = anticodon)) +
  geom_smooth(color="black", method="lm", se=T) +
  geom_point(aes(color=cset), alpha=0.6) +
  stat_cor() +
  geom_text_repel(aes(label=cond)) +
  theme_classic()
(p1 | p2)/(p3 | p4)
ggsave("figure_panels/CAA_scatters.pdf", width = 8, height = 5.5)

# GO ontology
gotab = select(org.Sc.sgd.db, rownames(codfreq), columns=c("GOALL","ORF"), keytype="ORF")
mapgo = go2term(unique(gotab$GO))
gotab$goterm = as.character(sapply(gotab$GO, function(x) mapgo[mapgo$go_id==x,"Term"]))
gotab_bp = unique(gotab[gotab$ONTOLOGY=="BP",c("goterm","ORF")])
gotab_cc = unique(gotab[gotab$ONTOLOGY=="CC",c("goterm","ORF")])
gotab_mf = unique(gotab[gotab$ONTOLOGY=="MF",c("goterm","ORF")])

keggtab <- getGeneKEGGLinks(species="sce")
mappaths = getKEGGPathwayNames(species="sce")
keggtab$pathname = sapply(keggtab$PathwayID, function(x) substr(mappaths[mappaths$PathwayID==x,"Description"],0, nchar(mappaths[mappaths$PathwayID==x,"Description"])-43))
keggtab = unique(keggtab[,c("pathname","GeneID")])

#### GO enrichment ####
gsea_results <- list()
for (codon in c("TCT","TCC","TCA","TCG","AGT","AGC","TTA","TTG","CTT","CTC","CTA","CTG")) {
  gene_list <- sort(codfreq[,codon], decreasing = T)
  gsea_results[[codon]] <- GSEA(gene_list, pvalueCutoff = 0.05, TERM2GENE = gotab_bp)
}

# Combine results for plotting
combined <- merge_result(gsea_results)
write.table(combined@compareClusterResult, "results/GSEA_codons_SerLeu.tsv", sep = "\t")

# Plot
combined@compareClusterResult$logp <- -log10(combined@compareClusterResult$p.adjust)
dotplot(combined,
        showCategory = 2,
        split = ".sign",
        title = "GSEA",
        by = "p.adjust",
        size = "logp",
        color="NES") + 
  scale_fill_gradientn(colors = biColors, name = "NES",
                       limits  = c(-3, 3),
                       values  = scales::rescale(c(-3, 1, 3), to = c(0, 1)),
                       oob = scales::squish) 
ggsave("figure_panels/GSEA_codons.pdf", width = 10, height = 10)
