# ==============================================================================
# 7-6_initiation_vs_elongation.R
# A-site codon occupancies, gene elongation-time estimates (geometric mean of codon
# occupancies) and their correlation with ribosomes per kb, compared with an
# initiation estimate from the 5'UTR MPRA of Akirtava et al. 2024.
#
# Manuscript: Fig. 2D-E; Extended Data Fig. 3E
# Requires  : outputs of 7-3 (RxT) and 7-5 (A-site occupancies, provided on Zenodo); 2-2 and 4-1 for the exploratory second part
# Inputs    : data/sacCer3.no_dubious_ORFs_coding.fa
#             results/RxT.tsv
#             data/20241210_RNA_RPF/{glucose,nitrogen}_asite_occupancy.tsv
#             data/20240529_RNA_RPF/phosphate_asite_occupancy.tsv
#             data/Akirtava_et_al_2024_Supp_Data_Rev1.xlsx
# Outputs   : figure_panels/Asite_occupancies.pdf
#             figure_panels/elongation_vs_inittiation.pdf
#             figure_panels/RxT_vs_elong_G0_G10_rep4.pdf
#             results/elongation_initiation_estimates.tsv
#             results/elongation_initiation_correlations.tsv
# Run from the repository root (all paths are relative to it).
# ==============================================================================

library(gplots)
library(ggplot2)
library(seqinr)
library(RColorBrewer)
library(reshape)
library(ggpubr)
library(ggrepel)
library(patchwork)
library(MetaboAnalystR)
library(limma)
library(clusterProfiler)
library(org.Sc.sgd.db)

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
gencodus = sort(colSums(codus))

# Load RxT
rxts = read.csv("results/RxT.tsv", sep = "\t")
rxts = sapply(colnames(rxts), function(x) log2(rxts[rownames(codus),x]/(rowSums(codus)*3/1000)))

# Codon pause scores
AsiteG = read.csv("data/20241210_RNA_RPF/glucose_asite_occupancy.tsv", sep="\t", row.names = 1)
AsiteN = read.csv("data/20241210_RNA_RPF/nitrogen_asite_occupancy.tsv", sep="\t", row.names = 1)
AsiteP = read.csv("data/20240529_RNA_RPF/phosphate_asite_occupancy.tsv", sep="\t", row.names = 1)
Asite = rbind(AsiteG, AsiteN, AsiteP)

# Plot avg
conds = sapply(rownames(Asite), function(x) substr(x,1,nchar(x)-1))
Asite_avg = t(sapply(unique(conds),function(x) colMeans(log2(Asite[conds==x,names(gencodus)]))))
Colors=rev(brewer.pal(11,"RdBu"))
Colors=colorRampPalette(Colors)(100)
pdf("figure_panels/Asite_occupancies.pdf", width = 10, height = 6)
heatmap.2(Asite_avg, trace="none", Rowv=F, Colv=F, col=Colors, margins = c(4,5),breaks=seq(-2,2,0.04),symkey=F)
dev.off()

# Compute gene occupancy based on the codus * codon occuopancies
geneocc = log2(apply(codus,1,function(x) apply(sapply(names(x), function(y) (Asite[,y]^(1/sum(x)))^x[y]),1,prod)))
rownames(geneocc) = rownames(Asite)
normgeneocc = t(apply(geneocc,1,scale))
colnames(normgeneocc) = rownames(codus)

# Codon optimality
codweights = gencodus/max(gencodus)
codopt = data.frame(scale(apply(codus,1,function(x) prod(sapply(names(x), function(y) (codweights[y]^(1/sum(x)))^x[y])))))
colnames(codopt) = "codopt"

# Initiation scores from Akirtava et al (2025)
utrs = readxl::read_excel("data/Akirtava_et_al_2024_Supp_Data_Rev1.xlsx", sheet = 2)
utrs$gene = do.call(rbind, sapply(utrs$`@UTR`, strsplit, "\\|"))[,1]

# Keep average for each gene, normalize by gene length
initscore = data.frame(scale(sapply(rownames(codus), function(x) mean(unlist(utrs[utrs$gene==x,"MeanYFP"]),na.rm=T))))
colnames(initscore) = "initscore"

# Kozak
features = readxl::read_excel("data/Akirtava_et_al_2024_Supp_Data_Rev1.xlsx", sheet = 4)
features$gene = do.call(rbind, sapply(features$Gene, strsplit, "\\|"))[,1]
kozakscore = data.frame(scale(sapply(rownames(codus), function(x) mean(unlist(features[features$gene==x,"LSMKozakStart"]),na.rm=T))))
colnames(kozakscore) = "kozakscore"
# Export elongation (per sample) and initiation (per gene) estimates for the Supplementary Table
write.table(cbind(t(geneocc), initscore[colnames(geneocc),, drop=F]), "results/elongation_initiation_estimates.tsv", sep="\t", quote = F)

# Compute initiation and elongation contribution per condition
fulldf = cbind(rxts[colnames(geneocc),], t(normgeneocc), initscore[colnames(geneocc),, drop=F], 
               kozakscore[colnames(geneocc),, drop=F], codopt[colnames(geneocc),, drop=F])
cordf = data.frame(row.names = colnames(rxts))
coefdf = data.frame(row.names = colnames(rxts))
for (c in colnames(rxts)){
  repnum = if(substr(c,nchar(c),nchar(c))=="4"){"2"}else{"1"}
  cond = strsplit(c,"_")[[1]][1]
  cordf[c,"short"] = cor(as.numeric(fulldf[,c]),as.numeric(fulldf[,sprintf("%s_%s%s",toupper(cond),"short",repnum)]), use = "pairwise", method = "pearson")
  cordf[c,"long"] = cor(as.numeric(fulldf[,c]),as.numeric(fulldf[,sprintf("%s_%s%s",toupper(cond),"long",repnum)]), use = "pairwise", method = "pearson")
  cordf[c,"init"] = cor(as.numeric(fulldf[,c]),as.numeric(fulldf[,"initscore"]), use = "pairwise", method = "pearson")
  cordf[c,"short_p"] = cor.test(as.numeric(fulldf[,c]),as.numeric(fulldf[,sprintf("%s_%s%s",toupper(cond),"short",repnum)]), use = "pairwise", method = "pearson")$p.value
  cordf[c,"long_p"] = cor.test(as.numeric(fulldf[,c]),as.numeric(fulldf[,sprintf("%s_%s%s",toupper(cond),"long",repnum)]), use = "pairwise", method = "pearson")$p.value
  cordf[c,"init_p"] = cor.test(as.numeric(fulldf[,c]),as.numeric(fulldf[,"initscore"]), use = "pairwise", method = "pearson")$p.value
  #cordf[c,"kozak"] = cor(as.numeric(fulldf[,c]),as.numeric(fulldf[,"kozakscore"]), use = "pairwise", method = "pearson")
  #cordf[c,"codopt"] = cor(as.numeric(fulldf[,c]),as.numeric(fulldf[,"codopt"]), use = "pairwise", method = "pearson")
  longlm = lm(sprintf("%s ~ initscore + %s_%s%s",c,toupper(cond),"long",repnum), fulldf)
  shortlm = lm(sprintf("%s ~ initscore + %s_%s%s",c,toupper(cond),"short",repnum), fulldf)
  cordf[c,"L:init+long"] = sqrt(summary(longlm)$r.squared)
  cordf[c,"S:init+short"] = sqrt(summary(shortlm)$r.squared)
  coefdf[c,"init_L"] = coef(longlm)[2]
  coefdf[c,"long_L"] = coef(longlm)[3]
  coefdf[c,"init_S"] = coef(shortlm)[2]
  coefdf[c,"short_S"] = coef(shortlm)[3]
}
write.table(cbind(cordf, coefdf[rownames(cordf),]), "results/elongation_initiation_correlations.tsv", sep="\t", quote = F)

# Plot cors
cordf$sample = rownames(cordf)
meltdf = melt(cordf[,!grepl("_",colnames(cordf))], id.vars = "sample")
meltdf$time = sapply(meltdf$sample, function(x) as.numeric(substr(strsplit(x,"_")[[1]][1],2,4)))
meltdf$cond = substr(meltdf$sample,1,1)
meltdf$pval = apply(meltdf,1,function(x) if(sprintf("%s_p",x["variable"]) %in% colnames(cordf)){unlist(cordf[x["sample"],sprintf("%s_p",x["variable"])])}else{NA})
meltdf["padj"] = p.adjust(meltdf$pval,method="bonferroni")
p1= ggplot(meltdf[!grepl(":",meltdf$variable),], aes(x=time, y=value, color=variable)) + 
  facet_grid("~cond") +
  geom_hline(yintercept=0, linewidth=0.3) +
  geom_point(aes(shape=padj<0.05)) + 
  stat_summary(geom="line", linewidth=1, fun = "mean") +
  scale_color_manual(values = c(init="grey",short="#d7614d",long="#4b7b6a")) +
  scale_shape_manual(values = c(1,16)) +
  ylab("R_Pearson") + theme_test()
p2= ggplot(meltdf, aes(x=time, y=value^2, color=variable, linetype=variable)) + 
  facet_grid("~cond") +
  geom_point() + 
  stat_summary(geom="line", linewidth=1, fun = "mean") +
  scale_color_manual(values = c(init="grey",short="#d7614d",long="#4b7b6a",`S:init+short`="#d7614d", `L:init+long`="#4b7b6a")) + 
  scale_linetype_manual(values = c(init="solid",short="solid",long="solid",`S:init+short`="dashed", `L:init+long`="dashed")) + 
  ylab("R2") + theme_bw()

# Plot coef
coefdf$sample = rownames(coefdf)
meltdf = melt(coefdf, id.vars = "sample")
meltdf$time = sapply(meltdf$sample, function(x) as.numeric(substr(strsplit(x,"_")[[1]][1],2,4)))
meltdf$cond = substr(meltdf$sample,1,1)
p3=ggplot(meltdf, aes(x=time, y=value, color=variable, linetype=variable)) + 
  facet_grid("~cond") +
  geom_point() + 
  stat_summary(geom="line", linewidth=1, fun = "mean") +
  scale_color_manual(values = c(short_S="#d7614d",long_L="#4b7b6a",init_S="#d7614d", init_L="#4b7b6a")) + 
  scale_linetype_manual(values = c(short_S="solid",long_L="solid",init_S="dashed", init_L="dashed")) + 
  ylab("Coefficient (MLR)") + theme_bw()
p1 / p2 / p3
ggsave("figure_panels/elongation_vs_inittiation.pdf", width = 10, height = 8)

plotdf = cbind(rxts[colnames(geneocc),], t(geneocc))
p1 = ggplot(plotdf, aes(x=G0_short2,y=g0_rep4)) + 
  geom_point(alpha=0.1, size=0.5, color="#ca403c", shape=16) + geom_smooth(method = "lm", color="black", se=T) + stat_cor() + theme_bw()
sum(rowSums(!is.na(plotdf[,c("G0_short2","g0_rep4")]))==2)
p2 = ggplot(plotdf, aes(x=G10_short2,y=g10_rep4)) + 
  geom_point(alpha=0.1, size=0.5, color="#ca403c", shape=16) + geom_smooth(method = "lm", color="black", se=T) + stat_cor() + theme_bw()
sum(rowSums(!is.na(plotdf[,c("G10_short2","g10_rep4")]))==2)
p1 | p2
ggsave("figure_panels/RxT_vs_elong_G0_G10_rep4.pdf", width = 6, height = 3)
