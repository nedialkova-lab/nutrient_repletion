# ==============================================================================
# 9-2_ribonet_seqlets.R
# Characterisation of high-confidence 5'UTR GAUY seqlets identified by RiboNet /
# TF-MoDISco: GO/KEGG enrichment and codon-usage enrichment (GSEA, Wilcoxon) of
# GAUY-containing genes.
#
# Manuscript: Results text (no GO/codon-usage enrichment of GAUY genes); Supplementary Table
# Inputs    : results/seqlets/GATC_seqlets.tsv
#             results/seqlets/allgenes.tsv
#             data/sacCer3.no_dubious_ORFs_coding.fa
# Outputs   : results/GOenrich_GATC.xlsx
#             results/codus_GATC.tsv
#             results/codus_GATC_stats.tsv
# Run from the repository root (all paths are relative to it).
# ==============================================================================

library(RColorBrewer)
library(eulerr)
library(dplyr)
library(clusterProfiler)
library(limma)
library(org.Sc.sgd.db)
library(reshape)
library(openxlsx)
library(ggplot2)
library(UniProt.ws)
library(seqinr)
library(ggpubr)

# Load data
gatc_seqlets = read.csv("results/seqlets/GATC_seqlets.tsv", sep="\t")
gatc_seqlets = gatc_seqlets[gatc_seqlets$maxattr>0.01,]
allgenes = read.csv("results/seqlets/allgenes.tsv", sep="\t", header = F)$V1

# GO ontology
gotab = UniProt.ws::select(org.Sc.sgd.db, allgenes, columns=c("GOALL","ORF"), keytype="ORF")
mapgo = go2term(unique(gotab$GO))
gotab$goterm = as.character(sapply(gotab$GO, function(x) mapgo[mapgo$go_id==x,"Term"]))
gotab_bp = unique(gotab[gotab$ONTOLOGY=="BP",c("goterm","ORF")])
gotab_cc = unique(gotab[gotab$ONTOLOGY=="CC",c("goterm","ORF")])
gotab_mf = unique(gotab[gotab$ONTOLOGY=="MF",c("goterm","ORF")])

keggtab <- getGeneKEGGLinks(species="sce")
mappaths = getKEGGPathwayNames(species="sce")
keggtab$pathname = sapply(keggtab$PathwayID, function(x) substr(mappaths[mappaths$PathwayID==x,"Description"],0, nchar(mappaths[mappaths$PathwayID==x,"Description"])-43))
keggtab = unique(keggtab[,c("pathname","GeneID")])

# Venn diagram
overlap = list(G = gatc_seqlets[gatc_seqlets$glucose!="","gene"],
               N = gatc_seqlets[gatc_seqlets$nitrogen!="","gene"],
               P = gatc_seqlets[gatc_seqlets$phosphate!="","gene"])
fit <- euler(overlap)
plot(fit, fills = c("#f4d7d7","#f4f4d7","#d7e6f4"), edges=T, legend=F, labels=T, quantities=T,main="GATC")

#### GO enrichment ####
genesGATC = gatc_seqlets$gene
res.go1 <- enricher(genesGATC, pvalueCutoff = 0.05, TERM2GENE = gotab_bp, universe = allgenes)
res.go2 <- enricher(genesGATC, pvalueCutoff = 0.05, TERM2GENE = gotab_cc, universe = allgenes)
res.go3 <- enricher(genesGATC, pvalueCutoff = 0.05, TERM2GENE = gotab_mf, universe = allgenes)
write.xlsx(list(BP=res.go1@result, CC=res.go2@result, MF=res.go3@result), "results/GOenrich_GATC.xlsx")

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

# Create GATC label column
codfreq = as.data.frame(codfreq)
codfreq$group <- ifelse(rownames(codfreq) %in% gatc_seqlets$gene, 
                        "GATC", 
                        "Non-GATC")
write.table(codfreq, "results/codus_GATC.tsv", sep = "\t")

# Melt to long format
codfreq_long <- melt(codfreq, 
                     id.vars    = "group",
                     variable.name = "codon",
                     value.name    = "frequency")
testdf = compare_means(value ~ group, codfreq_long, group.by = "variable",
                    method   = "wilcox.test", p.adjust.method = "fdr")
write.table(testdf, "results/codus_GATC_stats.tsv", sep = "\t")

ggplot(codfreq_long, aes(x = variable, 
                         y = value, 
                         fill = group)) +
  geom_boxplot(outlier.size = 0.5,
               outlier.alpha = 0.5) +
  scale_fill_manual(values = c("GATC"     = "#B2182B",
                               "Non-GATC" = "#2166AC")) +
  theme_bw() +   coord_flip() 
