# ==============================================================================
# 3-3_LiPMS_batchcorrection.R
# LiP-MS / trypsin-only (TrP) data extraction with LiPAnalyzeR and batch correction
# (limma::removeBatchEffect, PCR-strip batch) of LiP peptides, TrP peptides and
# TrP proteins.
#
# Manuscript: (pre-processing for Fig. 5 and Extended Data Fig. 6)
# Inputs    : data/20240507_lipms_upshifts_5199/Orbi5199_LiP+Ctrl_HybridLib_Report.tsv
#             data/20240507_lipms_upshifts_5199/240507_lipms.xlsx
# Outputs   : results/{LiPPep,TrPPep,TrPProt}.txt
#             results/proteomics_batchcorrected.tsv
# Run from the repository root (all paths are relative to it).
# ==============================================================================

library(LiPAnalyzeR)
library(openxlsx)
library(UniProt.ws)
library(reshape)
library(RColorBrewer)
library(limma)
library(ggplot2)
library(dplyr)

# Load data
hybriddata = read.csv("data/20240507_lipms_upshifts_5199/Orbi5199_LiP+Ctrl_HybridLib_Report.tsv", sep="\t")
sampletable = read.xlsx("data/20240507_lipms_upshifts_5199/240507_lipms.xlsx", rowNames = F)

# Separate lip and trp
trpdata = hybriddata[hybriddata$R.Replicate %in% 73:144,]
lipdata = hybriddata[hybriddata$R.Replicate %in% 1:72,]

# Annotate sample names
lipdata["Sample"] = apply(lipdata["R.Replicate"],1,function(x) sampletable[sampletable$LiP_MSnumber==x,"Sample.name"])
trpdata["Sample"] = apply(trpdata["R.Replicate"],1,function(x) sampletable[sampletable$TrP_MSnumber==x,"Sample.name"])
lipdata["Condition"] = apply(lipdata["R.Replicate"],1,function(x) sampletable[sampletable$LiP_MSnumber==x,"Condition"])
trpdata["Condition"] = apply(trpdata["R.Replicate"],1,function(x) sampletable[sampletable$TrP_MSnumber==x,"Condition"])

#### LiPAnalyzeR RUV ####

# Extract data
QuantityList <- extractSpectroData(spectroLiP = lipdata, spectroTrP = trpdata, sampleName = "Sample")

# Create annotation files
annotPepProt <- getPepProtAnnot(trpdata, lipdata, "PEP.StrippedSequence", 
                                "PEP.StrippedSequence", "PG.ProteinGroups", "PEP.DigestType....Trypsin.P..LysC.", 
                                "PEP.PeptidePosition", "PEP.AllOccurringProteinAccessions", 
                                "PEP.NrOfMissedCleavages", "PEP.IsProteotypic")
annotSample <- getSampleAnnot(reportOut = lipdata, sampleName = "Sample", sampleCondition = "Condition", contrastCoding=NULL)
annotSample["timepoint"] = as.numeric(sapply(rownames(annotSample),function(x) sampletable[sampletable$Sample.name==x,"Time.point"]))
annotSample["strip"] = as.factor(sapply(rownames(annotSample),function(x) sampletable[sampletable$Sample.name==x,"Strip"]))
annotSample["rep"] = as.factor(sapply(rownames(annotSample),function(x) paste0(sampletable[sampletable$Sample.name==x,c("Condition","Replicate")], collapse="")))

# Fetch gene names
flatprots = unlist(apply(annotPepProt["Protein"],1,function(x) unlist(strsplit(x, ";"))))
mapping = mapUniProt(query=flatprots, from = "UniProtKB_AC-ID", to = "Gene_Name")
annotPepProt["gene_name"] = apply(annotPepProt["Protein"], 1, function(x) paste0(mapping[mapping$From %in% unlist(strsplit(x,";")),"To"], collapse = ";"))

# Filter and log transform
QuantityList_p <- preprocessQuantityMatrix(quantityList = QuantityList,
                                         annotPP = annotPepProt,
                                         annotS = annotSample,
                                         filterTryptic = F,
                                         filterMinLogQuant=T,
                                         thresholdMinLogQuant=5,
                                         filterNA="none")

# Batch correct QuantityList
cond = rep(sprintf("%s%s", annotSample[colnames(QuantityList$LiPPep),"Condition"],annotSample[colnames(QuantityList$LiPPep),"timepoint"]),3)
datatype = c(rep("LiPPep",72),rep("TrPPep",72),rep("TrPProt",72))
design = model.matrix(~ cond + cond:datatype)
batch = c(sprintf("LiPPep_%s",annotSample[colnames(QuantityList$LiPPep),"strip"]),
          sprintf("TrPPep_%s",annotSample[colnames(QuantityList$TrPPep),"strip"]),
          sprintf("TrPProt_%s",annotSample[colnames(QuantityList$TrPProt),"strip"]))
QuantityList_c = removeBatchEffect(cbind(QuantityList_p$LiPPep,QuantityList_p$TrPPep,QuantityList_p$TrPProt), batch=batch, design=design)
LiPPep_c = QuantityList_c[,1:72]
TrPPep_c = QuantityList_c[,73:144]
TrPProt_c = QuantityList_c[,145:216]

# Save
write.table(LiPPep_c,"results/LiPPep.txt", sep="\t", quote = F)
write.table(TrPPep_c,"results/TrPPep.txt", sep="\t", quote = F)
write.table(TrPProt_c,"results/TrPProt.txt", sep="\t", quote = F)

# Corplot between replicates
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}
orderdf = LiPPep_c[,order(colnames(LiPPep_c))]
upper_tri <- melt(get_upper_tri(round(cor(orderdf, method = "spearman", use = "pairwise"),2)), na.rm = TRUE)
ggplot(data = upper_tri, aes(X2, X1, fill = value))+
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "#053061", high = "#67001F", mid = "white", midpoint = 0.95, limit = c(0.9,1), space = "Lab", name="Spearman\nR", na.value = "lightgrey") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 9, hjust = 1),
        axis.title.x=element_blank(), axis.title.y=element_blank())+
  coord_fixed()

# For proteins, save also at protein level
trpprot_dedup = as.data.frame(2^TrPProt_c)
trpprot_dedup$Protein = sapply(rownames(trpprot_dedup), function(x) annotPepProt[x,"Protein"])
trpprot_dedup = trpprot_dedup %>% group_by(Protein) %>% 
  summarise(across(everything(), \(x) mean(x, na.rm = T)), .groups = 'drop')  %>%
  as.data.frame()
write.table(trpprot_dedup, "results/proteomics_batchcorrected.tsv", sep="\t", quote = F, row.names = F)
