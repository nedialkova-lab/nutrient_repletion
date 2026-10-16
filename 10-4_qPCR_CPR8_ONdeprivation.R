# ==============================================================================
# 10-4_qPCR_CPR8_ONdeprivation.R
# RT-qPCR of Citrine relative to mCherry (ddCt) for the CPR8 5'UTR reporter strains.
#
# Manuscript: Extended Data Fig. 7K
# Inputs    : data/qPCR/20251020/20251020_CPR8_ONdeprivation.txt
# Outputs   : figure_panels/CPR8_log2fc_qPCR.pdf
#             results/qPCR_CPR8_dCt.tsv
# Run from the repository root (all paths are relative to it).
# ==============================================================================

library(ggplot2)
library(reshape)
library(ggpubr)
library(patchwork)

# Open data
data_file <- read.csv("data/qPCR/20251020/20251020_CPR8_ONdeprivation.txt", sep = "\t", skip = 1)
data_file = data_file[!(grepl("Sample",data_file$Name)|grepl("NRT",data_file$Name)|grepl("H2O",data_file$Name)),]
data_file$strain = sapply(data_file$Name, function(x) sprintf("YDN%s",strsplit(x,"_")[[1]][1]))
strainsmap = c(YDN269="HXT6-WT",YDN270="CPR8-duplicated",YDN271="CPR8-scrambled",YDN272="CPR8-WT",
               YDN273="HXT6-duplicated",YDN274="HXT6-scrambled1",YDN275="HXT6-scrambled2",
               YDN276="RTR1-ctrl",YDN277="TYW1-duplicated",YDN278="TYW1-scrambled",YDN279="TYW1-WT")
data_file$strain_def = sapply(data_file$strain, function(x) strsplit(strainsmap[x],"-")[[1]][2])
data_file$strain_def = factor(as.character(data_file$strain_def), levels=c("scrambled2","scrambled1","scrambled","WT","duplicated","ctrl"))
data_file$cond = sapply(data_file$Name, function(x) strsplit(x,"_")[[1]][2])
data_file$rep = sapply(data_file$cond, function(x) substr(x,nchar(x), nchar(x)))
data_file$cond = sapply(data_file$cond, function(x) strsplit(substr(x,1, 3),"-")[[1]][1])
data_file$cond = factor(as.character(data_file$cond), levels=c("Glc","G0","G1","G2"))
data_file$gene = sapply(data_file$Name, function(x) strsplit(x,"_")[[1]][3])

# Summary stats
meandf = cast(data_file, strain + strain_def + cond + rep ~ gene, value = "Cp", fun.aggregate = "mean", na.rm=T)
sddf = cast(data_file, strain + strain_def + cond + rep ~ gene, value = "Cp", fun.aggregate = "sd", na.rm=T)

# Compute dCt
dCtdf = cast(data_file, strain + strain_def + cond + rep ~ gene, value = "Cp", fun.aggregate = "mean", na.rm=T)
shift <- log2(20/5)  # = 2, take dilution of samples into account
dCtdf$dCt = dCtdf$Cit - dCtdf$mCh - shift
dCtdf$log2fc = -dCtdf$dCt

# ΔΔCt: normalize to time 0 within each strain
dCtdf$ddCt = apply(as.data.frame(dCtdf),1,function(x) as.numeric(x["dCt"]) - as.numeric(dCtdf[(dCtdf$cond=="G0")&(dCtdf$strain==x["strain"])&(dCtdf$rep==x["rep"]),"dCt"]))
dCtdf$log2fc_vs_g0 = -dCtdf$ddCt
write.table(dCtdf, "results/qPCR_CPR8_dCt.tsv", sep="\t", quote = F, row.names = F)

# Compare harvesting
ggplot(dCtdf, aes(x=cond, y=log2fc, color=strain_def)) +
  geom_point(position=position_jitterdodge(dodge.width=0.65,jitter.width=0)) +
  stat_summary(fun=mean, geom="point", size=7, shape="_", position=position_dodge(0.65)) +
  theme_bw()
ggsave("figure_panels/CPR8_log2fc_qPCR.pdf", width = 6, height = 3.2)
ggplot(dCtdf, aes(x=cond, y=log2fc_vs_g0, color=strain_def)) +
  geom_point(position=position_jitterdodge(dodge.width=0.65,jitter.width=0)) +
  stat_summary(fun=mean, geom="point", size=7, shape="_", position=position_dodge(0.65)) +
  theme_bw()
