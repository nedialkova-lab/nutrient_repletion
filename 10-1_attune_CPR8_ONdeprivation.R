# ==============================================================================
# 10-1_attune_CPR8_ONdeprivation.R
# Flow cytometry (Attune NxT) of the CPR8 5'UTR reporter strains: gating (debris,
# singlets), background subtraction and Citrine/mCherry ratios relative to the
# glucose-depleted condition. QC plots and the gated table are written into the
# data folder itself.
#
# Manuscript: Fig. 6G; Extended Data Fig. 7H
# Inputs    : data/Attune/20251009_CPR8_ONdeprivation/*.fcs
# Outputs   : data/Attune/20251009_CPR8_ONdeprivation/{*.png, CPR8depr_bg_corrected.tsv}
#             figure_panels/CPR8_log2ratio.pdf
#             figure_panels/CPR8_deltalog2ratio.pdf
#             results/flow_CPR8_ratio_percentiles.tsv
# Run from the repository root (all paths are relative to it).
# ==============================================================================

library(flowCore)
library(flowStats)
library(flowViz)
library(plyr)
library(tidyverse)
library(ggplot2)
library(ggrastr)
library(reshape)
library(ggpubr)
library(patchwork)
flowViz.par.set(theme = trellis.par.get(), reset = TRUE)

plots_on <- TRUE
data_dir <- "data/Attune/20251009_CPR8_ONdeprivation/"

#Import and gate flow data - run on all .fcs files in the working directory folder 
#Set variable names for flow data files - must correspond to .fcs file parameters
var.names <- c("time","FSC.A","SSC.A","Citrine.A","Citrine2.A", "mCherry.A",
               "FSC.H","SSC.H","Citrine.H","Citrine2.H", "mCherry.H",
               "FSC.W","SSC.W","Citrine.W","Citrine2.W", "mCherry.W")


#load the fcs files into a set
flowData <- read.flowSet(files = list.files(path = data_dir, pattern = "*.fcs"),
                         path = data_dir, transformation = F)
flowCore::colnames(flowData) <- var.names
sampleNames(flowData) <- sapply(sapply(list.files(path = data_dir, pattern = ".fcs"),
                                function(x){paste0(strsplit(x,"_")[[1]][5:length(strsplit(x,"_")[[1]])], collapse="_")}), function(x) substr(x,1,nchar(x)-4))

# gating: remove debris with low fsc and ssc
nondebrisfilter <- rectangleGate("FSC.A"=c(2e4, Inf), "SSC.A"=c(1e4, Inf))
#sqrcut = matrix(c(0,1e6,1e6,0,8e5,0),ncol=2,nrow=3)
#budfilter <- polygonGate(sqrcut)
nondebrisresults <- flowCore::filter(flowData, nondebrisfilter)
nondebris <- flowCore::Subset(flowData, nondebrisresults)
#normfilter <- norm2Filter("FSC.A", "SSC.A", scale.factor=1.5)
singletfilter <- norm2Filter("FSC.H", "FSC.A", scale.factor=3)
#gateresults <- flowCore::filter(nonbudcells, normfilter & singletfilter)
gateresults <- flowCore::filter(nondebris, singletfilter)
flowDataGated <- flowCore::Subset(nondebris, gateresults)

# diagnostic plots of flow data
if(plots_on == TRUE) {
  #plot data with gate
  png(sprintf("%s/%s",data_dir,"fsc-ssc-nondebrisgate.png"), width = 20, height = 20, units = "in", res = 300 )
  print(xyplot(`SSC.A` ~ `FSC.A`, data = flowData, smooth = FALSE, filter = nondebrisfilter, ylim = c(0,5e5), xlim = c(0,1e6)))
  dev.off()
  # plot gate data only
  png(sprintf("%s/%s",data_dir,"fsc-ssc-gated-only.png"), width = 20, height = 20, units = "in", res = 300 )
  print(xyplot(`SSC.A` ~ `FSC.A`, data = nondebris, smooth = FALSE, ylim = c(0,5e5), xlim = c(0,1e6)))
  dev.off()
  #plot singlets data
  png(sprintf("%s/%s",data_dir,"fscH-fscA-singletgate.png"), width = 20, height = 20, units = "in", res = 300 )
  print(xyplot(`FSC.H` ~ `FSC.A`, data = nondebris, smooth = FALSE, filter = singletfilter, ylim = c(0,8e5), xlim = c(0,1e6)))
  dev.off()
  # plot singlets data gate only
  png(sprintf("%s/%s",data_dir,"fscH-fscA-gated-only.png"), width = 20, height = 20, units = "in", res = 300 )
  print(xyplot(`FSC.H` ~ `FSC.A`, data = flowDataGated, smooth = FALSE, ylim = c(0,8e5), xlim = c(0,1e6)))
  dev.off()
  # plot gated fluorescence data only
  png(sprintf("%s/%s",data_dir,"fluor-gated-only.png"), width = 20, height = 20, units = "in", res = 300 )
  print(xyplot(`Citrine.H` ~ `mCherry.H`, data = flowDataGated, smooth = FALSE, ylim = c(0,30000), xlim = c(0,12000)))
  dev.off()
  png(sprintf("%s/%s",data_dir,"fluor-all.png"), width = 20, height = 20, units = "in", res = 300 )
  print(xyplot(`Citrine.H` ~ `mCherry.H`, data = flowData, smooth = FALSE, ylim = c(0,30000), xlim = c(0,12000)))
  dev.off()
}

# convert the flowSet data to matrices then dataframes and combine them into one big dataframe
gated_data <- fsApply(flowDataGated, exprs, simplify = F, use.exprs = F)
gated_data <- lapply(gated_data, data.frame)
gated_data <- ldply(gated_data)
table(gated_data$.id)

#find WT background
citbackground <- gated_data %>%
  dplyr::filter(.id == 'YSBN6_Glc') %>%
  select(Citrine.H) %>%
  summarise(median = median(Citrine.H))

mchbackground <- gated_data %>%
  dplyr::filter(.id == 'YSBN6_Glc') %>%
  select(mCherry.H) %>%
  summarise(median = median(mCherry.H))

#background correction
gated_data$Citrine.cor <- gated_data$Citrine.H - citbackground$median
gated_data$mCherry.cor <- gated_data$mCherry.H - mchbackground$median

#if fluorescence is <1, set it to 1
gated_data[gated_data$Citrine.cor<1,"Citrine.cor"] = 1
gated_data[gated_data$mCherry.cor<1,"mCherry.cor"] = 1

# ratio of citrine to mCherry for each data point
gated_data$ratio <- with(gated_data, Citrine.cor / mCherry.cor)
gated_data$ratio.H <- with(gated_data, Citrine.H / mCherry.H)
write_delim(gated_data, sprintf("%s/%s",data_dir,"CPR8depr_bg_corrected.tsv"), delim = "\t")

# Compute percentiles
strainsmap = c(YDN269="HXT6-WT",YDN270="CPR8-duplicated",YDN271="CPR8-scrambled",YDN272="CPR8-WT",
               YDN273="HXT6-duplicated",YDN274="HXT6-scrambled1",YDN275="HXT6-scrambled2",
               YDN276="RTR1-ctrl",YDN277="TYW1-duplicated",YDN278="TYW1-scrambled",YDN279="TYW1-WT")
percdf = as.data.frame(sapply(unique(gated_data$.id), function(x) apply(gated_data[gated_data$.id==x,c("Citrine.cor","mCherry.cor","ratio")],2,function(y) quantile(y, probs=seq(0,1,0.0001)))))
percdf = percdf[,!(colnames(percdf) %in% c("YSBN6_Glc","YSBN6_G0"))] # remove YSBN6
percdf$perc = rep(seq(0,1,0.0001),3)
percdf$fluor = c(rep("Citrine.cor",10001),rep("mCherry.cor",10001),rep("ratio",10001))
percdf = melt(percdf,id.vars=c("perc","fluor"))
percdf$strain = sapply(as.character(percdf$variable), function(x) strsplit(x,"_")[[1]][1])
percdf$strain_def = sapply(percdf$strain, function(x) strsplit(strainsmap[x],"-")[[1]][2])
percdf$strain_def = factor(as.character(percdf$strain_def), levels=c("scrambled2","scrambled1","scrambled","WT","duplicated","ctrl"))
percdf$rep = sapply(as.character(percdf$variable), function(x) strsplit(x,"_")[[1]][3])
percdf$time = sapply(as.character(percdf$variable), function(x) substr(strsplit(x,"_")[[1]][2],1,3))
percdf$time = factor(as.character(percdf$time), levels=c("Glc","G0","G1","G2"))

# compute fc
percdf = percdf[order(percdf$strain),]
orderedstrains = unique(percdf[order(percdf$strain),"strain"])
percdf["g0"] = c(sapply(orderedstrains, function(x) rep(percdf[(percdf$time=="G0")&(percdf$strain==x),"value"],4)))
percdf["delta_log2_ratio"] = log2(percdf$value/percdf$g0)
write.table(percdf[percdf$fluor=="ratio",], "results/flow_CPR8_ratio_percentiles.tsv", sep="\t", quote = F, row.names = F)

# Plot
ggplot(percdf[(percdf$fluor=="mCherry.cor"),], aes(x=strain_def, y=value, fill=time)) +
  facet_grid(~ rep, scales="free_x") +
  geom_point_rast(position=position_jitterdodge(jitter.width=0.2), alpha=0.2, size=0.5, color="grey") +
  geom_boxplot(outlier.shape = NA, alpha=0.5) +
  coord_cartesian(ylim = c(0, 8000)) +
  theme_bw()
ggplot(percdf[(percdf$fluor=="Citrine.cor"),], aes(x=strain_def, y=value, fill=time)) +
  facet_grid(~ rep, scales="free_x") +
  geom_point_rast(position=position_jitterdodge(jitter.width=0.2), alpha=0.2, size=0.5, color="grey") +
  geom_boxplot(outlier.shape = NA, alpha=0.5) +
  coord_cartesian(ylim = c(0, 8000)) +
  theme_bw()
ggplot(percdf[(percdf$fluor=="ratio"),], aes(x=strain_def, y=log2(value), fill=time)) +
  facet_wrap(~ rep, scales="free_x", ncol=1) +
  geom_point_rast(position=position_jitterdodge(jitter.width=0.2), alpha=0.2, size=0.5, color="grey") +
  geom_boxplot(outlier.shape = NA, alpha=0.5) +
  coord_cartesian(ylim = c(-3.5, 1.5)) +
  theme_bw()
ggsave("figure_panels/CPR8_log2ratio.pdf", width = 6, height = 6)
ggplot(percdf[(percdf$fluor=="ratio"),], aes(x=time, y=delta_log2_ratio, fill=strain_def)) +
  facet_wrap(~ rep, scales="free_x", ncol=1) +
  geom_point_rast(position=position_jitterdodge(jitter.width=0.2), alpha=0.2, size=0.5, color="grey") +
  geom_boxplot(outlier.shape = NA, alpha=0.5) +
  coord_cartesian(ylim = c(-0.8, 2.5)) +
  theme_bw()
ggsave("figure_panels/CPR8_deltalog2ratio.pdf", width = 6, height = 6)
my_comparisons <- list(c("duplicated", "WT"),  c("duplicated","scrambled"), c("WT", "scrambled"), c("WT", "scrambled1"), c("duplicated","scrambled1"),
                       c("WT", "scrambled2"), c("duplicated","scrambled2"), c("scrambled1","scrambled2"))
compare_means(delta_log2_ratio ~ strain_def, percdf[(percdf$fluor=="ratio")&(percdf$time!="G0"),], method = "wilcox.test", alternative = "greater", comparisons = my_comparisons, group.by = c("time","rep"), paired = T)
