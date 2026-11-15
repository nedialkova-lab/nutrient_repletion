# ==============================================================================
# 10-2_attune_TYW1_Ctrls_ONdeprivation.R
# Flow cytometry of the TYW1 5'UTR reporters (replicates 1-3), the synonymous
# fast/slow Citrine controls (RTR1 5'UTR) and HXT6 reporters. Same gating as 10-1.
# Merges the former 10-2 (replicates 1-2 + controls, 20251216) and 10-3
# (replicate 3 + HXT6, 20260116) scripts. Each run is gated and background-
# corrected with its own YSBN6 sample; the bootstrap statistics then pool the
# three TYW1 biological replicates into a single panel.
#
# Manuscript: Extended Data Fig. 7G, 7J
# Inputs    : data/Attune/20251216_TYW1_Ctrls/*.fcs
#             data/Attune/20260116_HXT6_TYW1/*.fcs
# Outputs   : data/Attune/20251216_TYW1_Ctrls/{*.png, TYW1_Ctrls_depr_bg_corrected.tsv}
#             data/Attune/20260116_HXT6_TYW1/{*.png, HXT6_TYW1_depr_bg_corrected.tsv}
#             figure_panels/RTR1_log2ratio.pdf
#             figure_panels/TYW1_{log2ratio,deltalog2ratio}_{r12,r3}.pdf
#             figure_panels/TYW1_bootstrap_{responses,contrasts}.pdf
#             results/flow_TYW1_Ctrls_r12_ratio_percentiles.tsv
#             results/flow_TYW1_r3_ratio_percentiles.tsv
#             results/flow_TYW1_bootstrap_{sample_medians,responses,contrasts,replicate_summary}.tsv
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
library(patchwork)
flowViz.par.set(theme = trellis.par.get(), reset = TRUE)

plots_on <- TRUE
data_dir_r12 <- "data/Attune/20251216_TYW1_Ctrls/"
data_dir_r3 <- "data/Attune/20260116_HXT6_TYW1/"

#Set variable names for flow data files - must correspond to .fcs file parameters
var.names <- c("time","FSC.A","SSC.A","Citrine.A","Citrine2.A", "mCherry.A",
               "FSC.H","SSC.H","Citrine.H","Citrine2.H", "mCherry.H",
               "FSC.W","SSC.W","Citrine.W","Citrine2.W", "mCherry.W")

strainsmap = c(YDN269="HXT6-WT",YDN270="CPR8-duplicated",YDN271="CPR8-scrambled",YDN272="CPR8-WT",
               YDN273="HXT6-duplicated",YDN274="HXT6-scrambled1",YDN275="HXT6-scrambled2",
               YDN276="RTR1-ctrl",YDN277="TYW1-duplicated",YDN278="TYW1-scrambled",YDN279="TYW1-WT",
               YDN280="HXT6-WT", YDN282="RTR1-citmax",YDN283="RTR1-cit0",YDN284="TYW1-duplicated")

# ------------------------------------------------------------------------------
# Import, gate and background-correct one Attune run (all .fcs files in data_dir)
# ------------------------------------------------------------------------------
gate_and_correct <- function(data_dir, out_name) {
  #load the fcs files into a set
  flowData <- read.flowSet(files = list.files(path = data_dir, pattern = "*.fcs"),
                           path = data_dir, transformation = F)
  flowCore::colnames(flowData) <- var.names
  sampleNames(flowData) <- sapply(sapply(list.files(path = data_dir, pattern = ".fcs"),
                                  function(x){paste0(strsplit(x,"_")[[1]][5:length(strsplit(x,"_")[[1]])], collapse="_")}), function(x) substr(x,1,nchar(x)-4))

  # gating: remove debris with low fsc and ssc
  nondebrisfilter <- rectangleGate("FSC.A"=c(2e4, Inf), "SSC.A"=c(1e4, Inf))
  nondebrisresults <- flowCore::filter(flowData, nondebrisfilter)
  nondebris <- flowCore::Subset(flowData, nondebrisresults)
  singletfilter <- norm2Filter("FSC.H", "FSC.A", scale.factor=3)
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
  print(table(gated_data$.id))

  #find WT background
  citbackground <- gated_data %>%
    dplyr::filter(.id == 'YSBN6') %>%
    select(Citrine.H) %>%
    summarise(median = median(Citrine.H))

  mchbackground <- gated_data %>%
    dplyr::filter(.id == 'YSBN6') %>%
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
  write_delim(gated_data, sprintf("%s/%s",data_dir,out_name), delim = "\t")
  gated_data
}

# ------------------------------------------------------------------------------
# Percentiles of Citrine.cor, mCherry.cor and ratio per sample, and fold change vs G0
# ------------------------------------------------------------------------------
compute_percentiles <- function(gated_data) {
  percdf = as.data.frame(sapply(unique(gated_data$.id), function(x) apply(gated_data[gated_data$.id==x,c("Citrine.cor","mCherry.cor","ratio")],2,function(y) quantile(y, probs=seq(0,1,0.0001)))))
  percdf = percdf[,colnames(percdf)!="YSBN6"] # remove YSBN6
  percdf$perc = rep(seq(0,1,0.0001),3)
  percdf$fluor = c(rep("Citrine.cor",10001),rep("mCherry.cor",10001),rep("ratio",10001))
  percdf = melt(percdf,id.vars=c("perc","fluor"))
  percdf$strain = sapply(as.character(percdf$variable), function(x) strsplit(x,"_")[[1]][1])
  percdf$strain_def = sapply(percdf$strain, function(x) strsplit(strainsmap[x],"-")[[1]][2])
  percdf$gene = sapply(percdf$strain, function(x) strsplit(strainsmap[x],"-")[[1]][1])
  percdf$strain_def = factor(as.character(percdf$strain_def), levels=c("scrambled2","scrambled1","scrambled","WT","duplicated","cit0","citmax","ctrl"))
  percdf$rep = sapply(as.character(percdf$variable), function(x) strsplit(x,"_")[[1]][3])
  percdf$time = sapply(as.character(percdf$variable), function(x) substr(strsplit(x,"_")[[1]][2],1,3))
  percdf$time = factor(as.character(percdf$time), levels=c("Glc","G0","G1","G2"))

  # compute fc
  percdf = percdf[order(percdf$strain),]
  orderedstrains = unique(percdf[order(percdf$strain),"strain"])
  percdf["g0"] = c(unlist(sapply(orderedstrains, function(x) rep(percdf[(percdf$time=="G0")&(percdf$strain==x),"value"],4))))  # c(unlist()) handles both matrix (equal reps per strain) and list output
  percdf["delta_log2_ratio"] = log2(percdf$value/percdf$g0)
  percdf
}

# ==============================================================================
# Replicates 1-2 + RTR1 Citrine controls (20251216)
# ==============================================================================
gated_data_r12 <- gate_and_correct(data_dir_r12, "TYW1_Ctrls_depr_bg_corrected.tsv")
percdf <- compute_percentiles(gated_data_r12)
write.table(percdf[percdf$fluor=="ratio",], "results/flow_TYW1_Ctrls_r12_ratio_percentiles.tsv", sep="\t", quote = F, row.names = F)

# Plot
ggplot(percdf[(percdf$fluor=="mCherry.cor"),], aes(x=strain_def, y=value, fill=time)) +
  facet_grid(~ gene*rep, scales="free_x") +
  geom_point_rast(position=position_jitterdodge(jitter.width=0.2), alpha=0.2, size=0.5, color="grey") +
  geom_boxplot(outlier.shape = NA, alpha=0.5) +
  coord_cartesian(ylim = c(0, 8000)) +
  theme_bw()
ggplot(percdf[(percdf$fluor=="Citrine.cor"),], aes(x=strain_def, y=value, fill=time)) +
  facet_grid(~ gene*rep, scales="free_x") +
  geom_point_rast(position=position_jitterdodge(jitter.width=0.2), alpha=0.2, size=0.5, color="grey") +
  geom_boxplot(outlier.shape = NA, alpha=0.5) +
  coord_cartesian(ylim = c(0, 16000)) +
  theme_bw()

# RTR1
ggplot(percdf[(percdf$fluor=="ratio")&(percdf$gene=="RTR1"),], aes(x=strain_def, y=log2(value), fill=time)) +
  facet_wrap(~ rep, scales="free_x", ncol=1) +
  geom_point_rast(position=position_jitterdodge(jitter.width=0.2), alpha=0.2, size=0.5, color="grey") +
  geom_boxplot(outlier.shape = NA, alpha=0.5) +
  coord_cartesian(ylim = c(-7, 3.0)) +
  theme_bw()
ggsave("figure_panels/RTR1_log2ratio.pdf", width = 6, height = 6)
# TYW1
ggplot(percdf[(percdf$fluor=="ratio")&(percdf$gene=="TYW1"),], aes(x=strain_def, y=log2(value), fill=time)) +
  facet_wrap(~ rep, scales="free_x", ncol=1) +
  geom_point_rast(position=position_jitterdodge(jitter.width=0.2), alpha=0.2, size=0.5, color="grey") +
  geom_boxplot(outlier.shape = NA, alpha=0.5) +
  coord_cartesian(ylim = c(-1.0, 3.0)) +
  theme_bw()
ggsave("figure_panels/TYW1_log2ratio_r12.pdf", width = 6, height = 6)
ggplot(percdf[(percdf$fluor=="ratio")&(percdf$gene=="TYW1"),], aes(x=time, y=delta_log2_ratio, fill=strain_def)) +
  facet_wrap(~ rep, scales="free_x", ncol=1) +
  geom_point_rast(position=position_jitterdodge(jitter.width=0.2), alpha=0.2, size=0.5, color="grey") +
  geom_boxplot(outlier.shape = NA, alpha=0.5) +
  coord_cartesian(ylim = c(-0.5, 2.0)) +
  theme_bw()
ggsave("figure_panels/TYW1_deltalog2ratio_r12.pdf", width = 6, height = 6)

# ==============================================================================
# Replicate 3 + HXT6 (20260116)
# ==============================================================================
gated_data_r3 <- gate_and_correct(data_dir_r3, "HXT6_TYW1_depr_bg_corrected.tsv")
percdf <- compute_percentiles(gated_data_r3)
write.table(percdf[percdf$fluor=="ratio",], "results/flow_TYW1_r3_ratio_percentiles.tsv", sep="\t", quote = F, row.names = F)

# Plot
ggplot(percdf[(percdf$fluor=="mCherry.cor"),], aes(x=strain_def, y=value, fill=time)) +
  facet_grid(~ gene*rep, scales="free_x") +
  geom_point_rast(position=position_jitterdodge(jitter.width=0.2), alpha=0.2, size=0.5, color="grey") +
  geom_boxplot(outlier.shape = NA, alpha=0.5) +
  coord_cartesian(ylim = c(0, 8000)) +
  theme_bw()
ggplot(percdf[(percdf$fluor=="Citrine.cor"),], aes(x=strain_def, y=value, fill=time)) +
  facet_grid(~ gene*rep, scales="free_x") +
  geom_point_rast(position=position_jitterdodge(jitter.width=0.2), alpha=0.2, size=0.5, color="grey") +
  geom_boxplot(outlier.shape = NA, alpha=0.5) +
  coord_cartesian(ylim = c(0, 16000)) +
  theme_bw()

# TYW1
ggplot(percdf[(percdf$fluor=="ratio")&(percdf$gene=="TYW1"),], aes(x=strain_def, y=log2(value), fill=time)) +
  facet_wrap(~ rep, scales="free_x", ncol=1) +
  geom_point_rast(position=position_jitterdodge(jitter.width=0.2), alpha=0.2, size=0.5, color="grey") +
  geom_boxplot(outlier.shape = NA, alpha=0.5) +
  coord_cartesian(ylim = c(-1.0, 3.0)) +
  theme_bw()
ggsave("figure_panels/TYW1_log2ratio_r3.pdf", width = 6, height = 3.2)

# TYW1
ggplot(percdf[(percdf$fluor=="ratio")&(percdf$gene=="TYW1"),], aes(x=time, y=delta_log2_ratio, fill=strain_def)) +
  facet_wrap(~ rep, scales="free_x", ncol=1) +
  geom_point_rast(position=position_jitterdodge(jitter.width=0.2), alpha=0.2, size=0.5, color="grey") +
  geom_boxplot(outlier.shape = NA, alpha=0.5) +
  coord_cartesian(ylim = c(-0.5, 1.5)) +
  theme_bw()
ggsave("figure_panels/TYW1_deltalog2ratio_r3.pdf", width = 6, height = 3.2)

# ==============================================================================
# Bootstrap reanalysis (replaces the former percentile-level Wilcoxon tests)
# Per-culture medians of the single-cell log2(ratio), ratio = Citrine.cor/mCherry.cor
# (background-corrected per run, as saved in the *_bg_corrected.tsv tables), with
# within-culture cell-bootstrap 95% CIs. No quantile binning, no cell/percentile-
# level p-values. CIs are cell-sampling uncertainty within cultures, NOT between
# biological replicates. Replicates 1-2 and 3 come from different runs and are
# pooled here as three biological replicates in one panel.
# ==============================================================================
set.seed(20261115)
B <- 1000L  # bootstrap draws per sample

bs_gene <- "TYW1"
bs_strains <- c(YDN284 = "duplicated", YDN278 = "scrambled", YDN279 = "WT")
bs_variants <- c("scrambled", "WT", "duplicated")
bs_conditions <- c("Glc", "G0", "G1", "G2")  # same x axis as the original figure

gated_data <- dplyr::bind_rows(gated_data_r12, gated_data_r3)

# Parse sample IDs (e.g., YDN284_G1_rep3); use corrected single-cell ratio.
cells <- gated_data %>%
  dplyr::filter(grepl("^YDN(278|279|284)_(Glc|G0|G1|G2)_rep[0-9]+$", .id),
                is.finite(ratio), ratio > 0) %>%
  dplyr::transmute(
    sample_id = .id,
    strain = sub("_.*", "", .id),
    condition = sub("^[^_]+_([^_]+)_.*$", "\\1", .id),
    replicate = sub("^.*_(rep[0-9]+)$", "\\1", .id),
    variant = unname(bs_strains[strain]),
    log_ratio = log2(ratio)
  )
if (nrow(cells) == 0) stop("No TYW1 cell data found. Check sample IDs and gated_data$ratio.")
# replicate labels must be unique across the two runs (rep1-2 vs rep3)
stopifnot(length(intersect(cells$sample_id[cells$sample_id %in% gated_data_r12$.id],
                           cells$sample_id[cells$sample_id %in% gated_data_r3$.id])) == 0)

# Every replicate needs all conditions for each of the 3 variants.
coverage <- cells %>% dplyr::distinct(sample_id, replicate, variant, condition)
needed <- expand.grid(replicate = sort(unique(cells$replicate)),
                      variant = bs_variants,
                      condition = bs_conditions,
                      stringsAsFactors = FALSE)
missing <- dplyr::anti_join(needed, coverage, by = c("replicate", "variant", "condition"))
if (nrow(missing)) {
  print(missing)
  stop("Missing matched samples: do not calculate contrasts until these are resolved.")
}

# Resampling is performed WITHIN each independent cell sample, not across cultures.
# The b-th draw of a given G0 sample is reused in all contrasts containing that baseline.
# This preserves shared-baseline dependence in the bootstrap contrast distributions.
bootstrap_median <- function(x, B) {
  n <- length(x)
  if (n < 2L) stop("Too few cells in a sample")
  vapply(seq_len(B), function(b) median(x[sample.int(n, n, replace = TRUE)]), numeric(1))
}

samples <- cells %>%
  dplyr::group_by(replicate, variant, condition) %>%
  dplyr::summarise(n_cells = dplyr::n(),
                   estimate = median(log_ratio),
                   draws = list(bootstrap_median(log_ratio, B)),
                   .groups = "drop")
stopifnot(nrow(samples) == nrow(needed))

# Give each sample an unambiguous key; list-column contains B bootstrap draws.
keys <- paste(samples$replicate, samples$variant, samples$condition, sep = "|")
point <- setNames(samples$estimate, keys)
boot <- setNames(samples$draws, keys)
key <- function(rep, variant, condition) paste(rep, variant, condition, sep = "|")
ci <- function(draws) unname(quantile(draws, c(.025, .975), names = FALSE))

# 1. Baseline (G0, depleted) normalized response for each variant within each replicate.
# 2. Difference-in-differences of WT / duplicated relative to scrambled.
# G0 is the embedded reference (0 by construction); its x slot is kept in the plots
# to match the original figure, but its points are not drawn.
response_rows <- list()
contrast_rows <- list()
for (rep in sort(unique(samples$replicate))) {
  for (time in bs_conditions) {
    response_est <- list()
    response_boot <- list()
    for (v in bs_variants) {
      kt <- key(rep, v, time)
      k0 <- key(rep, v, "G0")
      response_est[[v]] <- unname(point[[kt]] - point[[k0]])
      response_boot[[v]] <- boot[[kt]] - boot[[k0]]
      bounds <- ci(response_boot[[v]])
      response_rows[[length(response_rows) + 1L]] <- data.frame(
        gene = bs_gene, replicate = rep, condition = time, variant = v,
        response_log2 = response_est[[v]],
        lower95 = bounds[1], upper95 = bounds[2],
        n_cells_time = samples$n_cells[match(kt, keys)],
        n_cells_G0 = samples$n_cells[match(k0, keys)]
      )
    }
    for (v in c("WT", "duplicated")) {
      difference <- response_est[[v]] - response_est[["scrambled"]]
      # Same bootstrap draws for scrambled baseline are reused across both comparisons.
      draws <- response_boot[[v]] - response_boot[["scrambled"]]
      bounds <- ci(draws)
      contrast_rows[[length(contrast_rows) + 1L]] <- data.frame(
        gene = bs_gene, replicate = rep, condition = time,
        contrast = paste(v, "vs scrambled"),
        effect_log2 = difference,
        lower95 = bounds[1], upper95 = bounds[2],
        bootstrap_ci_excludes_zero = (bounds[1] > 0 || bounds[2] < 0)
      )
    }
  }
}

responses <- dplyr::bind_rows(response_rows) %>%
  dplyr::mutate(variant = factor(variant, levels = bs_variants),
                condition = factor(condition, levels = bs_conditions))
contrasts <- dplyr::bind_rows(contrast_rows) %>%
  dplyr::mutate(contrast = factor(contrast, levels = c("WT vs scrambled", "duplicated vs scrambled")),
                condition = factor(condition, levels = bs_conditions))

# Descriptive summary across biological replicates: NO between-replicate p-value.
effect_summary <- contrasts %>%
  dplyr::filter(condition != "G0") %>%
  dplyr::group_by(gene, condition, contrast) %>%
  dplyr::summarise(n_biological_replicates = dplyr::n(),
                   mean_effect_log2 = mean(effect_log2),
                   min_effect_log2 = min(effect_log2),
                   max_effect_log2 = max(effect_log2),
                   positive_replicates = sum(effect_log2 > 0), .groups = "drop")

write.table(samples %>% dplyr::select(-draws), sprintf("results/flow_%s_bootstrap_sample_medians.tsv", bs_gene), sep = "\t", quote = F, row.names = F)
write.table(responses, sprintf("results/flow_%s_bootstrap_responses.tsv", bs_gene), sep = "\t", quote = F, row.names = F)
write.table(contrasts, sprintf("results/flow_%s_bootstrap_contrasts.tsv", bs_gene), sep = "\t", quote = F, row.names = F)
write.table(effect_summary, sprintf("results/flow_%s_bootstrap_replicate_summary.tsv", bs_gene), sep = "\t", quote = F, row.names = F)
print(responses)
print(contrasts)
print(effect_summary)

# Single-panel plots: same x axis (Glc, G0, G1, G2) and construct colors as the
# original boxplots (scrambled = muted green, WT = ochre, duplicated = lavender).
# Points are biological replicates (shape), whiskers are cell-bootstrap 95% CIs.
variant_colors <- c(scrambled = "#78A99D", WT = "#DBB579", duplicated = "#B99AC8")
bs_reps <- sort(unique(responses$replicate))
replicate_shapes <- setNames(c(16, 17, 15, 18, 8)[seq_along(bs_reps)], bs_reps)
rep_offset <- if (length(bs_reps) > 1) setNames(seq(-0.05, 0.05, length.out = length(bs_reps)), bs_reps) else setNames(0, bs_reps)
# manual nested dodge: condition -> variant -> replicate
add_x <- function(df, var_levels, var_halfwidth) {
  var_offset <- setNames(seq(-var_halfwidth, var_halfwidth, length.out = length(var_levels)), var_levels)
  df$x <- as.numeric(df$condition) + var_offset[as.character(df$variant)] + rep_offset[as.character(df$replicate)]
  df
}
bs_theme <- theme_classic(base_size = 12) +
  theme(panel.grid.major.y = element_line(colour = "grey92", linewidth = .3),
        axis.text.x = element_text(size = 11),
        axis.title.y = element_text(margin = margin(r = 10)),
        legend.position = "bottom", legend.box = "vertical",
        legend.margin = margin(0, 0, 0, 0), legend.spacing.y = grid::unit(2, "pt"),
        plot.caption = element_text(size = 8.5, colour = "grey35", hjust = 0, margin = margin(t = 9)),
        plot.margin = margin(10, 12, 7, 10))
bs_caption <- sprintf("Points: individual biological replicates (n = %d).\nError bars: 95%% cell-bootstrap CIs within each replicate (B = %d).",
                      length(bs_reps), B)

# Plot 1: response of each construct vs depleted (G0), analogous to the original boxplots.
p_response <- ggplot(add_x(dplyr::filter(responses, condition != "G0"), bs_variants, 0.28),
                     aes(x = x, y = response_log2, colour = variant, shape = replicate)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey65", linewidth = .35) +
  geom_errorbar(aes(ymin = lower95, ymax = upper95), width = .05, linewidth = .6) +
  geom_point(size = 2.6) +
  scale_x_continuous(breaks = seq_along(bs_conditions), labels = bs_conditions,
                     expand = expansion(add = .45)) +
  scale_colour_manual(values = variant_colors) +
  scale_shape_manual(values = replicate_shapes) +
  labs(x = NULL, colour = "Construct", shape = "Biological replicate",
       y = expression(Delta~"median "*log[2]*"(Citrine/mCherry) vs depleted"),
       title = bs_gene, caption = bs_caption) +
  bs_theme
print(p_response)
ggsave(sprintf("figure_panels/%s_bootstrap_responses.pdf", bs_gene), p_response,
       width = 5.5, height = 4.6, device = "pdf")

# Plot 2: difference-in-differences relative to scrambled
# (dashed zero line = scrambled).
contrasts_plot <- contrasts %>%
  dplyr::mutate(variant = factor(ifelse(as.character(contrast) == "WT vs scrambled", "WT", "duplicated"),
                                 levels = c("WT", "duplicated")))
p_effect <- ggplot(add_x(dplyr::filter(contrasts_plot, condition != "G0"), c("WT", "duplicated"), 0.17),
                   aes(x = x, y = effect_log2, colour = variant, shape = replicate)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = variant_colors[["scrambled"]], linewidth = .5) +
  geom_errorbar(aes(ymin = lower95, ymax = upper95), width = .05, linewidth = .6) +
  geom_point(size = 2.6) +
  scale_x_continuous(breaks = seq_along(bs_conditions), labels = bs_conditions,
                     expand = expansion(add = .45)) +
  scale_y_continuous(expand = expansion(mult = c(.12, .12))) +
  scale_colour_manual(values = variant_colors[c("WT", "duplicated")]) +
  scale_shape_manual(values = replicate_shapes) +
  labs(x = NULL, colour = "Construct (vs scrambled)", shape = "Biological replicate",
       y = expression(Delta*Delta~"median "*log[2]*"(Citrine/mCherry)"),
       title = bs_gene, caption = bs_caption) +
  bs_theme
print(p_effect)
ggsave(sprintf("figure_panels/%s_bootstrap_contrasts.pdf", bs_gene), p_effect,
       width = 5.5, height = 4.6, device = "pdf")
