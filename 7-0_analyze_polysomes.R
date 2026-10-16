# ==============================================================================
# 7-0_analyze_polysomes.R
# Polysome profiles (A260 along 10-45% sucrose gradients) and percentage of
# polysomes (area under the curve) along the repletions.
#
# Manuscript: Fig. 1E-F; Extended Data Fig. 2D
# Inputs    : data/polysomes/*_polysomes.csv
# Outputs   : figure_panels/polyprof_{Glc,Nit,Pho}.pdf
#             figure_panels/polyperc_{Glc,Nit,Pho}.pdf
# Run from the repository root (all paths are relative to it).
# ==============================================================================

library(ggplot2)
library(DescTools)
library(RColorBrewer)

# Load data GLUCOSE
g0 = read.csv("data/polysomes/240313_R3_G0_polysomes.csv", skip = 42); g0$cond = "g0"; g0$rna = 112
g1 = read.csv("data/polysomes/240313_R3_G1_polysomes.csv", skip = 42); g1$cond = "g1"; g1$rna = 64
g2 = read.csv("data/polysomes/240313_R3_G2_polysomes.csv", skip = 42); g2$cond = "g2"; g2$rna = 71
g4 = read.csv("data/polysomes/240313_R3_G4_polysomes.csv", skip = 42); g4$cond = "g4"; g4$rna = 113
g10 = read.csv("data/polysomes/240313_R3_G10_polysomes.csv", skip = 42); g10$cond = "g10"; g10$rna = 87
g20 = read.csv("data/polysomes/240313_R3_G20_polysomes.csv", skip = 42); g20$cond = "g20"; g20$rna = 118
mergedf = rbind(g0,g1,g2,g4,g10,g20)
mergedf$normabs = mergedf$Absorbance*100/mergedf$rna
mergedf$cond = factor(as.character(mergedf$cond), levels=c("g0","g1","g2","g4","g10","g20"))

p = ggplot(mergedf, aes(Position.mm., normabs, color = cond)) +
  geom_rect(data=data.frame(xmin=22, xmax=31, ymin=-Inf, ymax=Inf),
            aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
            alpha=0.4, fill="lightgrey", inherit.aes=FALSE) +
  geom_rect(data=data.frame(xmin=38, xmax=70, ymin=-Inf, ymax=Inf),
            aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
            alpha=0.4, fill="lightgrey", inherit.aes=FALSE) +
  geom_line() +
  xlim(10, 80) + ylim(0, 3.6) +
  scale_color_manual(values=c("#B2182B","#D6604D","#F4A582","#92C5DE","#4393C3","#2166AC")) +
  theme_bw()
ggsave("figure_panels/polyprof_Glc.pdf", p, width = 4, height = 2.2)

ratiodf = data.frame(cond=factor(c("g0","g1","g2","g4","g10","g20"), levels = c("g0","g1","g2","g4","g10","g20")))
ratiodf$poly = sapply(ratiodf$cond, function(x) AUC(mergedf[mergedf$cond==x,"Position.mm."],mergedf[mergedf$cond==x,"normabs"], from=38, to=70))
ratiodf$mono = sapply(ratiodf$cond, function(x) AUC(mergedf[mergedf$cond==x,"Position.mm."],mergedf[mergedf$cond==x,"normabs"], from=22, to=31))
ratiodf$poly2mono = ratiodf$poly/ratiodf$mono
ratiodf$mono2poly = ratiodf$mono/ratiodf$poly
ratiodf$total = sapply(ratiodf$cond, function(x) AUC(mergedf[mergedf$cond==x,"Position.mm."],mergedf[mergedf$cond==x,"normabs"], from=22, to=70))
ratiodf$polyperc = ratiodf$poly*100/ratiodf$total
p = ggplot(ratiodf, aes(cond,polyperc,fill=cond)) + 
  geom_bar(stat="identity") +
  ylim(0,80) +
  scale_fill_manual(values=c("#B2182B","#D6604D","#F4A582","#92C5DE","#4393C3","#2166AC")) +
  theme_bw()
ggsave("figure_panels/polyperc_Glc.pdf", p, width = 3.2, height = 2.2)


# Load data NITROGEN
n0 = read.csv("data/polysomes/240306_R3_N0_polysomes.csv", skip = 42); n0$cond = "n0"; n0$rna = 114
n1 = read.csv("data/polysomes/240306_R3_N1_polysomes.csv", skip = 42); n1$cond = "n1"; n1$rna = 75
n4 = read.csv("data/polysomes/240306_R3_N4_polysomes.csv", skip = 42); n4$cond = "n4"; n4$rna = 66
n20 = read.csv("data/polysomes/240306_R3_N20_polysomes.csv", skip = 42); n20$cond = "n20"; n20$rna = 56

mergedf = rbind(n0,n1,n4,n20)
mergedf$normabs = mergedf$Absorbance*100/mergedf$rna
mergedf$cond = factor(as.character(mergedf$cond), levels=c("n0","n1","n4","n20"))

p = ggplot(mergedf, aes(Position.mm.,normabs, color=cond)) +
  geom_rect(data=data.frame(xmin=22, xmax=31, ymin=-Inf, ymax=Inf),
            aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
            alpha=0.4, fill="lightgrey", inherit.aes=FALSE) +
  geom_rect(data=data.frame(xmin=38, xmax=70, ymin=-Inf, ymax=Inf),
            aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
            alpha=0.4, fill="lightgrey", inherit.aes=FALSE) +
  geom_line() +
  xlim(10,80) + ylim(0,1.7) +
  scale_color_manual(values=c("#B2182B","#D6604D","#92C5DE","#2166AC")) +
  theme_bw()
ggsave("figure_panels/polyprof_Nit.pdf", p, width = 4, height = 2.2)

ratiodf = data.frame(cond=factor(c("n0","n1","n4","n20"), levels = c("n0","n1","n4","n20")))
ratiodf$poly = sapply(ratiodf$cond, function(x) AUC(mergedf[mergedf$cond==x,"Position.mm."],mergedf[mergedf$cond==x,"normabs"], from=38, to=70))
ratiodf$mono = sapply(ratiodf$cond, function(x) AUC(mergedf[mergedf$cond==x,"Position.mm."],mergedf[mergedf$cond==x,"normabs"], from=22, to=31))
ratiodf$poly2mono = ratiodf$poly/ratiodf$mono
ratiodf$mono2poly = ratiodf$mono/ratiodf$poly
ratiodf$total = sapply(ratiodf$cond, function(x) AUC(mergedf[mergedf$cond==x,"Position.mm."],mergedf[mergedf$cond==x,"normabs"], from=22, to=70))
ratiodf$polyperc = ratiodf$poly*100/ratiodf$total
p = ggplot(ratiodf, aes(cond,polyperc,fill=cond)) + 
  geom_bar(stat="identity") +
  ylim(0,80) +
  scale_fill_manual(values=c("#B2182B","#D6604D","#92C5DE","#2166AC")) +
  theme_bw()
ggsave("figure_panels/polyperc_Nit.pdf", p, width = 2.6, height = 2.2)


# Load data PHOSPHATE
# No need to calibrate per RNA since same amount was loaded to each. Align profiles based on monosome peak
p0 = read.csv("data/polysomes/240514_R3_P0_polysomes.csv", skip = 42); p0$cond = "p0"
p1 = read.csv("data/polysomes/240514_R3_P1_polysomes.csv", skip = 42); p1$cond = "p1"
p1$Position.mm. = p1$Position.mm. - 1
p2 = read.csv("data/polysomes/240514_R3_P2_polysomes.csv", skip = 42); p2$cond = "p2"
p4 = read.csv("data/polysomes/240514_R3_P4_polysomes.csv", skip = 42); p4$cond = "p4"
p4$Position.mm. = p4$Position.mm. - 1
p7 = read.csv("data/polysomes/240514_R3_P7_polysomes.csv", skip = 42); p7$cond = "p7"
p7$Position.mm. = p7$Position.mm. - 1
p10 = read.csv("data/polysomes/240514_R3_P10_polysomes.csv", skip = 42); p10$cond = "p10"
p10$Position.mm. = p10$Position.mm. - 1
mergedf = rbind(p0,p1,p2,p4,p7,p10)
mergedf$normabs = mergedf$Absorbance
mergedf$cond = factor(as.character(mergedf$cond), levels=c("p0","p1","p2","p4","p7","p10"))

p = ggplot(mergedf, aes(Position.mm.,normabs, color=cond)) +
  geom_rect(data=data.frame(xmin=22, xmax=31, ymin=-Inf, ymax=Inf),
            aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
            alpha=0.4, fill="lightgrey", inherit.aes=FALSE) +
  geom_rect(data=data.frame(xmin=38, xmax=70, ymin=-Inf, ymax=Inf),
            aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
            alpha=0.4, fill="lightgrey", inherit.aes=FALSE) +
  geom_line() +
  xlim(10,80) + ylim(0,0.6) +
  scale_color_manual(values=c("#B2182B","#D6604D","#F4A582","#92C5DE","#4393C3","#2166AC")) +
  theme_bw()
ggsave("figure_panels/polyprof_Pho.pdf", p, width = 4, height = 2.2)

ratiodf = data.frame(cond=factor(c("p0","p1","p2","p4","p7","p10"), levels = c("p0","p1","p2","p4","p7","p10")))
ratiodf$poly = sapply(ratiodf$cond, function(x) AUC(mergedf[mergedf$cond==x,"Position.mm."],mergedf[mergedf$cond==x,"normabs"], from=38, to=70))
ratiodf$mono = sapply(ratiodf$cond, function(x) AUC(mergedf[mergedf$cond==x,"Position.mm."],mergedf[mergedf$cond==x,"normabs"], from=22, to=31))
ratiodf$poly2mono = ratiodf$poly/ratiodf$mono
ratiodf$mono2poly = ratiodf$mono/ratiodf$poly
ratiodf$total = sapply(ratiodf$cond, function(x) AUC(mergedf[mergedf$cond==x,"Position.mm."],mergedf[mergedf$cond==x,"normabs"], from=22, to=70))
ratiodf$polyperc = ratiodf$poly*100/ratiodf$total
p = ggplot(ratiodf, aes(cond,polyperc,fill=cond)) + 
  geom_bar(stat="identity") +
  ylim(0,80) +
  scale_fill_manual(values=c("#B2182B","#D6604D","#F4A582","#92C5DE","#4393C3","#2166AC")) +
  theme_bw()
ggsave("figure_panels/polyperc_Pho.pdf", p, width = 3.2, height = 2.2)

