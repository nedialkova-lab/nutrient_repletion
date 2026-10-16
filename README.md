# Repletions

Code to reproduce the analyses and figures of:

> Hernandez-Alias X, Kaya S, Doychev K, Schroeder F, Hermann B, Madran B, Majidian S, Gagneur J, Nedialkova DD.
> **Metabolic state rapidly reprograms translation through initiation and elongation control.**

## Quick start

```bash
# 1. Clone the repository
git clone https://github.com/nedialkova-lab/nutrient_repletions.git
cd nutrient_repletions

# 2. Download the source data from Zenodo into the repository root

# 3. Unzip it at the repository root: this creates data/, results/seqlets/ and source_data/
unzip repletions_data.zip

# 4. Run the scripts from the repository root, in numerical order (see below)
Rscript 1-1_explore_proteomics.R
```

All paths in the code are relative to the repository root, so every script (and notebook) must be run with the
repository root as working directory (e.g. open the folder as an RStudio project, or `setwd()` to it). Outputs are
written to `results/` (analysis tables), `source_data/` (figure source data), `figure_panels/` (PDF panels) and
`plots/` (interactive HTML). The figures were assembled from the PDF panels.

## Repository structure

```
nutrient_repletions/
├── <N-M>_<analysis>.{R/py/ipynb}   # analysis scripts, run in numerical order
├── data/          # input data, from Zenodo (not tracked)
├── results/       # analysis results (results/seqlets/ comes from Zenodo)
├── source_data/   # figure source data (source_data/LiPMS_vs_Metab.tsv comes from Zenodo)
├── figure_panels/ # PDF figure panels
└── plots/         # interactive HTML plots
```

Supplementary Tables are assembled from the files written to `results/` and `source_data/`.

Each script starts with a header describing its purpose, the manuscript panels it produces, its inputs, its outputs
and the scripts it depends on.

## Scripts and figures

| Script | Analysis | Figures | Depends on |
|---|---|---|---|
| `1-1_explore_proteomics.R` | Total proteome of depleted cultures vs SD | ED Fig. 1 | |
| `2-1_targeted_metabolomics_4973.R` | Targeted amino acid metabolomics | ED Fig. 5G-H | |
| `2-2_unfiltered_metabolomics_4973.R` | Untargeted metabolomics: normalisation, batch correction, differential analysis | Fig. 1C-D; ED Fig. 2C | |
| `2-3_clustering_metabolomics_4973.R` | Clustering and mummichog pathway enrichment of metabolites | Fig. 1B; ED Fig. 2A-B | 2-2 |
| `3-1_explore_ctrl_lipms.R` | Exploratory: total proteome vs LiP-MS TrP controls | – | |
| `3-2_explore_LiPMS_5199.R` | Exploratory: LiP-MS quality control (PCA, replicates, batch) | – | |
| `3-3_LiPMS_batchcorrection.R` | LiP-MS / TrP extraction and batch correction | – | |
| `3-4_abundances_timecourse_5199.R` | Differential protein abundance | Fig. 5B | 3-3 |
| `3-5_clustering_abundances_5199.R` | Clustering / GO of protein abundance | Fig. 5C | 3-4 |
| `4-1_analyze_LiPMS_5199.R` | LiP-MS structural residuals and differential analysis | Fig. 5B; ED Fig. 6D | 3-3, 3-4 |
| `4-2_clustering_LiPMS_5199.R` | Clustering / GO of LiP-MS peptides | Fig. 5C; ED Fig. 6A-C | 4-1 |
| `5-1_LiPMS_vs_Metab.R` | Peptide–metabolite correlations (and significant subset) | Fig. 5D, F, H-L; ED Fig. 6E-I | 2-2, 4-1 |
| `5-2_structure_mapping_ribosome.ipynb` | Significant peptides on the 80S ribosome (PDB 8Z70) | Fig. 5E | 4-1, 5-1 |
| `5-2_structure_mapping_other.ipynb` | Significant peptides on MES1, SES1, CDC60, STM1 (AlphaFold) and TEF1/2 (PDB 1F60) | Fig. 5H, K; ED Fig. 6H | 4-1, 5-1 |
| `6-1_explore_RNAseq.R` | RNA-seq differential expression | Fig. 2A | |
| `6-2_clustering_RNAseq.R` | Clustering / GO of RNA-seq | Sup. Table | 6-1 |
| `7-0_analyze_polysomes.R` | Polysome profiles | Fig. 1E-F; ED Fig. 2D | |
| `7-1_explore_riboseq.R` | Ribosome profiling differential analysis | Fig. 2A | 6-1 |
| `7-2_clustering_riboseq.R` | Clustering / GO of ribosome footprints | Sup. Table | 7-1 |
| `7-3_spikein_normalization.R` | Spike-in calibration of ribosomes per transcript (RxT) | Fig. 2B; ED Fig. 3A | |
| `7-4_clustering_RxT.R` | Clustering / GO of RxT | Fig. 2C; ED Fig. 3B-C | 7-3 |
| `7-5_compute_asite_occupancy.py` | A-site codon occupancies from RiboNet pre-processing outputs | ED Fig. 3E | RiboNet |
| `7-6_initiation_vs_elongation.R` | Codon occupancies, elongation and initiation estimates | Fig. 2D-E; ED Fig. 3E | 7-3, 7-5 |
| `8-1_explore_tRNAseq.R` | tRNA abundance and aminoacylation | Fig. 4C; ED Fig. 5C | |
| `8-2_tRNA_vs_codon_usage.R` | Codon occupancy vs tRNA correlations, codon GSEA | Fig. 4D-E; ED Fig. 5D-F | 7-3, 7-5, 8-1 |
| `9-1_ribonet_candidates.R` | Exploratory: RxT/RPF/RNA time courses of candidate genes | – | 6-1, 7-1, 7-3 |
| `9-2_ribonet_seqlets.R` | GO and codon-usage enrichment of GAUY genes | Sup. Table | RiboNet |
| `10-1_attune_CPR8_ONdeprivation.R` | Flow cytometry, CPR8 5'UTR reporters | Fig. 6G; ED Fig. 7H | |
| `10-2_attune_TYW1_Ctrls_ONdeprivation.R` | Flow cytometry, TYW1 reporters and Citrine controls (rep. 1-2) | ED Fig. 7G, J | |
| `10-3_attune_TYW1_Ctrls_ONdeprivation.R` | Flow cytometry, TYW1 reporters (rep. 3) | ED Fig. 7J | |
| `10-4_qPCR_CPR8_ONdeprivation.R` | RT-qPCR of reporters | ED Fig. 7K | |

## Software requirements

Analyses were run with R (Bioconductor) and Python 3. Versions of the main tools are given in the Methods.

R packages (CRAN): `cluster`, `DescTools`, `dendextend`, `dplyr`, `eulerr`, `factoextra`, `fdrtool`, `GeneNet`,
`ggplot2`, `ggpubr`, `ggrastr`, `ggrepel`, `gplots`, `heatmaply`, `htmlwidgets`, `openxlsx`, `patchwork`, `plyr`,
`RColorBrewer`, `readxl`, `reshape`, `reshape2`, `seqinr`, `stringr`, `tidyverse`.

R packages (Bioconductor): `clusterProfiler`, `edgeR`, `flowCore`, `flowStats`, `flowViz`, `ggtree`, `graph`,
`limma`, `MicrobiomeProfiler`, `msigdbr`, `org.Sc.sgd.db`, `Rgraphviz`, `RUVSeq`, `UniProt.ws`.

R packages (GitHub): `MetaboAnalystR` (xia-lab/MetaboAnalystR, v4.3.0), `LiPAnalyzeR` (Nagel et al., Mol. Cell. Proteomics 2025).

Python packages: `numpy`, `pandas`, `pyarrow`, `requests`, `py3Dmol`,
`freesasa`, `gemmi`, `scipy`, `seaborn`, `matplotlib`, `plotnine`.

## Contact

Xavier Hernandez-Alias (xhernandez@ub.edu), Julien Gagneur (gagneur@in.tum.de),
Danny D. Nedialkova (danny.nedialkova@maxperutzlabs.ac.at).
