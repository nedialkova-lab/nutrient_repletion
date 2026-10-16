# ==============================================================================
# 7-5_compute_asite_occupancy.py
# A-site codon occupancies per sample (short 20-23 nt and long 27-32 nt footprints)
# from the RiboNet pre-processing outputs of each repletion dataset.
#
# Manuscript: Extended Data Fig. 3E; inputs of 7-6 (Fig. 2D-E)
# Requires  : RiboNet pre-processing ('ribonet preprocessing') of each dataset
# Inputs    : <ribonet_preprocessing>/x_inputs.parquet, y_outputs.parquet, sample_file.txt
# Outputs   : data/20241210_RNA_RPF/glucose_asite_occupancy.tsv
#             data/20241210_RNA_RPF/nitrogen_asite_occupancy.tsv
#             data/20240529_RNA_RPF/phosphate_asite_occupancy.tsv
#             (provided on Zenodo), e.g.:
#   python 7-5_compute_asite_occupancy.py <glucose_preprocessing>/x_inputs.parquet \
#       <glucose_preprocessing>/y_outputs.parquet data/20241210_RNA_RPF/glucose_asite_occupancy.tsv \
#       --sample_file <glucose_preprocessing>/sample_file.txt --trim 15
# ==============================================================================

"""
Compute A-site codon occupancy (pause scores) from ribosome profiling parquet files.

For each gene/condition/replicate:
  1. Extract CDS counts, trimming the first and last ``trim`` nucleotides.
  2. Normalise each position by the per-(gene, condition, replicate) mean count.
  3. Average normalised pause scores by A-site codon and individual sample.

A-site assignment (nucleotide resolution):
  CDS-relative position p  ->  A-site codon index = (p + 3) // 3

Output: TSV with individual samples (sample_id) as rows and codons as columns.
When a sample_file is provided, replicate indices are mapped to their sample_ids;
otherwise rows are labelled as "<track>_rep<n>".

Usage:
    python 7-5_compute_asite_occupancy.py x_inputs.parquet y_outputs.parquet output.tsv
    python 7-5_compute_asite_occupancy.py x_inputs.parquet y_outputs.parquet output.tsv \\
        --sample_file preprocessing/sample_file.txt --trim 15
"""

import argparse
import sys
import numpy as np
import pandas as pd

STOP_CODONS = {"TAA", "TAG", "TGA"}


def compute_asite_occupancy(
    x_inputs_path: str,
    y_outputs_path: str,
    sample_file: "pd.DataFrame | str | None" = None,
    trim: int = 15,
) -> pd.DataFrame:
    """Compute per-sample A-site codon occupancy from ribosome profiling data.

    Parameters
    ----------
    x_inputs_path : str
        Path to x_inputs.parquet produced by RiboPreprocess.
    y_outputs_path : str
        Path to y_outputs.parquet produced by RiboPreprocess.
    sample_file : pd.DataFrame, str, or None
        Sample metadata (DataFrame or path to TSV) with at least columns
        ``sample_id`` and ``track``.  When provided, rows in the output are
        labelled with ``sample_id``; otherwise with ``<track>_rep<n>``.
    trim : int
        Nucleotides to exclude from each end of the CDS (default 15).

    Returns
    -------
    pd.DataFrame
        Table with sample rows and codon columns containing mean pause scores.
        Stop codons and any codons with fewer than 3 bases are excluded.
        Columns are sorted alphabetically.
    """
    x_c = pd.read_parquet(x_inputs_path)
    y_c = pd.read_parquet(y_outputs_path)

    # Resolve sample_file → track-to-sample_id mapping
    if isinstance(sample_file, str):
        sample_file = pd.read_csv(sample_file, sep="\t")

    # Build mapping: track -> ordered list of sample_ids
    if sample_file is not None and "track" in sample_file.columns:
        track_to_samples: dict = {}
        for track, grp in sample_file.groupby("track", sort=False):
            track_to_samples[track] = grp["sample_id"].tolist()
    else:
        track_to_samples = {}

    cds_y = y_c["cds"]
    conditions = list(cds_y.columns)  # tracks

    # Accumulate per (sample_label, A-site codon)
    pause_sum: dict = {}
    pause_count: dict = {}

    genes_in_y = set(cds_y.index)
    n_skipped_short = 0
    n_skipped_missing = 0

    for gene in x_c.index:
        if gene not in genes_in_y:
            n_skipped_missing += 1
            continue

        cds_len = int(x_c.loc[gene, "cds_len"])
        cds_seq = x_c.loc[gene, "cds"][:cds_len]

        start = trim
        end = cds_len - trim

        if end <= start:
            n_skipped_short += 1
            continue

        codons = [cds_seq[i : i + 3] for i in range(0, cds_len, 3)]

        # A-site codon for each CDS-relative position in [start, end)
        positions = np.arange(start, end)
        a_codon_idx = (positions + 3) // 3
        asite_arr = np.array(
            [
                codons[idx] if idx < len(codons) and len(codons[idx]) == 3 else ""
                for idx in a_codon_idx
            ]
        )
        valid_codon_mask = asite_arr != ""

        for cond in conditions:
            rep_list = cds_y.loc[gene, cond]
            if rep_list is None:
                continue

            try:
                rep_mat = np.stack([np.asarray(r, dtype=float) for r in rep_list])
            except Exception:
                continue

            trimmed = rep_mat[:, start:end]          # (n_reps, n_pos)
            rep_means = np.nanmean(trimmed, axis=1)  # (n_reps,)

            sample_ids = track_to_samples.get(cond)

            for rep_idx in range(rep_mat.shape[0]):
                rep_mean = rep_means[rep_idx]
                if rep_mean <= 0 or np.isnan(rep_mean):
                    continue

                pause = trimmed[rep_idx] / rep_mean  # (n_pos,)

                if sample_ids is not None and rep_idx < len(sample_ids):
                    label = sample_ids[rep_idx]
                else:
                    label = f"{cond}_rep{rep_idx}"

                include = valid_codon_mask & ~np.isnan(pause)
                for p_val, a_codon in zip(pause[include], asite_arr[include]):
                    key = (label, a_codon)
                    if key in pause_sum:
                        pause_sum[key] += p_val
                        pause_count[key] += 1
                    else:
                        pause_sum[key] = p_val
                        pause_count[key] = 1

    if n_skipped_missing > 0:
        print(
            f"[warn] {n_skipped_missing} gene(s) in x_inputs not found in y_outputs — skipped.",
            file=sys.stderr,
        )
    if n_skipped_short > 0:
        print(
            f"[warn] {n_skipped_short} gene(s) with CDS shorter than 2×{trim} nt — skipped.",
            file=sys.stderr,
        )

    if not pause_sum:
        raise RuntimeError("No valid data found — check that the parquet files are compatible.")

    rows = [
        {
            "sample": label,
            "Asite": asite,
            "mean_pause": pause_sum[(label, asite)] / pause_count[(label, asite)],
        }
        for label, asite in pause_sum
    ]
    result = pd.DataFrame(rows).pivot(index="sample", columns="Asite", values="mean_pause")
    result.columns.name = None
    result.index.name = "cond"

    stopcod = [c for c in STOP_CODONS if c in result.columns]
    result = result.drop(columns=stopcod, errors="ignore")
    result = result[sorted(result.columns)]

    # Sort samples by readlength first, then by condition
    if sample_file is not None and "sample_id" in sample_file.columns:
        order_cols = [c for c in ("readlength", "condition") if c in sample_file.columns]
        order_df = sample_file[["sample_id"] + order_cols].drop_duplicates("sample_id")
        order_df = order_df[order_df["sample_id"].isin(result.index)]
        if order_cols:
            order_df = order_df.sort_values(order_cols)
        ordered = order_df["sample_id"].tolist()
        extra = [s for s in result.index if s not in set(ordered)]
        result = result.loc[ordered + extra]

    return result


def plot_asite_occupancy(result: pd.DataFrame) -> "matplotlib.figure.Figure":
    """Plot a heatmap of log2 A-site codon pause scores.

    Parameters
    ----------
    result : pd.DataFrame
        Output of :func:`compute_asite_occupancy` (samples × codons).

    Returns
    -------
    matplotlib.figure.Figure
    """
    import matplotlib.pyplot as plt
    import seaborn as sns
    from matplotlib.colors import LinearSegmentedColormap

    color_map = LinearSegmentedColormap.from_list("", ["#762a83", "#f7f7f7", "#1b7837"])

    v = 2
    fig, ax = plt.subplots(
        figsize=(15, max(0.4 * result.shape[0], 3)),
    )
    sns.heatmap(
        np.log2(result),
        cmap=color_map,
        center=0.0,
        vmin=-v,
        vmax=v,
        ax=ax,
        yticklabels=True,
    )
    ax.set_title("A-site codon occupancy (log2 pause score)")
    ax.set_xlabel("A-site codon")
    ax.set_ylabel("Sample")
    fig.tight_layout()
    return fig


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Compute per-sample A-site codon occupancy from ribosome profiling parquet files."
    )
    parser.add_argument("x_inputs", help="Path to x_inputs.parquet")
    parser.add_argument("y_outputs", help="Path to y_outputs.parquet")
    parser.add_argument("output", help="Path to output TSV file")
    parser.add_argument(
        "--sample_file",
        default=None,
        help="Path to sample_file.txt (TSV with sample_id and track columns). "
             "When provided, rows are labelled with sample_id instead of track_repN.",
    )
    parser.add_argument(
        "--trim",
        type=int,
        default=15,
        help="Nucleotides to exclude from each end of the CDS (default: 15)",
    )
    parser.add_argument(
        "--plot",
        default=None,
        metavar="PDF",
        help="If provided, save a heatmap of log2 pause scores to this PDF path.",
    )
    args = parser.parse_args()

    print(f"Loading {args.x_inputs} ...", file=sys.stderr)
    print(f"Loading {args.y_outputs} ...", file=sys.stderr)

    result = compute_asite_occupancy(
        args.x_inputs,
        args.y_outputs,
        sample_file=args.sample_file,
        trim=args.trim,
    )

    result.to_csv(args.output, sep="\t")
    print(
        f"Saved {result.shape[0]} samples × {result.shape[1]} codons → {args.output}",
        file=sys.stderr,
    )

    if args.plot:
        import matplotlib
        matplotlib.use("Agg")
        fig = plot_asite_occupancy(result)
        fig.savefig(args.plot)
        print(f"Heatmap saved → {args.plot}", file=sys.stderr)


if __name__ == "__main__":
    main()
