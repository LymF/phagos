#!/usr/bin/env python3
"""
Circular genome map — publication-quality figure for bacteriophage genomes.

Tracks (outside -> inside):
  1. Genome scale (tick marks; 1 kb minor, 5 kb major with labels)
  2. CDS forward strand  — colored by PHROG functional category
  3. CDS reverse strand  — same color scheme
  4. Feature markers     — tRNA and CRISPR arrays
  5. GC content          — deviation from genome mean (red = above, blue = below)
  6. GC skew             — (G-C)/(G+C) in sliding windows (green/purple)
  7. Coverage            — read depth (optional; requires samtools depth TSV or BAM)

Required input files:
  GBK   — Phold/Pharokka GenBank output.
          CDS features must have: locus_tag, function, translation qualifiers.

  FASTA — Single-record FASTA of the phage genome.

  DEPO  — DepoScope output CSV (optional).
          Expected columns: phage_ID, gene_ID, gene_sequence, protein_sequence, scores_DepoScore

  DEPTH — (optional) samtools depth -a output: 3 columns, tab-separated,
          no header: chrom  pos  depth

Usage:
  python3 genome_map.py --gbk phold_out/phold.gbk \\
                        --fasta genome.fasta \\
                        --depo deposcope/phage_genes.csv \\
                        --name "Phage X1" \\
                        --out genome_map.pdf

  # Add coverage layer
  python3 genome_map.py ... --depth coverage.tsv

Dependencies:
  pip install pycirclize matplotlib biopython numpy
"""

import argparse
from pathlib import Path

import matplotlib as mpl
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
from Bio import SeqIO
from pycirclize import Circos
from pycirclize.parser import Genbank

# Default paths (override with CLI args)
GBK_PATH     = "phold_out/phold.gbk"
DEPO_CSV     = "deposcope/phage_genes.csv"
GENOME_FASTA = "X1.fasta"
DEPTH_TSV    = "coverage.tsv"
OUT_PDF      = "genome_map.pdf"

# Color scheme — Okabe-Ito + Wong palette (colorblind-safe, Nature Methods)
CATEGORY_COLORS = {
    "head and packaging":                          "#0072B2",
    "tail":                                        "#009E73",
    "DNA, RNA and nucleotide metabolism":          "#E69F00",
    "lysis":                                       "#D55E00",
    "connector":                                   "#CC79A7",
    "transcription regulation":                    "#56B4E9",
    "moron, auxiliary metabolic gene and host takeover": "#F0E442",
    "other":                                       "#999999",
    "unknown function":                            "#DDDDDD",
}
DEPO_COLOR     = "#FF1F5B"
TRNA_COLOR     = "#9932CC"
CRISPR_COLOR   = "#FF8C00"
FORWARD_ALPHA  = 1.0
REVERSE_ALPHA  = 0.85

# Figure style
mpl.rcParams.update({
    "font.family":        "sans-serif",
    "font.sans-serif":    ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size":          7,
    "axes.titlesize":     8,
    "axes.labelsize":     7,
    "xtick.labelsize":    6,
    "ytick.labelsize":    6,
    "legend.fontsize":    7,
    "pdf.fonttype":       42,
    "ps.fonttype":        42,
    "figure.dpi":         300,
    "savefig.dpi":        300,
    "savefig.bbox":       "tight",
})


def gc_content_windows(seq: str, window: int = 500, step: int = 100):
    """Return (midpoint, gc_fraction) pairs for sliding windows."""
    positions, values = [], []
    n = len(seq)
    for start in range(0, n - window + 1, step):
        sub = seq[start : start + window].upper()
        gc  = (sub.count("G") + sub.count("C")) / len(sub)
        positions.append(start + window // 2)
        values.append(gc)
    positions.append(n)
    values.append(values[0])
    return np.array(positions), np.array(values)


def gc_skew_windows(seq: str, window: int = 500, step: int = 100):
    """Return (midpoint, skew) pairs where skew = (G-C)/(G+C)."""
    positions, values = [], []
    n = len(seq)
    for start in range(0, n - window + 1, step):
        sub = seq[start : start + window].upper()
        g, c = sub.count("G"), sub.count("C")
        skew = (g - c) / (g + c) if (g + c) > 0 else 0.0
        positions.append(start + window // 2)
        values.append(skew)
    positions.append(n)
    values.append(values[0])
    return np.array(positions), np.array(values)


def load_depth(depth_path: str, genome_size: int, bin_size: int = 200):
    """Read samtools depth -a output and return binned (midpoint, mean_depth) arrays."""
    bins   = np.zeros(genome_size // bin_size + 1)
    counts = np.zeros_like(bins)
    with open(depth_path) as fh:
        for line in fh:
            parts = line.split()
            if len(parts) < 3:
                continue
            pos, dep = int(parts[1]) - 1, float(parts[2])
            idx = pos // bin_size
            if idx < len(bins):
                bins[idx]   += dep
                counts[idx] += 1
    means = np.where(counts > 0, bins / counts, 0)
    mids  = np.arange(len(means)) * bin_size + bin_size // 2
    mids  = mids[mids < genome_size]
    means = means[:len(mids)]
    return mids.astype(float), means


def depth_from_bam(bam_path: str, out_tsv: str = "coverage.tsv"):
    """Run samtools depth and return the output path."""
    import subprocess
    cmd = f"samtools depth -a -@ 4 {bam_path} > {out_tsv}"
    subprocess.run(cmd, shell=True, check=True)
    return out_tsv


def main():
    parser = argparse.ArgumentParser(
        description="Circular genome map for bacteriophage genomes",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--gbk",     default=GBK_PATH,     help="Phold/Pharokka GenBank file")
    parser.add_argument("--fasta",   default=GENOME_FASTA, help="Genome FASTA")
    parser.add_argument("--depo",    default=DEPO_CSV,     help="DepoScope predictions CSV (optional)")
    parser.add_argument("--depth",   default=DEPTH_TSV,    help="samtools depth -a TSV (optional)")
    parser.add_argument("--bam",     default=None,         help="BAM file (auto-runs samtools depth)")
    parser.add_argument("--name",    default=None,         help="Phage name for title")
    parser.add_argument("--out",     default=OUT_PDF,      help="Output PDF path")
    parser.add_argument("--labels",  action="store_true",  default=True,
                        help="Show CDS numbers and product labels (default: True)")
    parser.add_argument("--no-labels", dest="labels", action="store_false",
                        help="Hide CDS labels")
    parser.add_argument("--label-unknown", action="store_true", default=False,
                        help="Also label unknown/hypothetical proteins")
    parser.add_argument("--label-min-size", type=int, default=0,
                        help="Minimum CDS size (bp) to show label (default: 0)")
    args = parser.parse_args()

    if args.depth and not Path(args.depth).exists():
        args.depth = None

    print("\n[1/6] Loading annotation...")
    gbk    = Genbank(args.gbk)
    genome = next(SeqIO.parse(args.fasta, "fasta"))
    seq    = str(genome.seq)
    gsize  = len(seq)
    print(f"      Genome: {gsize:,} bp  |  {len(list(gbk.extract_features('CDS')))} CDS")

    # Load DepoScore data if provided
    depo_scores: dict[str, float] = {}
    depo_label:  dict[str, str]   = {}
    if args.depo and Path(args.depo).exists():
        print("[2/6] Loading DepoScope scores...")
        seq_to_score:   dict[str, float] = {}
        seq_to_gene_id: dict[str, str]   = {}
        with open(args.depo) as fh:
            next(fh)
            for line in fh:
                parts = line.strip().split(",")
                if len(parts) >= 5:
                    gene_id  = parts[1].strip()
                    prot_seq = parts[3].strip()
                    score    = float(parts[4])
                    seq_to_score[prot_seq]   = score
                    seq_to_gene_id[prot_seq] = gene_id

        for feat in SeqIO.read(args.gbk, "genbank").features:
            if feat.type != "CDS":
                continue
            translation = feat.qualifiers.get("translation", [""])[0].strip()
            if translation in seq_to_score:
                locus = feat.qualifiers.get("locus_tag",
                        feat.qualifiers.get("gene", ["?"]))[0]
                depo_scores[locus] = seq_to_score[translation]
                depo_label[locus]  = seq_to_gene_id[translation].split("_")[-1]

        matched = [(k, v) for k, v in depo_scores.items() if v >= 0.5]
        print(f"      {len(seq_to_score)} proteins -> {len(depo_scores)} matched")
        for locus, score in matched:
            print(f"      DepoScore hit: {locus} ({depo_label[locus]})  score={score:.4f}")

    # Collect CDS list for numbering
    cds_ordered = []
    for feat in SeqIO.read(args.gbk, "genbank").features:
        if feat.type == "CDS":
            cds_ordered.append(feat)
    cds_ordered.sort(key=lambda f: int(f.location.start))
    locus_to_num = {}
    for idx, feat in enumerate(cds_ordered, 1):
        locus = feat.qualifiers.get("locus_tag",
                feat.qualifiers.get("gene", [f"CDS{idx}"]))[0]
        locus_to_num[locus] = idx

    # tRNA and CRISPR positions
    trna_positions:   list[tuple[int, int]] = []
    crispr_positions: list[tuple[int, int]] = []
    for feat in SeqIO.read(args.gbk, "genbank").features:
        s = int(feat.location.start)
        e = int(feat.location.end)
        if feat.type == "tRNA":
            trna_positions.append((s, e))
        elif feat.type in ("repeat_region", "CRISPR"):
            note = " ".join(feat.qualifiers.get("note", []) +
                            feat.qualifiers.get("rpt_type", []))
            if "CRISPR" in note.upper() or feat.type == "CRISPR":
                crispr_positions.append((s, e))
    print(f"      tRNA: {len(trna_positions)}  |  CRISPR: {len(crispr_positions)}")

    # Coverage
    depth_pos = depth_val = None
    if args.bam:
        print("[2b] Computing coverage from BAM...")
        depth_file  = depth_from_bam(args.bam)
        depth_pos, depth_val = load_depth(depth_file, gsize)
    elif args.depth:
        print("[2b] Loading coverage depth file...")
        depth_pos, depth_val = load_depth(args.depth, gsize)
    else:
        print("     No coverage data provided.")

    # GC content and skew
    print("[3/6] Computing GC content and skew...")
    gc_pos,   gc_val   = gc_content_windows(seq, window=500, step=50)
    skew_pos, skew_val = gc_skew_windows(seq,    window=500, step=50)
    gc_mean            = np.mean(gc_val[:-1])

    print("[4/6] Building Circos layout...")
    R = {
        "scale_outer":  93,
        "fwd_outer":    92,  "fwd_inner":  86,
        "rev_outer":    85,  "rev_inner":  79,
        "feat_outer":   78,  "feat_inner": 74,
        "gc_outer":     73,  "gc_inner":   65,
        "skew_outer":   64,  "skew_inner": 56,
        "cov_outer":    55,  "cov_inner":  47,
    }

    circos = Circos(sectors={genome.id: gsize})
    phage_name = args.name or Path(args.fasta).stem
    circos.text(f"{phage_name}\n{gsize // 1000} kb", size=9, weight="bold", r=15)

    sector = circos.sectors[0]

    # Scale track (clock style, inner)
    scale = sector.add_track((42, 47))
    scale.axis(fc="none", ec="none")

    ruler_1k_pos = list(range(0, gsize, 1000))
    scale.xticks(ruler_1k_pos, labels=None, tick_length=1.5, outer=False,
                line_kws=dict(ec="#888888", lw=0.4))

    ruler_5k_pos = list(range(0, gsize, 5000))
    ruler_5k_labels = [f"{p // 1000}kb" if p > 0 else "" for p in ruler_5k_pos]
    scale.xticks(ruler_5k_pos, labels=ruler_5k_labels, tick_length=2.5, outer=False,
                label_size=4, label_orientation="vertical", label_margin=1,
                line_kws=dict(ec="#333333", lw=0.6))

    # CDS tracks
    fwd_track = sector.add_track((R["fwd_inner"], R["fwd_outer"]))
    rev_track = sector.add_track((R["rev_inner"], R["rev_outer"]))
    fwd_track.axis(fc="none", ec="none")
    rev_track.axis(fc="none", ec="none")

    label_track = sector.add_track((R["scale_outer"] + 1, R["scale_outer"] + 2))
    label_track.axis(fc="none", ec="none")

    UNKNOWN_FUNCS = {"unknown function", "hypothetical protein", ""}

    for feat in gbk.extract_features("CDS"):
        start   = int(feat.location.start)
        end     = int(feat.location.end)
        strand  = feat.location.strand
        quals   = feat.qualifiers

        fn      = quals.get("function", ["unknown function"])[0]
        color   = CATEGORY_COLORS.get(fn, CATEGORY_COLORS["unknown function"])

        dscore  = depo_scores.get(quals.get("locus_tag", [""])[0], 0.0)
        if dscore >= 1.0:
            color = DEPO_COLOR

        gene_deg = (end - start) / gsize * 360
        head_deg = min(1.5, gene_deg * 0.25)
        head_deg = max(head_deg, 0.3)

        locus = quals.get("locus_tag", quals.get("gene", ["?"]))[0]
        cds_num = locus_to_num.get(locus, "?")
        mid = (start + end) / 2

        if strand == 1:
            fwd_track.arrow(start, end, head_length=head_deg, shaft_ratio=0.5,
                            fc=color, ec="none", alpha=FORWARD_ALPHA)
            track_for_num = fwd_track
        else:
            rev_track.arrow(end, start, head_length=head_deg, shaft_ratio=0.5,
                            fc=color, ec="none", alpha=REVERSE_ALPHA)
            track_for_num = rev_track

        if gene_deg >= 3.5:
            track_for_num.text(str(cds_num), x=mid, size=4.5, color="white",
                              weight="bold", ha="center", va="center")

        # CDS labels with leader lines
        if args.labels:
            product = quals.get("product", ["hypothetical protein"])[0]
            fn_norm = fn.strip().lower()
            is_unk  = fn_norm in UNKNOWN_FUNCS or "hypothetical" in product.lower()

            size_ok = (end - start) >= args.label_min_size
            show    = size_ok and (not is_unk or args.label_unknown or dscore >= 1.0)

            if show:
                if len(product) > 45:
                    product = product[:42] + ".."
                label_txt = product
                if dscore >= 1.0:
                    label_txt = f"* {product}"

                label_track.annotate(
                    mid, label_txt,
                    min_r=R["scale_outer"] + 1,
                    max_r=R["scale_outer"] + 5,
                    label_size=4.5,
                    line_kws=dict(color=color, lw=0.5),
                    text_kws=dict(color="#222222", weight="normal",
                                  ha="left" if mid < gsize / 2 else "right"),
                )

    # Feature track
    feat_track = sector.add_track((R["feat_inner"], R["feat_outer"]))
    feat_track.axis(fc="none", ec="none")
    for s, e in trna_positions:
        feat_track.rect(s, e, fc=TRNA_COLOR, ec="none", alpha=0.9)
    for s, e in crispr_positions:
        feat_track.rect(s, e, fc=CRISPR_COLOR, ec="none", alpha=0.9)

    # GC content track
    gc_track = sector.add_track((R["gc_inner"], R["gc_outer"]))
    gc_track.axis(fc="#F8F8F8", ec="#CCCCCC", lw=0.5)

    gc_dev = gc_val - gc_mean
    gc_abs_max = np.abs(gc_dev[:-1]).max() or 1e-6

    r_mid_gc   = (R["gc_inner"] + R["gc_outer"]) / 2
    r_range_gc = (R["gc_outer"] - R["gc_inner"]) / 2
    for i in range(len(gc_pos) - 1):
        x0, x1 = gc_pos[i], gc_pos[i + 1]
        val     = gc_dev[i]
        if abs(val) < 1e-9:
            continue
        norm  = val / gc_abs_max
        r_val = r_mid_gc + norm * r_range_gc
        color_gc = "#D55E00" if val >= 0 else "#0072B2"
        lim = (r_mid_gc, r_val) if val >= 0 else (r_val, r_mid_gc)
        gc_track.rect(x0, x1, r_lim=lim, fc=color_gc, ec="none", alpha=0.7)

    gc_track.line(gc_pos, np.full_like(gc_pos, r_mid_gc), lw=0.4, color="#666666", ls="dashed")

    # GC skew track
    skew_track = sector.add_track((R["skew_inner"], R["skew_outer"]))
    skew_track.axis(fc="#F8F8F8", ec="#CCCCCC", lw=0.5)

    skew_abs_max = np.abs(skew_val[:-1]).max() or 1e-6
    r_mid_sk     = (R["skew_inner"] + R["skew_outer"]) / 2
    r_range_sk   = (R["skew_outer"] - R["skew_inner"]) / 2

    for i in range(len(skew_pos) - 1):
        x0, x1 = skew_pos[i], skew_pos[i + 1]
        val     = skew_val[i]
        if abs(val) < 1e-9:
            continue
        norm  = val / skew_abs_max
        r_val = r_mid_sk + norm * r_range_sk
        color_sk = "#009E73" if val >= 0 else "#CC79A7"
        lim = (r_mid_sk, r_val) if val >= 0 else (r_val, r_mid_sk)
        skew_track.rect(x0, x1, r_lim=lim, fc=color_sk, ec="none", alpha=0.7)

    skew_track.line(skew_pos, np.full_like(skew_pos, r_mid_sk), lw=0.4, color="#666666", ls="dashed")

    # Coverage track
    if depth_pos is not None and depth_val is not None:
        cov_track = sector.add_track((R["cov_inner"], R["cov_outer"]))
        cov_track.axis(fc="#F8F8F8", ec="#CCCCCC", lw=0.5)
        cov_max = depth_val.max() or 1

        for i in range(len(depth_pos) - 1):
            x0, x1 = depth_pos[i], depth_pos[i + 1]
            norm    = depth_val[i] / cov_max
            r_val   = R["cov_inner"] + norm * (R["cov_outer"] - R["cov_inner"])
            cov_track.rect(x0, x1, r_lim=(R["cov_inner"], r_val),
                          fc="#56B4E9", ec="none", alpha=0.8)

    print("[5/6] Rendering figure...")
    fig = circos.plotfig(figsize=(8, 11))
    fig.subplots_adjust(bottom=0.22, top=0.94, left=0.05, right=0.87)
    ax  = fig.axes[0]

    # Legend
    legend_items = []
    for cat, color in CATEGORY_COLORS.items():
        if cat == "unknown function":
            label = "Unknown function"
        elif cat == "moron, auxiliary metabolic gene and host takeover":
            label = "Moron / AMG"
        else:
            label = cat.capitalize()
        legend_items.append(mpatches.Patch(fc=color, ec="#AAAAAA", lw=0.3, label=label))
    legend_items += [
        mpatches.Patch(fc=DEPO_COLOR,   ec="none", label="Depolymerase (DepoScore = 1.0)"),
        mpatches.Patch(fc=TRNA_COLOR,   ec="none", label="tRNA"),
        mpatches.Patch(fc=CRISPR_COLOR, ec="none", label="CRISPR array"),
        mpatches.Patch(fc="#D55E00",    ec="none", label=f"GC > {gc_mean*100:.1f}% (mean)"),
        mpatches.Patch(fc="#0072B2",    ec="none", label=f"GC < {gc_mean*100:.1f}% (mean)"),
        mpatches.Patch(fc="#009E73",    ec="none", label="GC skew G>C"),
        mpatches.Patch(fc="#CC79A7",    ec="none", label="GC skew C>G"),
    ]

    fig.legend(
        handles=legend_items, loc="lower center", bbox_to_anchor=(0.46, 0.01),
        frameon=True, framealpha=0.95, edgecolor="#CCCCCC", fontsize=6,
        title="Functional category", title_fontsize=6.5,
        handlelength=1.2, handleheight=0.9, ncol=4,
    )

    # Track labels
    track_labels = [
        (R["fwd_outer"]  + R["fwd_inner"])  / 2, "CDS (+)",
        (R["rev_outer"]  + R["rev_inner"])  / 2, "CDS (-)",
        (R["feat_outer"] + R["feat_inner"]) / 2, "tRNA / CRISPR",
        (R["gc_outer"]   + R["gc_inner"])   / 2, "GC content",
        (R["skew_outer"] + R["skew_inner"]) / 2, "GC skew",
    ]
    if depth_pos is not None:
        track_labels += [(R["cov_outer"] + R["cov_inner"]) / 2, "Coverage"]

    for i in range(0, len(track_labels), 2):
        r_norm = track_labels[i] / 100
        label  = track_labels[i + 1]
        ax.text(1.04, 0.5 + (r_norm - 0.69) * 0.9, label, transform=ax.transAxes,
                va="center", ha="left", fontsize=6, color="#444444")

    fig.suptitle(f"{phage_name} — Genomic map", fontsize=9, fontweight="bold", y=0.97)

    print("[6/6] Saving output files...")
    out_pdf = Path(args.out)
    out_png = out_pdf.with_suffix(".png")

    fig.savefig(out_pdf, format="pdf", bbox_inches="tight")
    fig.savefig(out_png, format="png", dpi=600, bbox_inches="tight")

    print(f"\nDone:")
    print(f"  {out_pdf}  (vector PDF for submission)")
    print(f"  {out_png}  (300 DPI PNG for preview)")

    if depth_pos is None:
        fasta = args.fasta
        print(f"""
To add coverage track, map reads and run:
  bwa mem {fasta} reads_R1.fastq.gz reads_R2.fastq.gz | \\
      samtools sort -o reads_mapped.bam && samtools index reads_mapped.bam
  samtools depth -a reads_mapped.bam > coverage.tsv
  python3 genome_map.py --depth coverage.tsv
""")


if __name__ == "__main__":
    main()
