#!/usr/bin/env python3
"""
comparative_prep.py — Automates comparative genomics for phage characterization.

Given taxmyPHAGE output and the query genome, this script:
  1. Parses taxonomy from Summary_taxonomy.tsv (genus -> family -> order fallback)
  2. Downloads all reference genomes for that taxon via ICTV MSL + NCBI
  3. Runs FastANI all-vs-all (query + references)
  4. Calculates intergenomic similarity (VIRIDIC-equivalent via BLASTn)
  5. Selects the top-N most similar genomes for Clinker
  6. Downloads their GenBank files for Clinker input

Usage:
  python3 comparative_prep.py \\
    --taxmyphage results_B3/taxmy/Summary_taxonomy.tsv \\
    --genome     results_B3/B3.fasta \\
    --msl        ICTV_Master_Species_List_2025_MSL41.v1.xlsx \\
    --email      user@email.com \\
    --outdir     results_B3/comparative \\
    --threads    40 \\
    --top-clinker 6
"""

import argparse
import csv
import shutil
import subprocess
import sys
import textwrap
from pathlib import Path

import pandas as pd
from Bio import Entrez, SeqIO
from Bio.Blast import NCBIXML

sys.path.insert(0, str(Path(__file__).parent))
import ictv_ncbi_fetcher as fetcher

RANK_ORDER = ["Genus", "Family", "Order"]


# ─────────────────────────────────────────────────────────────────────────────
# Taxonomy helpers (shared logic with phylo_prep)
# ─────────────────────────────────────────────────────────────────────────────

def parse_taxmyphage(tsv_path: str) -> dict:
    df = pd.read_csv(tsv_path, sep="\t", dtype=str).fillna("")
    col_map = {}
    for col in df.columns:
        cl = col.strip().lower()
        if cl in ("genus", "predicted_genus"):         col_map["genus"]  = col
        elif cl in ("family", "predicted_family"):     col_map["family"] = col
        elif cl in ("order", "predicted_order"):       col_map["order"]  = col
        elif cl in ("class", "predicted_class"):       col_map["class_"] = col
    row = df.iloc[0]
    def get(k): c = col_map.get(k, ""); return row[c].strip() if c and c in row.index else ""
    return {"genus": get("genus"), "family": get("family"),
            "order": get("order"), "class_": get("class_")}


def choose_query_rank(tax: dict, msl_df: pd.DataFrame, min_species: int = 5):
    for rank in RANK_ORDER:
        key = rank.lower()
        value = tax.get(key, "").strip()
        if not value or "new" in value.lower() or value.lower() in ("", "unknown", "unclassified"):
            continue
        if rank not in msl_df.columns:
            continue
        n = msl_df[rank].str.lower().str.contains(value.lower(), regex=False, na=False).sum()
        if n >= min_species:
            return rank, value
        print(f"  [{rank}] '{value}' has only {n} species in MSL -> going up one rank")
    raise ValueError("Could not determine taxon from taxmyPHAGE output.")


# ─────────────────────────────────────────────────────────────────────────────
# Genome download
# ─────────────────────────────────────────────────────────────────────────────

def download_reference_genomes(species_list: list, species_meta: dict,
                                out_fasta: str, log: list) -> list:
    species_acc = fetcher.search_accessions(species_list, db="nuccore", seq_type="genome")
    return fetcher.fetch_genomes(species_acc, species_meta, out_fasta, log)


def download_genbank(accession: str, out_path: Path) -> bool:
    """Download a single GenBank record."""
    try:
        handle = fetcher.entrez_call(Entrez.efetch, db="nuccore", id=accession,
                                     rettype="gb", retmode="text")
        out_path.write_text(handle.read())
        handle.close()
        return True
    except Exception as exc:
        print(f"  [WARN] GBK download failed for {accession}: {exc}")
        return False


# ─────────────────────────────────────────────────────────────────────────────
# FastANI
# ─────────────────────────────────────────────────────────────────────────────

def run_fastani(genome_list_path: Path, out_tsv: Path, threads: int) -> pd.DataFrame:
    cmd = [
        "fastANI",
        "--ql", str(genome_list_path),
        "--rl", str(genome_list_path),
        "-t", str(threads),
        "-o", str(out_tsv),
    ]
    print(f"  $ {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"  [WARN] fastANI stderr: {result.stderr[:500]}")

    if not out_tsv.exists():
        return pd.DataFrame()

    df = pd.read_csv(out_tsv, sep="\t", header=None,
                     names=["query", "reference", "ani", "mapped", "total"])
    return df


def ani_vs_query(ani_df: pd.DataFrame, query_path: str) -> pd.DataFrame:
    """Filter FastANI results to only rows where the query is the phage genome."""
    q = str(Path(query_path).resolve())
    mask = ani_df["query"].apply(lambda x: Path(x).resolve() == Path(q))
    sub = ani_df[mask].copy()
    sub["ref_name"] = sub["reference"].apply(lambda x: Path(x).stem)
    return sub.sort_values("ani", ascending=False)


# ─────────────────────────────────────────────────────────────────────────────
# VIRIDIC-equivalent intergenomic similarity (BLASTn-based)
# ─────────────────────────────────────────────────────────────────────────────

def calc_intergenomic_similarity(query_fasta: str, ref_fasta: str,
                                  work_dir: Path, threads: int) -> float:
    """
    Approximate VIRIDIC intergenomic similarity between two genomes.
    Formula: sum(aligned_length * pident) / (2 * max_genome_length) * 100
    Replicates the VIRIDIC BLASTn approach.
    """
    db_path = work_dir / "blastdb"
    out_xml = work_dir / "blast_result.xml"

    try:
        subprocess.run(
            ["makeblastdb", "-in", ref_fasta, "-dbtype", "nucl",
             "-out", str(db_path)],
            capture_output=True, check=True
        )
        subprocess.run(
            ["blastn", "-query", query_fasta, "-db", str(db_path),
             "-out", str(out_xml), "-outfmt", "5",
             "-perc_identity", "70", "-num_threads", str(threads),
             "-task", "blastn"],
            capture_output=True, check=True
        )
    except (subprocess.CalledProcessError, FileNotFoundError):
        return 0.0

    query_len = sum(len(r.seq) for r in SeqIO.parse(query_fasta, "fasta"))
    ref_len   = sum(len(r.seq) for r in SeqIO.parse(ref_fasta,   "fasta"))
    max_len   = max(query_len, ref_len)
    if max_len == 0:
        return 0.0

    score = 0.0
    with open(out_xml) as fh:
        for record in NCBIXML.parse(fh):
            for aln in record.alignments:
                for hsp in aln.hsps:
                    score += hsp.align_length * (hsp.identities / hsp.align_length)

    return min((score / (2 * max_len)) * 100, 100.0)


def run_viridic_similarity(query_fasta: str, ref_fastas: list[str],
                            out_tsv: Path, threads: int) -> pd.DataFrame:
    """
    Calculate VIRIDIC-equivalent similarity of query vs each reference.
    """
    rows = []
    work_dir = out_tsv.parent / "viridic_tmp"
    work_dir.mkdir(exist_ok=True)

    total = len(ref_fastas)
    for i, ref in enumerate(ref_fastas, 1):
        ref_name = Path(ref).stem
        print(f"  [{i}/{total}] VIRIDIC sim -> {ref_name}", end=" ", flush=True)
        tmp = work_dir / ref_name
        tmp.mkdir(exist_ok=True)
        sim = calc_intergenomic_similarity(query_fasta, ref, tmp, threads)
        rows.append({"reference": ref, "name": ref_name, "viridic_sim": round(sim, 2)})
        print(f"-> {sim:.1f}%")

    df = pd.DataFrame(rows).sort_values("viridic_sim", ascending=False)
    df.to_csv(out_tsv, sep="\t", index=False)
    print(f"  [OK] VIRIDIC similarities -> {out_tsv}")
    return df


# ─────────────────────────────────────────────────────────────────────────────
# VIRIDIC (optional — if viridic is installed locally)
# ─────────────────────────────────────────────────────────────────────────────

def try_run_viridic(genome_list: list[str], out_dir: Path, threads: int) -> bool:
    """Try to run VIRIDIC locally. Returns True if successful."""
    for viridic_cmd in ("viridic.pl", "viridic", "VIRIDIC"):
        if shutil.which(viridic_cmd):
            cmd = [viridic_cmd, "-i", ",".join(genome_list),
                   "-o", str(out_dir), "-t", str(threads)]
            print(f"  Running VIRIDIC: {' '.join(cmd)}")
            subprocess.run(cmd)
            return True
    return False


# ─────────────────────────────────────────────────────────────────────────────
# Clinker GBK preparation
# ─────────────────────────────────────────────────────────────────────────────

def download_clinker_gbks(top_refs: pd.DataFrame, metadata_df: pd.DataFrame,
                           out_dir: Path) -> list[Path]:
    """
    Download GenBank files for the top-N references.
    Uses the accession from the download metadata when available.
    """
    out_dir.mkdir(exist_ok=True)
    gbks = []

    # Build name -> accession map from metadata
    acc_map = {}
    if "accession" in metadata_df.columns and "species" in metadata_df.columns:
        for _, row in metadata_df.iterrows():
            stem = Path(row.get("accession", "")).stem
            acc_map[row.get("species", "")] = row["accession"]
            acc_map[stem] = row["accession"]

    for _, row in top_refs.iterrows():
        name = row.get("name") or row.get("ref_name", "")
        acc  = acc_map.get(name, name)  # fallback: use name as accession

        out_path = out_dir / f"{name}.gbk"
        if out_path.exists():
            print(f"  [SKIP] {out_path.name} already exists")
            gbks.append(out_path)
            continue

        print(f"  Downloading GBK: {acc} ({name})")
        if download_genbank(acc, out_path):
            gbks.append(out_path)

    return gbks


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────

def main():
    p = argparse.ArgumentParser(
        prog="comparative_prep.py",
        description=textwrap.dedent(__doc__),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--taxmyphage",   required=True)
    p.add_argument("--genome",       required=True,  help="Query genome FASTA")
    p.add_argument("--gbk",          default=None,   help="Query GBK for Clinker (phold.gbk)")
    p.add_argument("--msl",          required=True)
    p.add_argument("--email",        required=True)
    p.add_argument("--api-key",      dest="api_key", default=None)
    p.add_argument("--outdir",       default="comparative")
    p.add_argument("--threads",      type=int, default=8)
    p.add_argument("--top-clinker",  dest="top_clinker", type=int, default=6,
                   help="Number of closest genomes to include in Clinker (default: 6)")
    p.add_argument("--min-species",  dest="min_species", type=int, default=5)
    p.add_argument("--skip-viridic", dest="skip_viridic", action="store_true",
                   help="Skip VIRIDIC-equivalent similarity calculation")
    args = p.parse_args()

    out_dir = Path(args.outdir)
    out_dir.mkdir(parents=True, exist_ok=True)

    print("\n" + "="*60)
    print("  Comparative Genomics Prep")
    print("="*60)

    # ── Setup Entrez ──────────────────────────────────────────────────────────
    Entrez.email = args.email
    Entrez.tool  = "comparative_prep"
    if args.api_key:
        Entrez.api_key = args.api_key
        fetcher.NCBI_DELAY = 0.12

    # ── 1. Taxonomy ───────────────────────────────────────────────────────────
    print(f"\n[1] Taxonomy: {args.taxmyphage}")
    tax = parse_taxmyphage(args.taxmyphage)
    print(f"  Order: {tax['order']} | Family: {tax['family']} | Genus: {tax['genus']}")

    # ── 2. MSL ────────────────────────────────────────────────────────────────
    print(f"\n[2] MSL: {args.msl}")
    msl_df = fetcher.load_msl(args.msl)
    query_rank, query_value = choose_query_rank(tax, msl_df, args.min_species)
    print(f"  -> Selected rank: {query_rank} '{query_value}'")

    df_filtered  = fetcher.filter_msl(msl_df, query_value, query_rank)
    fetcher.summarize_filtered(df_filtered)
    species_list = df_filtered["Species"].tolist()
    species_meta = df_filtered.set_index("Species").to_dict("index")

    # ── 3. Download reference genomes ─────────────────────────────────────────
    ref_fasta = out_dir / f"{query_value.replace(' ', '_')}_references.fna"
    meta_tsv  = out_dir / "metadata_genomes.tsv"
    log: list = []

    if ref_fasta.exists():
        print(f"\n[3] Reference genomes already downloaded: {ref_fasta}")
        metadata = list(csv.DictReader(open(meta_tsv), delimiter="\t")) if meta_tsv.exists() else []
    else:
        print(f"\n[3] Downloading reference genomes ({len(species_list)} species)...")
        metadata = download_reference_genomes(
            species_list, species_meta, str(ref_fasta), log
        )
        fetcher.write_tsv(metadata, str(meta_tsv), label="Genome metadata")
        fetcher.write_tsv(log, str(out_dir / "not_found_genomes.tsv"), label="Not-found log")

    # Split multi-FASTA into individual files for FastANI / VIRIDIC
    refs_dir = out_dir / "ref_fastas"
    refs_dir.mkdir(exist_ok=True)
    ref_fasta_files = []
    if ref_fasta.exists():
        for rec in SeqIO.parse(str(ref_fasta), "fasta"):
            out_f = refs_dir / f"{rec.id.replace('/', '_')}.fasta"
            if not out_f.exists():
                SeqIO.write(rec, str(out_f), "fasta")
            ref_fasta_files.append(str(out_f))
    print(f"  {len(ref_fasta_files)} reference genomes available")

    # ── 4. FastANI all-vs-all ─────────────────────────────────────────────────
    genome_list_path = out_dir / "genome_list.txt"
    all_fastas = ref_fasta_files + [str(Path(args.genome).resolve())]
    genome_list_path.write_text("\n".join(all_fastas) + "\n")

    fastani_tsv = out_dir / "fastani_results.tsv"
    if fastani_tsv.exists():
        print(f"\n[4] FastANI already computed: {fastani_tsv}")
        ani_df = pd.read_csv(fastani_tsv, sep="\t", header=None,
                             names=["query", "reference", "ani", "mapped", "total"])
    else:
        print(f"\n[4] FastANI all-vs-all ({len(all_fastas)} genomes)...")
        try:
            ani_df = run_fastani(genome_list_path, fastani_tsv, args.threads)
            print(f"  [OK] -> {fastani_tsv}")
        except FileNotFoundError:
            print("  [WARN] fastANI not found -- skipping")
            ani_df = pd.DataFrame()

    # ── 5. VIRIDIC-equivalent similarity ──────────────────────────────────────
    viridic_tsv = out_dir / "viridic_similarity.tsv"
    if args.skip_viridic:
        viridic_df = pd.DataFrame()
    elif viridic_tsv.exists():
        print(f"\n[5] VIRIDIC similarity already computed: {viridic_tsv}")
        viridic_df = pd.read_csv(viridic_tsv, sep="\t")
    else:
        print(f"\n[5] Computing intergenomic similarity (VIRIDIC-equivalent)...")
        if try_run_viridic(all_fastas, out_dir / "viridic_out", args.threads):
            viridic_df = pd.DataFrame()  # VIRIDIC writes its own outputs
        else:
            viridic_df = run_viridic_similarity(
                args.genome, ref_fasta_files, viridic_tsv, args.threads
            )

    # ── 6. Select top-N for Clinker ───────────────────────────────────────────
    print(f"\n[6] Selecting top-{args.top_clinker} genomes for Clinker...")

    if not viridic_df.empty and "viridic_sim" in viridic_df.columns:
        top_refs = viridic_df.head(args.top_clinker)
        sim_col  = "viridic_sim"
        print("  (using VIRIDIC similarity)")
    elif not ani_df.empty:
        top_refs = ani_vs_query(ani_df, args.genome).head(args.top_clinker)
        sim_col  = "ani"
        print("  (using FastANI -- VIRIDIC not available)")
    else:
        print("  [WARN] No similarity data available. Clinker will require manual genome selection.")
        top_refs = pd.DataFrame()
        sim_col  = None

    if not top_refs.empty:
        print(f"\n  Top-{args.top_clinker} closest references:")
        name_col = "name" if "name" in top_refs.columns else "ref_name"
        for _, row in top_refs.iterrows():
            sim = f"{row[sim_col]:.1f}%" if sim_col else ""
            print(f"    {row.get(name_col, '?'):<40} {sim}")

    # ── 7. Download GBKs for Clinker ──────────────────────────────────────────
    clinker_dir = out_dir / "clinker_genomes"
    clinker_dir.mkdir(exist_ok=True)

    # Copy query GBK — use genome stem as filename so Clinker labels it correctly
    if args.gbk and Path(args.gbk).exists():
        query_stem    = Path(args.genome).stem
        query_gbk_dst = clinker_dir / f"{query_stem}_query.gbk"
        shutil.copy(args.gbk, query_gbk_dst)
        print(f"\n[7] Query GBK copied -> {query_gbk_dst}")

    metadata_df = pd.DataFrame(metadata)
    gbk_files   = []
    if not top_refs.empty:
        print(f"  Downloading GBKs for {args.top_clinker} reference genomes...")
        gbk_files = download_clinker_gbks(top_refs, metadata_df, clinker_dir)

    # ── Summary ───────────────────────────────────────────────────────────────
    print(f"\n{'='*60}")
    print("  SUMMARY — Comparative Prep")
    print(f"{'='*60}")
    print(f"  Taxon queried    : {query_rank} '{query_value}'")
    print(f"  Species in MSL   : {len(species_list)}")
    print(f"  References DL'd  : {len(ref_fasta_files)}")
    print(f"  FastANI          : {fastani_tsv}")
    print(f"  VIRIDIC sim      : {viridic_tsv if viridic_tsv.exists() else 'not computed'}")
    print(f"  GBKs for Clinker : {clinker_dir}/ ({len(gbk_files)+1} files)")
    print()
    if gbk_files or (args.gbk and (clinker_dir / f"{Path(args.genome).stem}_query.gbk").exists()):
        all_gbks = list(clinker_dir.glob("*.gbk")) + list(clinker_dir.glob("*.gb"))
        print("  Next step — Clinker:")
        gbk_str = " ".join(str(g) for g in sorted(all_gbks))
        print(f"    clinker {gbk_str} -i 0.3 -o clinker_synteny.html -p clinker_synteny.svg")
    print(f"{'='*60}\n")


if __name__ == "__main__":
    main()
