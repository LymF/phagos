#!/usr/bin/env python3
"""
phylo_prep.py — Prepare marker-protein FASTAs for phage phylogenetic analysis.

Given taxmyPHAGE output and a phold GBK, this script:
  1. Parses the taxonomy (family, genus, order) from Summary_taxonomy.tsv
  2. Chooses the right ICTV taxon level to query (genus -> family -> order)
  3. Picks an outgroup from a distant clade via the MSL
  4. Downloads marker proteins for all species in the clade (ictv_ncbi_fetcher)
  5. Extracts the query phage's own proteins from phold.gbk and appends them

Output: one .faa per marker in --outdir, ready for MAFFT -> IQ-TREE2.

Usage:
  python3 phylo_prep.py \\
    --taxmyphage results_B3/taxmy/Summary_taxonomy.tsv \\
    --gbk        results_B3/phold/phold.gbk \\
    --msl        ../results-paper-X1/ICTV_Master_Species_List_2025_MSL41.v1.xlsx \\
    --markers    "terminase large subunit,major capsid protein" \\
    --name       "Ralstonia phage B3" \\
    --email      user@email.com \\
    --outdir     results_B3/phylo/markers
"""

import argparse
import sys
import textwrap
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Import sibling script — assumes both live in the same directory
sys.path.insert(0, str(Path(__file__).parent))
import ictv_ncbi_fetcher as fetcher

# ─────────────────────────────────────────────────────────────────────────────
# Taxonomy helpers
# ─────────────────────────────────────────────────────────────────────────────

RANK_ORDER = ["Genus", "Family", "Order", "Class"]

# For each rank, the rank one level up (for outgroup selection)
OUTGROUP_RANK = {
    "Genus":  "Family",
    "Family": "Order",
    "Order":  "Class",
    "Class":  "Class",
}

# Well-characterised reference genera to use as outgroup when nothing better
# is found in the MSL. These cover most dsDNA phage classes.
FALLBACK_OUTGROUPS = [
    "Tequatrovirus",    # T4-like, Myoviridae — almost universal for Caudoviricetes
    "Teseptimavirus",   # T7-like, Autographiviridae
    "Lambdavirus",      # lambda-like, Siphoviridae
]


def parse_taxmyphage(tsv_path: str) -> dict:
    """
    Read taxmyPHAGE Summary_taxonomy.tsv and return a dict with keys:
      genus, family, order, class_, species, ictv_status
    Handles both classified and 'new genus/species' outcomes.
    """
    df = pd.read_csv(tsv_path, sep="\t", dtype=str).fillna("")

    # taxmyPHAGE uses different column naming across versions — be tolerant
    col_map = {}
    for col in df.columns:
        cl = col.strip().lower()
        if cl in ("genus", "predicted_genus"):
            col_map["genus"] = col
        elif cl in ("family", "predicted_family"):
            col_map["family"] = col
        elif cl in ("order", "predicted_order"):
            col_map["order"] = col
        elif cl in ("class", "predicted_class"):
            col_map["class_"] = col
        elif "species" in cl and "ictv" not in cl:
            col_map.setdefault("species", col)
        elif "ictv" in cl and "status" in cl:
            col_map["ictv_status"] = col
        elif cl == "classification":
            col_map["ictv_status"] = col

    row = df.iloc[0]

    def get(key: str) -> str:
        c = col_map.get(key, "")
        return row[c].strip() if c and c in row.index else ""

    return {
        "genus":       get("genus"),
        "family":      get("family"),
        "order":       get("order"),
        "class_":      get("class_"),
        "species":     get("species"),
        "ictv_status": get("ictv_status"),
    }


def choose_query_rank(tax: dict, msl_df: pd.DataFrame,
                      min_species: int = 5) -> tuple[str, str]:
    """
    Walk RANK_ORDER top-down until we find a rank that contains >= min_species
    in the MSL.  Returns (rank_name, taxon_value).
    """
    for rank in RANK_ORDER:
        key = "class_" if rank == "Class" else rank.lower()
        value = tax.get(key, "").strip()
        if not value or value.lower() in ("", "unknown", "unclassified"):
            continue
        # Skip "new genus" / "new family" strings
        if "new" in value.lower():
            continue
        col = rank if rank in msl_df.columns else None
        if col is None:
            continue
        n = msl_df[col].str.lower().str.contains(
            value.lower(), regex=False, na=False
        ).sum()
        if n >= min_species:
            return rank, value
        print(f"  [{rank}] '{value}' has only {n} species in MSL -> going up one rank")
    raise ValueError("Could not determine a valid query taxon from taxmyPHAGE output.")


def choose_outgroup(query_rank: str, query_value: str,
                    msl_df: pd.DataFrame) -> tuple[str, str]:
    """
    Choose an outgroup genus from a clade one level above the query rank.
    Returns (outgroup_rank, outgroup_value).

    Strategy:
      - Find all entries in the same parent clade as the query
      - Exclude entries that share the query value
      - Return the most represented genus (most species in MSL) in that set
    """
    parent_rank = OUTGROUP_RANK.get(query_rank, "Order")
    child_rank  = "Genus"

    # Find the parent clade value for the query
    mask = msl_df.get(query_rank, pd.Series(dtype=str)).str.lower().str.contains(
        query_value.lower(), regex=False, na=False
    )
    parent_val = ""
    if parent_rank in msl_df.columns and mask.any():
        vals = msl_df.loc[mask, parent_rank].dropna().unique()
        vals = [v for v in vals if v.strip()]
        parent_val = vals[0] if vals else ""

    if parent_val and parent_rank in msl_df.columns:
        # All genera within the same parent clade but OUTSIDE query taxon
        parent_mask = msl_df[parent_rank].str.lower().str.contains(
            parent_val.lower(), regex=False, na=False
        )
        exclude_mask = msl_df.get(query_rank, pd.Series(dtype=str)).str.lower().str.contains(
            query_value.lower(), regex=False, na=False
        )
        candidate_df = msl_df[parent_mask & ~exclude_mask]

        if not candidate_df.empty and child_rank in candidate_df.columns:
            top_genus = (
                candidate_df[child_rank]
                .replace("", pd.NA)
                .dropna()
                .value_counts()
                .index[0]
            )
            print(f"  Outgroup: '{top_genus}' [{child_rank}] "
                  f"(from {parent_rank} '{parent_val}', outside {query_rank} '{query_value}')")
            return child_rank, top_genus

    # Fallback to hardcoded genera
    for genus in FALLBACK_OUTGROUPS:
        if genus.lower() not in query_value.lower():
            print(f"  Outgroup (fallback): '{genus}' [Genus]")
            return "Genus", genus

    raise ValueError("Could not determine a suitable outgroup.")


# ─────────────────────────────────────────────────────────────────────────────
# Protein extraction from phold GBK
# ─────────────────────────────────────────────────────────────────────────────

def extract_query_proteins(gbk_path: str,
                            markers: list[str],
                            alias_map: dict[str, list[str]],
                            phage_name: str) -> dict[str, SeqRecord]:
    """
    Search phold.gbk for each marker protein and return one SeqRecord per marker.
    Uses the same alias matching as ictv_ncbi_fetcher.
    """
    found: dict[str, SeqRecord] = {}

    for rec in SeqIO.parse(gbk_path, "genbank"):
        for feat in rec.features:
            if feat.type != "CDS":
                continue
            translation = feat.qualifiers.get("translation", [""])[0]
            if not translation:
                continue
            for marker in markers:
                if marker in found:
                    continue
                if fetcher._feature_matches(feat, alias_map[marker]):
                    locus   = feat.qualifiers.get("locus_tag", ["query"])[0]
                    product = feat.qualifiers.get("product", [marker])[0]
                    sr = SeqRecord(
                        seq=Seq(translation),
                        id=f"QUERY|{locus}",
                        description=f"{phage_name} | {product}",
                    )
                    found[marker] = sr
                    print(f"  [query] {marker}: found -> {locus} ({product})")

    for marker in markers:
        if marker not in found:
            print(f"  [query] {marker}: NOT FOUND in {gbk_path}")

    return found


def append_query_to_fastas(query_proteins: dict[str, SeqRecord],
                            marker_fastas: dict[str, Path]) -> None:
    """Append each query protein to the corresponding marker FASTA."""
    for marker, sr in query_proteins.items():
        fasta_path = marker_fastas.get(marker)
        if fasta_path and fasta_path.exists():
            with open(fasta_path, "a") as fh:
                SeqIO.write(sr, fh, "fasta")
            print(f"  [append] {marker} -> {fasta_path}")
        else:
            print(f"  [WARNING] No FASTA to append {marker} to")


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────

def main() -> None:
    p = argparse.ArgumentParser(
        prog="phylo_prep.py",
        description=textwrap.dedent(__doc__),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--taxmyphage", required=True,
                   help="Path to taxmyPHAGE Summary_taxonomy.tsv")
    p.add_argument("--gbk",        required=True,
                   help="Path to phold.gbk (query phage annotation)")
    p.add_argument("--msl",        default="../results-paper-X1/ICTV_Master_Species_List_2025_MSL41.v1.xlsx",
                   help="Path to ICTV MSL .xlsx")
    p.add_argument("--markers",    default="terminase large subunit,major capsid protein",
                   help="Comma-separated marker protein names")
    p.add_argument("--name",       default="Query phage",
                   help="Display name for the query phage (used in FASTA headers)")
    p.add_argument("--email",      required=True,
                   help="E-mail for NCBI Entrez")
    p.add_argument("--api-key",    dest="api_key", default=None,
                   help="NCBI API key (optional)")
    p.add_argument("--outdir",     default="phylo_markers",
                   help="Output directory for marker FASTAs")
    p.add_argument("--min-species", dest="min_species", type=int, default=5,
                   help="Minimum species in MSL to accept a rank (default: 5)")
    p.add_argument("--no-outgroup", dest="no_outgroup", action="store_true",
                   help="Skip outgroup retrieval")
    args = p.parse_args()

    print("\n" + "="*60)
    print("  Phylogenetic Marker Prep")
    print("="*60)

    out_dir = Path(args.outdir)
    out_dir.mkdir(parents=True, exist_ok=True)

    markers = [m.strip() for m in args.markers.split(",") if m.strip()]

    # ── 1. Parse taxmyPHAGE ──────────────────────────────────────────────────
    print(f"\n[1] Parsing taxmyPHAGE: {args.taxmyphage}")
    tax = parse_taxmyphage(args.taxmyphage)
    print(f"  Class:  {tax['class_'] or '—'}")
    print(f"  Order:  {tax['order']  or '—'}")
    print(f"  Family: {tax['family'] or '—'}")
    print(f"  Genus:  {tax['genus']  or '—'}")
    print(f"  Status: {tax['ictv_status'] or '—'}")

    # ── 2. Load MSL ───────────────────────────────────────────────────────────
    print(f"\n[2] Loading MSL: {args.msl}")
    msl_df = fetcher.load_msl(args.msl)

    # ── 3. Choose query rank ──────────────────────────────────────────────────
    print(f"\n[3] Choosing taxon level (min {args.min_species} species in MSL)...")
    query_rank, query_value = choose_query_rank(tax, msl_df, args.min_species)
    print(f"  -> Querying at {query_rank} level: '{query_value}'")

    # ── 4. Filter MSL -> species list ─────────────────────────────────────────
    df_filtered = fetcher.filter_msl(msl_df, query_value, query_rank)
    fetcher.summarize_filtered(df_filtered)
    species_list = df_filtered["Species"].tolist()
    species_meta = df_filtered.set_index("Species").to_dict("index")

    # ── 5. Choose outgroup ────────────────────────────────────────────────────
    outgroup_value = None
    outgroup_rank  = None
    if not args.no_outgroup:
        print(f"\n[5] Choosing outgroup...")
        outgroup_rank, outgroup_value = choose_outgroup(query_rank, query_value, msl_df)
        print(f"  -> Outgroup: {outgroup_rank} '{outgroup_value}'")

    # ── 6. Setup fetcher ──────────────────────────────────────────────────────
    from Bio import Entrez
    Entrez.email = args.email
    Entrez.tool  = "phylo_prep"
    if args.api_key:
        Entrez.api_key = args.api_key
        fetcher.NCBI_DELAY = 0.12

    alias_map = fetcher.build_alias_map(markers, {})
    print("\n  Marker aliases:")
    for m, terms in alias_map.items():
        print(f"    {m!r}: {terms}")

    # ── 7. Search + download marker proteins ──────────────────────────────────
    print(f"\n[7] Searching NCBI protein accessions for {len(species_list)} species...")
    species_acc = fetcher.search_accessions(species_list, db="nuccore", seq_type="gene")

    prefix = str(out_dir / query_value.replace(" ", "_"))
    log: list[dict] = []
    print(f"\n[7b] Downloading marker sequences...")
    metadata = fetcher.fetch_genes_or_proteins(
        species_acc, markers, alias_map, species_meta,
        db="nuccore", seq_type="gene",
        out_prefix=prefix, log=log,
    )

    # Paths to each marker FASTA
    marker_fastas = {
        m: Path(f"{prefix}_{m.replace(' ', '_')}.faa")
        for m in markers
    }

    # ── 8. Download outgroup ──────────────────────────────────────────────────
    if outgroup_value and not args.no_outgroup:
        print(f"\n[8] Downloading outgroup '{outgroup_value}'...")
        og_log: list[dict] = []

        # Find one representative species from the outgroup genus
        og_df     = fetcher.filter_msl(msl_df, outgroup_value, outgroup_rank)
        og_species = og_df["Species"].tolist()
        og_meta   = og_df.set_index("Species").to_dict("index")

        if og_species:
            og_acc  = fetcher.search_accessions(og_species[:5], db="nuccore", seq_type="gene")
            og_prefix = str(out_dir / f"outgroup_{outgroup_value.replace(' ', '_')}")
            fetcher.fetch_genes_or_proteins(
                og_acc, markers, alias_map, og_meta,
                db="nuccore", seq_type="gene",
                out_prefix=og_prefix, log=og_log,
            )
            # Append outgroup sequences to the main marker FASTAs
            for m in markers:
                og_fasta = Path(f"{og_prefix}_{m.replace(' ', '_')}.faa")
                main_fasta = marker_fastas[m]
                if og_fasta.exists() and main_fasta.exists():
                    with open(main_fasta, "a") as fh_main:
                        for rec in SeqIO.parse(str(og_fasta), "fasta"):
                            rec.description = f"OUTGROUP|{outgroup_value} | {rec.description}"
                            SeqIO.write(rec, fh_main, "fasta")
                    og_fasta.unlink()
                    print(f"  [OK] Outgroup appended -> {main_fasta}")
        else:
            print(f"  [WARNING] No species found in MSL for outgroup '{outgroup_value}'")

    # ── 9. Extract and append query phage proteins ────────────────────────────
    print(f"\n[9] Extracting {args.name} proteins from {args.gbk}...")
    query_proteins = extract_query_proteins(args.gbk, markers, alias_map, args.name)
    append_query_to_fastas(query_proteins, marker_fastas)

    # ── 10. Write logs ────────────────────────────────────────────────────────
    fetcher.write_tsv(metadata, str(out_dir / "metadata.tsv"),  label="Metadata")
    fetcher.write_tsv(log,      str(out_dir / "not_found.tsv"), label="Not-found log")

    # ── Summary ───────────────────────────────────────────────────────────────
    print(f"\n{'='*60}")
    print("  SUMMARY")
    print(f"{'='*60}")
    print(f"  Query taxon     : {query_rank} '{query_value}'")
    print(f"  Species queried : {len(species_list)}")
    outgroup_str = f"{outgroup_rank} '{outgroup_value}'" if outgroup_value else "none"
    print(f"  Outgroup        : {outgroup_str}")
    print(f"  Markers         : {', '.join(markers)}")
    print()
    for m in markers:
        fasta = marker_fastas[m]
        if fasta.exists():
            n = sum(1 for _ in SeqIO.parse(str(fasta), "fasta"))
            print(f"  {m}: {n} sequences -> {fasta}")
        else:
            print(f"  {m}: no sequences written")
    print()
    print("  Next steps:")
    for m in markers:
        fasta = marker_fastas[m]
        stem  = str(out_dir / m.replace(" ", "_"))
        print(f"    mafft --localpair --maxiterate 1000 --thread N {fasta} > {stem}_aligned.faa")
        print(f"    trimal -in {stem}_aligned.faa -out {stem}_trimmed.faa -automated1")
        print(f"    iqtree2 -s {stem}_trimmed.faa -m MFP -T N --alrt 1000 --bb 1000 -pre {stem}_tree")
    print(f"{'='*60}\n")


if __name__ == "__main__":
    main()
