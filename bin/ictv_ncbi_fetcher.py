"""
ictv_ncbi_fetcher.py — Interactive pipeline for retrieving virus sequences from the ICTV MSL.

Steps:
  1. Filter the MSL by taxon (any hierarchical rank)
  2. Retrieve sequences from NCBI:
       a) Complete genome (nucleotide)
       b) Specific genes (nucleotide) — tolerant to annotation naming variation
       c) Specific proteins — same
  3. Define outgroup (NCBI taxon or user-supplied FASTA)
  4. Output: FASTA files + metadata TSV + not-found log

Dependencies: biopython, pandas, openpyxl, requests
Interactive:  python3 ictv_ncbi_fetcher.py
CLI:          python3 ictv_ncbi_fetcher.py --help
"""

import argparse
import csv
import re
import shutil
import sys
import time
from pathlib import Path
from typing import Optional

import pandas as pd
import requests
from Bio import Entrez, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# ─────────────────────────────────────────────────────────────────────────────
# Defaults — overridden by CLI args or interactive prompts
# ─────────────────────────────────────────────────────────────────────────────
DEFAULT_MSL    = "ICTV_Master_Species_List_2025_MSL41.v1.xlsx"
MSL_SHEET      = "MSL"
MSL_URL        = "https://ictv.global/vmr/current"  # auto-download fallback

TAXON_RANKS    = [
    "Realm", "Subrealm", "Kingdom", "Subkingdom",
    "Phylum", "Subphylum", "Class", "Subclass",
    "Order", "Suborder", "Family", "Subfamily",
    "Genus", "Subgenus", "Species",
]

RETRY_ATTEMPTS = 5
RETRY_WAIT     = 3    # seconds; doubles on each attempt
NCBI_DELAY     = 0.4  # seconds between calls without API key (≤3 req/s)
                      # set to 0.12 when an API key is supplied (≤10 req/s)

# ─────────────────────────────────────────────────────────────────────────────
# Utilities
# ─────────────────────────────────────────────────────────────────────────────

def prompt(msg: str, default: str = "") -> str:
    """Prompt with an optional default value shown in brackets."""
    hint = f" [{default}]" if default else ""
    val  = input(f"{msg}{hint}: ").strip()
    return val or default


def prompt_choice(msg: str, options: list[str]) -> str:
    """Numbered selection from a list of options."""
    print(f"\n{msg}")
    for i, opt in enumerate(options, 1):
        print(f"  {i}. {opt}")
    while True:
        raw = input("Choice (number): ").strip()
        if raw.isdigit() and 1 <= int(raw) <= len(options):
            return options[int(raw) - 1]
        print("  Invalid option, try again.")


def entrez_call(func, *args, **kwargs):
    """Call an Entrez function with exponential-backoff retry."""
    wait = RETRY_WAIT
    for attempt in range(1, RETRY_ATTEMPTS + 1):
        try:
            result = func(*args, **kwargs)
            time.sleep(NCBI_DELAY)
            return result
        except Exception as exc:
            if attempt == RETRY_ATTEMPTS:
                raise
            print(f"    [retry {attempt}/{RETRY_ATTEMPTS}] {exc} — waiting {wait}s...")
            time.sleep(wait)
            wait *= 2


# ─────────────────────────────────────────────────────────────────────────────
# Step 1 — MSL
# ─────────────────────────────────────────────────────────────────────────────

def _detect_msl_sheet(path: Path) -> str:
    """
    Return the correct sheet name for the species table.
    Tries MSL_SHEET first, then falls back to case-insensitive matching,
    then to the sheet whose first row contains 'Species' and 'Genus'.
    """
    xl = pd.ExcelFile(path)
    sheets = xl.sheet_names

    # Exact match
    if MSL_SHEET in sheets:
        return MSL_SHEET

    # Case-insensitive match
    for s in sheets:
        if s.strip().lower() == MSL_SHEET.lower():
            return s

    # Heuristic: find the sheet that looks like a taxonomy table
    for s in sheets:
        try:
            header = pd.read_excel(path, sheet_name=s, nrows=1, dtype=str)
            cols   = [c.strip() for c in header.columns]
            if "Species" in cols and "Genus" in cols:
                print(f"[MSL] Sheet '{MSL_SHEET}' not found; using '{s}' instead.")
                return s
        except Exception:
            continue

    raise ValueError(
        f"Could not find the MSL species sheet in {path}.\n"
        f"Available sheets: {sheets}"
    )


def load_msl(msl_path: str) -> pd.DataFrame:
    """Load the ICTV MSL. Downloads it automatically if the file is not found."""
    path = Path(msl_path)
    if not path.exists():
        print(f"[MSL] File not found: {path}")
        print(f"[MSL] Downloading from {MSL_URL}...")
        r = requests.get(MSL_URL, timeout=60)
        r.raise_for_status()
        path.write_bytes(r.content)
        print(f"[MSL] Saved to {path}")

    print(f"[MSL] Reading {path}...")
    sheet = _detect_msl_sheet(path)
    df    = pd.read_excel(path, sheet_name=sheet, dtype=str).fillna("")
    print(f"[MSL] {len(df):,} species loaded (sheet: '{sheet}').")
    return df


def filter_msl(df: pd.DataFrame, taxon: str, rank: Optional[str] = None) -> pd.DataFrame:
    """
    Filter the MSL for entries matching *taxon*.
    If *rank* is None, searches all taxonomic rank columns.
    Matching is case-insensitive and accepts partial strings.
    """
    taxon_lower = taxon.strip().lower()
    cols = [rank] if rank else TAXON_RANKS

    mask = pd.Series([False] * len(df), index=df.index)
    for col in cols:
        if col in df.columns:
            mask |= df[col].str.lower().str.contains(taxon_lower, regex=False, na=False)

    return df[mask].copy()


def summarize_filtered(df: pd.DataFrame) -> None:
    """Print taxon counts by rank for the filtered dataset."""
    print(f"\n{'─'*55}")
    print(f"  {len(df)} species found")
    print(f"{'─'*55}")
    for rank in ["Order", "Family", "Subfamily", "Genus"]:
        if rank in df.columns:
            counts = df[rank].replace("", pd.NA).dropna().value_counts()
            if not counts.empty:
                top = ", ".join(f"{k}({v})" for k, v in counts.head(8).items())
                print(f"  {rank:12s}: {top}")
    print(f"{'─'*55}\n")


# ─────────────────────────────────────────────────────────────────────────────
# Step 2 — NCBI sequence retrieval
# ─────────────────────────────────────────────────────────────────────────────

# Built-in synonyms for common phage marker proteins.
# Extend at runtime via --extra-aliases (CLI) or the interactive prompt.
BUILT_IN_ALIASES: dict[str, list[str]] = {
    "terminase large subunit": [
        "terminase large", "large terminase", "TerL", "major terminase",
        "large subunit terminase", "ATP-dependent terminase",
    ],
    "terminase small subunit": [
        "terminase small", "small terminase", "TerS",
    ],
    "major capsid protein": [
        "MCP", "capsid protein", "head protein", "coat protein",
        "major head protein", "major coat protein",
    ],
    "tail fiber protein": [
        "tail fiber", "fiber protein", "tail spike", "receptor binding protein",
        "tail spike protein", "tailspike",
    ],
    "RNA polymerase": [
        "RNAP", "phage RNA polymerase", "DNA-dependent RNA polymerase",
    ],
    "endolysin": [
        "lysin", "lytic enzyme", "peptidoglycan hydrolase", "muralytic enzyme",
    ],
    "baseplate": [
        "baseplate protein", "base plate",
    ],
}


def build_alias_map(targets: list[str],
                    extra_aliases: dict[str, list[str]]) -> dict[str, list[str]]:
    """
    For each requested target, return its full set of search terms
    (the target itself + any built-in or user-supplied aliases).
    """
    combined = {**BUILT_IN_ALIASES, **extra_aliases}
    alias_map: dict[str, list[str]] = {}

    for t in targets:
        t_lower = t.strip().lower()
        terms   = [t]
        for canonical, syns in combined.items():
            all_variants = [canonical.lower()] + [s.lower() for s in syns]
            if t_lower in all_variants or any(t_lower in v for v in all_variants):
                terms = list({canonical} | set(syns) | {t})
                break
        alias_map[t] = terms

    return alias_map


def _feature_matches(feat, terms: list[str], use_regex: bool = False) -> bool:
    """Return True if any of *terms* appears in the feature's annotation fields."""
    text_fields: list[str] = []
    for q in ("product", "gene", "note", "function", "standard_name"):
        text_fields.extend(feat.qualifiers.get(q, []))
    combined = " ".join(text_fields).lower()

    for term in terms:
        t = term.lower()
        if use_regex:
            if re.search(t, combined):
                return True
        else:
            if t in combined:
                return True
    return False


def search_accessions(species_list: list[str],
                      db: str,
                      seq_type: str) -> dict[str, list[str]]:
    """
    Search NCBI for accession IDs for each species.
    db       : 'nuccore' for genomes/genes, 'protein' for proteins.
    seq_type : 'genome' | 'gene' | 'protein' (refines the query string).
    Returns  : {species_name: [accID, ...]}
    """
    results: dict[str, list[str]] = {}
    total = len(species_list)

    for i, species in enumerate(species_list, 1):
        base_query = f'"{species}"[Organism]'
        query = f'{base_query} AND "complete genome"[Title]' if seq_type == "genome" \
                else base_query

        print(f"  [{i}/{total}] {species[:55]:<55}", end=" ", flush=True)

        try:
            handle = entrez_call(Entrez.esearch, db=db, term=query,
                                 retmax=20, idtype="acc")
            record = Entrez.read(handle)
            handle.close()
            ids = record.get("IdList", [])
            results[species] = ids
            print(f"-> {len(ids)} acc")
        except Exception as exc:
            results[species] = []
            print(f"-> ERROR: {exc}")

    return results


def fetch_genomes(species_acc: dict[str, list[str]],
                  species_meta: dict[str, dict],
                  out_fasta: str,
                  log: list[dict]) -> list[dict]:
    """
    Download one complete genome FASTA per species.
    Prefers RefSeq (NC_*) accessions when available.
    Returns a list of metadata dicts.
    """
    metadata: list[dict] = []
    fetched = 0
    total   = sum(1 for v in species_acc.values() if v)

    with open(out_fasta, "w") as fh:
        for species, acc_ids in species_acc.items():
            if not acc_ids:
                log.append({"species": species, "type": "genome",
                            "status": "no_accession"})
                continue

            # Prefer RefSeq
            acc_ids = sorted(acc_ids, key=lambda x: (0 if x.startswith("NC_") else 1))
            acc     = acc_ids[0]
            fetched += 1
            print(f"  [{fetched}/{total}] {acc}  ({species[:40]})")

            try:
                handle = entrez_call(Entrez.efetch, db="nuccore", id=acc,
                                     rettype="fasta", retmode="text")
                rec = SeqIO.read(handle, "fasta")
                handle.close()

                rec.id          = acc
                rec.description = f"{species} | {rec.description}"
                SeqIO.write(rec, fh, "fasta")

                meta = species_meta.get(species, {})
                metadata.append({
                    "accession": acc,
                    "species":   species,
                    **{k: meta.get(k, "") for k in ["Family", "Genus", "Order", "Genome"]},
                    "length":    len(rec.seq),
                    "type":      "genome",
                })
            except Exception as exc:
                log.append({"species": species, "accession": acc,
                            "type": "genome", "status": str(exc)})

    return metadata


def fetch_genes_or_proteins(species_acc: dict[str, list[str]],
                             targets: list[str],
                             alias_map: dict[str, list[str]],
                             species_meta: dict[str, dict],
                             db: str,
                             seq_type: str,
                             out_prefix: str,
                             log: list[dict],
                             use_regex: bool = False) -> list[dict]:
    """
    Download sequences for each species and extract entries matching the requested
    targets. Writes one FASTA file per target.

    - seq_type == 'gene'    : fetches nuccore GenBank, extracts CDS nucleotide sequences.
    - seq_type == 'protein' : fetches protein GenPept records directly (db='protein').
      Accessions returned by search_accessions for proteins are protein IDs, so they
      must be fetched from db='protein', NOT from nuccore.
    """
    metadata: list[dict] = []

    fasta_handles = {
        t: open(f"{out_prefix}_{t.replace(' ', '_')}.faa", "w")
        for t in targets
    }

    total = len(species_acc)
    for i, (species, acc_ids) in enumerate(species_acc.items(), 1):
        if not acc_ids:
            for t in targets:
                log.append({"species": species, "target": t,
                            "type": seq_type, "status": "no_accession"})
            continue

        print(f"  [{i}/{total}] {species[:50]}", end=" ", flush=True)

        try:
            if seq_type == "protein":
                # ── Protein mode ──────────────────────────────────────────────
                # Batch-fetch all protein accessions for this species at once.
                batch_ids = acc_ids[:20]
                handle    = entrez_call(Entrez.efetch, db="protein",
                                        id=",".join(batch_ids),
                                        rettype="gp", retmode="text")
                prot_records = list(SeqIO.parse(handle, "genbank"))
                handle.close()

                found_targets: set[str] = set()
                for prec in prot_records:
                    # Build a searchable text from description + feature qualifiers
                    search_text = prec.description.lower()
                    for feat in prec.features:
                        for q in ("product", "gene", "note", "function"):
                            search_text += " " + " ".join(feat.qualifiers.get(q, []))
                    search_text = search_text.lower()

                    for target in targets:
                        terms = alias_map[target]
                        matched = any(
                            (re.search(t.lower(), search_text) if use_regex
                             else t.lower() in search_text)
                            for t in terms
                        )
                        if not matched:
                            continue

                        seq_str = str(prec.seq)
                        product = prec.description[:100]

                        rec_out = SeqRecord(
                            seq=Seq(seq_str),
                            id=prec.id,
                            description=f"{species} | {product}",
                        )
                        SeqIO.write(rec_out, fasta_handles[target], "fasta")

                        meta = species_meta.get(species, {})
                        metadata.append({
                            "accession":     prec.id,
                            "locus":         prec.id,
                            "species":       species,
                            **{k: meta.get(k, "") for k in
                               ["Family", "Genus", "Order", "Genome"]},
                            "length":        len(seq_str),
                            "target":        target,
                            "product_found": product,
                            "type":          seq_type,
                        })
                        found_targets.add(target)

                hits = len(found_targets)
                print(f"-> {hits}/{len(targets)} targets found "
                      f"({len(prot_records)} proteins scanned)")

                for t in targets:
                    if t not in found_targets:
                        log.append({"species": species, "accession": batch_ids[0],
                                    "target": t, "type": seq_type,
                                    "status": "feature_not_found"})

            else:
                # ── Gene mode ─────────────────────────────────────────────────
                # Fetch the full nucleotide GenBank and extract CDS features.
                acc_ids_sorted = sorted(acc_ids,
                                        key=lambda x: (0 if x.startswith("NC_") else 1))
                acc    = acc_ids_sorted[0]
                handle = entrez_call(Entrez.efetch, db="nuccore", id=acc,
                                     rettype="gb", retmode="text")
                record = SeqIO.read(handle, "genbank")
                handle.close()

                found_targets = set()
                for feat in record.features:
                    if feat.type not in ("CDS", "gene", "mat_peptide"):
                        continue
                    for target in targets:
                        if not _feature_matches(feat, alias_map[target], use_regex):
                            continue

                        # Extract protein translation, not nucleotide — .faa output
                        # is used for protein-based phylogenetics (MAFFT + IQ-TREE2)
                        seq_str = feat.qualifiers.get("translation", [""])[0]
                        if not seq_str:
                            continue  # skip pseudogenes / features without translation
                        locus   = feat.qualifiers.get("locus_tag",
                                  feat.qualifiers.get("gene", [acc]))[0]
                        product = feat.qualifiers.get("product", [target])[0]

                        rec_out = SeqRecord(
                            seq=Seq(seq_str),
                            id=f"{acc}|{locus}",
                            description=f"{species} | {product}",
                        )
                        SeqIO.write(rec_out, fasta_handles[target], "fasta")

                        meta = species_meta.get(species, {})
                        metadata.append({
                            "accession":     acc,
                            "locus":         locus,
                            "species":       species,
                            **{k: meta.get(k, "") for k in
                               ["Family", "Genus", "Order", "Genome"]},
                            "length":        len(seq_str),
                            "target":        target,
                            "product_found": product,
                            "type":          seq_type,
                        })
                        found_targets.add(target)

                print(f"-> {len(found_targets)}/{len(targets)} targets found  [{acc}]")

                for t in targets:
                    if t not in found_targets:
                        log.append({"species": species, "accession": acc,
                                    "target": t, "type": seq_type,
                                    "status": "feature_not_found"})

        except Exception as exc:
            print(f"-> ERROR: {exc}")
            for t in targets:
                log.append({"species": species, "accession": acc_ids[0] if acc_ids else "?",
                            "target": t, "type": seq_type, "status": str(exc)})

    for fh in fasta_handles.values():
        fh.close()

    return metadata


# ─────────────────────────────────────────────────────────────────────────────
# Step 3 — Outgroup
# ─────────────────────────────────────────────────────────────────────────────

def fetch_outgroup(taxon: str,
                   seq_type: str,
                   targets: list[str],
                   alias_map: dict[str, list[str]],
                   out_fasta: str,
                   log: list[dict]) -> None:
    """
    Retrieve one representative sequence for the outgroup taxon from NCBI.
    Prefers RefSeq reference genomes.
    """
    db    = "nuccore"
    query = f'"{taxon}"[Organism] AND "complete genome"[Title]'
    print(f"  Searching outgroup: {query}")

    try:
        handle = entrez_call(Entrez.esearch, db=db, term=query, retmax=5)
        record = Entrez.read(handle)
        handle.close()
        ids = record.get("IdList", [])

        if not ids:
            print(f"  [WARNING] No results for outgroup '{taxon}'")
            log.append({"species": taxon, "type": "outgroup", "status": "no_result"})
            return

        ids = sorted(ids, key=lambda x: (0 if "NC_" in x else 1))
        acc = ids[0]
        print(f"  Downloading outgroup: {acc}")

        if seq_type == "genome":
            handle = entrez_call(Entrez.efetch, db=db, id=acc,
                                 rettype="fasta", retmode="text")
            with open(out_fasta, "w") as fh:
                fh.write(handle.read())
            handle.close()
        else:
            handle = entrez_call(Entrez.efetch, db="nuccore", id=acc,
                                 rettype="gb", retmode="text")
            gbk_rec = SeqIO.read(handle, "genbank")
            handle.close()
            with open(out_fasta, "w") as fh:
                for feat in gbk_rec.features:
                    if feat.type not in ("CDS", "gene"):
                        continue
                    for target in targets:
                        if _feature_matches(feat, alias_map[target]):
                            seq_str = (str(feat.extract(gbk_rec.seq))
                                       if seq_type == "gene"
                                       else feat.qualifiers.get("translation", [""])[0])
                            if seq_str:
                                fh.write(f">{acc}|outgroup|{taxon}\n{seq_str}\n")
                                break

        print(f"  [OK] Outgroup saved to {out_fasta}")

    except Exception as exc:
        print(f"  [ERROR] Outgroup: {exc}")
        log.append({"species": taxon, "type": "outgroup", "status": str(exc)})


# ─────────────────────────────────────────────────────────────────────────────
# Step 4 — Output helpers
# ─────────────────────────────────────────────────────────────────────────────

def write_tsv(rows: list[dict], path: str, label: str = "rows") -> None:
    if not rows:
        return
    fieldnames = list({k for r in rows for k in r})
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames,
                                delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    print(f"[OK] {label} -> {path}  ({len(rows)} rows)")


# ─────────────────────────────────────────────────────────────────────────────
# Interactive mode
# ─────────────────────────────────────────────────────────────────────────────

def run(args) -> None:
    """Unified entry point for both interactive and CLI modes."""
    print("\n" + "═"*60)
    print("  ICTV → NCBI Sequence Fetcher")
    print("═"*60)

    # ── NCBI credentials ─────────────────────────────────────────
    email   = args.email   or prompt("NCBI Entrez e-mail")
    api_key = args.api_key or prompt("NCBI API key (optional, press Enter to skip)", "")
    Entrez.email = email
    Entrez.tool  = "ictv_ncbi_fetcher"
    if api_key:
        Entrez.api_key = api_key
        global NCBI_DELAY
        NCBI_DELAY = 0.12

    # ── MSL ──────────────────────────────────────────────────────
    if args.msl:
        msl_path = args.msl
    elif getattr(args, "_interactive", False):
        msl_path = prompt("MSL file path (.xlsx)", DEFAULT_MSL)
    else:
        msl_path = DEFAULT_MSL
    df_msl = load_msl(msl_path)

    # ── Step 1: taxon filter ──────────────────────────────────────
    taxon = args.taxon or prompt("Taxon of interest (e.g. Pradovirus, Autographiviridae)")

    rank = args.rank
    if not rank:
        choice = prompt_choice(
            "Which rank to search? (select 'All ranks' to search everywhere)",
            ["All ranks"] + TAXON_RANKS,
        )
        rank = None if choice == "All ranks" else choice

    df = filter_msl(df_msl, taxon, rank)

    if df.empty:
        print(f"\n[ERROR] No species found for '{taxon}'.")
        sys.exit(1)

    summarize_filtered(df)

    if prompt(f"Continue with {len(df)} species? (y/n)", "y").lower() != "y":
        sys.exit(0)

    species_list = df["Species"].tolist()
    species_meta = df.set_index("Species").to_dict("index")

    # ── Step 2: sequence type ─────────────────────────────────────
    seq_type = args.seq_type or prompt_choice(
        "Type of data to retrieve:",
        ["genome", "gene", "protein"],
    )

    targets:       list[str]        = []
    alias_map:     dict[str, list]  = {}
    extra_aliases: dict[str, list]  = {}
    use_regex                       = False

    if seq_type in ("gene", "protein"):
        raw = args.targets or prompt(
            "Sequences of interest (comma-separated)\n"
            "  e.g.: terminase large subunit, major capsid protein"
        )
        targets = [t.strip() for t in raw.split(",") if t.strip()]

        extra_raw = args.extra_aliases or prompt(
            "Extra aliases as 'target:alias1;alias2,...' (Enter to skip)", ""
        )
        if extra_raw:
            for pair in extra_raw.split(","):
                if ":" in pair:
                    key, vals = pair.split(":", 1)
                    extra_aliases[key.strip()] = [v.strip() for v in vals.split(";")]

        alias_map = build_alias_map(targets, extra_aliases)
        use_regex = prompt("Use regex for matching? (y/n)", "n").lower() == "y"

        print("\nSearch terms per target:")
        for t, terms in alias_map.items():
            print(f"  {t!r:40s} -> {terms}")

    # ── Output directory ──────────────────────────────────────────
    out_dir = Path(args.outdir or prompt("Output directory",
                                         f"ncbi_{taxon.replace(' ', '_')}"))
    out_dir.mkdir(parents=True, exist_ok=True)

    log:          list[dict] = []
    all_metadata: list[dict] = []

    # ── Search accessions ─────────────────────────────────────────
    print(f"\n{'─'*55}")
    print(f"  Searching NCBI accessions for {len(species_list)} species...")
    print(f"{'─'*55}")
    db          = "nuccore" if seq_type in ("genome", "gene") else "protein"
    species_acc = search_accessions(species_list, db, seq_type)

    # ── Download sequences ────────────────────────────────────────
    print(f"\n{'─'*55}")
    print("  Downloading sequences...")
    print(f"{'─'*55}")

    if seq_type == "genome":
        out_fasta = str(out_dir / f"{taxon.replace(' ', '_')}_genomes.fna")
        all_metadata.extend(
            fetch_genomes(species_acc, species_meta, out_fasta, log)
        )
        print(f"[OK] Genomes -> {out_fasta}")

    else:
        prefix = str(out_dir / taxon.replace(" ", "_"))
        meta   = fetch_genes_or_proteins(
            species_acc, targets, alias_map, species_meta,
            db, seq_type, prefix, log, use_regex,
        )
        all_metadata.extend(meta)
        for t in targets:
            n = sum(1 for m in meta if m.get("target") == t)
            print(f"[OK] {t!r}: {n} sequences -> {prefix}_{t.replace(' ', '_')}.faa")

    # ── Step 3: outgroup ──────────────────────────────────────────
    print(f"\n{'─'*55}")
    og_choice = prompt_choice(
        "Outgroup:",
        ["NCBI taxon (automatic)", "User-supplied FASTA", "No outgroup"],
    )

    if og_choice == "NCBI taxon (automatic)":
        og_taxon = args.outgroup or prompt("Outgroup taxon name")
        fetch_outgroup(og_taxon, seq_type, targets, alias_map,
                       str(out_dir / "outgroup.fna"), log)

    elif og_choice == "User-supplied FASTA":
        og_src = prompt("Path to outgroup FASTA")
        shutil.copy(og_src, out_dir / "outgroup.fasta")
        print(f"[OK] Outgroup copied to {out_dir / 'outgroup.fasta'}")

    # ── Write metadata and log ────────────────────────────────────
    write_tsv(all_metadata, str(out_dir / "metadata.tsv"),   label="Metadata")
    write_tsv(log,          str(out_dir / "not_found.tsv"),  label="Not-found log")

    # ── Summary ───────────────────────────────────────────────────
    found = sum(1 for v in species_acc.values() if v)
    print(f"\n{'═'*60}")
    print("  SUMMARY")
    print(f"{'═'*60}")
    print(f"  Species in taxon          : {len(species_list)}")
    print(f"  Species with accession    : {found}")
    print(f"  Sequences downloaded      : {len(all_metadata)}")
    print(f"  Not-found / error entries : {len(log)}")
    print(f"  Output directory          : {out_dir}/")
    print(f"{'═'*60}\n")


# ─────────────────────────────────────────────────────────────────────────────
# CLI
# ─────────────────────────────────────────────────────────────────────────────

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="ictv_ncbi_fetcher.py",
        description="Retrieve virus sequences from NCBI based on the ICTV MSL.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Interactive mode (recommended for first use)
  python3 ictv_ncbi_fetcher.py

  # Complete genomes for Pradovirus
  python3 ictv_ncbi_fetcher.py \\
      --taxon Pradovirus --rank Genus --seq-type genome \\
      --email user@email.com --outdir out_pradovirus

  # Marker proteins for Autographiviridae with an outgroup
  python3 ictv_ncbi_fetcher.py \\
      --taxon Autographiviridae --rank Family --seq-type protein \\
      --targets "terminase large subunit,major capsid protein,tail fiber protein" \\
      --outgroup T7virus \\
      --email user@email.com --outdir out_autographiviridae
        """,
    )
    p.add_argument("--msl",           default=None,
                   help=f"Path to MSL .xlsx file (default: {DEFAULT_MSL})")
    p.add_argument("--taxon",         default=None,
                   help="Taxon of interest")
    p.add_argument("--rank",          default=None, choices=TAXON_RANKS,
                   help="Taxonomic rank to search (default: all ranks)")
    p.add_argument("--seq-type",      default=None, dest="seq_type",
                   choices=["genome", "gene", "protein"],
                   help="Type of sequence to retrieve")
    p.add_argument("--targets",       default=None,
                   help="Comma-separated gene/protein names (required for gene/protein mode)")
    p.add_argument("--extra-aliases", default=None, dest="extra_aliases",
                   help="Extra aliases: 'target:alias1;alias2,...' (comma-separated pairs)")
    p.add_argument("--outgroup",      default=None,
                   help="Outgroup taxon name (retrieved automatically from NCBI)")
    p.add_argument("--outdir",        default=None,
                   help="Output directory")
    p.add_argument("--email",         default=None,
                   help="E-mail for NCBI Entrez (required)")
    p.add_argument("--api-key",       default=None, dest="api_key",
                   help="NCBI API key (optional; increases rate limit to 10 req/s)")
    return p


def main() -> None:
    parser = build_parser()

    # No arguments → interactive mode
    if len(sys.argv) == 1:
        args = parser.parse_args([])
        args._interactive = True
        run(args)
        return

    args = parser.parse_args()
    args._interactive = False

    if not args.email:
        parser.error("--email is required in CLI mode")
    if not args.taxon:
        parser.error("--taxon is required in CLI mode")
    if not args.seq_type:
        parser.error("--seq-type is required in CLI mode")
    if args.seq_type in ("gene", "protein") and not args.targets:
        parser.error("--targets is required when --seq-type is gene or protein")

    run(args)


if __name__ == "__main__":
    main()
