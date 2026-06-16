#!/usr/bin/env python3
"""
clinker_phrogs_colors.py — Recolor Clinker synteny HTML by PHROGS functional category.

The Clinker HTML contains an embedded JSON dataset with gene groups (homologous clusters).
This script:
  1. Parses the embedded JSON from the Clinker HTML
  2. Reads the source GBK files to extract PHROGS functional categories per gene
  3. For each group, assigns the color of the majority PHROGS category
  4. Injects a fixed legend overlay into the HTML
  5. Writes the modified HTML

Usage:
  python3 clinker_phrogs_colors.py \\
      --clinker_html figures/clinker_synteny.html \\
      --gbk_dir     comparative/clinker_genomes/ \\
      --output      figures/clinker_phrogs.html
"""

import argparse
import json
import re
from collections import Counter
from pathlib import Path

from Bio import SeqIO

# Colorblind-safe palette (Okabe-Ito + Wong, matching genome_map.py)
PHROGS_COLORS: dict[str, str] = {
    "head and packaging":                                "#0072B2",
    "tail":                                              "#009E73",
    "DNA, RNA and nucleotide metabolism":                "#E69F00",
    "lysis":                                             "#D55E00",
    "connector":                                         "#CC79A7",
    "transcription regulation":                          "#56B4E9",
    "moron, auxiliary metabolic gene and host takeover": "#F0E442",
    "other":                                             "#999999",
    "unknown function":                                  "#DDDDDD",
}

_LEGEND_ITEM = (
    '<div style="display:flex;align-items:center;margin:2px 0;">'
    '<span style="display:inline-block;width:14px;height:14px;'
    'background:{color};border-radius:2px;border:1px solid #aaa;'
    'margin-right:7px;flex-shrink:0;"></span>'
    '<span style="font-size:11px;">{label}</span></div>'
)

_LEGEND_LABELS = {
    "moron, auxiliary metabolic gene and host takeover": "Moron / AMG / host takeover",
    "DNA, RNA and nucleotide metabolism": "DNA/RNA metabolism",
}


def _legend_html() -> str:
    items = "".join(
        _LEGEND_ITEM.format(
            color=color,
            label=_LEGEND_LABELS.get(cat, cat.capitalize()),
        )
        for cat, color in PHROGS_COLORS.items()
    )
    return (
        '<div id="phrogs-legend" style="'
        "position:fixed;bottom:20px;right:20px;z-index:9999;"
        "background:rgba(255,255,255,0.95);border:1px solid #ccc;"
        "border-radius:6px;padding:10px 14px;"
        "font-family:Arial,sans-serif;"
        'box-shadow:0 2px 8px rgba(0,0,0,0.15);max-width:270px;">'
        '<div style="font-weight:bold;margin-bottom:6px;font-size:12px;">'
        "PHROG functional category</div>"
        f"{items}</div>"
    )


# ─────────────────────────────────────────────────────────────────────────────
# GBK → function map
# ─────────────────────────────────────────────────────────────────────────────

def load_gbk_functions(gbk_dir: str) -> dict[str, str]:
    """
    Read all GBK files and return { locus_tag | product_lower -> phrog_category }.
    Pharokka/Phold store the PHROGS category in the 'function' qualifier.
    """
    func_map: dict[str, str] = {}
    gbk_path = Path(gbk_dir)

    for gbk_file in sorted(gbk_path.glob("*.gb*")):
        try:
            for rec in SeqIO.parse(str(gbk_file), "genbank"):
                for feat in rec.features:
                    if feat.type != "CDS":
                        continue
                    fn = feat.qualifiers.get("function", ["unknown function"])[0].strip().lower()
                    if not fn:
                        fn = "unknown function"
                    for qual in ("locus_tag", "gene"):
                        for key in feat.qualifiers.get(qual, []):
                            if key.strip():
                                func_map[key.strip()] = fn
                    for prod in feat.qualifiers.get("product", []):
                        if prod.strip():
                            func_map[prod.strip().lower()] = fn
        except Exception as exc:
            print(f"  [WARN] Could not parse {gbk_file.name}: {exc}")

    print(f"  [gbk] {len(func_map)} gene→function entries loaded from {gbk_dir}")
    return func_map


# ─────────────────────────────────────────────────────────────────────────────
# Clinker HTML parsing
# ─────────────────────────────────────────────────────────────────────────────

def _find_json_object(html: str, search_from: int) -> tuple[int, int]:
    """
    Starting from search_from, locate the next '{' and return (start, end)
    indices of the complete JSON object using bracket counting.
    Correctly handles nested objects and quoted strings.
    """
    start = html.find('{', search_from)
    if start == -1:
        raise ValueError("No JSON object found after position %d" % search_from)

    depth     = 0
    in_string = False
    escape    = False

    for i in range(start, len(html)):
        ch = html[i]
        if escape:
            escape = False
            continue
        if ch == '\\' and in_string:
            escape = True
            continue
        if ch == '"':
            in_string = not in_string
            continue
        if in_string:
            continue
        if ch == '{':
            depth += 1
        elif ch == '}':
            depth -= 1
            if depth == 0:
                return start, i + 1

    raise ValueError("Unmatched '{' in HTML — JSON object is not closed")


def extract_clinker_json(html: str) -> tuple[str, dict]:
    """
    Extract the embedded JSON data from a Clinker HTML file.
    Returns (prefix_plus_json_string, parsed_dict) where the first element
    is the exact substring in *html* that will be replaced later.

    Uses bracket counting instead of regex to correctly handle deeply nested
    Clinker JSON objects (regex with .*? stops at the first '}' found inside
    a nested object and breaks json.loads).
    """
    keyword_patterns = [
        r'(?:const|var|let)\s+data\s*=\s*',
        r'plot\(',
    ]
    for kw_pat in keyword_patterns:
        m = re.search(kw_pat, html)
        if not m:
            continue
        try:
            json_start, json_end = _find_json_object(html, m.end())
            data = json.loads(html[json_start:json_end])
            # original_str spans from keyword start to end of JSON object
            original_str = html[m.start():json_end]
            return original_str, data
        except (ValueError, json.JSONDecodeError):
            continue

    raise ValueError(
        "Could not find embedded JSON data in Clinker HTML. "
        "Ensure this is a valid clinker output file."
    )


# ─────────────────────────────────────────────────────────────────────────────
# Gene UID → PHROGS function mapping
# ─────────────────────────────────────────────────────────────────────────────

def _normalize_category(raw: str) -> str:
    """Map a raw function string to a known PHROGS_COLORS key."""
    raw = raw.strip().lower()
    for cat in PHROGS_COLORS:
        if cat in raw or raw in cat:
            return cat
    return "unknown function"


def build_uid_func_map(data: dict, func_map: dict[str, str]) -> dict[str, str]:
    """
    For every gene in the Clinker JSON, resolve its PHROGS category.
    Uses the gene label (which Clinker sets to locus_tag or product) as lookup key.
    """
    uid_to_func: dict[str, str] = {}

    for locus in data.get("loci", []):
        for gene in locus.get("genes", []):
            uid   = gene.get("uid", "")
            label = gene.get("label", "").strip()

            fn = (
                func_map.get(label)
                or func_map.get(label.lower())
                or next(
                    (v for k, v in func_map.items()
                     if k and label and (k in label.lower() or label.lower() in k)),
                    "unknown function",
                )
            )
            uid_to_func[uid] = _normalize_category(fn)

    return uid_to_func


# ─────────────────────────────────────────────────────────────────────────────
# Recolor groups
# ─────────────────────────────────────────────────────────────────────────────

def _dominant_category(uids: list[str], uid_to_func: dict[str, str]) -> str:
    funcs = [uid_to_func.get(u, "unknown function") for u in uids]
    known = [f for f in funcs if f != "unknown function"]
    if known:
        return Counter(known).most_common(1)[0][0]
    return "unknown function"


def recolor_groups(data: dict, uid_to_func: dict[str, str]) -> int:
    """
    Assign PHROGS colors to each homologous group in-place.
    Returns the number of groups recolored.
    """
    count = 0
    for group in data.get("groups", []):
        uids: list[str] = []

        # Clinker >= 0.0.26 stores per-locus gene lists inside group
        genes_field = group.get("genes", {})
        if isinstance(genes_field, dict):
            for gene_list in genes_field.values():
                if isinstance(gene_list, list):
                    uids.extend(gene_list)
        elif isinstance(genes_field, list):
            uids.extend(genes_field)

        # Older format: group["members"] is a list of {gene: uid} dicts
        for member in group.get("members", []):
            if isinstance(member, dict):
                uids.append(member.get("gene", ""))
            elif isinstance(member, str):
                uids.append(member)

        if uids:
            cat   = _dominant_category(uids, uid_to_func)
            color = PHROGS_COLORS.get(cat, PHROGS_COLORS["unknown function"])
            group["colour"] = color
            count += 1

    return count


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(
        prog="clinker_phrogs_colors.py",
        description="Recolor Clinker HTML output by PHROGS functional category.",
    )
    parser.add_argument("--clinker_html", required=True,
                        help="Input Clinker HTML file (clinker_synteny.html)")
    parser.add_argument("--gbk_dir",     default=None,
                        help="Directory with the GBK files used as Clinker input "
                             "(enables accurate PHROGS category lookup)")
    parser.add_argument("--output",      required=True,
                        help="Output recolored HTML file")
    args = parser.parse_args()

    print("\n" + "="*60)
    print("  Clinker PHROGS Recoloring")
    print("="*60)

    html_path = Path(args.clinker_html)
    if not html_path.exists():
        raise FileNotFoundError(f"Clinker HTML not found: {html_path}")

    html = html_path.read_text(encoding="utf-8")
    print(f"  Input : {html_path}  ({len(html):,} bytes)")

    # Load gene→function from GBK files
    func_map: dict[str, str] = {}
    if args.gbk_dir and Path(args.gbk_dir).exists():
        func_map = load_gbk_functions(args.gbk_dir)
    else:
        print("  [WARN] --gbk_dir not provided; colors based on gene labels only")

    # Extract and parse embedded JSON
    print("\n  Extracting Clinker JSON data...")
    original_str, data = extract_clinker_json(html)

    # Map gene UIDs to PHROGS categories
    uid_to_func = build_uid_func_map(data, func_map)
    n_known = sum(1 for v in uid_to_func.values() if v != "unknown function")
    print(f"  Gene assignments: {n_known}/{len(uid_to_func)} with known PHROGS category")

    # Recolor
    print("\n  Recoloring gene groups...")
    n_recolored = recolor_groups(data, uid_to_func)
    print(f"  {n_recolored} groups recolored")

    # Rebuild HTML: replace the JSON object inside original_str with new JSON
    new_json      = json.dumps(data, separators=(",", ":"))
    json_offset   = original_str.index('{')
    new_matched   = original_str[:json_offset] + new_json
    html_out      = html.replace(original_str, new_matched, 1)

    # Inject legend before </body>
    legend = _legend_html()
    if "</body>" in html_out:
        html_out = html_out.replace("</body>", f"{legend}\n</body>", 1)
    else:
        html_out += legend

    # Write output
    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(html_out, encoding="utf-8")

    print(f"\n  Output: {out_path}  ({len(html_out):,} bytes)")
    print("  Done.\n")


if __name__ == "__main__":
    main()
