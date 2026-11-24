#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Dict, List, Tuple

from Bio import SeqIO


def load_clusters(path: Path) -> List[Tuple[str, List[str]]]:
    clusters: List[Tuple[str, List[str]]] = []
    with path.open() as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            members = row["members"].split(";") if row.get("members") else []
            clusters.append((row["cluster_id"], members))
    return clusters


def load_core_tab(path: Path) -> Tuple[List[str], List[dict]]:
    with path.open() as f:
        reader = csv.DictReader(f, delimiter="\t")
        rows = list(reader)
        fields = reader.fieldnames or []
    fixed = {"CHR", "POS", "REF", "Reference", "LOCUS_TAG", "GENE", "PRODUCT", "EFFECT"}
    samples = [c for c in fields if c not in fixed]
    return samples, rows


def build_feature_map(genbank_path: Path) -> Dict[int, Dict[str, str]]:
    pos_to_annot: Dict[int, Dict[str, str]] = {}
    record = SeqIO.read(str(genbank_path), "genbank")
    preferred = {"CDS", "rRNA", "tRNA", "ncRNA"}
    secondary = {"gene"}

    def add_features(allowed: set[str], overwrite: bool = False) -> None:
        for feat in record.features:
            if feat.type not in allowed:
                continue
            parts = list(getattr(feat.location, "parts", [feat.location]))
            quals = feat.qualifiers
            gene = ";".join(quals.get("gene", []))
            locus = ";".join(quals.get("locus_tag", []))
            product = ";".join(quals.get("product", []))
            for part in parts:
                start = int(part.start) + 1
                end = int(part.end)
                for pos in range(start, end + 1):
                    if overwrite or pos not in pos_to_annot:
                        pos_to_annot[pos] = {
                            "locus_tag": locus,
                            "gene": gene,
                            "product": product,
                            "feature_type": feat.type,
                        }

    # First pass keeps richer functional annotations from CDS/functional RNA.
    add_features(preferred, overwrite=False)
    # Second pass fills gaps with generic gene features.
    add_features(secondary, overwrite=False)
    return pos_to_annot


def defining_sites(
    rows: List[dict],
    samples: List[str],
    members: List[str],
    min_outside_diff: float = 0.9,
) -> List[Tuple[str, int, str, str, float]]:
    member_set = [m for m in members if m in samples]
    outside = [s for s in samples if s not in member_set]
    out: List[Tuple[str, int, str, str, float]] = []
    if not member_set:
        return out

    missing = {"", "-", "N", "n"}
    for row in rows:
        member_alleles = [row[m] for m in member_set]
        if any(a in missing for a in member_alleles):
            continue
        if len(set(member_alleles)) != 1:
            continue
        cluster_allele = member_alleles[0]
        outside_alleles = [row[s] for s in outside if row[s] not in missing]
        if not outside_alleles:
            continue
        frac_diff = sum(1 for a in outside_alleles if a != cluster_allele) / len(outside_alleles)
        if frac_diff >= min_outside_diff:
            out.append(
                (
                    row["CHR"],
                    int(row["POS"]),
                    row["Reference"],
                    cluster_allele,
                    frac_diff,
                )
            )
    return out


def singleton_unique_count(rows: List[dict], samples: List[str], sample: str) -> int:
    missing = {"", "-", "N", "n"}
    unique = 0
    for row in rows:
        a = row[sample]
        if a in missing:
            continue
        others = [row[s] for s in samples if s != sample and row[s] not in missing]
        if others and all(o != a for o in others):
            unique += 1
    return unique


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Annotate cluster-defining SNP sites from Snippy core.tab using a GenBank reference."
    )
    parser.add_argument("--clusters", required=True, help="Path to cluster TSV (cluster_id,size,members)")
    parser.add_argument("--core-tab", required=True, help="Path to snippy core.tab")
    parser.add_argument("--genbank", required=True, help="Path to reference GenBank file")
    parser.add_argument("--out-prefix", required=True, help="Output prefix path")
    parser.add_argument("--min-outside-diff", type=float, default=0.9, help="Outside-difference threshold")
    args = parser.parse_args()

    clusters = load_clusters(Path(args.clusters))
    samples, rows = load_core_tab(Path(args.core_tab))
    annot = build_feature_map(Path(args.genbank))

    out_prefix = Path(args.out_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    tsv_path = out_prefix.with_suffix(".tsv")
    summary_path = out_prefix.with_name(out_prefix.name + "_summary.txt")

    # Detailed per-site table.
    with tsv_path.open("w") as w:
        w.write(
            "cluster_id\tcluster_size\tchr\tpos\tref\tcluster_allele\toutside_diff_frac\t"
            "locus_tag\tgene\tproduct\tfeature_type\n"
        )
        for cid, members in clusters:
            sites = defining_sites(rows, samples, members, args.min_outside_diff)
            for chr_, pos, ref, allele, frac in sites:
                a = annot.get(pos, {})
                w.write(
                    f"{cid}\t{len(members)}\t{chr_}\t{pos}\t{ref}\t{allele}\t{frac:.4f}\t"
                    f"{a.get('locus_tag','')}\t{a.get('gene','')}\t{a.get('product','')}\t"
                    f"{a.get('feature_type','')}\n"
                )

    # Human-readable summary.
    with summary_path.open("w") as w:
        w.write("Cluster mutation signature summary\n")
        w.write(f"min_outside_diff={args.min_outside_diff}\n\n")
        w.write("cluster_id\tsize\tdefining_sites\tunique_singleton_sites\n")
        for cid, members in clusters:
            sites = defining_sites(rows, samples, members, args.min_outside_diff)
            unique = singleton_unique_count(rows, samples, members[0]) if len(members) == 1 else ""
            w.write(f"{cid}\t{len(members)}\t{len(sites)}\t{unique}\n")

    print(f"Wrote: {tsv_path}")
    print(f"Wrote: {summary_path}")


if __name__ == "__main__":
    main()
