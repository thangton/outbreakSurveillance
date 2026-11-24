#!/usr/bin/env python3
from __future__ import annotations

import argparse
import math
from itertools import combinations
from pathlib import Path

from Bio import Phylo, SeqIO


def quantile_type7(values: list[int], q: float) -> float:
    vals = sorted(values)
    n = len(vals)
    if n == 1:
        return float(vals[0])
    h = (n - 1) * q + 1
    hf = math.floor(h)
    hc = math.ceil(h)
    if hf == hc:
        return float(vals[hf - 1])
    v1 = vals[hf - 1]
    v2 = vals[hc - 1]
    return float(v1 + (h - hf) * (v2 - v1))


def load_alignment(path: Path) -> dict[str, str]:
    seqs: dict[str, str] = {}
    for rec in SeqIO.parse(str(path), "fasta"):
        seqs[rec.id] = str(rec.seq).upper()
    if not seqs:
        raise ValueError(f"No sequences loaded from {path}")
    lengths = {len(v) for v in seqs.values()}
    if len(lengths) != 1:
        raise ValueError("Alignment contains unequal sequence lengths")
    return seqs


def pairwise_snp_table(seqs: dict[str, str]) -> list[tuple[str, str, int, int, float]]:
    valid = set("ACGT")
    rows: list[tuple[str, str, int, int, float]] = []
    names = list(seqs.keys())
    for a, b in combinations(names, 2):
        s1 = seqs[a]
        s2 = seqs[b]
        compared = 0
        snps = 0
        for x, y in zip(s1, s2):
            if x in valid and y in valid:
                compared += 1
                if x != y:
                    snps += 1
        p_dist = (snps / compared) if compared else float("nan")
        rows.append((a, b, snps, compared, p_dist))
    return rows


def build_components(
    names: list[str], rows: list[tuple[str, str, int, int, float]], threshold: int
) -> list[list[str]]:
    nonref = [n for n in names if n != "Reference"]
    idx = {n: i for i, n in enumerate(nonref)}
    parent = list(range(len(nonref)))

    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a: int, b: int) -> None:
        ra = find(a)
        rb = find(b)
        if ra != rb:
            parent[rb] = ra

    for a, b, snps, _, _ in rows:
        if a == "Reference" or b == "Reference":
            continue
        if snps <= threshold:
            union(idx[a], idx[b])

    comps: dict[int, list[str]] = {}
    for n in nonref:
        root = find(idx[n])
        comps.setdefault(root, []).append(n)
    return sorted((sorted(v) for v in comps.values()), key=lambda v: (-len(v), v))


def monophyly_checks(tree, components: list[list[str]]) -> list[tuple[int, bool, str, list[str]]]:
    tip_names = {t.name for t in tree.get_terminals() if t.name}
    out = []
    for comp in components:
        if len(comp) < 2:
            continue
        missing = [m for m in comp if m not in tip_names]
        if missing:
            out.append((len(comp), False, f"missing:{','.join(missing[:3])}", comp))
            continue
        anc = tree.common_ancestor(comp)
        descendants = {t.name for t in anc.get_terminals() if t.name}
        is_mono = descendants == set(comp)
        anc_name = anc.name if anc.name else "internal"
        out.append((len(comp), is_mono, anc_name, comp))
    return out


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Compute pairwise SNP distances from snippy core.aln and compare to tree clades."
    )
    parser.add_argument("--alignment", required=True, help="Path to snippy core.aln")
    parser.add_argument("--tree", required=True, help="Path to Newick tree")
    parser.add_argument("--out-dir", required=True, help="Output directory")
    args = parser.parse_args()

    aln_path = Path(args.alignment).resolve()
    tree_path = Path(args.tree).resolve()
    out_dir = Path(args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    seqs = load_alignment(aln_path)
    rows = pairwise_snp_table(seqs)
    names = list(seqs.keys())
    nonref_snps = [r[2] for r in rows if r[0] != "Reference" and r[1] != "Reference"]
    q1 = quantile_type7(nonref_snps, 0.25)
    q2 = quantile_type7(nonref_snps, 0.50)
    q3 = quantile_type7(nonref_snps, 0.75)
    threshold = int(math.floor(q1))
    components = build_components(names, rows, threshold)

    tree = Phylo.read(str(tree_path), "newick")
    checks = monophyly_checks(tree, components)

    pair_path = out_dir / "snippy_pairwise_snp.tsv"
    with pair_path.open("w") as w:
        w.write("sample_a\tsample_b\tsnp_count\tcompared_sites\tp_distance\n")
        for a, b, snps, compared, p_dist in sorted(rows, key=lambda r: (r[2], r[0], r[1])):
            w.write(f"{a}\t{b}\t{snps}\t{compared}\t{p_dist:.8f}\n")

    closest_path = out_dir / "snippy_closest_pairs_top30.tsv"
    with closest_path.open("w") as w:
        w.write("sample_a\tsample_b\tsnp_count\tcompared_sites\tp_distance\n")
        count = 0
        for a, b, snps, compared, p_dist in sorted(rows, key=lambda r: (r[2], r[0], r[1])):
            if a == "Reference" or b == "Reference":
                continue
            w.write(f"{a}\t{b}\t{snps}\t{compared}\t{p_dist:.8f}\n")
            count += 1
            if count >= 30:
                break

    cluster_path = out_dir / "snippy_clusters_q1.tsv"
    with cluster_path.open("w") as w:
        w.write("cluster_id\tsize\tmembers\n")
        for i, comp in enumerate(components, 1):
            w.write(f"C{i}\t{len(comp)}\t{';'.join(comp)}\n")

    summary_path = out_dir / "snippy_distance_summary.txt"
    with summary_path.open("w") as w:
        w.write(f"alignment: {aln_path}\n")
        w.write(f"tree: {tree_path}\n")
        w.write(f"samples_total: {len(names)} (including Reference)\n")
        w.write(f"pairwise_nonref_n: {len(nonref_snps)}\n")
        w.write(f"nonref_snp_min: {min(nonref_snps)}\n")
        w.write(f"nonref_snp_q1: {q1:.2f}\n")
        w.write(f"nonref_snp_median: {q2:.2f}\n")
        w.write(f"nonref_snp_q3: {q3:.2f}\n")
        w.write(f"nonref_snp_max: {max(nonref_snps)}\n")
        w.write(f"cluster_threshold_q1_floor: {threshold}\n")
        w.write(f"clusters_total: {len(components)}\n")
        w.write("\nmonophyly_check_for_multi_member_clusters:\n")
        for size, is_mono, anc, comp in checks:
            w.write(
                f"size={size}\tmonophyletic={is_mono}\tancestor={anc}\tmembers={';'.join(comp)}\n"
            )

    print(f"Wrote: {pair_path}")
    print(f"Wrote: {closest_path}")
    print(f"Wrote: {cluster_path}")
    print(f"Wrote: {summary_path}")


if __name__ == "__main__":
    main()
