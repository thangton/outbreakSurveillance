#!/usr/bin/env python3
"""Join per-stage outputs into one row per sample, plus an AMR-by-SNP-cluster table.

Any input that is missing or empty (stage skipped) just leaves its columns blank.
"""
import argparse
import csv
import re
from collections import Counter, defaultdict
from pathlib import Path


def rows(path, header=True, fieldnames=None):
    p = Path(path)
    if not p.is_file() or p.stat().st_size == 0:
        return []
    with p.open() as fh:
        lines = (l for l in fh if not l.startswith("##"))
        if header:
            return list(csv.DictReader(lines, delimiter="\t"))
        return list(csv.DictReader(lines, delimiter="\t", fieldnames=fieldnames))


def sid(name):
    """D18eb3e7d6.contigs.fa / D18eb3e7d6.contigs / path/D18eb3e7d6.faa -> D18eb3e7d6"""
    return re.sub(r"\.(contigs|faa|fa|fna|gff|emapper\.annotations|16S).*$", "", Path(name).name)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--indir", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--amr-by-cluster", required=True)
    a = ap.parse_args()
    d = Path(a.indir)
    table = defaultdict(dict)

    # D1 – assembly stats
    for r in rows(d / "quast" / "transposed_report.tsv"):
        t = table[sid(r["Assembly"])]
        t.update(contigs=r.get("# contigs"), total_length=r.get("Total length"),
                 N50=r.get("N50"), GC=r.get("GC (%)"))

    # D2 – gene / annotation counts
    for f in sorted((d / "gff").glob("*.gff")):
        table[sid(f.name)]["predicted_genes"] = sum(
            1 for l in f.open() if not l.startswith("#") and "\tCDS\t" in l)
    for f in sorted((d / "eggnog").glob("*.emapper.annotations")):
        table[sid(f.name)]["eggnog_annotated"] = sum(
            1 for l in f.open() if l.strip() and not l.startswith("#"))

    # D3 – typing / taxonomy / QC
    for r in rows(d / "mlst.tsv", header=False, fieldnames=["file", "scheme", "ST"]):
        table[sid(r["file"])].update(mlst_scheme=r["scheme"], ST=r["ST"])

    seen = set()
    for r in rows(d / "blast16s.tsv", header=False,
                  fieldnames=["qseqid", "evalue", "bitscore", "sscinames"]):
        if r["qseqid"] in seen:
            continue
        seen.add(r["qseqid"])
        table[r["qseqid"]].update(genus_16S=(r["sscinames"] or "").split(" ")[0],
                                  bitscore_16S=r["bitscore"])

    best = {}
    for r in rows(d / "skani.tsv"):
        q = sid(r["Query_file"])
        if q not in best or float(r["ANI"]) > float(best[q]["ANI"]):
            best[q] = r
    for q, r in best.items():
        words = r["Ref_name"].split()
        table[q].update(species_skani=" ".join(words[1:3]), ANI=r["ANI"],
                        skani_ref=words[0] if words else "")

    for r in rows(d / "checkm.tsv"):
        table[sid(r["Bin Id"])].update(completeness=r.get("Completeness"),
                                       contamination=r.get("Contamination"))

    for r in rows(d / "gunc.tsv"):
        table[sid(r["genome"])].update(gunc_CSS=r.get("clade_separation_score"),
                                       gunc_pass=r.get("pass.GUNC"))

    # D4 – SNP clusters + AMR
    cluster_of, cluster_size = {}, {}
    for r in rows(d / "clusters.tsv"):
        members = [m for m in r["members"].split(";") if m and m != "Reference"]
        for m in members:
            cluster_of[m] = r["cluster_id"]
        cluster_size[r["cluster_id"]] = len(members)
    for s, c in cluster_of.items():
        table[s]["snp_cluster"] = c

    amr_genes = defaultdict(list)
    for f in sorted((d / "amr").glob("*_amr.tsv")):
        s = f.name[: -len("_amr.tsv")]
        for r in rows(f):
            if r.get("Element type", "AMR") == "AMR":
                amr_genes[s].append(r["Gene symbol"])
        table[s]["amr_n"] = len(amr_genes[s])
        table[s]["amr_genes"] = ",".join(sorted(set(amr_genes[s])))

    cols = ["sample", "contigs", "total_length", "N50", "GC", "predicted_genes", "eggnog_annotated",
            "mlst_scheme", "ST", "genus_16S", "bitscore_16S", "species_skani", "ANI", "skani_ref",
            "completeness", "contamination", "gunc_CSS", "gunc_pass", "snp_cluster", "amr_n", "amr_genes"]
    with open(a.out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(cols)
        for s in sorted(table):
            w.writerow([s] + [table[s].get(c, "") for c in cols[1:]])

    # AMR determinant prevalence within each SNP cluster (D4 public-health summary)
    counts = Counter()
    for s, genes in amr_genes.items():
        for g in set(genes):
            counts[(cluster_of.get(s, "NA"), g)] += 1
    with open(a.amr_by_cluster, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["snp_cluster", "cluster_size", "gene", "n_samples"])
        for (c, g), n in sorted(counts.items(), key=lambda x: (x[0][0], -x[1], x[0][1])):
            w.writerow([c, cluster_size.get(c, ""), g, n])


if __name__ == "__main__":
    main()
