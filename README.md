# outbreakSurveillance

Nextflow DSL2 pipeline for bacterial outbreak surveillance from paired-end Illumina reads: assembly → annotation → typing/QC → AMR and SNP-based transmission clustering. Every step runs in a pinned container (Docker, Singularity or Apptainer), so the same command works on a laptop, a Linux workstation or an HPC cluster.

It merges the four parts of the BIOL 7210 Computational Genomics (Spring 2026) Team D project into one workflow. The defaults are the tools and parameters each group picked in its own benchmark, and they were applied to 34 *Neisseria gonorrhoeae* isolates.

```
reads ─► fastp ─► SPAdes|SKESA ─► contigs ─┬─► QUAST
           │                               ├─► Prodigal ─► eggNOG-mapper                    (D2)
           │                               │      └──────► GUNC                             (D3)
           │                               ├─► MLST ─► genus ─► NCBI refs ─► skani ─► species ─► CheckM
           │                               ├─► barrnap 16S ─► BLAST                         (D3)
           │                               ├─► AMRFinderPlus                                (D4)
           │                               └─► Parsnp ─────────────────┐
           └─► Snippy ─► snippy-core ─► IQ-TREE ─► SNP clusters ◄──────┘                    (D4)
                                                         │
                                   outbreak_summary.tsv + amr_by_cluster.tsv + MultiQC
```

## Quick start

```bash
# requirements: Nextflow >= 24.04, Java 17+, Docker (or Singularity/Apptainer)
git clone https://github.com/thangton/outbreakSurveillance.git && cd outbreakSurveillance

nextflow run main.nf -profile docker \
    --reads '/path/to/reads/*_R{1,2}_001.fastq.gz' \
    --outdir results
```

On an Apple Silicon Mac, add `arm` (`-profile docker,arm`). Biocontainers are amd64-only, so Docker runs them under emulation.

On PACE:

```bash
module load nextflow
nextflow run main.nf -profile pace --pace_account gts-<pi> --reads '...' --outdir results
```

Dry-run the wiring without any tools: `nextflow run main.nf -stub --reads '...'`.

## Input

You can give the reads in either of two ways:

- `--reads` is a glob over paired files. The sample ID is the part of the file name before the first `_`, so `D9d128a010_S01_L001_R1_001.fastq.gz` becomes `D9d128a010`.
- `--input` is a CSV with the columns `sample,fastq_1,fastq_2` (see `assets/samplesheet_example.csv`).

The comparative stage builds a tree, so it needs at least 4 samples.

## Parameters

| Parameter | Default | |
|---|---|---|
| `--assembler` | `spades` | `spades` (D1's recommendation) or `skesa` (faster; D4 used it) |
| `--eggnog_db` | – | eggNOG 5.0.2 data dir. eggNOG-mapper is skipped when this is unset |
| `--taxon` | from MLST | Genus used to download skani reference genomes |
| `--checkm_species` | from skani | Species for CheckM `taxonomy_wf`. Falls back to `lineage_wf --reduced_tree` |
| `--reference_gbk` | `assets/reference/CP016016.gb` | Reference for Snippy, Parsnp and SNP-site annotation |
| `--amrfinder_organism` | `Neisseria_gonorrhoeae` | Turns on AMRFinderPlus point-mutation screening |
| `--skip_checkm` / `--skip_gunc` / `--skip_comparative` | `false` | |
| `--db_cache` | `./databases` | Where downloaded databases are stored and reused |
| `--max_cpus` / `--max_memory` / `--max_time` | 16 / 64 GB / 48 h | Resource ceiling |

Tool arguments live in `conf/modules.config` as `ext.args`. For example, the fastp trimming D1 settled on is `-f 10 -l 50 -r --cut_right_window_size 10 --cut_right_mean_quality 30`.

## Databases

These download automatically on the first run into `--db_cache`, then get reused:

| Database | Size | Used by |
|---|---|---|
| NCBI 16S_ribosomal_RNA + taxdb | ~100 MB | BLAST genus check |
| NCBI reference genomes for the genus | varies | skani |
| CheckM data (2015-01-16) | ~275 MB | CheckM |
| GUNC proGenomes 2.1 | ~13 GB | GUNC |
| AMRFinderPlus | ~200 MB | AMRFinderPlus |

eggNOG is ~50 GB, so it isn't downloaded automatically:

```bash
mkdir eggnog_db && docker run --rm -v $PWD/eggnog_db:/db quay.io/biocontainers/eggnog-mapper:2.1.12--pyhdfd78af_0 \
    download_eggnog_data.py --data_dir /db -y
nextflow run main.nf ... --eggnog_db $PWD/eggnog_db
```

## Output

```
results/
├── outbreak_summary.tsv          one row per isolate: assembly stats, genes, ST, 16S genus, species/ANI,
│                                 completeness/contamination, GUNC, SNP cluster, AMR determinants
├── amr_by_cluster.tsv            AMR determinant prevalence within each SNP cluster
├── multiqc_report.html           fastp + QUAST
├── 01_assembly/                  contigs, fastp reports, QUAST
├── 02_annotation/                Prodigal gff/faa/fna, eggNOG annotations
├── 03_characterization/          mlst, barrnap, blast_16s, skani, checkm, gunc
├── 04_comparative/               amr, snippy, snippy_core, iqtree, parsnp,
│                                 snp_clusters (pairwise SNPs, clusters, defining sites, cladograms)
└── pipeline_info/                Nextflow timeline, report, trace, DAG
```

SNP clusters use D4's rule: two isolates are linked when their pairwise core-SNP distance is at or below the first quartile of all pairwise distances. Each cluster is then checked for monophyly on the IQ-TREE tree.

## Where each stage came from

| Stage | Origin | Changes from the course version |
|---|---|---|
| fastp → SPAdes → QUAST | D1 (6 trimmer × assembler combos benchmarked) | Python/conda controllers replaced by Nextflow processes |
| Prodigal → eggNOG-mapper | D2 (also tested Glimmer, InterProScan) | Same parameters (`-p meta -g 11`, eggNOG 5.0.2) |
| MLST, barrnap/BLAST, skani, CheckM, GUNC | D3 | Same genus → species → CheckM hand-off. GUNC now reuses D2's Prodigal proteins instead of re-calling genes |
| AMRFinderPlus, Snippy/IQ-TREE, Parsnp, clustering | D4 | D4's cluster and plotting scripts kept in `bin/` (path names generalized) |

## Notes

- The default reference and AMRFinder organism are set for *N. gonorrhoeae*. For another species, change `--reference_gbk` and `--amrfinder_organism`, and set `--taxon` if the MLST scheme name isn't the genus (for example, `ecoli_achtman_4`).
- Prodigal stays in D2's metagenome mode (`-p meta`). For single isolates, `ext.args = '-g 11'` (single mode) is the more usual choice.
- `ghcr.io/thangton/outbreak-surveillance-utils` holds the Python helpers and downloaders. `.github/workflows/docker.yml` builds it on push. Before it's published, build it locally with `docker build -t ghcr.io/thangton/outbreak-surveillance-utils:1.0.0 docker/`.
