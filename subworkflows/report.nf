process MULTIQC {
    label 'process_single'
    container 'quay.io/biocontainers/multiqc:1.25--pyhdfd78af_0'

    input:
    path inputs

    output:
    path 'multiqc_report.html'

    script:
    """
    multiqc -f .
    """

    stub:
    "touch multiqc_report.html"
}

process SUMMARY_REPORT {
    label 'process_single'
    container params.utils_container

    input:
    path quast,  stageAs: 'in/quast'
    path gff,    stageAs: 'in/gff/*'
    path eggnog, stageAs: 'in/eggnog/*'
    path mlst,   stageAs: 'in/mlst.tsv'
    path blast,  stageAs: 'in/blast16s.tsv'
    path skani,  stageAs: 'in/skani.tsv'
    path checkm, stageAs: 'in/checkm.tsv'
    path gunc,   stageAs: 'in/gunc.tsv'
    path amr,    stageAs: 'in/amr/*'
    path clusters, stageAs: 'in/clusters.tsv'

    output:
    path 'outbreak_summary.tsv'
    path 'amr_by_cluster.tsv'

    script:
    """
    build_report.py --indir in --out outbreak_summary.tsv --amr-by-cluster amr_by_cluster.tsv
    """

    stub:
    "touch outbreak_summary.tsv amr_by_cluster.tsv"
}
