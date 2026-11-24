#!/usr/bin/env nextflow
/*
 * outbreakSurveillance
 *   D1  read cleaning + de novo assembly      fastp -> SPAdes|SKESA -> QUAST
 *   D2  gene prediction + annotation          Prodigal -> eggNOG-mapper
 *   D3  typing, taxonomy, genome QC           MLST, 16S BLAST, skani, CheckM, GUNC
 *   D4  AMR + comparative genomics            AMRFinderPlus, Snippy -> IQ-TREE, Parsnp, SNP clusters
 */
nextflow.enable.dsl = 2

include { ASSEMBLY         } from './subworkflows/assembly'
include { ANNOTATION       } from './subworkflows/annotation'
include { CHARACTERIZATION } from './subworkflows/characterization'
include { COMPARATIVE      } from './subworkflows/comparative'
include { MULTIQC; SUMMARY_REPORT } from './subworkflows/report'

def NO_FILE = file("${projectDir}/assets/NO_FILE")

workflow {
    if (!params.input && !params.reads) error "Provide --input samplesheet.csv or --reads '<dir>/*_R{1,2}_001.fastq.gz'"
    if (!(params.assembler in ['spades', 'skesa'])) error "--assembler must be 'spades' or 'skesa'"

    // [ id, [R1, R2] ]
    ch_reads = params.input
        ? Channel.fromPath(params.input, checkIfExists: true)
            .splitCsv(header: true)
            .map { row -> tuple(row.sample, [file(row.fastq_1, checkIfExists: true), file(row.fastq_2, checkIfExists: true)]) }
        : Channel.fromFilePairs(params.reads, checkIfExists: true)
            .map { key, files -> tuple(key.tokenize('_')[0], files) }

    ASSEMBLY(ch_reads)
    ANNOTATION(ASSEMBLY.out.contigs)
    CHARACTERIZATION(ASSEMBLY.out.contigs, ANNOTATION.out.proteins)
    COMPARATIVE(ASSEMBLY.out.clean_reads, ASSEMBLY.out.contigs)

    MULTIQC(ASSEMBLY.out.fastp_json.mix(ASSEMBLY.out.quast).collect())

    SUMMARY_REPORT(
        ASSEMBLY.out.quast,
        ANNOTATION.out.gff.map { it[1] }.collect(),
        ANNOTATION.out.eggnog.map { it[1] }.collect().ifEmpty([NO_FILE]),
        CHARACTERIZATION.out.mlst,
        CHARACTERIZATION.out.blast16s.ifEmpty(NO_FILE),
        CHARACTERIZATION.out.skani.ifEmpty(NO_FILE),
        CHARACTERIZATION.out.checkm.ifEmpty(NO_FILE),
        CHARACTERIZATION.out.gunc.ifEmpty(NO_FILE),
        COMPARATIVE.out.amr.map { it[1] }.collect(),
        COMPARATIVE.out.clusters.ifEmpty(NO_FILE)
    )
}

workflow.onComplete {
    log.info workflow.success
        ? "Done. Summary: ${params.outdir}/outbreak_summary.tsv"
        : "Failed – see .nextflow.log"
}
