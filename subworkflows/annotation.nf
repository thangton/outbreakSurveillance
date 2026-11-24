// D2 – gene prediction + functional annotation
// Final D2 choice: Prodigal (-p meta -g 11) -> eggNOG-mapper (DIAMOND vs eggNOG 5.0.2)

process PRODIGAL {
    tag "$id"
    label 'process_single'
    container 'quay.io/biocontainers/prodigal:2.6.3--h7b50bb2_10'

    input:
    tuple val(id), path(contigs)

    output:
    tuple val(id), path("${id}.faa"), emit: faa
    tuple val(id), path("${id}.gff"), emit: gff
    tuple val(id), path("${id}.fna"), emit: fna

    script:
    """
    prodigal -i ${contigs} -f gff -o ${id}.gff -a ${id}.faa -d ${id}.fna \\
        ${task.ext.args ?: ''} > ${id}.prodigal.log 2>&1
    """

    stub:
    "touch ${id}.faa ${id}.gff ${id}.fna"
}

process EGGNOG_MAPPER {
    tag "$id"
    label 'process_medium'
    container 'quay.io/biocontainers/eggnog-mapper:2.1.12--pyhdfd78af_0'

    input:
    tuple val(id), path(faa)
    path db

    output:
    tuple val(id), path("${id}.emapper.annotations"), emit: annotations
    tuple val(id), path("${id}.emapper.seed_orthologs"), emit: orthologs

    script:
    """
    emapper.py -i ${faa} -o ${id} --data_dir ${db} --cpu ${task.cpus} ${task.ext.args ?: ''}
    """

    stub:
    "touch ${id}.emapper.annotations ${id}.emapper.seed_orthologs"
}

workflow ANNOTATION {
    take:
    contigs   // [ id, contigs.fa ]

    main:
    PRODIGAL(contigs)

    ch_eggnog = Channel.empty()
    if (params.eggnog_db) {
        EGGNOG_MAPPER(PRODIGAL.out.faa, file(params.eggnog_db, checkIfExists: true))
        ch_eggnog = EGGNOG_MAPPER.out.annotations
    } else {
        log.warn "--eggnog_db not set: skipping eggNOG-mapper (see README for the database download)"
    }

    emit:
    proteins = PRODIGAL.out.faa
    gff      = PRODIGAL.out.gff
    eggnog   = ch_eggnog
}
