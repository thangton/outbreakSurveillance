// D1 – read cleaning + de novo assembly
// Final D1 benchmark choice: fastp -> SPAdes (SKESA kept as the fast option used in D4)

process FASTP {
    tag "$id"
    label 'process_low'
    container 'quay.io/biocontainers/fastp:0.23.4--h5f740d0_0'

    input:
    tuple val(id), path(reads)

    output:
    tuple val(id), path("${id}_R{1,2}.trim.fq.gz"), emit: reads
    path "${id}.fastp.json",                          emit: json
    path "${id}.fastp.html",                          emit: html

    script:
    """
    fastp -i ${reads[0]} -I ${reads[1]} \\
        -o ${id}_R1.trim.fq.gz -O ${id}_R2.trim.fq.gz \\
        -w ${task.cpus} -j ${id}.fastp.json -h ${id}.fastp.html \\
        ${task.ext.args ?: ''}
    """

    stub:
    "touch ${id}_R1.trim.fq.gz ${id}_R2.trim.fq.gz ${id}.fastp.json ${id}.fastp.html"
}

process SPADES {
    tag "$id"
    label 'process_high'
    container 'quay.io/biocontainers/spades:3.15.5--h95f258a_1'

    input:
    tuple val(id), path(reads)

    output:
    tuple val(id), path("${id}.contigs.fa"), emit: contigs
    path "${id}.spades.log",                emit: log

    script:
    """
    spades.py -1 ${reads[0]} -2 ${reads[1]} \\
        -t ${task.cpus} -m ${task.memory.toGiga()} -o spades ${task.ext.args ?: ''}
    cp spades/contigs.fasta ${id}.contigs.fa
    cp spades/spades.log ${id}.spades.log
    """

    stub:
    "touch ${id}.contigs.fa ${id}.spades.log"
}

process SKESA {
    tag "$id"
    label 'process_medium'
    container 'quay.io/biocontainers/skesa:2.5.1--hdcf5f25_1'

    input:
    tuple val(id), path(reads)

    output:
    tuple val(id), path("${id}.contigs.fa"), emit: contigs
    path "${id}.skesa.log",                 emit: log

    script:
    """
    skesa --reads ${reads[0]},${reads[1]} --cores ${task.cpus} \\
        --memory ${task.memory.toGiga()} --contigs_out ${id}.contigs.fa \\
        ${task.ext.args ?: ''} 2> ${id}.skesa.log
    """

    stub:
    "touch ${id}.contigs.fa ${id}.skesa.log"
}

process QUAST {
    label 'process_medium'
    container 'quay.io/biocontainers/quast:5.2.0--py39pl5321h2add14b_1'

    input:
    path contigs

    output:
    path 'quast', emit: dir

    script:
    """
    quast.py -o quast -t ${task.cpus} ${contigs}
    """

    stub:
    "mkdir quast && touch quast/transposed_report.tsv"
}

workflow ASSEMBLY {
    take:
    reads   // [ id, [R1, R2] ]

    main:
    FASTP(reads)
    ch_asm  = params.assembler == 'skesa' ? SKESA(FASTP.out.reads) : SPADES(FASTP.out.reads)
    QUAST(ch_asm.contigs.map { it[1] }.collect())

    emit:
    clean_reads = FASTP.out.reads
    contigs     = ch_asm.contigs
    fastp_json  = FASTP.out.json
    quast       = QUAST.out.dir
}
