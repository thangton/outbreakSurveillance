// D3 – genotyping, taxonomy and genome QC
// MLST -> genus  ->  skani vs NCBI references -> species  ->  CheckM taxonomy_wf
// barrnap 16S -> BLAST (independent genus check); GUNC on Prodigal proteins from D2

include { DOWNLOAD_BLAST16S; DOWNLOAD_NCBI_REFS; DOWNLOAD_CHECKM_DB; DOWNLOAD_GUNC_DB } from './databases'

process MLST {
    label 'process_low'
    container 'quay.io/biocontainers/mlst:2.23.0--hdfd78af_0'

    input:
    path contigs

    output:
    path 'mlst.tsv', emit: tsv

    script:
    """
    mlst --threads ${task.cpus} ${contigs} > mlst.tsv
    """

    stub:
    "printf 'S1.contigs.fa\\tneisseria\\t7822\\n' > mlst.tsv"
}

process BARRNAP {
    tag "$id"
    label 'process_single'
    container 'quay.io/biocontainers/barrnap:0.9--hdfd78af_4'

    input:
    tuple val(id), path(contigs)

    output:
    path "${id}.16S.fa",   emit: rrna16s
    path "${id}.rrna.gff", emit: gff

    script:
    // keep the first 16S hit, renamed to the sample id
    """
    barrnap --kingdom bac --threads ${task.cpus} --outseq rrna.fa ${contigs} > ${id}.rrna.gff
    awk '/^>/ { if (p) exit; p = (\$0 ~ /16S_rRNA/); if (p) print ">${id}"; next } p' rrna.fa > ${id}.16S.fa
    """

    stub:
    "touch ${id}.16S.fa ${id}.rrna.gff"
}

process BLAST_16S {
    label 'process_low'
    container 'quay.io/biocontainers/blast:2.15.0--pl5321h6f7f691_1'

    input:
    path seqs
    path db

    output:
    path 'blast16s.tsv', emit: tsv

    script:
    """
    cat ${seqs} > query.fa
    if [ -s query.fa ]; then
        BLASTDB=${db} blastn -db 16S_ribosomal_RNA -query query.fa -num_threads ${task.cpus} \\
            -max_target_seqs 10 -outfmt '6 qseqid evalue bitscore sscinames' > blast16s.tsv
    else
        touch blast16s.tsv
    fi
    """

    stub:
    "touch blast16s.tsv"
}

process SKANI {
    label 'process_medium'
    container 'quay.io/biocontainers/skani:0.3.1--ha6fb395_0'

    input:
    path contigs
    path refs

    output:
    path 'skani.tsv', emit: tsv

    script:
    """
    ls ${contigs} > queries.txt
    ls ${refs}/*.fna > refs.txt
    skani dist -t ${task.cpus} --ql queries.txt --rl refs.txt -o skani.tsv
    """

    stub:
    """
    printf 'Ref_file\\tQuery_file\\tANI\\tAlign_fraction_ref\\tAlign_fraction_query\\tRef_name\\tQuery_name\\n' > skani.tsv
    printf 'r.fna\\tS1.contigs.fa\\t99.7\\t95\\t96\\tNZ_AP023069.1 Neisseria gonorrhoeae strain X\\tS1\\n' >> skani.tsv
    """
}

process CHECKM {
    label 'process_high'
    container 'quay.io/biocontainers/checkm-genome:1.2.2--pyhdfd78af_1'

    input:
    path contigs
    path db
    val species

    output:
    path 'checkm.tsv', emit: tsv

    script:
    // taxonomy_wf on the inferred species (as in D3); lineage_wf if the species is unknown to CheckM
    """
    export CHECKM_DATA_PATH=\$PWD/${db}
    mkdir bins
    for f in ${contigs}; do ln -s ../\$f bins/\${f%.contigs.fa}.fa; done

    if [ -n "${species}" ] && checkm taxonomy_wf -t ${task.cpus} -x fa species "${species}" bins checkm_out; then
        :
    else
        rm -rf checkm_out
        checkm lineage_wf --reduced_tree -t ${task.cpus} -x fa bins checkm_out
    fi
    checkm qa \$(ls checkm_out/*.ms | head -1) checkm_out -o 2 --tab_table -f checkm.tsv -t ${task.cpus}
    """

    stub:
    "touch checkm.tsv"
}

process GUNC {
    label 'process_high'
    container 'quay.io/biocontainers/gunc:1.0.6--pyhdfd78af_0'

    input:
    path faa
    path db

    output:
    path 'gunc.tsv', emit: tsv
    path 'gunc_out', emit: dir

    script:
    // reuse D2's Prodigal proteins; strip terminal '*' which GUNC/DIAMOND reject
    """
    mkdir genes
    for f in ${faa}; do sed 's/[*]\$//' \$f > genes/\$f; done
    gunc run --db_file \$(ls ${db}/*.dmnd | head -1) --input_dir genes --file_suffix .faa \\
        --gene_calls --threads ${task.cpus} --temp_dir . --out_dir gunc_out
    cp gunc_out/GUNC.*maxCSS_level.tsv gunc.tsv
    """

    stub:
    "mkdir gunc_out && touch gunc.tsv"
}

workflow CHARACTERIZATION {
    take:
    contigs    // [ id, contigs.fa ]
    proteins   // [ id, proteins.faa ]

    main:
    ch_all_contigs = contigs.map { it[1] }.collect()

    MLST(ch_all_contigs)

    // genus = most common MLST scheme (e.g. 'neisseria' -> 'Neisseria'), unless --taxon is given
    ch_genus = params.taxon
        ? Channel.value(params.taxon)
        : MLST.out.tsv.map { f ->
            def schemes = f.readLines().collect { it.tokenize('\t')[1] }.findAll { it && it != '-' }
            if (!schemes) error "MLST assigned no scheme; set --taxon <genus>"
            schemes.countBy { it }.max { it.value }.key.replaceAll(/[_0-9].*$/, '').capitalize()
        }

    BARRNAP(contigs)
    BLAST_16S(BARRNAP.out.rrna16s.collect(), DOWNLOAD_BLAST16S())

    SKANI(ch_all_contigs, DOWNLOAD_NCBI_REFS(ch_genus))

    // species = consensus of each sample's best skani hit, unless --checkm_species is given
    ch_species = params.checkm_species
        ? Channel.value(params.checkm_species)
        : SKANI.out.tsv.map { f ->
            def rows = f.readLines().drop(1).collect { it.split('\t') }
            def best = rows.groupBy { it[1] }.collect { q, rs -> rs.max { it[2] as double } }
            def sp   = best.collect { it[5].tokenize(' ') }.findAll { it.size() >= 3 }.collect { it[1..2].join(' ') }
            sp ? sp.countBy { it }.max { it.value }.key : ''
        }

    ch_checkm = Channel.empty()
    if (!params.skip_checkm) {
        CHECKM(ch_all_contigs, DOWNLOAD_CHECKM_DB(), ch_species)
        ch_checkm = CHECKM.out.tsv
    }

    ch_gunc = Channel.empty()
    if (!params.skip_gunc) {
        GUNC(proteins.map { it[1] }.collect(), DOWNLOAD_GUNC_DB())
        ch_gunc = GUNC.out.tsv
    }

    emit:
    mlst     = MLST.out.tsv
    blast16s = BLAST_16S.out.tsv
    skani    = SKANI.out.tsv
    checkm   = ch_checkm
    gunc     = ch_gunc
    genus    = ch_genus
    species  = ch_species
}
