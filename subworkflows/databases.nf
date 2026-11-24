// Reference prep + database downloads.
// Everything is written to --db_cache via storeDir, so each database downloads once and is reused.


process PREP_REFERENCE {
    label 'process_single'
    container params.utils_container

    input:
    path gbk

    output:
    path 'reference.fna'

    script:
    """
    python3 -c "from Bio import SeqIO; SeqIO.convert('${gbk}', 'genbank', 'reference.fna', 'fasta')"
    """

    stub:
    "touch reference.fna"
}

process DOWNLOAD_BLAST16S {
    label 'process_download'
    container params.utils_container
    storeDir params.db_cache

    output:
    path 'blast_16s'

    script:
    """
    mkdir blast_16s && cd blast_16s
    curl -fsSL https://ftp.ncbi.nlm.nih.gov/blast/db/16S_ribosomal_RNA.tar.gz | tar -xz
    curl -fsSL https://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz | tar -xz
    """

    stub:
    "mkdir blast_16s"
}

process DOWNLOAD_NCBI_REFS {
    tag "$genus"
    label 'process_download'
    container params.utils_container
    storeDir params.db_cache

    input:
    val genus

    output:
    path "ncbi_refs_${genus}"

    script:
    """
    datasets download genome taxon "${genus}" --reference \\
        --assembly-level complete,chromosome --include genome --filename refs.zip
    unzip -q refs.zip -d tmp
    mkdir ncbi_refs_${genus}
    find tmp -name '*.fna' -exec mv {} ncbi_refs_${genus}/ \\;
    rm -rf tmp refs.zip
    """

    stub:
    "mkdir ncbi_refs_${genus} && touch ncbi_refs_${genus}/ref.fna"
}

process DOWNLOAD_CHECKM_DB {
    label 'process_download'
    container params.utils_container
    storeDir params.db_cache

    output:
    path 'checkm_data'

    script:
    """
    mkdir checkm_data
    curl -fsSL https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz \\
        | tar -xz -C checkm_data
    """

    stub:
    "mkdir checkm_data"
}

process DOWNLOAD_GUNC_DB {
    label 'process_download'
    container 'quay.io/biocontainers/gunc:1.0.6--pyhdfd78af_0'
    storeDir params.db_cache

    output:
    path 'gunc_db'

    script:
    """
    mkdir gunc_db
    gunc download_db gunc_db -db progenomes
    """

    stub:
    "mkdir gunc_db && touch gunc_db/gunc_db_progenomes2.1.dmnd"
}

process AMRFINDER_UPDATE {
    label 'process_download'
    container 'quay.io/biocontainers/ncbi-amrfinderplus:3.12.8--h283d18e_0'
    storeDir params.db_cache

    output:
    path 'amrfinder_db'

    script:
    """
    amrfinder_update -d amrfinder_db
    """

    stub:
    "mkdir -p amrfinder_db/latest"
}
