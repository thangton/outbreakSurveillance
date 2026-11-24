// D4 – AMR + comparative genomics
// AMRFinderPlus on assemblies; Snippy (reads vs reference) -> snippy-core -> IQ-TREE -> SNP clusters;
// Parsnp core-genome tree from assemblies as an alignment-free cross-check

include { PREP_REFERENCE; AMRFINDER_UPDATE } from './databases'

process AMRFINDER {
    tag "$id"
    label 'process_low'
    container 'quay.io/biocontainers/ncbi-amrfinderplus:3.12.8--h283d18e_0'

    input:
    tuple val(id), path(contigs)
    path db

    output:
    tuple val(id), path("${id}_amr.tsv"), emit: tsv

    script:
    def organism = params.amrfinder_organism ? "-O ${params.amrfinder_organism}" : ''
    """
    amrfinder -n ${contigs} ${organism} -d ${db}/latest --threads ${task.cpus} \\
        --name ${id} -o ${id}_amr.tsv
    """

    stub:
    "touch ${id}_amr.tsv"
}

process SNIPPY {
    tag "$id"
    label 'process_medium'
    container 'quay.io/biocontainers/snippy:4.6.0--hdfd78af_4'

    input:
    tuple val(id), path(reads)
    path ref

    output:
    path "${id}", emit: dir

    script:
    """
    snippy --cpus ${task.cpus} --ram ${task.memory.toGiga()} --ref ${ref} \\
        --R1 ${reads[0]} --R2 ${reads[1]} --outdir ${id}
    """

    stub:
    "mkdir ${id}"
}

process SNIPPY_CORE {
    label 'process_low'
    container 'quay.io/biocontainers/snippy:4.6.0--hdfd78af_4'

    input:
    path dirs
    path ref

    output:
    path 'core.aln',      emit: aln
    path 'core.tab',      emit: tab
    path 'core.{full.aln,txt,vcf}'

    script:
    """
    snippy-core --ref ${ref} --prefix core ${dirs}
    """

    stub:
    "touch core.aln core.tab core.full.aln core.txt core.vcf"
}

process IQTREE {
    label 'process_medium'
    container 'quay.io/biocontainers/iqtree:2.3.6--hdbdd923_0'

    input:
    path aln

    output:
    path 'core.treefile', emit: tree
    path 'core.{iqtree,log}'

    script:
    """
    iqtree2 -s ${aln} -T ${task.cpus} --prefix core ${task.ext.args ?: ''}
    """

    stub:
    "touch core.treefile core.iqtree core.log"
}

process PARSNP {
    label 'process_medium'
    container 'quay.io/biocontainers/parsnp:2.0.5--hdcf5f25_0'

    input:
    path contigs
    path ref

    output:
    path 'parsnp.tree', emit: tree
    path 'parsnp.ggr'

    script:
    """
    mkdir genomes
    for f in ${contigs}; do ln -s ../\$f genomes/\${f%.contigs.fa}.fa; done
    parsnp -r ${ref} -d genomes -o out -p ${task.cpus} ${task.ext.args ?: ''}
    cp out/parsnp.tree out/parsnp.ggr .
    """

    stub:
    "touch parsnp.tree parsnp.ggr"
}

process SNP_CLUSTERS {
    label 'process_single'
    container params.utils_container

    input:
    path aln
    path tab
    path snippy_tree
    path parsnp_tree
    path gbk

    output:
    path 'snippy_clusters_q1.tsv', emit: clusters
    path '*'

    script:
    // D4 method: link samples whose pairwise core-SNP distance <= Q1 of all pairwise distances
    """
    analyze_snippy_pairwise_vs_tree.py --alignment ${aln} --tree ${snippy_tree} --out-dir .
    annotate_cluster_signatures.py --clusters snippy_clusters_q1.tsv --core-tab ${tab} \\
        --genbank ${gbk} --out-prefix snippy_cluster_defining_sites
    plot_cladogram.py --tree ${snippy_tree} --title 'Snippy core-SNP ML tree' --out-prefix snippy_cladogram
    plot_cladogram.py --tree ${parsnp_tree} --title 'Parsnp core-genome tree' --out-prefix parsnp_cladogram
    """

    stub:
    "touch snippy_clusters_q1.tsv"
}

workflow COMPARATIVE {
    take:
    reads     // [ id, [R1, R2] ]  cleaned
    contigs   // [ id, contigs.fa ]

    main:
    gbk = file(params.reference_gbk, checkIfExists: true)

    AMRFINDER(contigs, AMRFINDER_UPDATE())

    ch_clusters = Channel.empty()
    if (!params.skip_comparative) {
        ref_fna = PREP_REFERENCE(gbk)
        SNIPPY(reads, gbk)
        SNIPPY_CORE(SNIPPY.out.dir.collect(), gbk)
        IQTREE(SNIPPY_CORE.out.aln)
        PARSNP(contigs.map { it[1] }.collect(), ref_fna)
        SNP_CLUSTERS(SNIPPY_CORE.out.aln, SNIPPY_CORE.out.tab, IQTREE.out.tree, PARSNP.out.tree, gbk)
        ch_clusters = SNP_CLUSTERS.out.clusters
    }

    emit:
    amr      = AMRFINDER.out.tsv
    clusters = ch_clusters
}
