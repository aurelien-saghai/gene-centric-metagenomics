#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.path = null
params.restrict_read_length = null
params.nrfA_gpkg = "/opt/databases/graftM/nrfA_refpack_2021.gpkg"
params.nirK_gpkg = "/opt/databases/graftM/nirK_refpack_2021.gpkg"
params.nirS_gpkg = "/opt/databases/graftM/nirS_refpack_2021.gpkg"
params.nosZ_gpkg = "/opt/databases/graftM/nosZ_refpack_2021.gpkg"
params.outdir = "graftM"

params.help = false
params.genes = null

process UNZIP_FASTQ {

    tag "$sample_id"

    cpus 8

    input:
    tuple val(sample_id), path(read1)

    output:
    tuple val(sample_id), path("${sample_id}.fastq")

    script:
    """
    gzip -dc -- "${read1}" > "${sample_id}.fastq"
    """
}

process GRAFTM {

    tag "${sample_id}:${gene}"

    cpus 8
    maxForks 1

    publishDir { "${params.outdir}/${gene}" }, mode: 'copy',
        saveAs: { filename -> filename.tokenize('/').last() }

    input:
    tuple val(sample_id), path(fastq), val(gene), val(gpkg)

    output:
    path("${sample_id}/${sample_id}")
    path("${sample_id}_graftM.log")

    script:
    """
    graftM graft \
        --graftm_package "${gpkg}" \
        --forward "${fastq}" \
        --evalue 1e-5 \
        --resolve_placements \
        --search_method hmmsearch \
        --restrict_read_length "${params.restrict_read_length}" \
        --output_directory "${sample_id}" \
        --log "${sample_id}_graftM.log"

    """
}

workflow {

    if (params.help) {
        log.info """
        Usage:
          module load graftM/0.15.1
          module load gappa/0.8.0
          module load pplacer/1.1.19
          
          nextflow run graftM.nf --path path/to/metagenomes --restrict_read_length 150 --genes nrfA,nirK,nirS,nosZ

        Required arguments:
          --path                   Path to the folder with the metagenomes
          --restrict_read_length   Read length used by graftM
          --genes                  Comma-separated (no spaces), case-sensitive gene names

        Available genes:
          nrfA, nirK, nirS, nosZ
        """.stripIndent()
    } else {

    gene_packages = [
        nrfA: params.nrfA_gpkg,
        nirK: params.nirK_gpkg,
        nirS: params.nirS_gpkg,
        nosZ: params.nosZ_gpkg
    ]

    selected_genes = params.genes
        .split(',')
        .collect { it.trim() }
        .findAll { it }
        .unique()

    if (!selected_genes) {
        error 'Specify at least one gene with --genes, for example --genes nrfA'
    }

    unknown_genes = selected_genes.findAll { !gene_packages.containsKey(it) }
    if (unknown_genes) {
        error "Unknown gene(s): ${unknown_genes.join(', ')}. Available genes: ${gene_packages.keySet().join(', ')}"
    }

    reads_ch = Channel
        .fromPath("${params.path}/*R1.fastq.gz", checkIfExists: true)
        .map { read1 ->
            def sample_id = read1.name.replaceFirst(/_?R1\.fastq\.gz$/, '')
            tuple(sample_id, read1)
    }

    UNZIP_FASTQ(reads_ch)

    graftm_inputs = UNZIP_FASTQ.out.flatMap { sample_id, fastq ->
        selected_genes.collect { gene ->
            tuple(sample_id, fastq, gene, gene_packages[gene])
        }
    }

    GRAFTM(graftm_inputs)
    }
}