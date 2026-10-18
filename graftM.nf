#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// This workflow takes metagenomic reads and runs graftM classification for one or more
// marker genes (e.g. nrfA, nirK, nirS, nosZ). It is designed to process many samples
// in parallel while keeping each gene/sample combination as an independent task.

// Global parameters used by the workflow.
params.path = null
params.restrict_read_length = null
params.nrfA_gpkg = "databases/graftM_packages/nrfA_refpack_2021.gpkg"
params.nirK_gpkg = "databases/graftM_packages/nirK_refpack_2021.gpkg"
params.nirS_gpkg = "databases/graftM_packages/nirS_refpack_2021.gpkg"
params.nosZ_gpkg = "databases/graftM_packages/nosZ_refpack_2021.gpkg"
params.outdir = "graftM"
params.threads = 1
params.conda_env = "${System.getProperty('user.home')}/miniforge3/envs/graftM_workflow"

params.help = false
params.genes = null

// Decompress a gzipped read file to a plain FASTQ so that graftM can process it.
// Each task receives one sample's read file and writes a sample-specific FASTQ name.
process UNZIP_FASTQ {

    conda params.conda_env

    tag "$sample_id"

    cpus params.threads.toString().toInteger()

    input:
    tuple val(sample_id), path(read1)

    output:
    tuple val(sample_id), path("${sample_id}.fastq")

    script:
    """
    gzip -dc -- "${read1}" > "${sample_id}.fastq"
    """
}

// Run graftM for a single sample and a single gene.
// The publishDir section copies the output into a gene-specific folder in the final
// results directory while keeping the original file names simple.
process GRAFTM {

    conda params.conda_env

    tag "${sample_id}:${gene}"

    cpus params.threads.toString().toInteger()
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
        --threads "${task.cpus}" \
        --restrict_read_length "${params.restrict_read_length}" \
        --output_directory "${sample_id}" \
        --log "${sample_id}_graftM.log"

    """
}

workflow {

    // If the user asks for help, print the expected command syntax and exit without
    // running any analyses.
    if (params.help) {
        log.info """
        Usage:
        The specified conda environment must contain graftM and all its required programs.

          nextflow run graftM.nf \\
            --path path/to/metagenomes \\
            --restrict_read_length 150 \\
            --genes nrfA,nirK,nirS,nosZ \\
            --threads 8 \\
            --conda_env ~/miniforge3/envs/graftM_workflow \\
            -with-conda

        Required arguments:
          --path                   Path to the folder with the metagenomes
          --restrict_read_length   Read length used by graftM
          --genes                  Comma-separated (no spaces), case-sensitive gene names
          --threads                Number of CPU threads used per process (default: 1)
          --conda_env              Path to the conda environment containing required programs (default: ~/miniforge3/envs/graftM_workflow)
          -with-conda               Use conda to run the workflow

        Available genes:
          nrfA, nirK, nirS, nosZ
        """.stripIndent()
    } else {

    // Validate the thread count before launching any tasks.
    if (!(params.threads.toString() ==~ /[1-9][0-9]*/)) {
        error "The --threads value must be a positive integer; received: ${params.threads}"
    }
    def conda_env_path = file(params.conda_env.toString())
    if (!conda_env_path.toFile().isDirectory()) {
        error "Conda environment directory does not exist: ${conda_env_path}"
    }

    // Map each allowed gene name to the corresponding graftM reference package.
    gene_packages = [
        nrfA: params.nrfA_gpkg,
        nirK: params.nirK_gpkg,
        nirS: params.nirS_gpkg,
        nosZ: params.nosZ_gpkg
    ]

    // Parse the comma-separated --genes argument into a clean list of gene names.
    // This also removes blanks and duplicate entries before validation.
    selected_genes = params.genes
        .split(',')
        .collect { it.trim() }
        .findAll { it }
        .unique()

    if (!selected_genes) {
        error 'Specify at least one gene with --genes, for example --genes nrfA'
    }

    // Reject any gene names that are not supported by the workflow.
    unknown_genes = selected_genes.findAll { !gene_packages.containsKey(it) }
    if (unknown_genes) {
        error "Unknown gene(s): ${unknown_genes.join(', ')}. Available genes: ${gene_packages.keySet().join(', ')}"
    }

    // Discover all R1 reads in the input directory and derive a sample ID from the file name.
    // Example: sample_001_R1.fastq.gz becomes sample_001.
    reads_ch = Channel
        .fromPath("${params.path}/*R1.fastq.gz", checkIfExists: true)
        .map { read1 ->
            def sample_id = read1.name.replaceFirst(/_?R1\.fastq\.gz$/, '')
            tuple(sample_id, read1)
    }

    // Decompress each read and feed the resulting FASTQ into the graftM classification step.
    UNZIP_FASTQ(reads_ch)

    // For each sample and the selected genes, create one graftM task.
    // This expands the work into a matrix of sample x gene combinations.
    graftm_inputs = UNZIP_FASTQ.out.flatMap { sample_id, fastq ->
        selected_genes.collect { gene ->
            tuple(sample_id, fastq, gene, gene_packages[gene])
        }
    }

    // Launch the actual graftM analysis for each sample/gene combination.
    GRAFTM(graftm_inputs)
    }
}