#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// This workflow post-processes graftM output for a set of target genes.
// It starts from per-gene result directories, filters false positives using HMMER,
// and then splits the remaining placements into per-sample jplace files for downstream
// phylogenetic interpretation. The steps are intentionally staged so that each process
// validates intermediate outputs before moving on to the next stage.

params.results_dir = 'graftM'
params.outdir = 'graftM/graftM_processing'
params.conda_env = "${System.getProperty('user.home')}/miniforge3/envs/graftM_workflow"
params.genes = null
params.threads = 1
params.nrfA_clade_list = "${launchDir}/databases/clade-files/nrfA_ingroup_outgroup.txt"
params.nirK_clade_list = "${launchDir}/databases/clade-files/nirK_ingroup_only.txt"
params.nirS_clade_list = "${launchDir}/databases/clade-files/nirS_ingroup_only.txt"
params.nosZ_clade_list = "${launchDir}/databases/clade-files/nosZ_2021_clade_I_II_as_ingroup.txt"
params.nosZ_split_clade_list = "${launchDir}/databases/clade-files/nosZ_2021_clade_I_II.txt"
params.nrfA_decoy = "${launchDir}/databases/decoys/nrfA_K03385_decoys.hmm"
params.nirK_decoy = "${launchDir}/databases/decoys/nirK_K00368_decoys.hmm"
params.nirS_decoy = "${launchDir}/databases/decoys/nirS_K15864_decoys.hmm"
params.nosZ_decoy = "${launchDir}/databases/decoys/nosZ_K00376_decoys.hmm"
params.help = false

// Step 1: aggregate per-sample graftM placements across all samples for one gene,
// then keep only the relevant ingroup hits before writing a merged jplace file.
// This stage converts the raw graftM outputs into a cleaned set of candidate placements
// that can be checked against known gene-specific HMMs in the next step.
process PROCESS_GENE_PLACEMENTS {
    conda params.conda_env

    tag "$gene"
    cpus params.threads.toString().toInteger()

    input:
    tuple val(gene), path(gene_dir, stageAs: 'graftm_results'), path(clade_list, stageAs: 'clade_list.txt')

    output:
    tuple val(gene), path('step1'), emit: results

    script:
    """
    set -euo pipefail

    if [[ ! -s clade_list.txt ]]; then
        echo "Clade list is missing or empty for ${gene}: clade_list.txt" >&2
        exit 1
    fi

    shopt -s nullglob
    sample_dirs=(graftm_results/*/)
    if (( \${#sample_dirs[@]} == 0 )); then
        echo "No sample directories found in graftM results for ${gene}" >&2
        exit 1
    fi

    for sample_dir in "\${sample_dirs[@]}"; do
        sample_id=\$(basename "\${sample_dir%/}")
        if [[ "\$sample_id" == "step1" || "\$sample_id" == "step2" || "\$sample_id" == "step3" ]]; then
            continue
        fi
        placement=\$(find "\$sample_dir" -type f -name placements.jplace -print -quit)
        if [[ -z "\$placement" ]]; then
            echo "No placements.jplace found for sample \$sample_id (${gene})" >&2
            exit 1
        fi
        gappa edit accumulate \\
            --jplace-path "\$placement" \\
            --threshold 0.95 \\
            --out-dir step1/jplace_files \\
            --file-prefix "acc0.95_\${sample_id}_" \\
            --threads ${task.cpus}
    done

    cd step1
    accumulated=(jplace_files/acc0.95_*.jplace)
    if (( \${#accumulated[@]} == 0 )); then
        echo "GAPPA did not produce any accumulated jplace files for ${gene}" >&2
        exit 1
    fi

    gappa edit merge --jplace-path "\${accumulated[@]}" --out-dir merged_jplace --threads ${task.cpus}

    if [[ "${gene}" == "nrfA" ]]; then
        (
            cd merged_jplace
            guppy to_csv -o merge.csv --point-mass merge.jplace
            Rscript -e 'infile <- read.csv("merge.csv"); subfile <- subset(infile, !edge_num %in% c(0:146, 2353)); write.table(subfile["name"], file = "ingroup_merge.txt", sep = "\\t", row.names = FALSE, col.names = FALSE, quote = FALSE)'

            if [[ ! -s ingroup_merge.txt ]]; then
                echo "No ingroup query names were selected from nrfA merge.csv" >&2
                exit 1
            fi

            gappa edit filter --keep-names ingroup_merge.txt --jplace-path merge.jplace --threads ${task.cpus}
            if [[ ! -s filter.jplace ]]; then
                echo "GAPPA did not produce merged_jplace/filter.jplace for nrfA" >&2
                exit 1
            fi
            mv filter.jplace merge_filtered.jplace
        )
    fi

    if [[ "${gene}" == "nrfA" ]]; then
        if [[ ! -s merged_jplace/merge_filtered.jplace ]]; then
            echo "GAPPA did not produce merged_jplace/merge_filtered.jplace for nrfA" >&2
            exit 1
        fi
        guppy to_csv -o hits_acc0.95_before_hmmscan.csv --point-mass merged_jplace/merge_filtered.jplace
    else
        gappa edit extract --jplace-path merged_jplace/ --clade-list-file ../clade_list.txt --point-mass --threads ${task.cpus}

        if [[ ! -s samples/ingroup.jplace ]]; then
            echo "GAPPA did not produce samples/ingroup.jplace for ${gene}" >&2
            exit 1
        fi
        if [[ "${gene}" == "nosZ" ]]; then
            sed 's/-nan/0/g' samples/ingroup.jplace > samples/ingroup.jplace.tmp
            mv samples/ingroup.jplace.tmp samples/ingroup.jplace
        fi
        guppy to_csv -o hits_acc0.95_before_hmmscan.csv --point-mass samples/ingroup.jplace
    fi

    if [[ ! -s hits_acc0.95_before_hmmscan.csv ]]; then
        echo "Guppy did not produce hits_acc0.95_before_hmmscan.csv for ${gene}" >&2
        exit 1
    fi
    awk -F ',' 'NR > 1 { print \$2 }' hits_acc0.95_before_hmmscan.csv > hits_acc0.95_before_hmmscan_ids.txt

    """
}

// Step 2: remove likely false positives from the candidate hit set.
// The workflow gathers all ORFs from the graftM outputs, retains only ingroup IDs,
// scans them against a decoy HMM, and keeps only sequences that match the expected KO.
// The resulting jplace file is written out for subsequent splitting and reporting.
process FILTER_FALSE_POSITIVES {
    conda params.conda_env

    tag "$gene"
    cpus params.threads.toString().toInteger()

    publishDir path: "${params.outdir}/${gene}", mode: 'copy', saveAs: { filename ->
        filename.startsWith('publish_files/') ? filename.tokenize('/').last() : null
    }

    input:
    tuple val(gene), path(step1_dir, stageAs: 'step1'), path(gene_dir, stageAs: 'graftm_results'), path(decoy), path(decoy_h3f), path(decoy_h3i), path(decoy_h3m), path(decoy_h3p), val(target_ko)

    output:
    tuple val(gene), path('step2'), emit: results
    path('publish_files/*'), emit: published

    script:
    """
    set -euo pipefail

    mkdir -p step2

    if [[ ! -s step1/hits_acc0.95_before_hmmscan_ids.txt ]]; then
        echo "No ingroup query IDs found for ${gene}: step1/hits_acc0.95_before_hmmscan_ids.txt" >&2
        exit 1
    fi

    if [[ ! -s step1/merged_jplace/merge.jplace ]]; then
        echo "Merged jplace file is missing for ${gene}: step1/merged_jplace/merge.jplace" >&2
        exit 1
    fi

    find -L graftm_results -type f -name '*_orf.fa' -print > step2/orf_files.txt
    mapfile -t orf_files < step2/orf_files.txt
    if (( \${#orf_files[@]} == 0 )); then
        echo "No *_orf.fa files found in GraftM results for ${gene}" >&2
        exit 1
    fi
    cat "\${orf_files[@]}" > step2/hits_all.fa
    rm step2/orf_files.txt
    hits_all_count=\$(awk '/^>/ { count++ } END { print count + 0 }' step2/hits_all.fa)

    awk '
        {
            gsub(sprintf("%c", 13), "")
            if (\$1 != "") print \$1
        }
    ' step1/hits_acc0.95_before_hmmscan_ids.txt > step2/ingroup_ids.txt
    if [[ ! -s step2/ingroup_ids.txt ]]; then
        echo "No usable ingroup query IDs found for ${gene}" >&2
        exit 1
    fi

    awk '
        FNR == NR { ids[\$1] = 1; next }
        substr(\$0, 1, 1) == ">" {
            id = substr(\$0, 2)
            sub("[[:space:]].*", "", id)
            keep = id in ids
        }
        keep { print }
    ' step2/ingroup_ids.txt step2/hits_all.fa > step2/hits_acc0.95_before_hmmscan.fa
    if [[ ! -s step2/hits_acc0.95_before_hmmscan.fa ]]; then
        echo "No ingroup query sequences matched the GraftM ORF FASTA files for ${gene}" >&2
        exit 1
    fi
    ingroup_hits_count=\$(awk '/^>/ { count++ } END { print count + 0 }' step2/hits_acc0.95_before_hmmscan.fa)

    hmmscan --cpu ${task.cpus} -E 0.00001 \\
        --tblout step2/hits.hmmscan \\
        "${decoy}" \\
        step2/hits_acc0.95_before_hmmscan.fa

    awk -v target_ko="${target_ko}" '
        substr(\$0, 1, 1) == "#" { next }
        !seen[\$3]++ && \$1 == target_ko { print \$3 }
    ' step2/hits.hmmscan > step2/hits_true_ids.txt
    true_hits_count=\$(awk 'END { print NR + 0 }' step2/hits_true_ids.txt)

    printf 'gene\\thits_all_fa\\thits_acc0.95_before_hmmscan_fa\\thits_true_ids\\n%s\\t%s\\t%s\\t%s\\n' \\
        "${gene}" "\$hits_all_count" "\$ingroup_hits_count" "\$true_hits_count" \\
        > step2/hit_counts.tsv

    gappa edit filter \\
        --jplace-path step1/merged_jplace/merge.jplace \\
        --keep-names step2/hits_true_ids.txt \\
        --out-dir step2/acc0.95_hmmscan_filtered \\
        --file-prefix acc0.95_hmmscan \\
        --threads ${task.cpus}

    if [[ ! -s step2/acc0.95_hmmscan_filtered/acc0.95_hmmscanfilter.jplace ]]; then
        echo "GAPPA did not produce the HMM-filtered jplace file for ${gene}" >&2
        exit 1
    fi

    mkdir -p publish_files
    cp step2/acc0.95_hmmscan_filtered/acc0.95_hmmscanfilter.jplace \\
        publish_files/acc0.95_hmmscanfilter.jplace
    cp step2/hit_counts.tsv "publish_files/${gene}_hits_counts.tsv"
    """
}

// Step 3: split the filtered placements back into sample-specific jplace files.
// This allows each sample to be counted and interpreted independently, and also handles
// nosZ-specific clade splitting (Clade I vs Clade II) when requested.
process SPLIT_FILTERED_PLACEMENTS {
    conda params.conda_env

    tag "$gene"
    cpus params.threads.toString().toInteger()

    publishDir path: "${params.outdir}/${gene}", mode: 'copy', saveAs: { filename ->
        if (filename.startsWith('publish_files/') && filename.endsWith('_counts.tsv')) {
            filename.tokenize('/').last()
        } else if (filename.endsWith('.jplace')) {
            "acc0.95_filtered_jplace_split/${filename.tokenize('/').last()}"
        } else {
            null
        }
    }

    input:
    tuple val(gene), path(step2_dir, stageAs: 'step2'), path(gene_dir, stageAs: 'graftm_results'), path(split_clade_list, stageAs: 'nosZ_clade_list.txt')

    output:
    tuple val(gene), path('step3'), emit: results
    path('publish_files/*_counts.tsv'), emit: published
    path('publish_files/acc0.95_filtered_jplace_split/*.jplace'), optional: true, emit: split_files

    script:
    """
    set -euo pipefail

    mkdir -p step3/acc0.95_filtered_jplace_split

    if [[ ! -f step2/hits_true_ids.txt ]]; then
        echo "HMM-validated query ID file is missing for ${gene}: step2/hits_true_ids.txt" >&2
        exit 1
    fi
    if [[ ! -s step2/acc0.95_hmmscan_filtered/acc0.95_hmmscanfilter.jplace ]]; then
        echo "HMM-filtered jplace file is missing for ${gene}" >&2
        exit 1
    fi

    find -L graftm_results -type f -name '*_orf.fa' -print | sort > step3/all_fasta_paths.txt
    if [[ ! -s step3/all_fasta_paths.txt ]]; then
        echo "No *_orf.fa files found in GraftM results for ${gene}" >&2
        exit 1
    fi

    : > step3/hits_true_ids_sample_mapping.tsv
    : > step3/sample_ids.txt
    while IFS= read -r fasta; do
        sample_id=\$(basename "\$(dirname "\$fasta")")
        printf '%s\\n' "\$sample_id" >> step3/sample_ids.txt
        awk -v sample_id="\$sample_id" '
            BEGIN {
                while ((getline id < "step2/hits_true_ids.txt") > 0) ids[id] = 1
                close("step2/hits_true_ids.txt")
            }
            substr(\$0, 1, 1) == ">" {
                id = substr(\$0, 2)
                sub("[[:space:]].*", "", id)
                if (id in ids) print id "," sample_id
            }
        ' "\$fasta" >> step3/hits_true_ids_sample_mapping.tsv
    done < step3/all_fasta_paths.txt
    sort -u step3/sample_ids.txt -o step3/sample_ids.txt
    sort -u step3/hits_true_ids_sample_mapping.tsv -o step3/hits_true_ids_sample_mapping.tsv

    while IFS= read -r id; do
        [[ -z "\$id" ]] && continue
        match_count=\$(awk -F ',' -v id="\$id" '\$1 == id { count++ } END { print count + 0 }' step3/hits_true_ids_sample_mapping.tsv)
        if [[ "\$match_count" -ne 1 ]]; then
            echo "Expected query \$id to match exactly one sample FASTA for ${gene}; found \$match_count" >&2
            exit 1
        fi
    done < step2/hits_true_ids.txt

    mapping_count=\$(awk 'END { print NR + 0 }' step3/hits_true_ids_sample_mapping.tsv)
    filtered_jplace=step2/acc0.95_hmmscan_filtered/acc0.95_hmmscanfilter.jplace
    if [[ "${gene}" == "nosZ" ]]; then
        sed 's/-nan/0/g' "\$filtered_jplace" > step3/acc0.95_hmmscanfilter.jplace
        filtered_jplace=step3/acc0.95_hmmscanfilter.jplace
    fi

    if [[ "\$mapping_count" -gt 0 ]]; then
        gappa edit split \\
            --jplace-path "\$filtered_jplace" \\
            --split-file step3/hits_true_ids_sample_mapping.tsv \\
            --out-dir step3/acc0.95_filtered_jplace_split \\
            --threads ${task.cpus}

        shopt -s nullglob
        split_files=(step3/acc0.95_filtered_jplace_split/*.jplace)
        if (( \${#split_files[@]} == 0 )); then
            echo "GAPPA did not produce per-sample jplace files for ${gene}" >&2
            exit 1
        fi
    fi

    printf 'sample\\thit_count\\n' > step3/counts.tsv
    while IFS= read -r sample_id; do
        hit_count=\$(awk -v sample_id="\$sample_id" -F ',' '\$2 == sample_id { count++ } END { print count + 0 }' step3/hits_true_ids_sample_mapping.tsv)
        printf '%s\\t%s\\n' "\$sample_id" "\$hit_count" >> step3/counts.tsv
    done < step3/sample_ids.txt

    if [[ "${gene}" == "nosZ" ]]; then
        if [[ ! -s nosZ_clade_list.txt ]]; then
            echo "nosZ Clade I/II list is missing or empty" >&2
            exit 1
        fi

        printf 'jplace_file\\thit_count\\n' > step3/nosZ_I_counts.tsv
        printf 'jplace_file\\thit_count\\n' > step3/nosZ_II_counts.tsv
        if [[ "\$mapping_count" -gt 0 ]]; then
            gappa edit extract \\
                --jplace-path "\$filtered_jplace" \\
                --clade-list-file nosZ_clade_list.txt \\
                --point-mass \\
                --samples-out-dir step3/acc0.95_filtered_split \\
                --threads ${task.cpus}

            for clade in CladeI CladeII; do
                clade_jplace="step3/acc0.95_filtered_split/\${clade}.jplace"
                if [[ ! -s "\$clade_jplace" ]]; then
                    echo "GAPPA did not produce \$clade_jplace for nosZ" >&2
                    exit 1
                fi
                sed 's/-nan/0/g' "\$clade_jplace" > "\${clade_jplace}.tmp"
                mv "\${clade_jplace}.tmp" "\$clade_jplace"
                clade_name=\$(printf '%s' "\$clade" | tr '[:upper:]' '[:lower:]')
                clade_out_dir="step3/acc0.95_filtered_jplace_clade_\${clade_name#clade}_split"
                gappa edit split \\
                    --jplace-path "\$clade_jplace" \\
                    --split-file step3/hits_true_ids_sample_mapping.tsv \\
                    --out-dir "\$clade_out_dir" \\
                    --threads ${task.cpus}

                for split_file in "\$clade_out_dir"/*.jplace; do
                    [[ -f "\$split_file" ]] || continue
                    hit_count=\$(awk '/"n"[[:space:]]*:/ { count++ } END { print count + 0 }' "\$split_file")
                    printf '%s\\t%s\\n' "\$(basename "\$split_file")" "\$hit_count" \\
                        >> "step3/nosZ_\${clade_name#clade}_counts.tsv"
                done
            done
        fi
    fi

    mkdir -p publish_files
    cp step3/counts.tsv "publish_files/${gene}_counts.tsv"
    mkdir -p publish_files/acc0.95_filtered_jplace_split
    if [[ "\$mapping_count" -gt 0 ]]; then
        cp step3/acc0.95_filtered_jplace_split/*.jplace publish_files/acc0.95_filtered_jplace_split/
    fi
    """
}

workflow {
    // Show the command-line usage and required inputs before running any work.
    if (params.help) {
        log.info """
        The specified conda environment must contain all programs used by this workflow.

        Usage:
          nextflow run graftM_step1.nf \\
            --results_dir graftM \\
            --outdir graftM/graftM_processing \\
            --genes nrfA,nirK,nirS,nosZ \\
            --threads 8 \\
            --conda_env /path/to/graftM_workflow \\
            -with-conda

        Required inputs:
          --results_dir       Directory containing one graftM results directory per gene
          --outdir            Directory for flat per-gene processing outputs (default: graftM/graftM_processing)
          --genes             Comma-separated subset of nrfA,nirK,nirS,nosZ
          --threads           CPU threads used by each process (default: 1)
          --conda_env         Path to the conda environment containing all required programs (default: ~/miniforge3/envs/graftM_workflow)

        Gene-specific clade-list paths are set near the top of this script.
        Gene-specific decoy HMMs and their HMMER pressed indexes are also set there.
        nosZ Clade I/II splitting uses --nosZ_split_clade_list.

        Outputs:
          <outdir>/<gene>/acc0.95_hmmscanfilter.jplace
          <outdir>/<gene>/<gene>_hits_counts.tsv
          <outdir>/<gene>/<gene>_counts.tsv
          <outdir>/<gene>/acc0.95_filtered_jplace_split/*.jplace

        Gappa, guppy, HMMER and Rscript must be available in the specified environment.
        """.stripIndent()
    } else {
        // Validate the thread count and gene-specific configuration before submission.
        if (!(params.threads.toString() ==~ /[1-9][0-9]*/)) {
            error "The --threads value must be a positive integer; received: ${params.threads}"
        }
        if (!params.conda_env) {
            error 'Specify the path to the shared conda environment with --conda_env'
        }
        def conda_env_path = file(params.conda_env.toString())
        if (!conda_env_path.toFile().isDirectory()) {
            error "Conda environment directory does not exist: ${conda_env_path}"
        }

        // Each gene has its own clade definition, decoy HMM, and expected KO target.
        // These mappings are used to validate inputs and to pass the right files into
        // every downstream process.
        def clade_params = [
            nrfA: params.nrfA_clade_list,
            nirK: params.nirK_clade_list,
            nirS: params.nirS_clade_list,
            nosZ: params.nosZ_clade_list
        ]
        def decoy_params = [
            nrfA: params.nrfA_decoy,
            nirK: params.nirK_decoy,
            nirS: params.nirS_decoy,
            nosZ: params.nosZ_decoy
        ]
        def target_kos = [
            nrfA: 'K03385',
            nirK: 'K00368',
            nirS: 'K15864',
            nosZ: 'K00376'
        ]

        if (!params.genes) {
            error 'Specify one or more genes with --genes, for example --genes nrfA,nirK'
        }

        // Parse the user-supplied gene list into a clean set of supported names.
        def genes = params.genes
            .split(',')
            .collect { it.trim() }
            .findAll { it }
            .unique()

        if (!genes) {
            error 'Specify at least one gene with --genes, for example --genes nrfA'
        }

        def unknown_genes = genes.findAll { !clade_params.containsKey(it) }
        if (unknown_genes) {
            error "Unknown gene(s): ${unknown_genes.join(', ')}. Available genes: ${clade_params.keySet().join(', ')}"
        }

        // Prepare the input tuples for the first processing step.
        // Each tuple contains the gene name, the directory of graftM results, and the
        // gene-specific clade list used to select ingroup placements.
        def gene_inputs = genes.collect { gene ->
            def gene_path = file("${params.results_dir}/${gene}")
            def clade_path = clade_params[gene] ? file(clade_params[gene].toString()) : null
            def decoy_path = decoy_params[gene] ? file(decoy_params[gene].toString()) : null

            if (!gene_path.toFile().exists() || !gene_path.toFile().isDirectory()) {
                error "GraftM results directory does not exist: ${gene_path}"
            }
            if (!clade_path || !clade_path.toFile().exists() || !clade_path.toFile().isFile()) {
                error "Provide an existing --${gene}_clade_list file for gene ${gene}"
            }
            if (!decoy_path || !decoy_path.toFile().exists() || !decoy_path.toFile().isFile()) {
                error "Provide an existing --${gene}_decoy HMM file for gene ${gene}"
            }

            tuple(gene, gene_path, clade_path)
        }

        // Prepare the second-stage inputs, including the decoy HMM and its required
        // pressed-index files for hmmscan-based filtering.
        def gene_filter_inputs = genes.collect { gene ->
            def gene_path = file("${params.results_dir}/${gene}")
            def decoy_path = file(decoy_params[gene].toString())
            def decoy_indexes = ['.h3f', '.h3i', '.h3m', '.h3p'].collect { suffix ->
                file(decoy_path.toString() + suffix)
            }
            if (!decoy_indexes.every { it.toFile().exists() && it.toFile().isFile() }) {
                error "HMMER pressed index files are missing for decoy HMM ${decoy_path} (expected .h3f, .h3i, .h3m and .h3p files)"
            }
            tuple(
                gene,
                gene_path,
                decoy_path,
                decoy_indexes[0],
                decoy_indexes[1],
                decoy_indexes[2],
                decoy_indexes[3],
                target_kos[gene]
            )
        }

        // Run the first post-processing stage on each gene independently.
        def gene_inputs_ch = Channel.fromList(gene_inputs)
        PROCESS_GENE_PLACEMENTS(gene_inputs_ch)

        // Join the step-1 outputs with the filtering inputs so that each gene enters the
        // HMM-based false positive filter with the right reference files.
        def step2_inputs = PROCESS_GENE_PLACEMENTS.out.results
            .join(Channel.fromList(gene_filter_inputs))

        FILTER_FALSE_POSITIVES(step2_inputs)

        // Build the final inputs for sample splitting, using the nosZ-specific clade list
        // when needed and otherwise reusing the normal gene clade list.
        def step3_gene_inputs = genes.collect { gene ->
            def gene_path = file("${params.results_dir}/${gene}")
            def split_clade_path = gene == 'nosZ'
                ? file(params.nosZ_split_clade_list.toString())
                : file(clade_params[gene].toString())

            if (gene == 'nosZ' && (!split_clade_path.toFile().exists() || !split_clade_path.toFile().isFile())) {
                error "Provide an existing --nosZ_split_clade_list file for nosZ"
            }

            tuple(gene, gene_path, split_clade_path)
        }
        // Finally, split the cleaned placements into per-sample outputs and, for nosZ,
        // also split by Clade I/II classification.
        def step3_inputs = FILTER_FALSE_POSITIVES.out.results
            .join(Channel.fromList(step3_gene_inputs))

        SPLIT_FILTERED_PLACEMENTS(step3_inputs)
    }
}
