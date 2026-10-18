# gene-centric-metagenomics

Gene-centric metagenomics workflow for identifying and filtering target genes in metagenomic datasets using GraftM and downstream placement-based post-processing.

This repository contains two Nextflow workflows:

- `graftM.nf`: runs GraftM against metagenomic read sets for selected marker genes.
- `graftM_processing.nf`: post-processes GraftM placements to reduce false positives, filter candidate hits using HMMER, and produce per-sample and per-clade outputs.

The workflows are designed for gene-centric analyses of functional markers such as:

- `nrfA`
- `nirK`
- `nirS`
- `nosZ`

## Overview

The repository is intended for metagenomic analysis where the goal is to detect the presence and placement of specific functional genes across environmental samples. The first workflow uses GraftM packages to search reads for gene homologs and generate placement files, while the second workflow performs quality control and removes likely false positives based on predefined clade lists and decoy HMMs.

This is particularly useful when working with genes involved in nitrogen cycling or other functional marker genes, where gene tree placement can help distinguish true in-group hits from unrelated or misassigned sequences.

## Repository contents

- `graftM.nf` — main GraftM execution workflow
- `graftM_processing.nf` — downstream filtering and placement processing workflow
- `README.md` — project overview and usage notes

## Requirements

Before running the workflows, make sure the following tools are installed and available in your environment:

- Nextflow
- GraftM
- pplacer
- GAPPA
- HMMER
- R

The scripts assume the corresponding GraftM packages and supporting files are available on disk, and some defaults are defined directly in the workflow parameters.

## Quick start

### 1) Run GraftM on metagenomic reads

```bash
nextflow run graftM.nf \
  --path path/to/metagenomes \
  --restrict_read_length 150 \
  --genes nrfA,nirK,nirS,nosZ
```

This workflow looks for files matching `*R1.fastq.gz` under the supplied path, extracts them, and runs GraftM for each selected gene using the gene-specific package paths defined at the top of the script.

Notes:

- `--path` should point to a directory containing gzipped paired-end read files.
- `--genes` is a comma-separated list and is case-sensitive.
- The available genes are currently: `nrfA`, `nirK`, `nirS`, `nosZ`.

### 2) Post-process GraftM placements

```bash
nextflow run graftM_processing.nf \
  --results_dir graftM \
  --outdir graftM/graftM_processing \
  --genes nrfA,nirK,nirS,nosZ \
  --threads 8 \
  -with-conda
```

This workflow does the following:

- accumulates placement files across samples
- merges and filters jplace files
- extracts ingroup candidates
- runs HMMER-based filtering against decoy HMMs
- removes false positives
- produces per-sample and per-clade outputs

## Workflow details

### `graftM.nf`

This script runs GraftM on each metagenome sample for each selected target gene. The gene packages are set in the parameter block and include paths such as:

- `nrfA_gpkg`
- `nirK_gpkg`
- `nirS_gpkg`
- `nosZ_gpkg`

The workflow then publishes results under:

```text
graftM/<gene>/
```

### `graftM_processing.nf`

This second workflow is a filtering pipeline designed to refine GraftM hits. It:

- validates query IDs from placement results
- finds all ORF FASTA sequences from GraftM outputs
- creates candidate ingroup hits based on clade membership
- runs `hmmscan` against gene-specific decoy HMMs
- keeps only sequences with expected functional KO assignments
- writes filtered jplace outputs and count tables

The script expects gene-specific clade lists and decoy HMM files, for example:

- `databases/nrfA_ingroup_outgroup.txt`
- `databases/nirK_ingroup_only.txt`
- `databases/nirS_ingroup_only.txt`
- `databases/nosZ_2021_clade_I_II_as_ingroup.txt`
- `decoys/nrfA_K03385_decoys.hmm`

These paths are defined near the top of the script and should be adjusted to match your environment.

## Expected outputs

After the processing workflow, you can expect output directories like:

```text
graftM/graftM_processing/nrfA/
GraftM/graftM_processing/nirK/
... 
```

Typical outputs include:

- `acc0.95_hmmscanfilter.jplace`
- `<gene>_hits_counts.tsv`
- `<gene>_counts.tsv`
- `acc0.95_filtered_jplace_split/*.jplace`

These files summarize the gene hits that pass the post-processing filters and can be used for downstream interpretation or visualization.

## Notes

- The workflows are tuned for gene-centric marker analysis and are not a general metagenome assembly pipeline.
- Some paths (database files, GraftM packages, decoy HMMs, and pressed HMM indexes) are hard-coded and may need to be adapted to your local installation.
- For `nosZ`, the workflow includes special handling for clade I/II split analysis.

## Example command summary

```bash
# Step 1: generate GraftM placements
nextflow run graftM.nf --path /data/metagenomes --restrict_read_length 150 --genes nrfA,nirK,nirS,nosZ

# Step 2: filter and summarize placements
nextflow run graftM_processing.nf \
  --results_dir graftM \
  --outdir graftM/graftM_processing \
  --genes nrfA,nirK,nirS,nosZ \
  --threads 8
```

## License

This project does not currently declare a license in the repository root. If you plan to share or redistribute the workflow, consider adding an explicit license file such as MIT or BSD.

## Acknowledgements

This workflow uses community tools for gene-centric metagenomics, including GraftM, GAPPA, pplacer, and HMMER.
