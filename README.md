# gene-centric-metagenomics

Gene-centric metagenomics workflow for identifying and filtering target genes in metagenomic datasets using GraftM and downstream phylogenetic placement-based post-processing.

This repository contains two Nextflow workflows:

- `graftM.nf`: runs GraftM against metagenomic read sets for selected marker genes.
- `graftM_processing.nf`: post-processes GraftM placements to reduce false positives, filter candidate hits using HMMER, and produce per-sample and per-clade outputs.

The workflows are designed for gene-centric analyses of functional markers. For now it includes markers of functional guilds involved in inorganic nitrogen cycling:

- `nrfA` (nitrite ammonifiers)
- `nirK` (denitrifiers)
- `nirS` (denitrifiers)
- `nosZ` (nitrous oxide reducers)

## Overview

The repository is intended for metagenomic analysis where the goal is to detect the presence of specific functional genes across environmental samples using a phylogeny-based approach. The first workflow uses custom GraftM packages to search reads for gene homologs and generate placement files, while the second workflow performs quality control and removes likely false positives based on predefined clade lists and decoy HMMs. 

This is particularly useful when working with genes involved in nitrogen cycling or other functional marker genes, where gene tree placement can help distinguish true in-group hits from unrelated or misassigned sequences.

Related publications:
- https://doi.org/10.1038/s41467-023-44022-3
- https://doi.org/10.1038%2Fs41467-025-65319-5
- https://doi.org/10.1093%2Fismeco%2Fycag144
- https://doi.org/10.1093/ismeco/ycae020

## Repository contents

- `graftM.nf` — main GraftM execution workflow
- `graftM_processing.nf` — downstream filtering and placement processing workflow
- `README.md` — project overview and usage notes
- `databases` — folder containing GraftM packages, clade information and decoys files
- `environment.yml` — File to create a conda environment will all required programs. The path to the graftM_workflow environment can be specified with the --conda_env argument.

## Requirements

Before running the workflows, create a conda environment containing the following programs using the file `environment.yml`:
- GraftM (https://github.com/geronimp/graftM)
- pplacer (https://github.com/matsen/pplacer)
- GAPPA (https://github.com/lczech/gappa)
- HMMER (https://github.com/EddyRivasLab/hmmer)
- R

The scripts assume the corresponding GraftM packages and supporting files are available in a databases/ folder located in the same folder as the nextflow scripts. If needed, change the path in the parameter block. 

## Quick start

### 1) Run GraftM on metagenomic reads

```bash
nextflow run graftM.nf \
--path path/to/metagenomes \
--restrict_read_length 150 \
--genes nrfA,nirK,nirS,nosZ \
--threads 8 \
--conda_env /path/to/graftM_workflow \
-with-conda
```

This workflow looks for files matching `*R1.fastq.gz` under the supplied path, extracts them, and runs GraftM for each selected gene using the gene-specific package paths defined at the top of the script.

Notes:

- `--path` should point to a directory containing gzipped paired-end read files. GraftM is typically run on the forward reads only.
- `--genes` is a comma-separated list and is case-sensitive (no spaces).
- `--conda_env` should be adjusted to your local environment, e.g. ~/miniforge3/envs/graftM_workflow.
- The available genes are currently: `nrfA`, `nirK`, `nirS`, `nosZ`.

### 2) Post-process GraftM placements

```bash
nextflow run graftM_processing.nf \
  --results_dir graftM \
  --outdir graftM/graftM_processing \
  --genes nrfA,nirK,nirS,nosZ \
  --conda_env /path/to/graftM_workflow \
  --threads 8 \
  -with-conda
```

This workflow performs quality control, removes likely false positives based on predefined clade lists and decoy HMMs, and creates count tables. It uses GAPPA, Pplacer and HMMER.

## Workflow details

### `graftM.nf`

This script runs GraftM on each metagenome sample for each selected target gene. The paths to the gene packages are defined in the parameter block of the script and should be adjusted if needed (default: databases/graftM_packages/).

GraftM packages can be created with the create package command in graftM, see https://github.com/geronimp/graftM.

The workflow then publishes results under:

```text
graftM/<gene>/
```

### `graftM_processing.nf`

This second workflow is a filtering pipeline designed to refine GraftM hits. It:

- finds the most likely location in the phylogeny for each placement by accumulating the placement mass (likelihood weight ratio) of the placements of each read upwards (from tips to root), until the accumulated mass reaches the threshold. Used gappa edit accumulate and --threshold 0.95.
- creates candidate ingroup hits based on clade membership defined in the GraftM packages
- runs `hmmscan` against gene-specific decoy HMMs from KEGG
- keeps only sequences with expected functional KO assignments
- writes filtered jplace outputs and count tables

The script expects gene-specific clade lists and decoy HMM files, for example:

- `databases/clade-files/nrfA_ingroup_outgroup.txt`
- `databases/decoys/nrfA_K03385_decoys.hmm`

These paths are defined in the parameter block of the script and should be adjusted if needed.

## Expected outputs

After the processing workflow, you can expect output directories like:

```text
graftM/graftM_processing/nrfA/
graftM/graftM_processing/nirK/
```

Typical outputs include:

- `acc0.95_hmmscanfilter.jplace`
- `<gene>_hits_counts.tsv`
- `<gene>_counts.tsv`
- `acc0.95_filtered_jplace_split/*.jplace`

These files summarize the gene hits that pass the post-processing filters and can be used for downstream interpretation or visualization.

## Next step

Gappa (https://github.com/lczech/gappa/wiki) can be used for further analysis of jplace files.

## Notes

- The workflows are tuned for gene-centric marker analysis and are not a general metagenome assembly pipeline.
- Some paths (database files, GraftM packages, decoy HMMs, and pressed HMM indexes) are hard-coded and may need to be adapted to your local installation.
- For `nosZ`, the workflow includes special handling for clade I/II split analysis.

This workflow uses community tools for gene-centric metagenomics, including GraftM, GAPPA, pplacer, and HMMER.
