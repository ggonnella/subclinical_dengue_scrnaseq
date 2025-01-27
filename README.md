# Requirements

The analyses are perfomed using:
- R version 4.4.2
- 10X Cellranger v6 or v7
- 10X enclone v0.5.219
- mixcr v3.0.12
- Seurat v5, Platypus and several further R libraries (see ``R_LIBRARIES`` file)

# Organization

The repository is organized as follows:
- "cellranger" are the scripts/metadata used for running Cellranger
- "analysis-1" processing and visualization steps and results
- "helpers" helper scripts, containing functions common to multiple steps
- "scripts" command line scripts (e.g. to run all steps)

# Processing and Visualization Steps

Each step under "analysis-N" is classified as a "processing step" or
"visualization step".

Thereby:
- processing steps take the previous step data and generate
  new further processed data (e.g. filtering, normalization, clustering, etc.);
  for most steps (R), the data is saved in the ``rds`` format under
  ``results/vars``
- visualization steps generate plots (in ``results/plots``) and tables
  (in ``results/tables``) but do not change the underlying data structure.

# Steps naming

Each step is named as follows
``stepnumber.versionnumber.type.celltype.description``
matching the following regular expression:
``(\d\d\d)\.(\d)\.([PV])\.([BT])\.(.*)``

where:
- stepnumber is the incremental number of the step,
  starting from 001 for the T cells and 101 for the B cells analysis
- versionnumber is a version of the step, if the step is revised
- type is either P for processing or V for visualization
- celltype is either B for B cells or T for T cells
- description is a short textual description of the step

# Step organization

Each step is organized as follows:
- ``scripts`` contains the scripts (in most cases, R markdown) for the step
- ``results`` contains the results of the step
- ``rundir`` contains the output of the step (e.g. logs, intermediate files)

The scripts under scripts are symlinked into rundir and run there. The scripts
create the results directory and store the results there.

# Running the pipeline

## From reads to counts: Cellranger

To run Cellranger, change the ``cellranger/metadata`` files to match the paths
to the fastq files and the reference genome. Then run the
``cellranger/scripts/run_all.sh`` script.

In alternative, Cellranger counts from GEO can be put in the
``cellranger/results`` directory to run the pipeline from that point on.

## Upstream analysis (individual samples)

To run the upstream analysis, one case use the ``run_all_samples.sh`` scripts
in the 001/002 (T cells) and 101/102 (B cells) P and V steps.
Thereby it is necessary to specify the own paths to this repository
(prjpath parameter) and to R library containing Seurat 5 (libpath).

## Downstream analysis (merged samples)

To run the Rmd files, one can use the knit2html script under scripts.
Thereby it is necessary to specify the own paths to this repository
(prjpath parameter) and to R library containing Seurat 5 (libpath).

## Running the pipeline (upstream and downstream)

The ``run_[TB]_cell_analysis.sh`` script run all the steps for the
T and B analysis starting from the Cellranger output.
Thereby it is necessary to specify the own paths to this repository
(prjpath parameter) and to R library containing Seurat 5 (libpath).

## Collecting results

The results of the pipeline for the manuscript figures and supplementary
figures can be collected using the ``collect_results.sh`` script under
``analysis-1/results``.
