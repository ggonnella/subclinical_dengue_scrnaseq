This repository documents the single cell RNA-seq analyses performed
for the manuscript "Gonnella et al., 2025".

# Requirements

The analyses are perfomed using:
- 10X Cellranger
- Seurat
- Platypus
- some further R libraries

# Organization

The repository is organized as follows:
- "analysis-1" original computation, analyses steps and results
- "helpers" helper scripts, containing functions common to multiple steps

Note: if a revision will be necessary, "analysis-2" will be created for
further analyses steps and results.

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


