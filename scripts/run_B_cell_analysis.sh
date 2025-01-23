#!/usr/bin/env bash
#
# (c) 2024-2025 BAIA Unit, Institut Pasteur du Cambodge, Phnom Pehn, Cambodia
#

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
PRJ_DIR=$SCRIPT_DIR/..
ANALYSIS_DIR=$PRJ_DIR/analysis-1
STEPS_DIR=$ANALYSIS_DIR/steps

if [ "$1" == "" ]; then
  echo "Run all the steps of the B cell analysis" > /dev/stderr
  echo "" > /dev/stderr
  echo "Usage:" > /dev/stderr
  echo "  $0 <prefix>" > /dev/stderr
  echo "" > /dev/stderr
  echo "<prefix>: path to the project directory in the local system" > /dev/stderr
  echo "" > /dev/stderr
  exit 1
fi
PFX=$1

K2H=$SCRIPT_DIR/knit2html

for rmd in \
  106.0.P.B.merging.seurat/rundir/merging.Rmd \
  106.0.V.B.merged/rundir/post_merging.Rmd \
  \
  107.0.P.B.integration.harmony/rundir/harmony.Rmd \
  107.0.P.B.integration.harmony/rundir/post_processing.Rmd \
  107.0.V.B.integrated/rundir/post_harmony_analysis.Rmd \
  107.0.V.B.integrated/rundir/res_selection.Rmd \
  107.0.V.B.integrated/rundir/cmp_singleR_refs.Rmd \
  \
  108.0.P.B.split_plasmablasts/rundir/split.Rmd \
  108.0.V.B.splitted_plasmablasts/rundir/post_split_analysis.Rmd \
  \
  109.0.P.B.scMetabolism/rundir/run_scMetabolism.Rmd \
  109.0.V.B.metabolism_plots/rundir/scMetabolism_results.Rmd \
  \
  110.0.P.B.de_analysis/rundir/compute_de.Rmd \
  110.0.V.B.de_results/rundir/de_dotplots.Rmd \
  110.0.V.B.de_results/rundir/deseq_tables_and_vulcanos.Rmd \
  \
  111.0.P.B.gsea_analysis/rundir/compute_gsea.Rmd \
  111.0.V.B.gsea_results/rundir/gsea_go.Rmd \
  111.0.V.B.gsea_results/rundir/gsea_kegg.Rmd \
  \
  112.0.P.B.isolate_plasmablasts/rundir/isolate_plasmablasts.Rmd \
  112.0.V.B.isolated_plasmablasts/rundir/analysis_isolated_plasmablasts.Rmd \
  \
  113.0.P.B.BCR_analysis.platypus/rundir/platypus_create.Rmd \
  113.0.P.B.BCR_analysis.platypus/rundir/post_processing.Rmd \
  113.0.P.B.BCR_analysis.platypus/rundir/further_post_processing.Rmd \
  113.0.V.B.BCR_results/rundir/gene_usage.Rmd \
  113.0.V.B.BCR_results/rundir/shm_analysis.Rmd \
  113.0.V.B.BCR_results/rundir/cdr3_region.Rmd \
  113.0.V.B.BCR_results/rundir/clonotypes.Rmd \
  \
  114.0.P.B.enclone_preparation/rundir/output_barcode_table.Rmd \
  114.0.V.B.enclone_plots/rundir/create_enclone_plots.Rmd
do
  # if rundir does not exist, create it and link the scripts file
  rundir=$(dirname $STEPS_DIR/$rmd)
  if [ ! -d $rundir ]; then
    mkdir -p $rundir
    cd $rundir
    rmd_bn=$(basename $rmd)
    ln -s ../scripts/$rmd_bn .
    cd -
  fi
  # if html file does not exist, run the script
  rmd_pfx=$(basename $rmd .Rmd)
  if [ ! -f $rundir/$rmd_pfx.html ]; then
    $K2H $STEPS_DIR/$rmd pfx=$PFX
  fi
done
