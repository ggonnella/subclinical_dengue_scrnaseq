#!/usr/bin/env bash
#
# (c) 2024-2025 BAIA Unit, Institut Pasteur du Cambodge, Phnom Pehn, Cambodia
#

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
PRJ_DIR=$SCRIPT_DIR/..
ANALYSIS_DIR=$PRJ_DIR/analysis-1
STEPS_DIR=$ANALYSIS_DIR/steps

if [ "$2" == "" ]; then
  echo "Run all the steps of the T cell analysis" > /dev/stderr
  echo "" > /dev/stderr
  echo "Usage:" > /dev/stderr
  echo "  $0 <prjpath> <libpath>" > /dev/stderr
  echo "" > /dev/stderr
  echo "<prjpath>: path to the project directory in the local system" > /dev/stderr
  echo "" > /dev/stderr
  echo "<libpath>: path to the R library containing Seurat 5" > /dev/stderr
  echo "" > /dev/stderr
  exit 1
fi
PRJPATH=$1
LIBPATH=$2

K2H=$SCRIPT_DIR/knit2html

for rundir in \
  001.0.P.T.create_seurat_obj/rundir/ \
  001.0.V.T.before_filtering/ \
  \
  002.0.P.T.qc_filter/rundir/ \
  002.0.V.T.qc_filtered/
do
  PREVDIR=$PWD
  # if rundir does not exist, create it and link the scripts file
  if [ ! -d $rundir ]; then
    mkdir -p $rundir
    cd $rundir
    ln -s ../scripts/* .
    cd -
  fi
  cd $rundir
  ./run_all_samples.sh $PRJPATH $LIBPATH
  cd $PREVDIR
done

for rmd in \
  003.0.P.T.merging.seurat/rundir/merging.Rmd \
  003.0.V.T.merged/rundir/post_merging.Rmd \
  \
  004.0.P.T.integration.harmony/rundir/harmony.Rmd \
  004.0.P.T.integration.harmony/rundir/harmony_post_processing.Rmd \
  004.0.V.T.integrated/rundir/post_harmony_analysis.Rmd \
  004.0.V.T.integrated/rundir/res_selection.Rmd \
  004.0.V.T.integrated/rundir/cmp_singleR_refs.Rmd \
  \
  005.0.P.T.split_CD8EM/rundir/split.Rmd \
  005.0.V.T.splitted_CD8EM/rundir/post_split_analysis.Rmd \
  \
  006.0.P.T.scMetabolism/rundir/run_scMetabolism.Rmd \
  006.0.V.T.metabolism_plots/rundir/scMetabolism.Rmd \
  \
  007.0.P.T.compute_de/rundir/compute_de.Rmd \
  007.0.V.T.de_results/rundir/de_dotplots.Rmd \
  007.0.V.T.de_results/rundir/selected_genes_heatmap.Rmd \
  007.0.V.T.de_results/rundir/deseq_tables_and_vulcanos.Rmd \
  \
  008.0.P.T.compute_gsea/rundir/compute_gsea.Rmd \
  008.0.V.T.gsea_results/rundir/gsea_go.Rmd \
  008.0.V.T.gsea_results/rundir/gsea_kegg.Rmd \
  \
  009.0.P.T.isolate_CD8EM/rundir/isolate_effmemCD8.Rmd \
  009.0.V.T.isolated_CD8EM/rundir/analysis_isolated_CD8EM.Rmd \
  \
  010.0.P.T.TCR_analysis.platypus/rundir/platypus_create.Rmd \
  010.0.P.T.TCR_analysis.platypus/rundir/post_processing.Rmd \
  010.0.P.T.TCR_analysis.platypus/rundir/further_post_processing.Rmd \
  010.0.V.T.TCR_results/rundir/gene_usage.Rmd \
  010.0.V.T.TCR_results/rundir/cdr3_region.Rmd \
  010.0.V.T.TCR_results/rundir/clonotypes.Rmd
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
    $K2H $STEPS_DIR/$rmd prjpath=$PRJPATH libpath=$LIBPATH
  fi
done

