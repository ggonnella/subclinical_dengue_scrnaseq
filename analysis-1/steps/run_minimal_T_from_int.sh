#!/bin/bash

THIS_SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
prjpath=$(readlink -f $THIS_SCRIPT_DIR/../../)/

if [ "$1" == "" ]; then
  libpath="/home/giorgio/R-pkg-Seurat5/"
elif [ "$#" -eq 1 ]; then
  libpath=$1
else
  echo "Usage: $0 [library_path]"
  echo "  library_path: path to the R library containing Seurat"
  exit 1
fi

knit2html 004.0.P.T.integration.harmony/rundir/harmony_post_processing.Rmd $pr input_integrated_rdata=TRUE
knit2html 004.0.V.T.integrated/rundir/post_harmony_analysis.Rmd $pr
knit2html 005.0.P.T.split_CD8EM/rundir/split.Rmd $pr
knit2html 005.0.V.T.splitted_CD8EM/rundir/post_split_analysis.Rmd $pr
knit2html 007.0.P.T.compute_de/rundir/compute_de.Rmd $pr
knit2html 007.0.V.T.de_results/rundir/de_dotplots.Rmd $pr
knit2html 007.0.V.T.de_results/rundir/deseq_tables_and_vulcanos.Rmd $pr
knit2html 007.0.V.T.de_results/rundir/selected_genes_heatmap.Rmd $pr
knit2html 009.0.P.T.isolate_CD8EM/rundir/isolate_effmemCD8.Rmd $pr
knit2html 009.0.V.T.isolated_CD8EM/rundir/analysis_isolated_CD8EM.Rmd $pr &
knit2html 010.0.P.T.TCR_analysis.platypus/rundir/platypus_create.Rmd $pr
knit2html 010.0.P.T.TCR_analysis.platypus/rundir/post_processing.Rmd $pr
knit2html 010.0.P.T.TCR_analysis.platypus/rundir/further_post_processing.Rmd $pr
knit2html 010.0.V.T.TCR_results/rundir/gene_usage.Rmd $pr &
knit2html 010.0.V.T.TCR_results/rundir/cdr3_region.Rmd $pr &
knit2html 010.0.V.T.TCR_results/rundir/clonotypes.Rmd $pr &
wait
knit2html 005.0.P.T.split_CD8EM/rundir/split.Rmd
knit2html 005.0.V.T.splitted_CD8EM/rundir/post_split_analysis.Rmd
knit2html 007.0.P.T.compute_de/rundir/compute_de.Rmd
knit2html 007.0.V.T.de_results/rundir/de_dotplots.Rmd
knit2html 007.0.V.T.de_results/rundir/deseq_tables_and_vulcanos.Rmd
knit2html 007.0.V.T.de_results/rundir/selected_genes_heatmap.Rmd
knit2html 009.0.P.T.isolate_CD8EM/rundir/isolate_effmemCD8.Rmd
knit2html 009.0.V.T.isolated_CD8EM/rundir/analysis_isolated_CD8EM.Rmd &
knit2html 010.0.P.T.TCR_analysis.platypus/rundir/platypus_create.Rmd
knit2html 010.0.P.T.TCR_analysis.platypus/rundir/post_processing.Rmd
knit2html 010.0.P.T.TCR_analysis.platypus/rundir/further_post_processing.Rmd
knit2html 010.0.V.T.TCR_results/rundir/gene_usage.Rmd &
knit2html 010.0.V.T.TCR_results/rundir/cdr3_region.Rmd &
knit2html 010.0.V.T.TCR_results/rundir/clonotypes.Rmd &
wait
