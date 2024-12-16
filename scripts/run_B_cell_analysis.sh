#!/usr/bin/env bash

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
PRJ_DIR=$SCRIPT_DIR/..
ANALYSIS_DIR=$PRJ_DIR/analysis-1
STEPS_DIR=$ANALYSIS_DIR/steps

K2H=$SCRIPT_DIR/knit2html

$K2H $STEPS_DIR/106.0.P.B.merging.seurat/rundir/merging.Rmd
$K2H $STEPS_DIR/106.0.V.B.post_merging/rundir/post_merging_analysis.Rmd
$K2H $STEPS_DIR/106.0.V.B.post_merging/rundir/post_merging_res_selection.Rmd
