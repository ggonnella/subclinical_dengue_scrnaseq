#!/usr/bin/env bash

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
PRJ_DIR=$SCRIPT_DIR/..
ANALYSIS_DIR=$PRJ_DIR/analysis-1
STEPS_DIR=$ANALYSIS_DIR/steps

K2H=$SCRIPT_DIR/knit2html

$K2H $STEPS_DIR/006.0.P.T.merging.seurat/rundir/merging.Rmd
$K2H $STEPS_DIR/006.0.V.T.post_merging/rundir/post_merging.Rmd

$K2H $STEPS_DIR/007.0.V.T.post_integration.harmony/rundir/post_harmony_analysis.Rmd
$K2H $STEPS_DIR/007.0.P.T.integration.harmony/rundir/harmony.Rmd
$K2H $STEPS_DIR/007.0.P.T.integration.harmony/rundir/post_processing.Rmd

$K2H $STEPS_DIR/008.0.V.T.post_split_CD8EM/rundir/post_split_analysis.Rmd
$K2H $STEPS_DIR/008.0.P.T.split_CD8EM/rundir/split.Rmd
