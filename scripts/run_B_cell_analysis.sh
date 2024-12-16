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
  106.0.V.B.post_merging/rundir/post_merging_analysis.Rmd \
  106.0.V.B.post_merging/rundir/post_merging_res_selection.Rmd
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
