#!/bin/bash

THISSCRIPTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
knit2html=$THISSCRIPTDIR/../../../../scripts/knit2html

SAMPLES="sample_1_B sample_3_B sample_4_B"
SAMPLES="$SAMPLES sample_5_B sample_6_B sample_7_B sample_8_B"
SAMPLES="$SAMPLES sample_9_B sample_10_B sample_11_B sample_12_B"
SAMPLES="$SAMPLES sample_13_B sample_15_B sample_16_B"
n_threads=14

if [ "$2" != "" ]; then PSTR="prjpath=$1 libpath=$2"; else PSTR=""; fi

function run_sample { local SAMPLE=$1
  echo "Running sample $SAMPLE"
  ln -s -f before_filtering.Rmd before_filtering.$SAMPLE.Rmd
  $knit2html before_filtering.$SAMPLE.Rmd sample=$SAMPLE $PSTR
  rm -f before_filtering.$SAMPLE.Rmd
  echo "Finished sample $SAMPLE"
}
export -f run_sample

NSAMPLES=$(echo $SAMPLES | wc -w)
if command -v parallel >/dev/null 2>&1; then
  echo "Running $NSAMPLES samples in parallel using $n_threads threads"
  parallel --line-buffer -j $n_threads run_sample ::: $SAMPLES
else
  echo "Running $NSAMPLES samples in serial as GNU parallel is not installed"
  for SAMPLE in $SAMPLES; do run_sample $SAMPLE; done
fi

