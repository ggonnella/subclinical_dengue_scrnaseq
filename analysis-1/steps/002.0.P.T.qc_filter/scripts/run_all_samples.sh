#!/bin/bash

SAMPLES="sample_1_T sample_2_T sample_3_T sample_4_T"
SAMPLES="$SAMPLES sample_5_T sample_6_T sample_7_T sample_8_T"
SAMPLES="$SAMPLES sample_9_T sample_10_T sample_11_T sample_12_T"
SAMPLES="$SAMPLES sample_13_T sample_14_T "
n_threads=14

function run_sample { local SAMPLE=$1
  echo "Running sample $SAMPLE"
  ln -s -f filter_so.Rmd filter_so.$SAMPLE.Rmd
  knit2html filter_so.$SAMPLE.Rmd sample=$SAMPLE
  rm -f filter_so.$SAMPLE.Rmd
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

