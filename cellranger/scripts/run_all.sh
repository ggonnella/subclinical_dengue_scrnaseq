#!/bin/bash
T_SAMPLES="1 2 3 4 5 6 7 8 9 10 11 12 13 14"
B_SAMPLES="$T_SAMPLES 15 16"

for sample in $T_SAMPLES; do
  ./run_cellranger_multi.$sample.sh
  cellranger multi --id sample_${sample}_T --csv ../metadata/config.S${sample}T.csv
done
for sample in $B_SAMPLES; do
  ./run_cellranger_multi.$sample.sh
  cellranger multi --id sample_${sample}_B --csv ../metadata/config.S${sample}B.csv
done
