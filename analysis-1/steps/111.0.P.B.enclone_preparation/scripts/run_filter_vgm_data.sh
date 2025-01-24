#!/bin/bash

PARAMS_PFX="/srv/baia"
PRJ_ROOT="$PARAMS_PFX/inden"
ANALYSIS_ROOT="$PRJ_ROOT/analysis-1"
INPUT_STEP_DIR="$ANALYSIS_ROOT/steps/111.0.P.B.enclone_preparation"
INPUT_TABLES_DIR="$INPUT_STEP_DIR/results/tables"
if [ ! -d "$INPUT_TABLES_DIR" ]; then
  echo "Input tables directory does not exist: $INPUT_TABLES_DIR"
  exit 1
fi
for file in vgm.tsv; do
  if [ ! -f "${INPUT_TABLES_DIR}/$file" ]; then
    echo "VGM file not found: ${INPUT_TABLES_DIR}/$file"
    exit 1
  fi
done

THIS_STEP_DIR="$ANALYSIS_ROOT/steps/111.0.P.B.enclone_preparation"
if [ ! -d "$THIS_STEP_DIR" ]; then
  echo "Current step directory does not exist: $THIS_STEP_DIR"
  exit 1
fi
SCRIPTS_DIR="$THIS_STEP_DIR/scripts"
JSON_DIR="$THIS_STEP_DIR/results/json"
mkdir -p "$JSON_DIR"

METADATA_DIR="$THIS_STEP_DIR/metadata"
if [ ! -d "$METADATA_DIR" ]; then
  echo "Metadata directory not found: $METADATA_DIR"
  exit 1
fi
for file in config.tsv; do
  if [ ! -f "${METADATA_DIR}/$file" ]; then
    echo "Metadata file not found: ${METADATA_DIR}/$file"
    exit 1
  fi
done

echo "Creating contig annotations files"
OUTDIR=$JSON_DIR
${SCRIPTS_DIR}/filter_vgm_data.py "${INPUT_TABLES_DIR}/vgm.tsv" \
                     "${METADATA_DIR}/config.tsv" \
                     --verbose \
                     > "${OUTDIR}/all_contig_annotations.json"
echo "Output to file: ${OUTDIR}/all_contig_annotations.json"
