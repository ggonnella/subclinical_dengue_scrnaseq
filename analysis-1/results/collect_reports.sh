#!/bin/bash

THIS_SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
cd $THIS_SCRIPT_DIR

STEPSDIR=../steps

BRDIR=$THIS_SCRIPT_DIR/reports/B
TRDIR=$THIS_SCRIPT_DIR/reports/T
mkdir -p $BRDIR $TRDIR

echo "Collecting reports of T cell analysis in $TRDIR"
find $STEPSDIR/*T* -name "*.html" -exec cp '{}' $TRDIR \;
echo "Collecting reports of B cell analysis in $BRDIR"
find $STEPSDIR/*B* -name "*.html" -exec cp '{}' $BRDIR \;
