#!/bin/bash
if [ "$1" == "" ]; then echo "Usage $0 sfx"; exit 1; fi
for step in 10[123456789]* 110* 111*; do mv $step/results $step/results.$1; done
