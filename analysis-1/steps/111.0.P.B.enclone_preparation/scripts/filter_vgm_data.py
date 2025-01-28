#!/usr/bin/env python3
"""
Process and filter VGM data with annotations.

Usage:
  filter_vgm_data.py [options] <vgmdata> <config>

Arguments:
  <vgmdata>       TSV file containing the VGM data (first row as header).
  <config>        TSV file containing sample <tab> JSON all_contig_annotations.json path

Options:
  --celltype=<ct>      Filter by specific celltype (string).
  --condition=<cond>   Filter by specific condition (string).
  --prjpath=<path>     Paths in config TSV are relative to this directory (default: current directory).
  -v --verbose         Enable verbose output.
  -q --quiet           Suppress all output except errors.
  -h --help            Show this help message.
  -V --version         Show version.
"""

import sys
import json
import csv
from docopt import docopt
from loguru import logger
from collections import defaultdict

def set_logger(verbose, quiet):
    logger.remove()
    if quiet:
        logger.add(lambda _: None, level="DEBUG")
    elif verbose:
        logger.add(sys.stderr, level="INFO")
    else:
        logger.add(sys.stderr, level="WARNING")


def load_vgmdata(filename, celltype=None, condition=None):
    """Load and optionally filter the VGM data."""
    logger.info(f"Load VGM data from file {filename}")
    if celltype is not None:
        logger.info(f"=> filter by celltype {celltype}")
    if condition is not None:
        logger.info(f"=> filter by condition {condition}")
    with open(filename, "r") as file:
        reader = csv.DictReader(file, delimiter="\t")
        data = [row for row in reader]
        if celltype is not None:
            data = [row for row in data if row.get("celltype") == celltype]
        if condition is not None:
            data = [row for row in data if row.get("condition") == condition]
    return data


def collect_filtered_barcodes(vgm_data):
    logger.info("Collect filtered barcodes from VGM data")
    filtered_barcodes = defaultdict(set)
    for row in vgm_data:
        #
        # Barcodes are sometimes repeated in samples, thus they are saved in a
        # dictionary with the sample as key
        #
        filtered_barcodes[row["sample"]].add(row["barcode"])
    return filtered_barcodes


def load_annotations(config_file, prjpath):
    """Load and merge annotations from JSON files."""
    logger.info(f"Load annotations listed in configuration file {config_file}")
    sample_data = {}
    with open(config_file, "r") as file:
        for line in file:
            sample, json_file = line.strip().split("\t")
            logger.info(\
                f"Load annotations for sample {sample} from file {json_file}")
            with open(prjpath+"/"+json_file, "r") as file:
                sample_data[sample] = json.load(file)
    logger.info(f"Loaded annotations for {len(sample_data)} samples")
    return sample_data


def output_filtered_annotations(filtered_barcodes, annotations):
    """Filter annotations based on barcodes and sample and output."""

    logger.info("Output filtered annotations")
    # Beginning of the JSON array:
    sys.stdout.write("[\n")
    first_elem = True

    for sample, records in annotations.items():
        for record in records:
            if record.get("barcode") in filtered_barcodes[sample]:
                if not first_elem:
                    sys.stdout.write(",\n")
                else:
                    first_elem = False
                # remove outdated info field (clonotype ids)
                record.pop('info')
                sys.stdout.write(json.dumps(record))

    # End of the JSON array:
    sys.stdout.write("\n]\n")


def main(vgmdata, config_file, celltype, condition, prjpath):
    vgm_data = load_vgmdata(vgmdata, celltype, condition)
    filtered_barcodes = collect_filtered_barcodes(vgm_data)
    annotation_data = load_annotations(config_file, prjpath)
    output_filtered_annotations(filtered_barcodes, annotation_data)


if __name__ == "__main__":
    args = docopt(__doc__, version="1.0")
    set_logger(args["--verbose"], args["--quiet"])

    main(
        vgmdata=args["<vgmdata>"],
        config_file=args["<config>"],
        celltype=args["--celltype"],
        condition=args["--condition"],
        prjpath=args["--prjpath"],
    )
