#!/usr/bin/env python

"""
Parse the RepeatMasker `.tbl` summary files produced for the unaligned reads
of each sample (see scripts/unaligned_reads/*.sh) and combine them into a
single tidy CSV, one row per sample.
"""

import csv
import re
import pathlib
import logging
import argparse

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# label -> (regex, fields to extract)
# "n_elements" fields are None for labels that don't report an element count
# (bases masked / total interspersed repeats)
PATTERNS = {
    "n_sequences": (r"^sequences:\s+(\d+)", ("n_sequences",)),
    "total_length_bp": (r"^total length:\s+(\d+) bp", ("total_length_bp",)),
    "gc_pct": (r"^GC level:\s+([\d.]+) %", ("gc_pct",)),
    "bases_masked": (
        r"^bases masked:\s+(\d+) bp\s+\(\s*([\d.]+) %\)",
        ("bases_masked_bp", "bases_masked_pct"),
    ),
    "SINEs": (
        r"^SINEs:\s+(\d+)\s+(\d+) bp\s+([\d.]+) %",
        ("SINEs_n", "SINEs_bp", "SINEs_pct"),
    ),
    "LINEs": (
        r"^LINEs:\s+(\d+)\s+(\d+) bp\s+([\d.]+) %",
        ("LINEs_n", "LINEs_bp", "LINEs_pct"),
    ),
    "LTR": (
        r"^LTR elements:\s+(\d+)\s+(\d+) bp\s+([\d.]+) %",
        ("LTR_n", "LTR_bp", "LTR_pct"),
    ),
    "DNA": (
        r"^DNA elements:\s+(\d+)\s+(\d+) bp\s+([\d.]+) %",
        ("DNA_n", "DNA_bp", "DNA_pct"),
    ),
    "Unclassified": (
        r"^Unclassified:\s+(\d+)\s+(\d+) bp\s+([\d.]+) %",
        ("Unclassified_n", "Unclassified_bp", "Unclassified_pct"),
    ),
    "Total_interspersed": (
        r"^Total interspersed repeats:\s+(\d+) bp\s+([\d.]+) %",
        ("Total_interspersed_bp", "Total_interspersed_pct"),
    ),
    "Small_RNA": (
        r"^Small RNA:\s+(\d+)\s+(\d+) bp\s+([\d.]+) %",
        ("Small_RNA_n", "Small_RNA_bp", "Small_RNA_pct"),
    ),
    "Satellites": (
        r"^Satellites:\s+(\d+)\s+(\d+) bp\s+([\d.]+) %",
        ("Satellites_n", "Satellites_bp", "Satellites_pct"),
    ),
    "Simple_repeats": (
        r"^Simple repeats:\s+(\d+)\s+(\d+) bp\s+([\d.]+) %",
        ("Simple_repeats_n", "Simple_repeats_bp", "Simple_repeats_pct"),
    ),
    "Low_complexity": (
        r"^Low complexity:\s+(\d+)\s+(\d+) bp\s+([\d.]+) %",
        ("Low_complexity_n", "Low_complexity_bp", "Low_complexity_pct"),
    ),
}

# only the *_n, *_bp, *_pct columns feeding numeric conversion
NUMERIC_FIELDS = {
    field for _, fields in PATTERNS.values() for field in fields
}

FIELDNAMES = [
    "sample",
    "n_sequences", "total_length_bp", "gc_pct",
    "bases_masked_bp", "bases_masked_pct",
    "SINEs_n", "SINEs_bp", "SINEs_pct",
    "LINEs_n", "LINEs_bp", "LINEs_pct",
    "LTR_n", "LTR_bp", "LTR_pct",
    "DNA_n", "DNA_bp", "DNA_pct",
    "Unclassified_n", "Unclassified_bp", "Unclassified_pct",
    "Total_interspersed_bp", "Total_interspersed_pct",
    "Small_RNA_n", "Small_RNA_bp", "Small_RNA_pct",
    "Satellites_n", "Satellites_bp", "Satellites_pct",
    "Simple_repeats_n", "Simple_repeats_bp", "Simple_repeats_pct",
    "Low_complexity_n", "Low_complexity_bp", "Low_complexity_pct",
    # derived: low-complexity/simple-repeat share, our "complexity" proxy
    "low_complexity_share_pct",
]


def parse_tbl(tbl_file: pathlib.Path) -> dict:
    """
    Parse a single RepeatMasker `.tbl` file into a flat dict of stats.
    """

    sample = tbl_file.name.split(".fasta.tbl")[0]
    record = {"sample": sample}

    text = tbl_file.read_text()

    for _, (pattern, fields) in PATTERNS.items():
        match = re.search(pattern, text, flags=re.MULTILINE)

        if not match:
            logger.warning(f"Pattern for {fields} not found in {tbl_file}")
            continue

        for field, value in zip(fields, match.groups()):
            record[field] = float(value) if "." in value else int(value)

    record["low_complexity_share_pct"] = round(
        record.get("Simple_repeats_pct", 0.0) + record.get("Low_complexity_pct", 0.0), 4
    )

    return record


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Combine per-sample RepeatMasker .tbl summaries into one CSV."
    )
    parser.add_argument(
        "-i", "--input_folder", required=True, type=pathlib.Path,
        help="Folder containing the *.fasta.tbl RepeatMasker summary files"
    )
    parser.add_argument(
        "-o", "--output", required=True, type=pathlib.Path,
        help="Output CSV path"
    )
    args = parser.parse_args()

    if not args.input_folder.is_dir():
        logger.error(f"Input folder '{args.input_folder}' does not exist or is not a directory.")
        exit(1)

    tbl_files = sorted(args.input_folder.glob("*.fasta.tbl"))

    if not tbl_files:
        logger.error(f"No *.fasta.tbl files found in '{args.input_folder}'.")
        exit(1)

    records = []
    for tbl_file in tbl_files:
        logger.info(f"Parsing {tbl_file.name}")
        records.append(parse_tbl(tbl_file))

    args.output.parent.mkdir(parents=True, exist_ok=True)

    with open(args.output, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDNAMES)
        writer.writeheader()
        writer.writerows(records)

    logger.info(f"Wrote summary for {len(records)} samples to {args.output}")
