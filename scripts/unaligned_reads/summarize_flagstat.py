#!/usr/bin/env python

"""
Parse the `samtools flagstat` summaries produced for the permissive remap of
each sample's unaligned reads (see scripts/unaligned_reads/remap_unaligned_reads.sh)
and combine them into a single tidy CSV, one row per sample.
"""

import csv
import re
import pathlib
import logging
import argparse

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# label -> (regex, fields to extract)
PATTERNS = {
    "total": (r"^(\d+) \+ \d+ in total", ("n_total_records",)),
    "primary": (r"^(\d+) \+ \d+ primary$", ("n_input_reads",)),
    "secondary": (r"^(\d+) \+ \d+ secondary$", ("n_secondary",)),
    "supplementary": (r"^(\d+) \+ \d+ supplementary$", ("n_supplementary",)),
    "mapped": (
        r"^(\d+) \+ \d+ mapped \(([\d.]+)%",
        ("n_mapped_records", "mapped_records_pct"),
    ),
    "primary_mapped": (
        r"^(\d+) \+ \d+ primary mapped \(([\d.]+)%",
        ("n_primary_mapped", "primary_mapped_pct"),
    ),
}

FIELDNAMES = [
    "sample",
    "n_input_reads",       # unaligned reads fed into the permissive remap
    "n_primary_mapped",    # of those, how many now have a primary alignment
    "primary_mapped_pct",  # = n_primary_mapped / n_input_reads * 100
    "n_still_unmapped",    # reads that remain unaligned even with permissive settings
    "n_secondary",         # secondary (multi-mapping) alignment records
    "n_supplementary",     # supplementary (chimeric/split) alignment records
    "secondary_per_mapped_read",  # avg. secondary alignments per recovered read
    "n_total_records",     # primary + secondary + supplementary
    "n_mapped_records",    # mapped records among n_total_records (any type)
    "mapped_records_pct",
]


def parse_flagstat(flagstat_file: pathlib.Path) -> dict:
    """
    Parse a single `samtools flagstat` text file into a flat dict of stats.
    """

    sample = flagstat_file.name.split(".permissive.flagstat.txt")[0]
    record = {"sample": sample}

    text = flagstat_file.read_text()

    for _, (pattern, fields) in PATTERNS.items():
        match = re.search(pattern, text, flags=re.MULTILINE)

        if not match:
            logger.warning(f"Pattern for {fields} not found in {flagstat_file}")
            continue

        for field, value in zip(fields, match.groups()):
            record[field] = float(value) if "." in value else int(value)

    record["n_still_unmapped"] = record.get("n_input_reads", 0) - record.get("n_primary_mapped", 0)
    record["secondary_per_mapped_read"] = round(
        record.get("n_secondary", 0) / record["n_primary_mapped"], 4
    ) if record.get("n_primary_mapped") else 0.0

    return record


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Combine per-sample samtools flagstat summaries into one CSV."
    )
    parser.add_argument(
        "-i", "--input_folder", required=True, type=pathlib.Path,
        help="Folder containing the *.permissive.flagstat.txt files"
    )
    parser.add_argument(
        "-o", "--output", required=True, type=pathlib.Path,
        help="Output CSV path"
    )
    args = parser.parse_args()

    if not args.input_folder.is_dir():
        logger.error(f"Input folder '{args.input_folder}' does not exist or is not a directory.")
        exit(1)

    flagstat_files = sorted(args.input_folder.glob("*.permissive.flagstat.txt"))

    if not flagstat_files:
        logger.error(f"No *.permissive.flagstat.txt files found in '{args.input_folder}'.")
        exit(1)

    records = []
    for flagstat_file in flagstat_files:
        logger.info(f"Parsing {flagstat_file.name}")
        records.append(parse_flagstat(flagstat_file))

    args.output.parent.mkdir(parents=True, exist_ok=True)

    with open(args.output, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDNAMES)
        writer.writeheader()
        writer.writerows(records)

    logger.info(f"Wrote summary for {len(records)} samples to {args.output}")
