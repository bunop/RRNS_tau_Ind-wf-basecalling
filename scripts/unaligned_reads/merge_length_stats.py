#!/usr/bin/env python

"""
Merge the per-sample, per-group read-length CSVs produced by
scripts/unaligned_reads/length_analysis.sh into a single per-read CSV, and
print the raw read count per sample/group so the real proportions between
originally_mapped / recovered / still_unmapped are visible before any
subsampling is decided in the downstream Quarto report.

No statistics, filtering or subsampling happens here: this only concatenates
already-extracted rows and counts them.
"""

import csv
import gzip
import logging
import pathlib
import argparse
from collections import OrderedDict

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

SAMPLES = ["A19_jun", "A21_jun", "A25_jun", "N03_jun", "N07_jun", "N13_jun"]
GROUPS = ["originally_mapped", "recovered", "still_unmapped"]
FIELDNAMES = ["sample", "group", "read_id", "read_length", "n_secondary"]


def open_text(path: pathlib.Path, mode: str):
    """Open plain or gzip-compressed text transparently, based on the '.gz' suffix."""
    if path.suffix == ".gz":
        return gzip.open(path, mode + "t")
    return open(path, mode, newline="" if "w" in mode else None)


def main():
    parser = argparse.ArgumentParser(
        description="Merge per-sample/per-group read-length CSVs into one CSV and print group counts."
    )
    parser.add_argument(
        "-i", "--input_folder", required=True, type=pathlib.Path,
        help="Folder containing <sample>.<group>.csv.gz files (length_analysis.sh output)",
    )
    parser.add_argument(
        "-o", "--output", required=True, type=pathlib.Path,
        help="Output merged CSV path (write a '.gz' path to compress it directly)",
    )
    args = parser.parse_args()

    if not args.input_folder.is_dir():
        logger.error(f"Input folder '{args.input_folder}' does not exist or is not a directory.")
        raise SystemExit(1)

    counts = OrderedDict((sample, OrderedDict((g, 0) for g in GROUPS)) for sample in SAMPLES)

    args.output.parent.mkdir(parents=True, exist_ok=True)

    with open_text(args.output, "w") as out_handle:
        writer = csv.writer(out_handle)
        writer.writerow(FIELDNAMES)

        for sample in SAMPLES:
            for group in GROUPS:
                part = args.input_folder / f"{sample}.{group}.csv.gz"
                if not part.exists():
                    logger.warning(f"Missing {part}, skipping")
                    continue
                n = 0
                with open_text(part, "r") as in_handle:
                    for line in in_handle:
                        if not line.endswith("\n"):
                            line += "\n"
                        out_handle.write(line)
                        n += 1
                counts[sample][group] = n

    logger.info(f"Wrote merged CSV to {args.output}")

    header = f"{'sample':<10} {'originally_mapped':>18} {'recovered':>10} {'still_unmapped':>15} {'total':>10}"
    print(header)
    print("-" * len(header))
    for sample in SAMPLES:
        row = counts[sample]
        total = sum(row.values())
        print(f"{sample:<10} {row['originally_mapped']:>18} {row['recovered']:>10} {row['still_unmapped']:>15} {total:>10}")


if __name__ == "__main__":
    main()
