#!/usr/bin/env python

import csv
import gzip
import pathlib
import logging
import argparse

from dataclasses import dataclass, field

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

header = [
    "chrom",
    "chromStart",
    "chromEnd",
    "name",
    "score",
    "strand",
    "thickStart",
    "thickEnd",
    "color",
    "valid_coverage",
    "percent_modified",
    "count_modified",
    "count_canonical",
    "count_other_mod",
    "count_delete",
    "count_fail",
    "count_diff",
    "count_nocall"
]

@dataclass
class bedMethylRecord:
    chrom: str
    chromStart: int
    chromEnd: int
    name: str
    score: int
    strand: str
    thickStart: int
    thickEnd: int
    color: str
    valid_coverage: int
    percent_modified: float
    count_modified: int
    count_canonical: int
    count_other_mod: int
    count_delete: int
    count_fail: int
    count_diff: int
    count_nocall: int
    custom_score: int = field(init=False, default=0)

    def __post_init__(self):
        self.chromStart = int(self.chromStart)
        self.chromEnd = int(self.chromEnd)
        self.score = int(self.score)
        self.thickStart = int(self.thickStart)
        self.thickEnd = int(self.thickEnd)
        self.valid_coverage = int(self.valid_coverage)
        self.percent_modified = float(self.percent_modified)
        self.count_modified = int(self.count_modified)
        self.count_canonical = int(self.count_canonical)
        self.count_other_mod = int(self.count_other_mod)
        self.count_delete = int(self.count_delete)
        self.count_fail = int(self.count_fail)
        self.count_diff = int(self.count_diff)
        self.count_nocall = int(self.count_nocall)

        # custom fields
        self.custom_score = self.count_modified + self.count_canonical

    def validate(self) -> list[str]:
        errors = []

        if self.strand not in ["+", "-", '.']:
            errors.append(f"Invalid strand '{self.strand}'")

        return errors

    def to_cov(self) -> list:
        """
        Convert the bedMethyl record to a cov record.
        """

        # consider this record: bedMethylRecord(chrom='NC_037328.1', chromStart=35464,
        # chromEnd=35465, name='h', score=9, strand='-', thickStart=35464, thickEnd=35465,
        # color='255,0,0', valid_coverage=9, percent_modified=11.11, count_modified=1,
        # count_canonical=0, count_other_mod=8, count_delete=0, count_fail=0, count_diff=2,
        # count_nocall=0): modkit expect to have 8 different modified call from name,
        # and this percent modified is calculated as count_modified / valid_coverage * 100.
        # I cannot handle count_other_mod, count_delete, count_fail, count_diff, count_nocall
        # so I will recalculate percent_modified as
        # count_modified / (count_canonical + count_modified) * 100
        denominator = self.count_canonical + self.count_modified
        if denominator == 0:
            percent_modified = 0.0
        else:
            percent_modified = round((self.count_modified / denominator) * 100, 2)

        if percent_modified != self.percent_modified:
            logger.debug(
                f"Percent modified mismatch: modkit:{self.percent_modified} vs self:"
                f"{percent_modified} for {self.name} at {self.chrom}:{self.chromStart}-"
                f"{self.chromEnd}"
            )

        return [
            self.chrom,
            self.chromStart + 1,  # Convert to 1-based index
            self.chromEnd,
            percent_modified,
            self.count_modified,
            self.count_canonical
        ]


def open_bedMethyl_file(filename: pathlib.PosixPath, mode='rt'):
    """
    Open a BED file (possibly gzipped) for reading.
    """

    if filename.suffix == '.gz':
        handle = gzip.open(filename, mode)

    else:
        handle = open(filename, mode)

    reader = csv.DictReader(handle, fieldnames=header, delimiter='\t')

    for i, row in enumerate(reader):
        if (i+1) % 100000 == 0:
            logger.info(f"Processing record {i+1} in {filename}")

        record = bedMethylRecord(**row)

        record_errors = record.validate()
        if record_errors:
            raise ValueError(f"Record {i} in {filename} is invalid: {record_errors}")

        yield record

    logger.info(f"Processed {i+1} records in {filename}")

    handle.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Convert a bedMethyl file to a coverage file."
    )
    parser.add_argument(
        "-i", "--input_folder", required=True,
        help="Input folder containing bedMethyl files",
        type=pathlib.Path
    )
    parser.add_argument(
        "-o", "--output", required=True,
        help="Output folder for the coverage data",
        type=pathlib.Path
    )
    parser.add_argument(
        "-s", "--score", type=int,
        help="Minimum score for filtering records (default: %(default)s)",
    )
    parser.add_argument(
        "--custom_score", type=int, default=1,
        help=(
            "Custom score for filtering records (intended as count_modified + "
            "count_canonical - default: %(default)s)"
        )
    )
    parser.add_argument(
        "--force", action='store_true', default=False,
        help="Force overwrite of output files if they already exist (default: %(default)s)"
    )
    args = parser.parse_args()

    if not args.input_folder.is_dir():
        logger.error(
            f"Input folder '{args.input_folder}' does not exist or is not a directory.")
        exit(1)

    args.output.mkdir(parents=True, exist_ok=args.force)

    logger.info(f"Processing input folder: {args.input_folder}")

    for bed_file in args.input_folder.rglob("*.bed.gz"):
        logger.debug(f"Processing file: {bed_file.relative_to(args.input_folder)}")

        output_prefix = args.output / (bed_file.name.split('.bed')[0])

        handles = {}
        writers = {}

        for record in open_bedMethyl_file(bed_file):
            if args.score and record.score < args.score:
                continue

            if args.custom_score and record.custom_score < args.custom_score:
                continue

            if record.name not in handles:
                output_file = output_prefix.with_suffix(f".{record.name}.cov.gz")
                logger.info(f"Creating output file: {output_file}")
                handles[record.name] = gzip.open(output_file, mode='wt')
                writers[record.name] = csv.writer(handles[record.name], delimiter='\t')

            writer = writers[record.name]
            writer.writerow(record.to_cov())

        for handle in handles.values():
            handle.close()

    logger.info(f"Output written to: {args.output}")
    logger.info("Processing complete.")
