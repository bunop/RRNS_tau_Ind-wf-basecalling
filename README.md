
# Adaptive sampling RRNS Tau-Ind

## Update models

To update the models, run the following command:

```bash
singularity run ${NXF_SINGULARITY_CACHEDIR}/ontresearch-dorado-shae9327ad17e023b76e4d27cf287b6b9d3a271092b.img \
    dorado download --models-directory ${HOME}/Projects/dorado_models/
```

Then ensure that the `conf/custom.config` file has the correct path to the models:

```groovy
env {
    DRD_MODELS_PATH = '${HOME}/Projects/dorado_models'
}

singularity {
    runOptions = '-B ${HOME}/Projects/dorado_models'
}
```

> NOTE: the models are downloaded to the `${HOME}/Projects/dorado_models/` directory

## Collect data

Create a directory to collect the data, then make symlinks to the pod5 files:

```bash
mkdir data
cd data
ln -s /data/minknow/RRNS_tau_Ind/no_sample_id/20250521_1532_P2S-02962-A_PBE20707_21cee8e2/pod5
```

## Calling wf-basecalling pipeline

Called the latest `wf-basecalling` pipeline (`v1.5.2`) customize to support
different models with different context (https://github.com/bunop/wf-basecalling/tree/multiple_calling)
using:

```bash
nextflow run ~/Projects/wf-basecalling/ -profile singularity -resume \
    -c conf/custom-wf-basecalling.config -params-file conf/params-wf-basecalling-5mCG_5hmCG.json
```

> NOTE: is not possible to call `--duplex=true` and `--barcode_kit=SQK-NBD114-24`
> at the same time, *demultiplexing* should be done after basecalling.

### Calling wf-basecalling on Galileo cluster

To call the `wf-basecalling` pipeline on the Galileo cluster, use the following command:

```bash
scripts/launch-wf-basecalling-5mC-5hmC.sh
```

> NOTE: is not possible to call this pipeline since the max allowed time for CPU
> job is 8 hours, and those processes require more time to complete.
> TODO: test with `cnr-ibba/nf-dorado-calling` pipeline

## Post processing

### calling PyCOQC

Calculate quality control metrics using `pycoQC`:

```bash
singularity run $NXF_SINGULARITY_CACHEDIR/pip_pycoqc_setuptools_31d5a8754dcc1b68.sif \
    pycoQC -f output_wf-basecalling-5mCG_5hmCG/SAMPLE.summary.tsv.gz -o output_wf-basecalling-5mCG_5hmCG/SAMPLE.summary.html
```

### Join passed simplex and duplex reads

Called reads are divided into two groups: `simplex` and `duplex` for both passed
and failed reads:

```text
SAMPLE.fail.duplex.cram
SAMPLE.fail.simplex.cram
SAMPLE.pass.duplex.cram
SAMPLE.pass.simplex.cram
```

Let's join the `pass` reads into a single file:

```bash
cd output_wf-basecalling-5mCG_5hmCG
singularity run $NXF_SINGULARITY_CACHEDIR/depot.galaxyproject.org-singularity-samtools-1.21--h50ea8bc_0.img \
    samtools merge -o SAMPLE.pass.all.cram SAMPLE.pass.duplex.cram SAMPLE.pass.simplex.cram
singularity run $NXF_SINGULARITY_CACHEDIR/depot.galaxyproject.org-singularity-samtools-1.21--h50ea8bc_0.img \
    samtools index SAMPLE.pass.all.cram
cd ..
```

### Demultiplexing

Do demultiplexing using `dorado`:

```bash
singularity run $NXF_SINGULARITY_CACHEDIR/ontresearch-dorado-shae9327ad17e023b76e4d27cf287b6b9d3a271092b.img \
    dorado demux --kit-name SQK-NBD114-24 --threads 2 --verbose --output-dir demux-5mCG_5hmCG \
    --sample-sheet conf/samplesheet.csv output_wf-basecalling-5mCG_5hmCG/SAMPLE.pass.all.cram
```

## Call nf-core/methylong

Data are stored in `demux-5mCG_5hmCG` folder:

```bash
sbatch scripts/launch-methylong-5mCG_5hmCG-cpg.sh
sbatch scripts/launch-methylong-5mCG_5hmCG-traditional.sh
```

## Investigate unaligned reads

`nf-core/methylong` (v1.0.0) aligns with `minimap2 -y -Y -x lr:hq --secondary=no`
(`ONT_MINIMAP2_ALIGN`, confirmed from `work/*/.command.sh` of the actual runs):
the `lr:hq` preset is tuned for Q20+ reads and `--secondary=no` discards every
multi-mapping alignment outright, so a read whose best hit isn't clearly
unique never gets a chance at a secondary/lower-scoring placement. A
non-trivial fraction of reads per sample end up unaligned as a result.

We want to demonstrate two things about those unaligned reads, for each
sample of `test_methylong-5mCG_5hmCG-traditional`:

1. **What they are** — are they enriched for repetitive elements (SINEs,
   LINEs, satellites, low-complexity/simple repeats), which would explain why
   they fail to place uniquely under strict settings?
2. **Whether they're actually mappable** — if we relax the aligner (allow
   secondary/multi-mapping hits, lower the score threshold, use a
   higher-error-tolerant preset), how many of them do map, and how many
   remain genuinely unmappable? This tells us whether the pipeline's
   unaligned rate mostly reflects the strictness of its settings, or reads
   that don't belong to the reference at all.

The analysis produces two CSV reports in
`test_methylong-5mCG_5hmCG-traditional/unaligned_reads/`.

### 1. RepeatMasker profile of unaligned reads

For each sample, `scripts/unaligned_reads/unaligned_reads_<SAMPLE>.sh`
extracts the unmapped reads (`samtools view -f 4`) from the pipeline's BAM,
converts them to FASTA, and runs `RepeatMasker -species cow -xsmall` on them.
These are independent per-sample SLURM jobs, launched with:

```bash
sbatch scripts/unaligned_reads/unaligned_reads_A19_jun.sh
sbatch scripts/unaligned_reads/unaligned_reads_A21_jun.sh
sbatch scripts/unaligned_reads/unaligned_reads_A25_jun.sh
sbatch scripts/unaligned_reads/unaligned_reads_N03_jun.sh
sbatch scripts/unaligned_reads/unaligned_reads_N07_jun.sh
sbatch scripts/unaligned_reads/unaligned_reads_N13_jun.sh
```

Each job writes fastq/fasta under `unaligned_reads/{fastq,fasta}/` and a
`<SAMPLE>.fasta.tbl` RepeatMasker summary. Once all jobs have completed,
combine the per-sample `.tbl` files into a single tidy CSV:

```bash
python scripts/unaligned_reads/summarize_repeatmasker.py \
    -i test_methylong-5mCG_5hmCG-traditional/unaligned_reads/fasta \
    -o test_methylong-5mCG_5hmCG-traditional/unaligned_reads/repeatmasker_summary.csv
```

### 2. Recovery rate under a permissive remap

`scripts/unaligned_reads/remap_unaligned_reads.sh <SAMPLE>` re-extracts the
same unmapped reads (this time keeping the `MM`/`ML`/`MN` methylation tags,
needed so `minimap2 -y` can propagate them and `modkit` can still call
methylation downstream) and remaps them with more permissive settings than
the pipeline: `-x map-ont --secondary=yes -N 50 -p 0.5` instead of
`-x lr:hq --secondary=no`. It then runs `samtools flagstat` on the result.
Since the `#SBATCH --job-name` can't reference the sample argument, the job
name is overridden at submit time, one job per sample:

```bash
sbatch --job-name=remap_unaligned_A19_jun scripts/unaligned_reads/remap_unaligned_reads.sh A19_jun
sbatch --job-name=remap_unaligned_A21_jun scripts/unaligned_reads/remap_unaligned_reads.sh A21_jun
sbatch --job-name=remap_unaligned_A25_jun scripts/unaligned_reads/remap_unaligned_reads.sh A25_jun
sbatch --job-name=remap_unaligned_N03_jun scripts/unaligned_reads/remap_unaligned_reads.sh N03_jun
sbatch --job-name=remap_unaligned_N07_jun scripts/unaligned_reads/remap_unaligned_reads.sh N07_jun
sbatch --job-name=remap_unaligned_N13_jun scripts/unaligned_reads/remap_unaligned_reads.sh N13_jun
```

Each job writes to `unaligned_reads/remap/` and leaves a
`<SAMPLE>.permissive.flagstat.txt`. Once all jobs have completed, combine
them into a single tidy CSV:

```bash
python scripts/unaligned_reads/summarize_flagstat.py \
    -i test_methylong-5mCG_5hmCG-traditional/unaligned_reads/remap \
    -o test_methylong-5mCG_5hmCG-traditional/unaligned_reads/flagstat_summary.csv
```

`flagstat_summary.csv` reports, per sample, how many of the originally
unaligned reads got a primary alignment back (`n_primary_mapped` /
`primary_mapped_pct`) versus how many are still unmapped
(`n_still_unmapped`) — e.g. in the current run, 68-80% of "unaligned" reads
across samples do map once secondary alignments and a lower score threshold
are allowed, indicating the pipeline's strict settings are discarding a
substantial share of reads that land in multi-mapping/repetitive regions
rather than reads that are simply not present in the reference.

### 3. Read-length profile of the three read groups

RRNS reads are `MspI`-digested to 150-450 bp, so shorter reads have less
unique flanking context and may be more prone to ambiguous (secondary) or
failed alignment. To test that, `scripts/unaligned_reads/length_analysis.sh
<SAMPLE>` extracts, per sample, the per-read query length (and secondary
count where applicable) of three groups:

- `originally_mapped` — primary alignments from the pipeline's own BAM
  (`ont/<SAMPLE>/alignment/<SAMPLE>.bam`), used as-is, no remapping.
- `recovered` — primary alignments from the permissive remap BAM
  (`unaligned_reads/remap/<SAMPLE>.permissive.bam`, see above), i.e. reads
  unmapped in the original BAM that got a hit under the relaxed settings.
- `still_unmapped` — reads unmapped in both BAMs, with length read back from
  the original unaligned FASTQ.

These are per-sample SLURM jobs (same reasoning as above: the permissive
BAMs are 40-50 GB each, dominated by secondary alignments, so running all
six serially would be far slower than one job per sample):

```bash
sbatch --job-name=length_analysis_A19_jun scripts/unaligned_reads/length_analysis.sh A19_jun
sbatch --job-name=length_analysis_A21_jun scripts/unaligned_reads/length_analysis.sh A21_jun
sbatch --job-name=length_analysis_A25_jun scripts/unaligned_reads/length_analysis.sh A25_jun
sbatch --job-name=length_analysis_N03_jun scripts/unaligned_reads/length_analysis.sh N03_jun
sbatch --job-name=length_analysis_N07_jun scripts/unaligned_reads/length_analysis.sh N07_jun
sbatch --job-name=length_analysis_N13_jun scripts/unaligned_reads/length_analysis.sh N13_jun
```

Each job writes per-sample/per-group CSVs, gzipped directly, under
`stats-unaligned/length_analysis/per_sample/` (`<SAMPLE>.<group>.csv.gz`).
Once all jobs have completed, merge them into a single per-read CSV and
print the raw counts per sample/group:

```bash
python scripts/unaligned_reads/merge_length_stats.py \
    -i stats-unaligned/length_analysis/per_sample \
    -o stats-unaligned/length_analysis/read_length_stats.csv.gz
```

`read_length_stats.csv.gz` has one row per read (`sample`, `group`,
`read_id`, `read_length`, `n_secondary`), with no subsampling or statistical
filtering applied — that's left to the downstream Quarto report (both
pandas' `read_csv` and R's `readr` read `.gz` CSVs transparently). Pass a
plain `.csv` path to `-o` instead if an uncompressed output is preferred.
