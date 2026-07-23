#! /bin/env bash
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=48G
#SBATCH --output=logs/%x.log
#SBATCH --job-name=remap_unaligned
# #SBATCH directives are parsed before the script runs, so SAMPLE can't be
# interpolated here — override the job name at submit time (see Usage below).

# Remap reads that nf-core/methylong (v1.0.0) left unaligned, using more
# permissive minimap2 settings than the pipeline itself.
#
# The pipeline's ONT_MINIMAP2_ALIGN step runs:
#   minimap2 -y -Y -x lr:hq --secondary=no -t 12 <ref> - -a
# (confirmed from work/*/.command.sh of the actual runs)
# --secondary=no discards every multi-mapping alignment outright, so any
# read whose best hit isn't clearly unique never gets a chance at a
# secondary/lower-scoring placement. This script relaxes that:
#   --secondary=yes -N 50 -p 0.5   report up to 50 secondary alignments per
#                                   read, down to 50% of the primary's score
#                                   (minimap2 defaults: report off, ratio 0.8)
#   -x map-ont                     swap the lr:hq preset (tuned for Q20+
#                                   reads) for map-ont, which tolerates the
#                                   higher error rate / shorter or chimeric
#                                   reads expected among reads that failed
#                                   the stricter preset. Edit back to lr:hq
#                                   below if you want to keep the preset and
#                                   only change the secondary-mapping behavior.
#
# Usage: sbatch remap_unaligned_reads.sh <sample_name>
#   e.g. sbatch --job-name=remap_unaligned_N07_jun remap_unaligned_reads.sh N07_jun

set -euxo pipefail

SAMPLE=${1:?Usage: $0 <sample_name> (e.g. N07_jun)}

PROJECT_DIR=/home/core/RRNS_tau_Ind-wf-basecalling
RESULTS_DIR=${PROJECT_DIR}/test_methylong-5mCG_5hmCG-traditional
UNALIGNED_DIR=${RESULTS_DIR}/unaligned_reads
REMAP_DIR=${UNALIGNED_DIR}/remap
REFERENCE=${PROJECT_DIR}/genome/GCF_002263795.3_ARS-UCD2.0_genomic.fa
SINGULARITY_OPTS="--bind ${PROJECT_DIR} --bind /home/core/"
MINIMAP2_IMAGE=${NXF_SINGULARITY_CACHEDIR}/depot.galaxyproject.org-singularity-minimap2%3A2.28--he4a0461_3.img
SAMTOOLS_IMAGE=${NXF_SINGULARITY_CACHEDIR}/samtools-1.22.1--h96c455f_0.img
CPUS=${SLURM_CPUS_PER_TASK:-12}

source_bam="${RESULTS_DIR}/ont/${SAMPLE}/alignment/${SAMPLE}.bam"
tagged_fastq="${REMAP_DIR}/${SAMPLE}.unaligned.tagged.fastq.gz"
output_bam="${REMAP_DIR}/${SAMPLE}.permissive.bam"

mkdir -p "${REMAP_DIR}"

# The existing unaligned_reads/fastq/${SAMPLE}.fastq.gz was produced with a
# plain `samtools fastq -0 out -` (no -T), so it does NOT carry the MM/ML/MN
# methylation tags. Re-extract with tags kept, so minimap2 -y can propagate
# them into the new BAM and modkit can still call methylation downstream.
if [ ! -f "${tagged_fastq}" ]; then
    singularity exec ${SINGULARITY_OPTS} "${SAMTOOLS_IMAGE}" \
        samtools view -f 4 "${source_bam}" --threads "${CPUS}" | \
        singularity exec ${SINGULARITY_OPTS} "${SAMTOOLS_IMAGE}" \
        samtools fastq -T "MM,ML,MN" -0 "${tagged_fastq}" -
    echo "Created tagged FASTQ: ${tagged_fastq}"
else
    echo "Skipped tagged FASTQ (already exists): ${tagged_fastq}"
fi

# Permissive remap: keep secondary/multi-mapping alignments and weaker hits
# instead of discarding them like the pipeline does.
if [ ! -f "${output_bam}" ]; then
    singularity exec ${SINGULARITY_OPTS} "${MINIMAP2_IMAGE}" \
        minimap2 -y -Y -x map-ont --secondary=yes -N 50 -p 0.5 \
        -t "${CPUS}" \
        "${REFERENCE}" \
        "${tagged_fastq}" \
        -a | \
    singularity exec ${SINGULARITY_OPTS} "${SAMTOOLS_IMAGE}" \
        samtools sort -@ "${CPUS}" -o "${output_bam}##idx##${output_bam}.bai" --write-index -
    echo "Created remapped BAM: ${output_bam}"
else
    echo "Skipped remap (already exists): ${output_bam}"
fi

# Quick sanity check: how many previously-unaligned reads now map, and how
# many of those are secondary/multi-mapping placements.
singularity exec ${SINGULARITY_OPTS} "${SAMTOOLS_IMAGE}" \
    samtools flagstat "${output_bam}" | tee "${REMAP_DIR}/${SAMPLE}.permissive.flagstat.txt"
