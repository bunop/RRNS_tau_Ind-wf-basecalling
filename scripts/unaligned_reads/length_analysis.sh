#! /bin/env bash
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=logs/%x.log
#SBATCH --job-name=length_analysis
# #SBATCH directives are parsed before the script runs, so SAMPLE can't be
# interpolated here — override the job name at submit time (see Usage below).

# Per-read query length (and secondary-alignment count) extraction for the
# three RRNS read groups used to test the "short MspI fragments map worse"
# hypothesis (low mapping efficiency correlates with the 150-450 bp RRNS
# fragment size, since shorter reads have less unique flanking context):
#
#   A originally_mapped - primary alignments from the pipeline's own BAM
#     (ont/<sample>/alignment/<sample>.bam, minimap2 -x lr:hq --secondary=no).
#     No secondary alignments exist here by construction.
#   B recovered         - primary alignments from the permissive remap BAM
#     (unaligned_reads/remap/<sample>.permissive.bam, minimap2 -x map-ont
#     --secondary=yes -N 50 -p 0.5; see remap_unaligned_reads.sh), i.e. reads
#     that were unmapped in A but got a primary hit here. n_secondary counts
#     the secondary (flag 0x100) records sharing that read's QNAME.
#   C still_unmapped     - reads unmapped in both BAMs; length is read back
#     from the original unaligned FASTQ (unaligned_reads/fastq/<sample>.fastq.gz)
#     rather than from the permissive BAM's unmapped records.
#
# No subsampling, no length/quality filtering: every read in each group is
# kept as-is. Downstream statistics/subsampling belong in the Quarto report,
# not here.
#
# Usage: sbatch length_analysis.sh <sample_name>
#   e.g. sbatch --job-name=length_analysis_N07_jun length_analysis.sh N07_jun
#
# Writes one CSV per group per sample (no header) under
# stats-unaligned/length_analysis/per_sample/; merge_length_stats.py then
# combines all of them into the final single CSV. Splitting like this lets
# the six samples run as separate SLURM jobs in parallel instead of one
# serial script churning through ~50 GB permissive BAMs x6 back to back.

set -euxo pipefail

SAMPLE=${1:?Usage: $0 <sample_name> (e.g. N07_jun)}

PROJECT_DIR=/home/core/RRNS_tau_Ind-wf-basecalling
RESULTS_DIR=${PROJECT_DIR}/test_methylong-5mCG_5hmCG-traditional
UNALIGNED_DIR=${RESULTS_DIR}/unaligned_reads
REMAP_DIR=${UNALIGNED_DIR}/remap
OUT_DIR=${PROJECT_DIR}/stats-unaligned/length_analysis/per_sample
SINGULARITY_OPTS="--bind ${PROJECT_DIR} --bind /home/core/"
SAMTOOLS_IMAGE=${NXF_SINGULARITY_CACHEDIR}/samtools-1.22.1--h96c455f_0.img
CPUS=${SLURM_CPUS_PER_TASK:-8}

original_bam="${RESULTS_DIR}/ont/${SAMPLE}/alignment/${SAMPLE}.bam"
permissive_bam="${REMAP_DIR}/${SAMPLE}.permissive.bam"
unaligned_fastq="${UNALIGNED_DIR}/fastq/${SAMPLE}.fastq.gz"

mkdir -p "${OUT_DIR}"

group_a_csv="${OUT_DIR}/${SAMPLE}.originally_mapped.csv.gz"
group_b_csv="${OUT_DIR}/${SAMPLE}.recovered.csv.gz"
group_c_csv="${OUT_DIR}/${SAMPLE}.still_unmapped.csv.gz"
unmapped_ids="${OUT_DIR}/${SAMPLE}.still_unmapped_ids.txt"

# Outputs are gzipped directly (stats-unaligned/ is kept compressed), so
# clean up both the compressed outputs and the plain-text working file
# before regenerating them.
rm -f "${group_a_csv}" "${group_b_csv}" "${group_c_csv}" "${unmapped_ids}"

# --- Group A: originally_mapped ---------------------------------------------
# -F 0x904 drops unmapped (0x4), secondary (0x100) and supplementary (0x800)
# records in a single native samtools filter, leaving one primary row per
# mapped read. Read straight off the pipeline's own BAM, untouched.
singularity exec ${SINGULARITY_OPTS} "${SAMTOOLS_IMAGE}" \
    samtools view -@ "${CPUS}" -F 0x904 "${original_bam}" | \
    awk -F'\t' -v OFS=',' -v sample="${SAMPLE}" \
        '{ print sample, "originally_mapped", $1, length($10), "NA" }' | \
    gzip -c > "${group_a_csv}"

# --- Group B (recovered) + still-unmapped QNAMEs ----------------------------
# One pass over the whole permissive BAM classifies every record by FLAG
# (unmapped / secondary / primary mapped). The BAM is coordinate-sorted, not
# queryname-sorted, so a read's secondary alignments aren't adjacent to its
# primary record; primary lengths and secondary counts are therefore
# accumulated in memory and only joined together at end-of-file. Doing this
# in one pass (rather than one samtools view per flag combination) matters
# here because these BAMs are tens of GB, dominated by secondary alignments.
singularity exec ${SINGULARITY_OPTS} "${SAMTOOLS_IMAGE}" \
    samtools view -@ "${CPUS}" "${permissive_bam}" | \
    awk -F'\t' -v OFS=',' -v sample="${SAMPLE}" -v unmapped_out="${unmapped_ids}" '
        {
            flag = $2
            is_unmapped      = int(flag / 4)    % 2
            is_secondary     = int(flag / 256)  % 2
            is_supplementary = int(flag / 2048) % 2
            if (is_unmapped) {
                print $1 >> unmapped_out
            } else if (is_secondary) {
                sec_count[$1]++
            } else if (!is_supplementary) {
                primary_len[$1] = length($10)
            }
        }
        END {
            for (qname in primary_len) {
                sc = (qname in sec_count) ? sec_count[qname] : 0
                print sample, "recovered", qname, primary_len[qname], sc
            }
        }
    ' | gzip -c > "${group_b_csv}"

# --- Group C: still_unmapped, length recovered from the original FASTQ -----
zcat "${unaligned_fastq}" | awk -v sample="${SAMPLE}" -v ids_file="${unmapped_ids}" -v OFS=',' '
    BEGIN {
        while ((getline line < ids_file) > 0) { want[line] = 1 }
        close(ids_file)
    }
    NR % 4 == 1 { id = substr($1, 2); keep = (id in want); next }
    NR % 4 == 2 { if (keep) print sample, "still_unmapped", id, length($0), "NA"; next }
' | gzip -c > "${group_c_csv}"

# unmapped_ids was only a working file to bridge the BAM pass and the FASTQ
# pass above; drop it so it doesn't linger uncompressed in stats-unaligned/.
rm -f "${unmapped_ids}"

n_a=$(zcat "${group_a_csv}" | wc -l)
n_b=$(zcat "${group_b_csv}" | wc -l)
n_c=$(zcat "${group_c_csv}" | wc -l)

echo "${SAMPLE} originally_mapped: ${n_a}"
echo "${SAMPLE} recovered:         ${n_b}"
echo "${SAMPLE} still_unmapped:    ${n_c}"
