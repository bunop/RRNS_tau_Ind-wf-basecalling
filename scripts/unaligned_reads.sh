#! /bin/env bash
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --output=%x.log
#SBATCH --job-name=unaligned_reads

set -euxo pipefail

PROJECT_DIR=/home/core/RRNS_tau_Ind-wf-basecalling
RESULTS_DIR=${PROJECT_DIR}/test_methylong-5mCG_5hmCG-traditional
UNALIGNED_DIR=${RESULTS_DIR}/unaligned_reads
SINGULARITY_OPTS="--bind ${PROJECT_DIR}:${PROJECT_DIR}"
SAMTOOLS_IMAGE=${NXF_SINGULARITY_CACHEDIR}/samtools-1.22.1--h96c455f_0.img
SEQKIT_IMAGE=${NXF_SINGULARITY_CACHEDIR}/community.wave.seqera.io-library-seqkit-2.13.0--05c0a96bf9fb2751.img
REPEATMASKER_IMAGE=${NXF_SINGULARITY_CACHEDIR}/depot.galaxyproject.org-singularity-repeatmasker-4.1.2.p1--pl5321hdfd78af_1.img
CPUS=${SLURM_CPUS_PER_TASK:-4}

# search for aligned files
ALIGNED_FILES=$(find ${RESULTS_DIR} -type f -name "*.bam" -not -path "*/repair/*")

# create fastq output directory
mkdir -p "${UNALIGNED_DIR}/fastq"
mkdir -p "${UNALIGNED_DIR}/fasta"

# process each BAM file
for bam_file in ${ALIGNED_FILES}; do
    base_name=$(basename "${bam_file}" .bam)
    output_fastq="${UNALIGNED_DIR}/fastq/${base_name}.fastq.gz"
    output_fasta="${UNALIGNED_DIR}/fasta/${base_name}.fasta.gz"
    output_masked="${UNALIGNED_DIR}/fasta/${base_name}.fasta.masked.gz"

    # transform BAM to FASTQ for unaligned reads
    if [ ! -f "${output_fastq}" ]; then
        singularity exec ${SINGULARITY_OPTS} "${SAMTOOLS_IMAGE}" \
            samtools view -f 4 "${bam_file}" --threads $CPUS | \
            singularity exec ${SINGULARITY_OPTS} "${SAMTOOLS_IMAGE}" \
            samtools fastq -0 "${output_fastq}" -
        echo "Created FASTQ: ${output_fastq}"
    else
        echo "Skipped FASTQ (already exists): ${output_fastq}"
    fi

    # transform BAM to FASTA for unaligned reads
    if [ ! -f "${output_fasta}" ]; then
        singularity exec ${SINGULARITY_OPTS} "${SEQKIT_IMAGE}" \
            seqkit fq2fa --threads $CPUS "${output_fastq}" -o "${output_fasta}" --line-width 60
        echo "Created FASTA: ${output_fasta}"
    else
        echo "Skipped FASTA (already exists): ${output_fasta}"
    fi

    # run RepeatMasker on FASTA sequences
    if [ ! -f "${output_masked}" ]; then
        rm_dir="${UNALIGNED_DIR}/fasta"
        rm_basename="${base_name}.fasta"
        id_map="${rm_dir}/${base_name}.id_map.tsv"

        # RepeatMasker rejects sequence ids longer than 50 characters (e.g. Dorado
        # split-read ids "<uuid>;<parent_uuid>"), so replace headers with short
        # sequential ids and keep the original id in a mapping table.
        zcat "${output_fasta}" | awk -v mapfile="${id_map}" '
            /^>/ {
                n++
                newid = "seq_" n
                print newid "\t" substr($0, 2) >> mapfile
                print ">" newid
                next
            }
            { print }
        ' > "${rm_dir}/${rm_basename}"

        singularity exec ${SINGULARITY_OPTS} "${REPEATMASKER_IMAGE}" \
            bash -c "cd ${rm_dir} && RepeatMasker -species cow -pa $CPUS -xsmall -gff ${rm_basename}"

        # need to pack the RM fasta file
        gzip "${rm_dir}/${rm_basename}.masked"
        gzip "${id_map}"

        # remove the intermediate renamed-id FASTA, output files are already in place
        rm -f "${rm_dir}/${rm_basename}"
        echo "Created RepeatMasker output: ${output_masked}"
    else
        echo "Skipped RepeatMasker (already exists): ${output_masked}"
    fi
done
