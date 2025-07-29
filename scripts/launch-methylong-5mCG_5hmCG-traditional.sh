#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=4-00:00:00
#SBATCH --mem=16G
#SBATCH --error=methylong-5mCG_5hmCG-traditional.err
#SBATCH --output=methylong-5mCG_5hmCG-traditional.out
#SBATCH --job-name=methylong-5mCG_5hmCG-traditional
#SBATCH --account=IscrC_BioGPUPX    # account name
#SBATCH --partition=g100_usr_prod   # partition name (see https://docs.hpc.cineca.it/hpc/galileo.html#job-managing-and-slurm-partitions)
#SBATCH --qos=g100_qos_lprod        # quality of service (see https://docs.hpc.cineca.it/hpc/galileo.html#job-managing-and-slurm-partitions)

# set the path of institution-specific configuration files
# (required since we're working offline)
export CUSTOM_CONFIG_BASE=${WORK}/nf-configs

# mind to the pipeline version (required)
nextflow run nf-core/methylong -r 1.0.0 \
    --custom_config_base ${CUSTOM_CONFIG_BASE} \
    -config ${CUSTOM_CONFIG_BASE}/nfcore_custom.config \
    -config conf/custom-methylong-traditional.config \
    -profile ibba,galileo -resume -params-file conf/params-methylong-5mCG_5hmCG-traditional.json
