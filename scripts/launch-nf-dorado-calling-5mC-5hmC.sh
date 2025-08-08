#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=4-00:00:00
#SBATCH --mem=16G
#SBATCH --output=dorado-calling-5mC-5hmC.log
#SBATCH --job-name=dorado-calling-5mC-5hmC
#SBATCH --account=IscrC_BioGPUPX    # account name
#SBATCH --partition=g100_usr_prod   # partition name (see https://docs.hpc.cineca.it/hpc/galileo.html#job-managing-and-slurm-partitions)
#SBATCH --qos=g100_qos_lprod        # quality of service (see https://docs.hpc.cineca.it/hpc/galileo.html#job-managing-and-slurm-partitions)

# set the path of institution-specific configuration files
# (required since we're working offline)
export CUSTOM_CONFIG_BASE=${WORK}/nf-configs

# mind to the pipeline version (required)
nextflow run cnr-ibba/nf-dorado-calling -r dev \
    --custom_config_base ${CUSTOM_CONFIG_BASE} \
    -config ${CUSTOM_CONFIG_BASE}/nfcore_custom.config \
    -config conf/custom-nf-dorado-calling.config \
    -profile ibba,galileo -resume -params-file conf/params-nf-dorado-calling.json
