#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=16:00:00
#SBATCH --mem=16G
#SBATCH --output=bedMethyl2cov.log
#SBATCH --job-name=bedMethyl2cov

# Check if conda or micromamba is available and activate the environment accordingly
if command -v conda &> /dev/null; then
    eval "$(conda shell.bash hook)"
    conda activate python3.12
elif command -v micromamba &> /dev/null; then
    eval "$(micromamba shell hook -s bash)"
    micromamba activate python3.12
else
    echo "Neither conda nor micromamba found in PATH."
    exit 1
fi

# Execute the bedMethyl2cov script with the specified input and output directories
# CpG 5mCG and 5hmCG
if [ ! -d "output_methylong-5mCG_5hmCG-cpg/cov/" ]; then
    scripts/bedMethyl2cov.py \
        -i output_methylong-5mCG_5hmCG-cpg/ont/ \
        -o output_methylong-5mCG_5hmCG-cpg/cov/ \
        &> 5mCG_5hmCG-cpg-cov.log
fi

# traditional 5mCG and 5hmCG
if [ ! -d "output_methylong-5mCG_5hmCG-traditional/cov/" ]; then
    scripts/bedMethyl2cov.py \
        -i output_methylong-5mCG_5hmCG-traditional/ont/ \
        -o output_methylong-5mCG_5hmCG-traditional/cov/ \
        &> 5mCG_5hmCG-traditional-cov.log
fi

echo "bedMethyl2cov processing completed successfully."
