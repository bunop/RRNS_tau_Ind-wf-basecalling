
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
