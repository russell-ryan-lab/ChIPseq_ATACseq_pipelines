# Table of Contents

* [About](#about)
* [Before Starting](#before-starting)
  * [Managing Software with Conda](#managing-software-with-conda)
  * [Installing the Pipelines](#installing-the-pipelines)
  * [Genome Reference Requirements](#genome-reference-requirements)
  * [Lab-Specific Customizations](#lab-specific-customizations)
* [Quick-Start Example](#quick-start-example-Processing-raw-reads-from-an-ATAC-seq-experiment)
* [Additional Information](#additional-information)
  * [Fastq Inputs](#fastq-inputs)
  * [Examples](#examples)



## About

This repository contains snakefiles, configuration files, and a config generator script for processing ATACseq and ChIPseq experiments.

These snakefiles are designed so that input fastqs and their associated wildcards are explicitly defined in the config, rather than implicitly determined by filename. This provides the flexibility allowing input files to be named in any fashion, and stored in any location. The easiest way to process a new experiment after receiving the data from the UMich sequencing core is to use the config generator script. The config generator script can be modified as needed to facilitate loading of differently named/stored inputs, or configs could even be manually edited.

### Before Starting

#### Managing Software with Conda

The pipelines manage software dependencies using [conda](https://docs.conda.io/en/latest/miniconda.html). To install miniconda to your home directory on Great Lakes, follow the instructions below; they are based on the [miniconda installation instructions](https://docs.conda.io/en/latest/miniconda.html).

    # Download for 64-bit linux to your home directory:
    cd ~
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

    # Run the installer (choosing default options is fine)
    bash Miniconda3-latest-Linux-x86_64.sh

After installation, create a conda environment containing snakemake and pandas by running:

    conda create -n snakemake -c conda-forge -c bioconda snakemake=5.10.0 pandas=0.24.2

To use snakemake, activate the environment by running:

    conda activate snakemake

To ensure that snakemake is installed in the environment, run:

    snakemake --version

The result should be `5.10.0`, the latest snakemake release (as of 2/20/2020).

Finally, deactivate the environment with:

    conda deactivate

#### Installing the Pipelines

    #Navigate to where you want this pipeline to be located
    cd /nfs/turbo/path-rjhryan-turbo/lab-members/Travis/

    git clone https://github.com/russell-ryan-lab/ChIPseq_ATACseq_pipelines

#### Genome Reference Requirements

Prior to using the pipeline, certain genome reference files are required. In particular, BWA indices, chromosome sizes, and blacklist regions are required for both the ATAC-seq and ChIP-seq pipelines. The ATAC-seq pipeline also requires a TSS file. Each is described in more detail below.

##### BWA Indices

Most iGenomes references ([link](https://support.illumina.com/sequencing/sequencing_software/igenome.html)) come with pre-built BWA indices. For example, in the hg19 download from iGenomes, the BWA index is located in `Homo_sapiens/UCSC/hg19/Sequence/BWAIndex`.

To build a BWA index for a new genome or genome build see the [bwa manual](http://bio-bwa.sourceforge.net/bwa.shtml). Briefly,

    bwa index ref.fa

##### Chromosome Sizes

Chromosome size files should be tab-delimited text with two columns, the first being the chromosome name (e.g. chr1) and the second being the chromosome length in bp (e.g. 249250621). These files can be downloaded using the `fetchChromSizes` tool ([Linux download](https://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/fetchChromSizes)). Alternatively, iGenomes references often contain this file in their folder hierarchy. The chromosome naming convention (i.e. chr1 or 1) must match those of the reference genome used.

An example:

    chr1	249250621
    chr2	243199373
    chr3	198022430
    chr4	191154276
    chr5	180915260
    ...

##### Blacklist Regions

Blacklist regions should be in the form of a [BED](https://genome.ucsc.edu/FAQ/FAQformat.html#format1) file. Version 2 blacklist regions can be downloaded from the Boyle Lab ([link](https://github.com/Boyle-Lab/Blacklist/tree/master/lists)). Version 1 blacklist regions can be downloaded from the Kundaje Lab ([link](https://sites.google.com/site/anshulkundaje/projects/blacklists)).

The blacklist filtering step requires that the BED file used contains mutually disjoint regions. To ensure this is the case, one can use [`bedtools merge`](https://bedtools.readthedocs.io/en/latest/content/tools/merge.html) to create a final BED file for use in the pipeline.

An example:

    chr1	564449	570371
    chr1	724136	727043
    chr1	825006	825115
    chr1	2583334	2634374
    chr1	4363064	4363242
    chr1	5725866	5736651
    chr1	16839923	16841396
    chr1	38077347	38077423
    chr1	91852785	91853147
    ...

##### TSS File

TSS BED files are used to calculate TSS enrichment scores in the `ataqv` quality asseessment. The `ataqv` respository contains ready-made TSS files for hg19, mm9, and rn5 ([link](https://github.com/ParkerLab/ataqv/tree/master/data/tss)), but also includes the code used to generate them as an example for other references.

An example:

    chr1	11874	11874
    chr1	69091	69091
    chr1	321084	321084
    chr1	321146	321146
    chr1	322037	322037
    ...

#### Lab-specific customizations

After installation, certain configuration files require editing for the pipeline to function. These include:

- `config/cluster_config.yaml`: The `account` attribute should be changed to the appropriate Great Lakes account.
- `config/ATAC_general.yaml`: The paths to bwa indices, chromosome sizes, TSSs, and blacklists should be changed.
- `config/ChIP_histone_general.yaml`: The paths to bwa indices, chromosome sizes, and blacklists should be changed.
- `config/ChIP_TF_general.yaml`: The paths to bwa indices, chromosome sizes, and blacklists should be changed.
- `config/ChIP_TF_se_general.yaml`: The paths to bwa indices, chromosome sizes, and blacklists should be changed.

### Quick-start example: Processing raw reads from an ATAC-seq experiment

First, create the directory where you'd like to save the results of the pipeline:

    mkdir -p /path/to/results

Next, create `tmp/` and `logs/` directories for the pipeline to use:

    cd /path/to/results

    mkdir -p logs
    mkdir -p tmp

Next, generate a pipeline configuration file using the `scripts/config_creator.py` script. The following should be specified:

1. General configuration details in the form of a partial configuration file (for example, `config/ATAC_general.yaml`).
2. A comma-separated file with one line per sample giving the sample ID, the human-readable sample name, the genome for alignment, and the basepath to directory containing the fastqs for the sample (for example, `data/atac_test_data/atac_test_samplesheet.csv`, and below).
3. A path indicating where files generated by the pipeline should go.
4. A temporary directory to use during certain

    |lib|sample|genome|basepath|
    |---|------|------|--------|
    |1|s1_control|mm10|data/atac_test_data/Sample_1|
    |2|s1_treatment|mm10|data/atac_test_data/Sample_2|
    |3|s2_treatment|mm10|data/atac_test_data/Sample_3|
    |4|s3_treatment|mm10|data/atac_test_data/Sample_4|

Along with this quick example, there are several more configuration files which are included in the example subdirectory. These may be useful for tutorial purposes, and general (partial) configs can be re-used for analyses with similar parameter requirements.

    # Activate the snakemake environment from above
    conda activate snakemake

    #Use the config creator script
    /path/to/repository/scripts/config_creator.py \
        --general_input config/ATAC_general.yaml \
        --per_lib_input data/atac_test_data/atac_test_samplesheet.csv \
        --results_dir /path/to/results \
        --temp_dir /path/to/results/tmp \
    > /path/to/results/test_config.yaml

This will generate the file `/path/to/results/test_config.yaml`.

It may be instructive to open the example general input, per-lib input, and the resulting configuration file to understand what the config creator script has done. Using modified versions of the example CSVs and/or general configs in this same manner should allow users to create pipeline-ready configuration files for their experiments.

    # Assuming the snakemake environment is activated

    # Running a dry-run (-n flag)
    snakemake -n --snakefile /path/to/repository/ATACseq.smk --configfile /path/to/results/test_config.yaml

    # If the dry-run succeeds, then proceed to run the ATACseq pipeline on the cluster

    # First start a persistent session with screen or tmux
    screen -S atac_test

    # Launch the pipeline
    snakemake -p \
        --snakefile /path/to/repository/ATACseq.smk \
        --configfile /path/to/results/test_config.yaml \
        --use-conda \
        --latency-wait 60 \
        --cluster-config /path/to/repository/config/cluster_config.json \
        --cluster 'sbatch \
            --job-name={cluster.name} \
            --account={cluster.account} \
            --partition={cluster.partition} \
            --nodes={cluster.nodes} \
            --ntasks-per-node={cluster.ntask} \
            --mem={cluster.memory} \
            --time={cluster.time} \
            --output=logs/%x-%j.out'

Results will be located where they were specified in the configuration - in this example, they are located in `/path/to/results`. These include bams (aligned, filtered), called peaks, display files, ataqv results, and cluster logs.

### Additional Information

#### Fastq inputs

In the script `config_creator.py`, one of the inputs is a CSV file which has a row for each sample, and among other things, a column containing a directory path which contains the FASTQs for a given sample. This configuration offers great flexibility and control. Various scenarios are exemplified here where we can take various actions to prepare inputs so that they are pipeline-ready.

##### A note on filename restrictions

One requirement for the input fastq files is that the read number must be discernable based on the filename (.e.g `_R1.fastq.gz`, `_1.fastq`, etc.), so that they can be appropriately placed as read 1 or read 2 in the configuration file. Since all three scenarios outlined here ultimately rely on the same config_creator script, they will all have this same restriction.

If your fastq files are single-end and do not explicitly have `_1` or `_R1` in the filenames, use the option `--no_capture` in `config_creator.py` to treat all fastqs in a subdirectory as read 1.


##### Fastq files already in separate directories for each sample

Example: Sequencing results from UMich Advanced Genomics Core

    umagc_data/
    ├── atac_test_samplesheet.csv
    ├── Sample_1
    │   ├── 1_L001_R1.fastq.gz
    │   └── 1_L001_R2.fastq.gz
    ├── Sample_2
    │   ├── 2_L001_R1.fastq.gz
    │   ├── 2_L001_R2.fastq.gz
    │   ├── 2_L002_R1.fastq.gz
    │   └── 2_L002_R2.fastq.gz
    ├── Sample_3
    │   ├── 3_L001_R1.fastq.gz
    │   ├── 3_L001_R2.fastq.gz
    │   ├── 3_L002_R1.fastq.gz
    │   ├── 3_L002_R2.fastq.gz
    │   ├── 3_L003_R1.fastq.gz
    │   └── 3_L003_R2.fastq.gz
    └── Sample_4
        ├── 4_L001_R1.fastq.gz
        └── 4_L001_R2.fastq.gz

This is the simplest case, where no extra steps need to be performed. The directories of fastq files can be fed to `config_creator.py` via the sample info CSV as exemplified in the quickstart example. The config creator script will gather the paths to the input fastqs to all of a sample's read1 or read2 inputs. The resulting information in the configuration file will eventually be used in the pipeline to regard these inputs as a sample.

##### Fastq files for all samples in a single directory, (somewhat) similarly named

Example: Sequencing runs downloaded from NCBI Sequence Read Archive

    sra_data/
    ├── SraRunTable.txt
    ├── SRR5799398_1.fastq
    ├── SRR5799398_2.fastq
    ├── SRR5799401_1.fastq
    ├── SRR5799401_2.fastq
    ├── SRR5799402_1.fastq
    ├── SRR5799402_2.fastq
    ├── SRR5799447_1.fastq
    └── SRR5799447_2.fastq

This is another common case, where the fastqs from multiple samples reside within a single directory. For this scenario, an additional step is taken to place the fastq files into separate subdirectories on a per-sample basis; we've created the prepare_fastq_inputs.py script to simplify this process. After the fastqs are placed in subdirectories, they can be fed to `config_creator.py` via the sample info CSV as exemplified in the main readme
In this case, an additional requirement is that filenames must be similar enough that they can be parsed uniformly to identify files to be grouped together as a sample. With this script, a table is supplied which maps sample names to the unique portion of the fastq filenames (the portion that remains after uniformly parsing them). In case of SRA files, we can group BioSamples which are comprised of multiple sequencing runs (multiple SRR accessions). In a more general use-case, we could group samples that have fastqs performed over multiple lanes of sequencing.

Note: prepare_fastq_inputs.py hardlinks the files into the proper locations, instead of moving them. This preserves the original flat structure in the given fastq directory in case that is desired, but also introduces in that location the nested structure needed for `config_creator.py`. These links do not consume any extra storage space.

##### Fastq files in heterogeneous directory structure or with filenames too disparate to parse uniformly

Example: Some samples from a collaborator, some samples internal, etc.

If your input samples are not in one of the two structures above, it is still possible to make them pipeline-ready.
Once again, the read number information must be discernible from the filename [as noted above](#a-note-on-filename-restrictions). However, that's the only restriction. The key objective is to get each sample's fastq files into their own directory. If you are able to manually create directories and move the appropriate fastqs into those locations, then that is the extent of the work to be done. After that, `config_creator.py` can be run exactly as described in the first scenario above, by supplying the appropriate fastq directory for each sample in the CSV.

#### Examples

1. [Example with SRA data](doc/Example_running_SE_ChIPseq_from_SRA.md) - running single-end ChIPseq through se pipeline
2. Example with in-house ChIPseq data - running paired-end ChIPseq reads through pe pipeline
