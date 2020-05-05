## Snakemake pipelines for the Ryan Lab

This repo contains snakefiles, configuration files, and a config generator script for processing ATACseq and ChIPseq experiments.

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
2. A comma-separated file with one line per sample giving the sample ID, the human-readable sample name, the genome for alignment, and the basepath to directory containing the fastqs for the sample (for example, `data/atac_test_data/atac_test_samplesheet.csv`).
3. A path indicating where files generated by the pipeline should go.
4. A temporary directory to use during certain

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

### Continuing with ATACseq pipeline

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

### Results

Results will be located where they were specified in the configuration - in this example, they are located in `/path/to/results`. These include bams (aligned, filtered), called peaks, display files, ataqv results, and cluster logs.

### Further reading and examples

1. [Example with SRA data](doc/Example_running_SE_ChIPseq_from_SRA.md) - running single-end ChIPseq through se pipeline
2. Example with in-house ChIPseq data - running paired-end ChIPseq reads through pe pipeline
