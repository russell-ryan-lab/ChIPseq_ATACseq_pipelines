## Snakemake pipelines for the Ryan Lab

This repo contains snakefiles, configuration files, and a config generator script for processing ATACseq and ChIPseq experiments

These snakefiles are designed so that input fastqs and their associated wildcards are explicitly defined in the config, rather than implicitly determined by filename. This provides the flexibility allowing input files to be named in any fashion, and stored in any location. The easiest way to process a new experiment after receiving the data from the UMich sequencing core is to use the config generator script. The config generator script can be modified as needed to facilitate loading of differently named/stored inputs, or configs could even be manually edited

### Before Starting

The pipelines manage software dependencies using [conda](https://docs.conda.io/en/latest/miniconda.html). In order to get started, you'll need a conda environment with snakemake and pandas installed. If you don't have conda installed, you can install miniconda to your home directory on GreatLakes by following the instructions below, based on the [miniconda installation instructions](https://docs.conda.io/en/latest/miniconda.html):

    # Download for 64-bit linux to your home directory:
    cd ~
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
    # Run the installer (choosing default options is fine)
    bash Miniconda3-latest-Linux-x86_64.sh


After that is set up, you can create a snakemake environment by running:

    conda create -n snakemake -c bioconda snakemake=5.10.0 pandas=0.24.2

Then when you want you use snakemake, you can activate that environment by running:

    conda activate snakemake

At that point, you should be able to run snakemake. Test it by running:

    snakemake --version

The result should be `5.10.0`, the latest snakemake release (as of 2/20/2020)


You can deactivate the environment with:

    conda deactivate


### Lab-specific customizations

The file cluster_config.yaml contains information that high-performance computing clusters need for submitting jobs. One field that other labs using these pipelines must modify is the "account" section of cluster_config.json. This should be changed from `rjhryan` to the lab's own HPC account.

The general configuration files within `example/` contain paths to references (.e.g bwa indices, blacklist regions, etc.) which are managed by individuals. These paths should also be updated with new locations that are accessible to your lab.

TODO: add links to web locations with e.g. blacklists.

### Installing

    #Navigate to where you want this pipeline to be located
    cd /nfs/turbo/path-rjhryan-turbo/lab-members/Travis/
    git clone https://github.com/russell-ryan-lab/ChIPseq_ATACseq_pipelines

### Quick-start example: Processing the results of ATACseq experiment

First, generate a config file. The ngs_rawdata_config_creator script will automatically generate a snakemake-ready configuration file by combining:

1. General configuration details from a partial config file
2. Parameters specified for each separate library in a .csv file.

Along with this quick example, there are several more configuration files which are included in the example subdirectory. These may be useful for tutorial purposes, and general (partial) configs can be re-used for analyses with similar parameter requirements.

    #Use the config creator script
    ./scripts/ngs_rawdata_config_creator.py --general_input example/ATAC_general.json \
    --per_lib_input example/ATAC_2369_abridged_sample_info.csv \
    --results_dir /scratch/rjhryan_root/rjhryan/trsaari/example_ATAC_2369/ \
    --temp_dir /scratch/rjhryan_root/rjhryan/trsaari/tmp/ \
      > example/ATAC_2369_config.json

This will generate the file example/ATAC_2369_config.json

It may be instructive to open the example general input, per-lib input, and the resulting configuration file to understand what the config creator script has done. Using modified versions of the example CSVs and/or general configs in this same manner should allow users to create pipeline-ready configuration files for their experiments.

### Continuing with ATACseq pipeline

    #Running a dry-run (-n flag)
    snakemake -n --snakefile Snakefile_ATACseq --configfile example/ATAC_2369_config.json

    #If the dry-run succeeds, then proceed to run the ATACseq pipeline on the cluster
    #First start a persistent session with screen or tmux
    screen
    #Then launch the pipeline
    snakemake -p --use-conda --snakefile Snakefile_ATACseq --configfile example/ATAC_2369_config.json \
    --latency-wait 60 --cluster-config cluster_config.json \
    --cluster 'sbatch --job-name={cluster.name} --account={cluster.account} --partition={cluster.partition} --nodes={cluster.nodes} --ntasks-per-node={cluster.ntask} --mem={cluster.memory} --time={cluster.time} --output=%x-%j.out'

### Results

Results will be located where they were specified in the config - in this example, they are located in `/scratch/rjhryan_root/rjhryan/trsaari/example_ATAC_2369/`. These include bams (aligned, filtered), called peaks, display files, ataqv results, and cluster logs.

### Further reading and examples

1. [Example with SRA data](doc/Example_running_SE_ChIPseq_from_SRA.md) - running single-end ChIPseq through se pipeline
2. Example with in-house ChIPseq data - running paired-end ChIPseq reads through pe pipeline
