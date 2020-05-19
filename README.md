# Table of Contents

* [About](#about)
* [Before Starting](#before-starting)
  * [Managing Software with Conda](#managing-software-with-conda)
  * [Installing the Pipelines](#installing-the-pipelines)
    * [Lab-Specific Customizations](#lab-specific-customizations)
  * [Genome Reference Requirements](#genome-reference-requirements)
* [Quick-Start Example](#quick-start-example-Processing-raw-reads-from-an-ATAC-seq-experiment)
* [Additional Information](#additional-information)
  * [Fastq Inputs](#fastq-inputs)
  * [Examples](#examples)

## About

This repository contains snakefiles, configuration files, and a config generator script for processing ATACseq and ChIPseq experiments.

The ATACseq pipeline takes raw fastq files as input, performs alignment, filtering, and QC, calls peaks using MACS2, and produces bigwig files for viewing.

The ChIPseq pipeline has two portions. The first portion is similar to the ATAC-seq pipeline, in that it performs alignment, filtering, and produces bigwig files (though the tools and procedures sometimes differ from ATACseq). The second portion takes filtered bam files, and performs peak-calling and motif enrichment steps. These two portions can be run separately. See below for details.


### Before Starting

#### Managing Software with Conda

The pipelines manage software dependencies using [conda](https://docs.conda.io/en/latest/miniconda.html). To install miniconda to your home directory on Great Lakes, follow the instructions below; they are based on the [miniconda installation instructions](https://docs.conda.io/en/latest/miniconda.html).

    # Download for 64-bit linux to your home directory:
    cd ~
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

    # Run the installer (choosing default options is fine)
    bash Miniconda3-latest-Linux-x86_64.sh

    # You may need to source the .bashrc file to initialize conda after it is first installed
    source ~/.bashrc

After installation, create a parent conda environment for running the pipeline with:

    conda create -n atac_chip_pipeline -c conda-forge -c bioconda snakemake-minimal=5.10.0 pandas=0.24.2 pytest=5.4.1 nose=1.3.7

To use the pipeline, activate the environment by running:

    conda activate atac_chip_pipeline

Check that snakemake was installed by running:

    snakemake --version

The result should be `5.10.0`, the latest snakemake release (as of 2/20/2020).

Finally, deactivate the environment with:

    conda deactivate

#### Installing the Pipelines

    #Navigate to where you want this pipeline to be located
    cd /nfs/turbo/<lab-turbo>/pipelines

    git clone https://github.com/russell-ryan-lab/ChIPseq_ATACseq_pipelines

##### Setup Homer
Peak calling in the ChIP-seq pipeline is done with [Homer](http://homer.ucsd.edu/homer/). Homer requires some genome reference information of its own for peak calling and motif finding. This requires some additional setup, outlined below. The resulting Homer references will be installed in the `atac_chip_pipeline` environment. This is a peculiarity of installing Homer with Conda. We define the environment variable $HOMER_DIR in the following example to flexibly handle different conda installation locations. If miniconda was installed in your home directory as outlined [above](#managing-software-with-conda), then `$HOMER_DIR` will be `~/miniconda3/envs/atac_chip_pipeline/share/homer-4.10-0`.

    # Activate the pipeline atac_chip_pipeline
    conda activate atac_chip_pipeline

    # Define the location of the homer directory (This is relative to the conda installation)
    HOMER_DIR=$(dirname $(dirname $(which findMotifsGenome.pl)))/share/homer-4.10-0

    # List all the genomes
    perl $HOMER_DIR/configureHomer.pl -list

    # Install hg19 (install other genomes in the same manner)
    perl $HOMER_DIR/configureHomer.pl -keepScript -install hg19

##### Lab-specific customizations

After installation, the following paths need to be changed to be lab-specific:

- In `data/atac_test_data/atac_test_samplesheet.csv`, the `/nfs/med-bfx-activeprojects/Ryan_rjhryan_CU1/ChIPseq_ATACseq_pipelines` part of each basepath should be changed to point to the installation path of the pipeline.
- In each of the `config/*_general.yaml` files, the paths in the `blacklist`, `chrom_sizes`, and `tss` sections should be modified as follows. Paths beginning with `/nfs/med-bfx-activeprojects/Ryan_rjhryan_CU1/ChIPseq_ATACseq_pipelines` should be changed to point to the installation path of the pipeline.
- In each of the `config/*_general.yaml` files, the `bwa_index` section should point to a valid BWA index created by the user. See the [section](#bwa-indices) below for how to create a BWA index.
- In `config/cluster_config.yaml`, the `account` attribute should be changed to the appropriate Great Lakes account.
- In `tests/atac_test.sh`, lines 1-3 should be changed to point to the installation path of the pipeline.

##### Quick Test of the installation

If all the above are configured correctly, it should be possible to run a simple shell script to test the installation of the pipelines. These will create a test configuration file and call the pipeline with included test data.

To test the ATACseq pipeline, simply run:

    tests/atac_test.sh

At the end of the run (~15 min) snakemake should indicate that it is finished by printing:

    54 of 54 steps (100%) done

#### Genome Reference Requirements

Prior to using the pipeline, BWA indices, chromosome sizes, and blacklist regions are required for both the ATAC-seq and ChIP-seq pipelines, and the ATAC-seq pipeline requires a TSS file. Each is described in detail below, and they are included in `data/references` for hg19, hg38, and mm10. The sources are described in [`data/references/sources.md`](data/references/sources.md).

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
    ...

##### Blacklist Regions

Blacklist regions should be in the form of a [BED](https://genome.ucsc.edu/FAQ/FAQformat.html#format1) file. Version 2 blacklist regions can be downloaded from the Boyle Lab ([link](https://github.com/Boyle-Lab/Blacklist/tree/master/lists)). Version 1 blacklist regions can be downloaded from the Kundaje Lab ([link](https://sites.google.com/site/anshulkundaje/projects/blacklists)).

The blacklist filtering step requires that the BED file used contains mutually disjoint regions. To ensure this is the case, one can use [`bedtools merge`](https://bedtools.readthedocs.io/en/latest/content/tools/merge.html) to create a final BED file for use in the pipeline.

Note: If you are using an organism that does not have a blacklist available, it is recommended to create a dummy placeholder blacklist and include its path in general config files (e.g. `config/ATAC_general.yaml`). To do this, you can use the following command:

    echo 'chr1\t1\t2' > /path/to/dummy_blacklist.bed #Note: use appropriate chromosome identifiers for your organism.

An example:

    chr1	564449	570371
    chr1	724136	727043
    chr1	825006	825115
    ...

##### TSS File

TSS BED files are used to calculate TSS enrichment scores in the `ataqv` quality asseessment. The `ataqv` respository contains ready-made TSS files for hg19, mm9, and rn5 ([link](https://github.com/ParkerLab/ataqv/tree/master/data/tss)), but also includes the code used to generate them as an example for other references.

An example:

    chr1	11874	11874
    chr1	69091	69091
    chr1	321084	321084
    ...

#### Input files

Fastq inputs must have filenames in the form `_R1.fastq.gz` or `_R2.fastq.gz`, and all fastqs for a given sample must be contained in their own directories. These directories are listed in the sample sheet that is given to config_creator.py. [Here is a link to some helpful notes and tools](#fastq-inputs) available if your input files are not yet 'pipeline-ready'.

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
    conda activate atac_chip_pipeline

    #Use the config creator script
    /path/to/repository/scripts/config_creator.py \
        --general_input config/ATAC_general.yaml \
        --per_lib_input data/atac_test_data/atac_test_samplesheet.csv \
        --results_dir /path/to/results \
        --temp_dir /path/to/results/tmp \
    > /path/to/results/test_config.yaml

This will generate the file `/path/to/results/test_config.yaml`.

It may be instructive to open the example general input, per-lib input, and the resulting configuration file to understand what the config creator script has done. Using modified versions of the example CSVs and/or general configs in this same manner should allow users to create pipeline-ready configuration files for their experiments.

    # Assuming the atac_chip_pipeline environment is activated

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
        --cluster-status /path/to/repository/scripts/slurm_status.py \
        --cluster-config /path/to/repository/config/cluster_config.json \
        --cluster 'sbatch \
            --job-name={cluster.name} \
            --account={cluster.account} \
            --partition={cluster.partition} \
            --nodes={cluster.nodes} \
            --ntasks-per-node={cluster.ntask} \
            --mem={cluster.memory} \
            --time={cluster.time} \
            --parsable \
            --output=logs/%x-%j.out'

Results will be located where they were specified in the configuration - in this example, they are located in `/path/to/results`. These include bams (aligned, filtered), called peaks, display files, ataqv results, and cluster logs.

### Additional Information

#### Fastq inputs

In the script `config_creator.py`, one of the inputs is a CSV file which has a row for each sample, and among other things, a column containing a directory path which contains the FASTQs for a given sample. This configuration offers great flexibility and control. Various scenarios are exemplified here where we can take various actions to prepare inputs so that they are pipeline-ready.

##### Notes on filename restrictions

Fastq files must be gzipped, and must have the extension `.fastq.gz` for the config_creator script and the pipeline to work.

To gzip all fastqs in a directory recursively, you can use the following one-liner. This will work in all cases - if fastqs are all in a single directory, if they are located within subdirectories, or a mixture of the two.

    # Assuming fastq filenames end in ".fastq" here: If they end in ".fq", change the argument to -name "*.fq"
    find /path/to/fastq_dir -type f -name "*.fastq" -exec gzip {} \;

If your file extensions are not `.fastq.gz`, for example if they are `.fq.gz`, you can rename them all with the following command. Once again, this operates recursively in the fastq directory.

    find /path/to/fastq_dir -type f -name "*.fq.gz" -exec rename .fq.gz .fastq.gz {} \;


Another requirement for the input fastq files is that the read number must be discernable based on the filename in the form `_R1.fastq.gz` or `_R2.fastq.gz`, so that they can be appropriately placed as read 1 or read 2 in the configuration file.

[TODO: flesh this out]



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
