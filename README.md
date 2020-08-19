# Table of Contents

* [About](#about)
* [Before Starting](#before-starting)
    * [Managing Software with Conda](#managing-software-with-conda)
    * [Installing the Pipelines](#installing-the-pipelines)
        * [Lab-Specific Customizations](#lab-specific-customizations)
    * [Genome Reference Requirements](#genome-reference-requirements)
* [Quick-Start ATAC](#quick-start-example-processing-raw-reads-from-an-atac-seq-experiment)
* [Quick-Start ChIP](#quick-start-example-processing-raw-reads-from-chip-seq-experiment)
* [Quick-Start Homer Only](#quick-start-example-homer-only)
* [Additional Information](#additional-information)
    * [Genome Reference Files](#genome-reference-files)
    * [Notes on Control Samples for ChIPseq](#notes-on-control-samples-for-chipseq)
    * [Notes on conditional outputs for ChIPseq](#notes-on-conditional-outputs-for-chipseq)
    * [Testing additional peak parameters](#testing-additional-peak-parameters)
    * [Tuning cluster resource requirements](#tuning-cluster-resource-requirements)
    * [Fastq Inputs](#fastq-inputs)

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

After installation, create a parent conda environment for running the pipeline with the command below. Note that we're setting a prefix pointing to an `/nfs/turbo` location; This will place this single environment into that location, which may be separate from your other conda environments.
This due to the Homer package, which will require a large amount of space. Command:

    conda create \
        --prefix /nfs/turbo/<lab-turbo>/envs/atac_chip_pipeline \
        --channel conda-forge --channel bioconda \
        snakemake-minimal=5.10.0 pandas=0.24.2 homer=4.10 samtools=1.9 pytest=5.4.1 nose=1.3.7

To use the pipeline, activate the environment by running:

    conda activate /nfs/turbo/<lab-turbo>/envs/atac_chip_pipeline

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
Peak calling in the ChIP-seq pipeline is done with [Homer](http://homer.ucsd.edu/homer/). Homer requires some additional setup of genome references, outlined below. The resulting Homer references will be installed in the `atac_chip_pipeline` environment. It is a peculiarity of installing Homer with Conda that this step must be done separately. We define the environment variable $HOMER_DIR in the following example to flexibly handle different conda environment locations.

    # Activate the pipeline atac_chip_pipeline
    conda activate /nfs/turbo/<lab-turbo>/envs/atac_chip_pipeline

    # Define the location of the homer directory (This is relative to the conda installation)
    HOMER_DIR=$(dirname $(dirname $(which findMotifsGenome.pl)))/share/homer-4.10-0

    # List all the genomes
    perl $HOMER_DIR/configureHomer.pl -list

    # Install hg19 (install other genomes in the same manner)
    perl $HOMER_DIR/configureHomer.pl -keepScript -install hg19

##### Lab-specific Customizations

After installation, run the following script to adjust the default base paths in various files. For the `--hpc_account` flag, users can see which accounts they have access to by calling `my_accounts` on Great Lakes.

    cd /path/to/repository
    python setup.py --hpc_account account_name --env_path /nfs/turbo/<lab-turbo>/envs/atac_chip_pipeline

Finally, change the following items

- In each of the `yaml` files in `config/`, change BWA indices to appropriate paths.

Note, if the pipeline folder is moved after running `setup.py`, the internal paths will no longer be valid. Instead of moving the folder, users should follow the steps above in the new desired location.

##### Quick Test of the Installation

If all the above are configured correctly, it should be possible to run a simple shell script to test the installation of the pipelines. These will create a test configuration file and call the pipeline with included test data.

To test the ATACseq pipeline, simply run:

    user_tests/atac_test.sh

To test the ChIPseq pipeline, run:

    user_tests/chip_test.sh

The tests can take anywhere from 15 to 60 minutes.

#### Genome Reference Requirements

The following files are required for the pipelines to run properly:

- BWA indices
- Blacklist regions
- TSS locations (ATAC-seq only)

All but the BWA indices are included in the repository, and are located in `data/references`, and the sources are described in [`data/references/sources.md`](data/references/sources.md). We describe how to download BWA indices below. For more detail about blacklist regions, and TSS locations see the section on [Genome Reference Files](#genome-reference-files).

##### BWA Indices

Most iGenomes references ([link](https://support.illumina.com/sequencing/sequencing_software/igenome.html)) come with pre-built BWA indices. For example, in the hg19 download from iGenomes, the BWA index is located in `Homo_sapiens/UCSC/hg19/Sequence/BWAIndex`.

To build a BWA index for genomes not downloaded from iGenomes, see the [bwa manual](http://bio-bwa.sourceforge.net/bwa.shtml). Briefly,

    bwa index ref.fa

#### Input Files

Fastq inputs must have filenames in the form `_R1.fastq.gz`, `_R2.fastq.gz`, `_R1_001.fastq.gz`, `_R2_001.fastq.gz`, etc. This is so the config creator script can determine read number based on the filename. Specifically, the filenames should each contain an `_R1.`, `_R1_`, `_R2.`, or `_R2_` element.

Additionally, all fastqs for a given sample must be contained in their own directories. These directories are listed in the sample sheet that is given to `config_creator.py`.

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

[See the section with helpful notes and tools](#fastq-inputs) for situations where your input files are not yet 'pipeline-ready'.

### Quick-start Example: Processing Raw Reads from an ATAC-seq Experiment

We create some bash environment variables to point to the conda environment, repository directory, and results directory.

    conda_env_dir=/path/to/atac_chip_pipeline
    repo_dir=/path/to/ChIPseq_ATACseq_pipelines
    results_dir=/path/to/atac_results

Create the directory where you'd like to save the results of the pipeline:

    mkdir -p ${results_dir}
    cd ${results_dir}

Next, generate a pipeline configuration file using the `scripts/config_creator.py` script. The following need to be specified:

1. General configuration details in the form of a partial configuration file (for example, `config/ATAC_general.yaml`).
2. A comma-separated file with one line per sample giving the sample ID, the human-readable sample name, the genome for alignment, and the basepath to directory containing the fastqs for the sample (for example, `data/atac_test_data/atac_test_samplesheet.csv`, and below).
3. A path indicating where files generated by the pipeline should go.
4. A temporary directory to use during certain rules.

    |lib|sample|genome|basepath|
    |---|------|------|--------|
    |1|s1_control|mm10|data/atac_test_data/Sample_1|
    |2|s1_treatment|mm10|data/atac_test_data/Sample_2|
    |3|s2_treatment|mm10|data/atac_test_data/Sample_3|
    |4|s3_treatment|mm10|data/atac_test_data/Sample_4|

The configuration file is then created by running:

    # Activate the snakemake environment from above
    conda activate ${conda_env_dir}

    #Use the config creator script
    ${repo_dir}/scripts/config_creator.py \
        --general_input ${repo_dir}/config/ATAC_general.yaml \
        --per_lib_input ${repo_dir}/data/atac_test_data/atac_test_samplesheet.csv \
        --results_dir ${results_dir} \
        --temp_dir ${results_dir}/tmp \
    > ${results_dir}/test_atac_config.yaml

This will generate the file `${results_dir}/test_atac_config.yaml`.

It may be instructive to open the example general input, per-lib input, and the resulting configuration file to understand what the config creator script has done.

    # Assuming the atac_chip_pipeline environment is activated

    # Running a dry-run (-n flag)
    snakemake -n --snakefile ${repo_dir}/ATACseq.smk --configfile ${results_dir}/test_atac_config.yaml

    # If the dry-run succeeds, then proceed to run the ATACseq pipeline on the cluster

    # First start a persistent session with screen or tmux
    screen -S atac_test

    # Launch the pipeline
    snakemake -p \
        --snakefile ${repo_dir}/ATACseq.smk \
        --configfile ${results_dir}/test_atac_config.yaml \
        --use-conda \
        --latency-wait 60 \
        --jobs 100 \
        --cluster-status ${repo_dir}/scripts/slurm_status.py \
        --cluster-config ${repo_dir}/config/cluster_config.yaml \
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

Results will be located where they were specified in the configuration - in this example, they are located in `/path/to/atac_results`. These include bams (aligned, filtered), called peaks, display files, ataqv results, and cluster logs.

### Quick-start Example: Processing Raw Reads from ChIP-seq Experiment

Running ChIP-seq experiments from raw reads is similar to ATAC-seq with some minor changes.

1. The samplesheet passed to the `config_creator.py` script. In particular, the samplesheet should look like:

    |lib|sample|genome|input|homer_fmg_genome|basepath|
    |---|------|------|-----|----------------|--------|
    |GM12878_NKRF|GM12878_NKRF|hg19|GM12878_Input|hg19r|data/sra_chip_test_data/GM12878_NKRF/|
    |GM12878_Input|GM12878_Input|hg19|||data/sra_chip_test_data/GM12878_Input/|

    Note the addition of the `input` and `homer_fmg_genome` columns. The `input` column indicates which library `sample` is to be used as the control sample when finding peaks with Homer. The `homer_fmg_genome` indicates the Homer coded genome reference that is to be used when finding motifs, if that rule is applicable. For more additional information about the `input` column see [notes on control samples for ChIPseq](#notes-on-control-samples-for-chipseq).

2. The `--general_input` flag to the `config_creator.py` function has four options for ChIPseq, depending on whether the data is paired- or single-end and whether the data is a histone or transcription factor (TF). The options are:

    - `config/ChIP_histone_general_pe.yaml`,
    - `config/ChIP_histone_general_se.yaml`,
    - `config/ChIP_TF_general_pe.yaml`, or
    - `config/ChIP_TF_general_se.yaml`

We create some bash environment variables to point to the conda environment, repository directory, and results directory.

    conda_env_dir=/path/to/atac_chip_pipeline
    repo_dir=/path/to/ChIPseq_ATACseq_pipelines
    results_dir=/path/to/chip_results

The configuration file is then created by running the following, similar to ATACseq:

    # Activate the snakemake environment from above
    conda activate ${conda_env_dir}

    #Use the config creator script
    ${repo_dir}/scripts/config_creator.py \
        --general_input ${repo_dir}/config/ChIP_histone_pe.yaml \
        --per_lib_input ${repo_dir}/data/sra_chip_test_data/sra_chip_samplesheet.csv \
        --results_dir ${results_dir} \
        --temp_dir ${results_dir}/tmp \
    > ${results_dir}/test_chip_config.yaml

This will generate the file `/path/to/results/test_chip_config.yaml`. The pipeline can now be executed. One difference from ATACseq is that there is a paired-end pipeline (`ChIPseq_pe.smk`) and a single-end pipeline (`ChIPseq_se.smk`). Be sure to select the correct one.

    # Assuming the atac_chip_pipeline environment is activated

    # Running a dry-run (-n flag)
    snakemake -n --snakefile ${repo_dir}/ChIPseq_pe.smk --configfile ${results_dir}/test_chip_config.yaml

    # If the dry-run succeeds, then proceed to run the ChIPseq pipeline on the cluster

    # First start a persistent session with screen or tmux
    screen -S chip_test

    # Launch the pipeline
    snakemake -p \
        --snakefile ${repo_dir}/ChIPseq_pe.smk \
        --configfile ${results_dir}/test_chip_config.yaml \
        --use-conda \
        --latency-wait 60 \
        --jobs 100 \
        --cluster-status ${repo_dir}/scripts/slurm_status.py \
        --cluster-config ${repo_dir}/config/cluster_config.yaml \
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

Results will be located where they were specified in the configuration - in this example, they are located in `/path/to/chip_results`. These include bams (aligned, filtered), called peaks, found motifs, display files, and cluster logs. It is possible to alter the flow of the pipeline to exclude peak calling (should a different peak caller be preferred) or motif finding. See notes on conditional outputs for ChIPseq, [below](#notes-on-conditional-outputs-for-chipseq). It is also possible to explore the effect of different peak calling parameters with Homer in one pipeline run, see the section on testing additional homer findPeaks parameters, [below](#testing-additional-homer-findpeaks-parameters).

### Quick-Start Example: Homer Only

It is possible to start with filtered bam files and perform only the downstream steps makeTagDirectory, findPeaks, blacklist filtering, and findMotifsGenome. A constraint for this use-case is that the input bam files must be named and structured in a way that is compatible with the pipeline at large.
The filtered bams are assumed to have only properly aligned reads or read pairs, and removed PCR/optical duplicates.

Input directory structure:

    ${result_dir}/pruned/
    ├── OCILY_H3K27ac.pruned.bam
    └── OCILY_input.pruned.bam

The filtered bam filenames should be `sample.pruned.bam`, and should be located in a directory called `pruned`. The `pruned` directory should be a subdirectory of the results path given to the config_creator script.

The config_creator script has a flag `--homer_only` which will not perform any fastq-finding steps, and will instead produce a configuration file that has only the keys needed for that pipeline.

    ${repo_dir}/scripts/config_creator.py \
        --general_input ${repo_dir}/config/ChIP_histone_general_se.yaml \
        --per_lib_input ${repo_dir}/data/sra_mixed_chip_test_data/sra_chip_histone_all.csv \
        --results_dir ${result_dir} \
        --temp_dir ${result_dir}/tmp \
        --homer_only \
        > ${result_dir}/config_histone_homer.yaml

Note: The homer_only config won't end up containing any se or pe params, so either can be used to provide the general input information.

### Additional Information

#### Genome Reference Files

While blacklist and TSS files for hg19, hg38, and mm10 are included in the `data/references` folder, we describe the form they should take below should the need arise for custom versions of these files. Remember that the pipeline, by default, uses those in `data/references` so care should be taken to make sure general or project configurations are using the desired files.

##### Blacklist Regions

Blacklist regions should be in the form of a [BED](https://genome.ucsc.edu/FAQ/FAQformat.html#format1) file. Version 1 and 2 blacklist regions are included in the repository, but can be downloaded from the Boyle Lab ([link](https://github.com/Boyle-Lab/Blacklist/tree/master/lists)).

If non-standard blacklist regions are required, ensure that the regions are mutually disjoint. This can be done with [`bedtools merge`](https://bedtools.readthedocs.io/en/latest/content/tools/merge.html).

Note: If you are using an organism that does not have a blacklist available, create a dummy placeholder blacklist and include its path in the general config files (e.g. `config/ATAC_general.yaml`). To do this, you can use the following command:

    echo 'chr1\t1\t2' > /path/to/dummy_blacklist.bed # Note: use appropriate chromosome identifiers for your organism.

An example blacklist bedfile:

    chr1    0       750100  High Signal Region
    chr1    814500  845200  High Signal Region
    chr1    2052400 2056000 High Signal Region
    ...

##### TSS File

TSS BED files are used to calculate TSS enrichment scores in the `ataqv` quality asseessment. The `ataqv` respository contains ready-made TSS files for hg19, mm9, and rn5 ([link](https://github.com/ParkerLab/ataqv/tree/master/data/tss)), but also includes the code used to generate them as an example for other references. We have used this code to generate hg38 and mm10 TSSs.

An example TSS bedfile:

    chr1	11874	11874
    chr1	69091	69091
    chr1	321084	321084
    ...

#### Notes on control samples for ChIPseq

The `input` column in the ChIPseq samplesheet allows users to use the same control for more than one pulldown library, or to use matched controls and pulldowns. Simply use the appropriate library `sample` entry.

#### Notes on conditional outputs for ChIPseq

The ChIPseq pipeline is flexible in the sense that including or excluding elements of the samplesheet can determine what steps of the pipeline are run. Two examples for the example samplesheet below:

|lib|sample|genome|input|homer_fmg_genome|basepath|
|---|------|------|-----|----------------|--------|
|GM12878_NKRF|GM12878_NKRF|hg19|GM12878_Input|hg19r|data/sra_chip_test_data/GM12878_NKRF/|
|GM12878_Input|GM12878_Input|hg19|||data/sra_chip_test_data/GM12878_Input/|

1. If the `input` column is left blank for all samples, peaks will not be called.
2. If the `homer_fmg_genome` column is left blank all samples, motifs will not be found.

Alternatively, the resulting configuration file from `config_creator.py` could be altered to achieve the same result:

    sample_homer_fmg_genome:
        GM12878_NKRF: hg19r
    sample_input:
        GM12878_NKRF: GM12878_Input

1. Remove the `sample_input` block.
2. Remove the `sample_homer_fmg_genome` attribute.

#### Testing additional peak parameters

It may be useful to try different parameters when running macs2 in the ATAC-seq pipeline or Homer findPeaks in the ChIP-seq pipeline. To achieve this, additional sets of parameters can be specified in the configuration file, and they will be placed in separate subfolders adjacent to the default results. For example, in the ATAC-seq case:

    macs2_params:
        default_params: "-f BAMPE --keep-dup all"
        new_params: "-f BAMPE --keep-dup all --qvalue 0.1"

And in the ChIP-seq case:

    homer_findPeaks_params:
        default_params: "-style histone"
        new_params: "-style histone -minDist 150 -fdr 0.01"

Will create two sets of results and place them in their corresponding subfolders. Downstream steps will also be placed in subfolders. In particular, for each parameter set there will be an ATAQV report (in the ATAC-seq case) and Homer report (in the ChIP-seq case).

In the ChIP-seq case, the peaks folder will have structure as below. The ATAC-seq case is similar.

    /path/to/results/peaks
    ├── default_params
    │   ├── OCILY_H3K27ac.all.bed
    │   ├── OCILY_H3K27ac.all.hpeaks
    │   ├── OCILY_H3K27ac.BLfiltered.bed
    │   └── OCILY_H3K27ac.BLfiltered.hpeaks
    └── new_params
        ├── OCILY_H3K27ac.all.bed
        ├── OCILY_H3K27ac.all.hpeaks
        ├── OCILY_H3K27ac.BLfiltered.bed
        └── OCILY_H3K27ac.BLfiltered.hpeaks

A warning, when selecting parameters for either macs2 or Homer findPeaks, care must be taken to avoid parameters which change the actual files that are output. For example, the `--broad` flag in macs2 will cause the pipeline to fail because the output will be `.broadPeak` instead of the expected `.narrowPeak`.

#### Tuning cluster resource requirements

The cluster resources enumerated in `config/cluster_config.yaml` are set with values which should handle normal datasets by default. If the data you're using are much larger, it may make sense to increase these values if your jobs are killed for exceeding the walltime. Setting longer times will not incur extra cost if the jobs complete before their scheduled time, but it may potentially cause jobs to wait longer in the queue.

#### Fastq Inputs

In the `config_creator.py` script, one of the inputs is a CSV file which has a row for each sample, and among other things, a column containing a directory path which contains the FASTQs for a given sample. This configuration offers great flexibility and control. Various scenarios are exemplified here where we can take various actions to prepare inputs so that they are pipeline-ready.

##### Notes on Filename Restrictions

Fastq files must be gzipped, and must have the extension `.fastq.gz` for the config_creator script and the pipeline to work.

To gzip all fastqs in a directory recursively, you can use the following one-liner. This will work in all cases - if fastqs are all in a single directory, if they are located within subdirectories, or a mixture of the two.

    # Assuming fastq filenames end in ".fastq" here: If they end in ".fq", change the argument to -name "*.fq"
    find /path/to/fastq_dir -type f -name "*.fastq" -exec gzip {} \;

If your file extensions are not `.fastq.gz`, for example if they are `.fq.gz`, you can rename them all with the following command. Once again, this operates recursively in the fastq directory.

    find /path/to/fastq_dir -type f -name "*.fq.gz" -exec rename .fq.gz .fastq.gz {} \;

Another requirement for the input fastq files is that the read number must be discernable based on the filename in the form `_R1.fastq.gz` or `_R2.fastq.gz`, so that they can be appropriately placed as read 1 or read 2 in the configuration file. See the [input Files](#input-files) section.

##### Fastq Files Already in Separate Directories for Each Sample

This is the simplest case, where no extra steps need to be performed. The directories of fastq files can be fed to `config_creator.py` via the sample info CSV as exemplified in the quickstart example. The config creator script will gather the paths to the input fastqs to all of a sample's read1 or read2 inputs. The resulting information in the configuration file will eventually be used in the pipeline to regard these inputs as a sample.

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

##### Fastq Files for All Samples in a Single Directory

Example: Sequencing runs downloaded from NCBI Sequence Read Archive (SRA)

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

This is another common case, where the fastqs from multiple samples reside within a single directory. For this scenario, an additional step is taken to place the fastq files into separate subdirectories on a per-sample basis; we've created the prepare_fastq_inputs.py script to simplify this process. After the fastqs are placed in subdirectories, they can be fed to `config_creator.py` via the sample info CSV as exemplified in the examples above.

In this case, an additional requirement is that filenames must be similar enough that they can be parsed uniformly to identify files to be grouped together as a sample. With this script, a table is supplied which maps sample names to the unique portion of the fastq filenames (the portion that remains after uniformly parsing them). In case of SRA files, we can group BioSamples which are comprised of multiple sequencing runs (multiple SRR accessions). In a more general use-case, we could group samples that have fastqs performed over multiple lanes of sequencing.

Note: prepare_fastq_inputs.py hardlinks the files into the proper locations, instead of moving them. This preserves the original flat structure in the given fastq directory in case that is desired, but also introduces in that location the nested structure needed for `config_creator.py`. These links do not consume any extra storage space.

###### Example downloading and preparing SRA fastqs:

    cd SRA_data
    # Download them (Just download 10M spots for this example with -N & -X)
    for i in $(seq 6730191 1 6730200) ; do echo $i ; fastq-dump -F -N 10000 -X 10010000 --split-files SRR${i} ; done
    # Run prepare_fastq_inputs.py - SraRunTable
    cd ../
    ${repo_dir}/scripts/prepare_fastq_inputs.py -d SRA_data/ -m SraRunTable.txt --add_R

Note: When downloading fastqs from SRA using `fastq-dump`, supplying the `-F` flag preserves headers, preventing the use of SRA's read names. SRA read names with suffixes `.1` and `.2` cause bwa sampe to throw errors.

##### Fastq Files with Heterogeneous Directory Structure or Filenames

Example: Some samples from a collaborator, some samples internal, etc.

If your input samples are not in one of the two structures above, it is still possible to make them pipeline-ready.

Once again, the read number information must be discernible from the filename [as noted above](#notes-on-filename-restrictions). However, that's the only restriction. The key objective is to get each sample's fastq files into their own directory. If you are able to manually create directories and move the appropriate fastqs into those locations, then that is the extent of the work to be done. After that, `config_creator.py` can be run exactly as described in the first scenario above, by supplying the appropriate fastq directory for each sample in the CSV.
