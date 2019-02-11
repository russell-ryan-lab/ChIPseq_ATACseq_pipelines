## Snakemake NGS pipelines for the Ryan Lab

This repo contains snakefiles, configuration files, and a config generator script for processing NGS experiments

These snakefiles are designed so that input fastqs and their associated wildcards are explicitly defined in the config, rather than implicitly determined by filename. This provides the flexibility allowing input files to be named in any fashion, and stored in any location. The easiest way to process a new experiment after receiving the data from the UMich sequencing core is to use the config generator script. The config generator script can be modified as needed to facilitate loading of differently named/stored inputs, or configs could even be manually edited.

#### Example: Processing the results of ATACseq experiment

##### First, generate a config file:

    ./ngs_rawdata_config_creator.py --general_input example/config_general.json --per_lib_input example/ATAC_2369_sample_info.csv --results_dir /scratch/rjhryan_fluxod/trsaari/example_ATAC_2369/ --temp_dir /scratch/rjhryan_fluxod/trsaari/tmp/ > example/ATAC_2369_config.json

This will generate the config ATAC_2369_config.json

Next, run the fastqc snakefile - this combines the fastq inputs for each sample, merges them, and runs fastqc on the merged sample fastqs.

    #Run fastqc snakefile on the flux cluster (To run on cluster, make sure shell.prefix is uncommented in Snakefile):
    snakemake -j 20 --snakefile Snakefile_fastqc --configfile example/ATAC_2369_config.json --latency-wait 60 --cluster-config flux_config_ATACseq.json --cluster "qsub -N {cluster.name} -A {cluster.account} -q {cluster.queue} -l nodes={cluster.nodes}:ppn={cluster.ntask} -l mem={cluster.memory} -l walltime={cluster.time} {cluster.env}"

##### Inspect the results of fastqc:

* If read-trimming is needed, perform as necessary and generate a new config for the trimmed input fastq files
* If not, continue on with the ATACseq snakefile

##### Continuing with ATACseq pipeline

This pipeline has built-in flexibility regarding the aligners used to create the requisite merged, aligned BAM files. Currently BWA aln and BWA mem are implemented and can be included in the ATACseq snakefile. From there, duplicates are marked, bams are pruned, bigwigs are generated, peaks are called using MACS2, called peaks are filtered against blacklist regions, and ATAQV is run to gather and summarize quality metrics.

    #Running ATACseq snakefile on the flux cluster
    snakemake -j 20 --snakefile Snakefile_ATACseq --configfile example/ATAC_2369_config.json --latency-wait 60 --cluster-config flux_config_ATACseq.json --cluster "qsub -N {cluster.name} -A {cluster.account} -q {cluster.queue} -l nodes={cluster.nodes}:ppn={cluster.ntask} -l mem={cluster.memory} -l walltime={cluster.time} {cluster.env}"


