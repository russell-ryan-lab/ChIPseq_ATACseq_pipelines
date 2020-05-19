# Example: Running single-end ChIPseq data from SRA through the pipeline

## Download files from SRA, prepare them

    #Load the toolkit and navigate to where we want to store the raw fastqs
    module load Bioinformatics
    module load sratoolkit
    cd /scratch/rjhryan_root/rjhryan/trsaari/SRA_fastqs/Basso_2018_mef2B/
    #Download them (Just download 10M spots for this example with -N & -X)
    for i in $(seq 6730191 1 6730200) ; do echo $i ; fastq-dump -N 10000 -X 10010000 --split-files SRR${i} ; done

The ngs_rawdata_config_creator script expects sample fastq's to have the format {Sample}_L{lanenum}_R{readnum}.fastq.gz
Although this script has flexibility to allow for other naming formats, for simplicity in this example we'll just restructure them into the expected format.
It also expects each sample's fastqs to be contained within a subdirectory. The following will accomplish both of the above.

    #These are single-end data, with all of the samples in this directory. For this example, we can move these into subdirs & rename them like so
    for i in *_1.fastq ; do SAMPLENAME=$(basename $i "_1.fastq") ; echo $SAMPLENAME ; mkdir $SAMPLENAME && mv $i ${SAMPLENAME}/${SAMPLENAME}_L001_R1.fastq ; done
    #Now gzip them
    find . -type f -name "*.fastq" -exec gzip {} \;

## Put together a sample sheet .csv file

See `example/Basso_2018_mef2B_samplesheet.csv` for an example

## Create a config

    ngs_rawdata_config_creator.py -g example/ChIP_TF_se_general.json \
      -p example/Basso_2018_mef2B_samplesheet.csv \
      -r /scratch/rjhryan_root/rjhryan/trsaari/Basso_2018_mef2B \
      -t /scratch/rjhryan_root/rjhryan/trsaari/tmp \
        > example/Basso_2018_mef2B_config.json

## Try running a dry-run

    snakemake --snakefile Snakefile_alignment_bwa_aln_se --configfile example/Basso_2018_mef2B_config.json -n -p

## Run the pipeline

    #Use gnu screen or tmux to create a persistent session, then call snakemake:
    screen
    snakemake -j 100 --snakefile Snakefile_ChIPseq_se --configfile example/Basso_2018_mef2B_config.json \
      --latency-wait 60 --cluster-config cluster_config.json \
      --cluster 'sbatch --job-name={cluster.name} --account={cluster.account} --partition={cluster.partition} --nodes={cluster.nodes} --ntasks-per-node={cluster.ntask} --mem={cluster.memory} --time={cluster.time}'
