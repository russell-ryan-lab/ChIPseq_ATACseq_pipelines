# Example: Running single-end ChIPseq data from SRA through the pipeline

## Download files from SRA, prepare them

    #Load the toolkit and navigate to where we want to store the raw fastqs
    module load sratoolkit
    cd /scratch/rjhryan_fluxod/trsaari/SRA_fastqs/Basso_2018_mef2B/
    #Download them
    for i in $(seq 6730191 1 6730200) ; do echo $i ; fastq-dump --split-files SRR${i} ; done
    #
    #These are single-end data, with all of the samples in this directory. Move these into subdirs like so
    for i in *_1.fastq ; do SAMPLENAME=$(basename $i "_1.fastq") ; echo $SAMPLENAME ; mkdir $SAMPLENAME && mv $i ${SAMPLENAME}/${SAMPLENAME}_1.fastq ; done

## Put together a sample sheet .csv file

See `example/Basso_2018_mef2B_samplesheet.csv` for an example

## Create a config

    ngs_rawdata_config_creator.py -g example/ChIP_TF_se_general.json \
      -p example/Basso_2018_mef2B_samplesheet.csv \
      -r /scratch/rjhryan_fluxod/trsaari/Basso_2018_mef2B \
      --file_glob "*.fastq" \
      --capture_regex ".*_(\d)\.fastq"\
      --simulate_single_lane \
        > example/Basso_2018_mef2B_config.json

      #--file_glob and --capture_regex are set because our files are in non-default form:
      # SRA files form - Sample_1.fastq, default is Sample_R1_L001.fastq.gz
      #--simulate_single_lane - Used if fastqs are not split by sequencing lane, assign SE or PE fastq(s) to single lane

## Run the pipeline

    #Run as you normally would - use gnu screen or tmux to create a persistent session, then call snakemake:
    snakemake -j 100 --snakefile Snakefile_ChIPseq_se --configfile example/Basso_2018_mef2B_config.json \
      --latency-wait 60 --cluster-config flux_config.json \
      --cluster "qsub -N {cluster.name} -A {cluster.account} -q {cluster.queue} -l nodes={cluster.nodes}:ppn={cluster.ntask} -l mem={cluster.memory} -l walltime={cluster.time} {cluster.env}"
