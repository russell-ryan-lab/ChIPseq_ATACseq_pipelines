# Ran chip_mixed_test.sh with full datasets downloaded from SRA

base_dir=/nfs/med-bfx-activeprojects/Ryan_rjhryan_CU1
test_dir=${base_dir}/chip_test_mixed
repo_dir=${base_dir}/ChIPseq_ATACseq_pipelines

mkdir -p ${test_dir}/tmp
mkdir -p ${test_dir}/logs

cd ${test_dir}

conda activate atac_chip_pipeline

# Create the configs for the input and IP, run paired-end and single-end, respectively
${repo_dir}/scripts/config_creator.py \
    --general_input ${repo_dir}/config/ChIP_histone_general_se.yaml \
    --per_lib_input ${base_dir}/sra_chip_test_data2/sra_chip_histone_input.csv \
    --results_dir ${test_dir} \
    --temp_dir ${test_dir}/tmp \
    > ${test_dir}/config_histone_input.yaml

${repo_dir}/scripts/config_creator.py \
    --general_input ${repo_dir}/config/ChIP_histone_general_pe.yaml \
    --per_lib_input ${base_dir}/sra_chip_test_data2/sra_chip_histone_IP.csv \
    --results_dir ${test_dir} \
    --temp_dir ${test_dir}/tmp \
    > ${test_dir}/config_histone_IP.yaml

# Run input through SE pipeline with SE params, and IP through PE pipeline with PE params. These can run simultaneously.
# On GL (screen 0)
snakemake -p --snakefile ${repo_dir}/ChIPseq_se.smk --configfile ${test_dir}/config_histone_input.yaml \
    --latency-wait 60 --jobs 144 --cluster-status ${repo_dir}/scripts/slurm_status.py --cluster-config ${repo_dir}/config/cluster_config.yaml --use-conda \
    --cluster 'sbatch --job-name={cluster.name} --account={cluster.account} --partition={cluster.partition} --nodes={cluster.nodes} --ntasks-per-node={cluster.ntask} --mem={cluster.memory} --time={cluster.time} --output=logs/%x-%j.out --parsable'

# On GL (screen 1)
snakemake -p --snakefile ${repo_dir}/ChIPseq_pe.smk --configfile ${test_dir}/config_histone_IP.yaml \
    --latency-wait 60 --jobs 144 --cluster-status ${repo_dir}/scripts/slurm_status.py --cluster-config ${repo_dir}/config/cluster_config.yaml --use-conda \
    --cluster 'sbatch --job-name={cluster.name} --account={cluster.account} --partition={cluster.partition} --nodes={cluster.nodes} --ntasks-per-node={cluster.ntask} --mem={cluster.memory} --time={cluster.time} --output=logs/%x-%j.out --parsable'


# Create the config that will be used for homer_only workflow
${repo_dir}/scripts/config_creator.py \
    --general_input ${repo_dir}/config/ChIP_histone_general_se.yaml \
    --per_lib_input ${base_dir}/sra_chip_test_data2/sra_chip_histone.csv \
    --results_dir ${test_dir} \
    --temp_dir ${test_dir}/tmp \
    --homer_only \
    > ${test_dir}/config_histone_homer.yaml

snakemake -p --snakefile ${repo_dir}/homer_only.smk --configfile ${test_dir}/config_histone_homer.yaml \
    --latency-wait 60 --jobs 144 --cluster-status ${repo_dir}/scripts/slurm_status.py --cluster-config ${repo_dir}/config/cluster_config.yaml --use-conda \
    --cluster 'sbatch --job-name={cluster.name} --account={cluster.account} --partition={cluster.partition} --nodes={cluster.nodes} --ntasks-per-node={cluster.ntask} --mem={cluster.memory} --time={cluster.time} --output=logs/%x-%j.out --parsable'

# Browsed results and chose regions on 3 separate chromosomes with H3K27ac peaks.

# To subset:

### Load required software
module purge
module load singularity/3.5.2
singularity shell -B /nfs/med-bfx-activeprojects /nfs/turbo/epicore-active/common/singularity/seq_utilities.simg

### Establish some paths

project_dir=/nfs/med-bfx-activeprojects/Ryan_rjhryan_CU1
repo_dir=${project_dir}/ChIPseq_ATACseq_pipelines
data_dir=${repo_dir}/data/sra_mixed_chip_test_data
orig_dir=${project_dir}/chip_test_mixed

cd $orig_dir

### Find reads
# NOTES
# column 7 of the result of bedtools intersect has the readnames
# a /1 or /2 is introduced into the readname and has to be cut out
find ${orig_dir}/aligned -name "*mrkdup.bam" | sort | xargs -I{} basename {} ".mrkdup.bam" | \
    xargs -n 1 -P 2 -I{} sh -c '\
        bedtools intersect \
            -a ${3}/sra_mixed_chip_test_data_regions.bed \
            -b "${2}/aligned/${1}.mrkdup.bam" \
            -sorted \
            -g /nfs/turbo/epicore-active/genomes/Homo_sapiens/UCSC/hg19/hg19.chrom.sizes.fa.sort \
            -wb \
        | cut -f 7 | cut -d '/' -f 1 | sort -T . | uniq > "$1.read_names"' -- {} ${orig_dir} ${data_dir}

# Create new fastqs
find ${orig_dir} -name "*read_names" | sort | xargs -I{} basename {} ".read_names" | \
    xargs -n 1 -P 2 -I{} sh -c 'seqtk subseq "${2}/concat_reads/${1}_R1.fastq.gz" "${2}/${1}.read_names" > "${2}/${1}_R1.fastq"' -- {} ${orig_dir}

find ${orig_dir} -name "*read_names" | sort | xargs -I{} basename {} ".read_names" | \
    xargs -n 1 -P 2 -I{} sh -c 'seqtk subseq "${2}/concat_reads/${1}_R2.fastq.gz" "${2}/${1}.read_names" > "${2}/${1}_R2.fastq"' -- {} ${orig_dir}

# Zip them
gzip *fastq

# Delete .read_name
rm *read_names

# Copy to data_dir
cp *fastq.gz ${data_dir}/

###
### Cleanup data_dir of full sra data and replace with subsetted data

cd ${data_dir}

mkdir -p ${project_dir}/full_sra_mixed_chip_test_data
mv SRR*fastq.gz ${project_dir}/full_sra_mixed_chip_test_data
mv SAMN* ${project_dir}/full_sra_mixed_chip_test_data

# Exit singularity
exit

###
### Create tree structure

conda activate atac_chip_pipeline
cd ${repo_dir}

python ${repo_dir}/scripts/prepare_fastq_inputs.py --fastq_dir data/sra_chip_test_data --metadata_csv data/sra_chip_test_data/sra_flat_to_tree.txt --group_col BioSample --indiv_col Run

###
### Remove flat files

cd ${data_dir}
rm *fastq.gz
