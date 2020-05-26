###
### Establish some paths to begin

project_dir=/nfs/med-bfx-activeprojects/Ryan_rjhryan_CU1
repo_dir=${project_dir}/ChIPseq_ATACseq_pipelines
data_dir=${repo_dir}/data/sra_chip_test_data

###
### Download data from SRA

# Download SraRunTable.txt and SRR_Acc_list.txt from SRA Run Selector

module load Bioinformatics
module load sratoolkit

mkdir -p ${data_dir}
cd ${data_dir}

cat SRR_Acc_List.txt | xargs -L 1 -P 2 -I {} sh -c "fastq-dump --outdir ./ --split-files {}"

###
### Run prepare_fastq_inputs.py to form tree-structure

conda activate atac_chip_pipeline

python ${repo_dir}/scripts/prepare_fastq_inputs.py --fastq_dir . --metadata_csv SraRunTable.txt --group_col BioSample --indiv_col Run --add_R

###
### Run full-data through entire pipeline in TF/histone and PE/SE modes (4 total)
### Goal is to make sure we pick regions so any run configuration will work in testing

# Create directories
results_histone_se_dir=${project_dir}/full_sra_histone_se
results_histone_pe_dir=${project_dir}/full_sra_histone_pe
results_tf_se_dir=${project_dir}/full_sra_tf_se
results_tf_pe_dir=${project_dir}/full_sra_tf_pe

for dir in `echo ${results_histone_se_dir} ${results_histone_pe_dir} ${results_tf_se_dir} ${results_tf_pe_dir}`
do
    mkdir -p ${dir}
done


# Create config files in appropriate directories
${repo_dir}/scripts/config_creator.py \
    --general_input ${repo_dir}/config/ChIP_histone_general_se.yaml \
    --per_lib_input ${repo_dir}/data/sra_chip_test_data/full_sra_chip_samplesheet.csv \
    --results_dir ${results_histone_se_dir} \
    --temp_dir ${results_histone_se_dir}/tmp \
    > ${results_histone_se_dir}/config_histone_se.yaml

${repo_dir}/scripts/config_creator.py \
    --general_input ${repo_dir}/config/ChIP_histone_general_pe.yaml \
    --per_lib_input ${repo_dir}/data/sra_chip_test_data/full_sra_chip_samplesheet.csv \
    --results_dir ${results_histone_pe_dir} \
    --temp_dir ${results_histone_pe_dir}/tmp \
    > ${results_histone_pe_dir}/config_histone_pe.yaml

${repo_dir}/scripts/config_creator.py \
    --general_input ${repo_dir}/config/ChIP_TF_general_se.yaml \
    --per_lib_input ${repo_dir}/data/sra_chip_test_data/full_sra_chip_samplesheet.csv \
    --results_dir ${results_tf_se_dir} \
    --temp_dir ${results_tf_se_dir}/tmp \
    > ${results_tf_se_dir}/config_tf_se.yaml

${repo_dir}/scripts/config_creator.py \
    --general_input ${repo_dir}/config/ChIP_TF_general_pe.yaml \
    --per_lib_input ${repo_dir}/data/sra_chip_test_data/full_sra_chip_samplesheet.csv \
    --results_dir ${results_tf_pe_dir} \
    --temp_dir ${results_tf_pe_dir}/tmp \
    > ${results_tf_pe_dir}/config_tf_pe.yaml


# Full runs on GL, this may take a while, best to create a screen for each and run in parallel
project_dir=/nfs/med-bfx-activeprojects/Ryan_rjhryan_CU1
repo_dir=${project_dir}/ChIPseq_ATACseq_pipelines
data_dir=${repo_dir}/data/sra_chip_test_data
results_histone_se_dir=${project_dir}/full_sra_histone_se
results_histone_pe_dir=${project_dir}/full_sra_histone_pe
results_tf_se_dir=${project_dir}/full_sra_tf_se
results_tf_pe_dir=${project_dir}/full_sra_tf_pe

cd ${results_histone_se_dir}
snakemake -p --snakefile ${repo_dir}/ChIPseq_se.smk --configfile config_histone_se.yaml \
    --latency-wait 60 --jobs 144 --use-conda --notemp \
    --cluster-config ${repo_dir}/config/cluster_config.yaml --cluster-status ${repo_dir}/scripts/slurm_status.py \
    --cluster '\
        sbatch --job-name={cluster.name} --account={cluster.account} --partition={cluster.partition} --nodes={cluster.nodes} \
        --ntasks-per-node={cluster.ntask} --mem={cluster.memory} --time={cluster.time} --parsable --output=logs/%x-%j.out'

cd ${results_histone_pe_dir}
snakemake -p --snakefile ${repo_dir}/ChIPseq_pe.smk --configfile config_histone_pe.yaml \
    --latency-wait 60 --jobs 144 --use-conda --notemp --conda-prefix ${results_histone_se_dir}/.snakemake/conda \
    --cluster-config ${repo_dir}/config/cluster_config.yaml --cluster-status ${repo_dir}/scripts/slurm_status.py \
    --cluster '\
        sbatch --job-name={cluster.name} --account={cluster.account} --partition={cluster.partition} --nodes={cluster.nodes} \
        --ntasks-per-node={cluster.ntask} --mem={cluster.memory} --time={cluster.time} --parsable --output=logs/%x-%j.out'

cd ${results_tf_se_dir}
snakemake -p --snakefile ${repo_dir}/ChIPseq_se.smk --configfile config_tf_se.yaml \
    --latency-wait 60 --jobs 144 --use-conda --notemp --conda-prefix ${results_histone_se_dir}/.snakemake/conda \
    --cluster-config ${repo_dir}/config/cluster_config.yaml --cluster-status ${repo_dir}/scripts/slurm_status.py \
    --cluster '\
        sbatch --job-name={cluster.name} --account={cluster.account} --partition={cluster.partition} --nodes={cluster.nodes} \
        --ntasks-per-node={cluster.ntask} --mem={cluster.memory} --time={cluster.time} --parsable --output=logs/%x-%j.out'

cd ${results_tf_pe_dir}
snakemake -p --snakefile ${repo_dir}/ChIPseq_pe.smk --configfile config_tf_pe.yaml \
    --latency-wait 60 --jobs 144 --use-conda --notemp --conda-prefix ${results_histone_se_dir}/.snakemake/conda \
    --cluster-config ${repo_dir}/config/cluster_config.yaml --cluster-status ${repo_dir}/scripts/slurm_status.py \
    --cluster '\
        sbatch --job-name={cluster.name} --account={cluster.account} --partition={cluster.partition} --nodes={cluster.nodes} \
        --ntasks-per-node={cluster.ntask} --mem={cluster.memory} --time={cluster.time} --parsable --output=logs/%x-%j.out'

###
### Subsample data based on data/sra_chip_test_data/sra_chip_test_data_regions.bed
### These regions contain peaks common to all the modes of running the data above
### We only need to subset them from one of the datasets, choose TF PE

module purge
module load singularity/3.5.2
singularity shell -B /nfs/med-bfx-activeprojects /nfs/turbo/epicore-active/common/singularity/seq_utilities.simg

project_dir=/nfs/med-bfx-activeprojects/Ryan_rjhryan_CU1
repo_dir=${project_dir}/ChIPseq_ATACseq_pipelines
data_dir=${repo_dir}/data/sra_chip_test_data
orig_dir=${project_dir}/full_sra_tf_pe

cd ${orig_dir}

# Find reads
# NOTES
# column 7 of the result of bedtools intersect has the readnames
# a /1 or /2 is introduced into the readname and has to be cut out
find ${orig_dir}/aligned -name "*mrkdup.bam" | sort | xargs -I{} basename {} ".mrkdup.bam" | \
    xargs -n 1 -P 2 -I{} sh -c '\
        bedtools intersect \
            -a ${3}/sra_chip_test_data_regions.bed \
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

mkdir -p ${project_dir}/full_sra_chip_test_data
mv SRR*fastq.gz ${project_dir}/full_sra_chip_test_data
mv SAMN* ${project_dir}/full_sra_chip_test_data

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
