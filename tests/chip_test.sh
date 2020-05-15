base_dir=/nfs/med-bfx-activeprojects/Ryan_rjhryan_CU1
test_dir=${base_dir}/chip_test_run
repo_dir=${base_dir}/ChIPseq_ATACseq_pipelines

mkdir -p ${test_dir}/tmp
mkdir -p ${test_dir}/logs

cd ${test_dir}

conda activate atac_chip_pipeline

${repo_dir}/scripts/config_creator.py \
    --general_input ${repo_dir}/config/ChIP_histone_general_se.yaml \
    --per_lib_input ${repo_dir}/data/sra_chip_test_data/sra_chip_histone.csv \
    --results_dir ${test_dir} \
    --temp_dir ${test_dir}/tmp \
    > ${test_dir}/config_histone.yaml

${repo_dir}/scripts/config_creator.py \
    --general_input ${repo_dir}/config/ChIP_TF_general_se.yaml \
    --per_lib_input ${repo_dir}/data/sra_chip_test_data/sra_chip_TF.csv \
    --results_dir ${test_dir} \
    --temp_dir ${test_dir}/tmp \
    > ${test_dir}/config_TF.yaml

# On GL
snakemake -p --snakefile ${repo_dir}/ChIPseq_se.smk --configfile ${test_dir}/config_histone.yaml \
    --latency-wait 60 --jobs 144 --cluster-config ${repo_dir}/config/cluster_config.yaml --use-conda \
    --cluster 'sbatch --job-name={cluster.name} --account={cluster.account} --partition={cluster.partition} --nodes={cluster.nodes} --ntasks-per-node={cluster.ntask} --mem={cluster.memory} --time={cluster.time} --output=logs/%x-%j.out'

snakemake -p --snakefile ${repo_dir}/ChIPseq_se.smk --configfile ${test_dir}/config_TF.yaml \
    --latency-wait 60 --jobs 144 --cluster-config ${repo_dir}/config/cluster_config.yaml --use-conda \
    --cluster 'sbatch --job-name={cluster.name} --account={cluster.account} --partition={cluster.partition} --nodes={cluster.nodes} --ntasks-per-node={cluster.ntask} --mem={cluster.memory} --time={cluster.time} --output=logs/%x-%j.out'
