base_dir=/nfs/med-bfx-activeprojects/Ryan_rjhryan_CU1
test_dir=${base_dir}/atac_test_run2
repo_dir=${base_dir}/ChIPseq_ATACseq_pipelines

mkdir -p ${test_dir}/tmp
mkdir -p ${test_dir}/logs

cd ${test_dir}

conda activate snakemake

${repo_dir}/scripts/config_creator.py \
    --general_input ${repo_dir}/config/ATAC_general.yaml \
    --per_lib_input ${repo_dir}/data/atac_test_data/atac_test_samplesheet.csv \
    --results_dir ${test_dir} \
    --temp_dir ${test_dir}/tmp \
    > ${test_dir}/config.yaml

# On GL
snakemake -p --snakefile ${repo_dir}/ATACseq.smk --configfile ${test_dir}/config.yaml \
    --latency-wait 60 --jobs 144 --cluster-config ${repo_dir}/config/cluster_config.yaml --use-conda \
    --cluster 'sbatch --job-name={cluster.name} --account={cluster.account} --partition={cluster.partition} --nodes={cluster.nodes} --ntasks-per-node={cluster.ntask} --mem={cluster.memory} --time={cluster.time} --output=logs/%x-%j.out'
