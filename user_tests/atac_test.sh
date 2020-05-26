conda_env_dir=/nfs/turbo/<lab-turbo>/envs/atac_chip_pipeline

repo_dir=/nfs/med-bfx-activeprojects/Ryan_rjhryan_CU1/ChIPseq_ATACseq_pipelines
results_dir=${repo_dir}/atac_test_run

mkdir -p ${results_dir}
cd ${results_dir}

conda activate ${conda_env_dir}

${repo_dir}/scripts/config_creator.py \
    --general_input ${repo_dir}/config/ATAC_general.yaml \
    --per_lib_input ${repo_dir}/data/atac_test_data/atac_test_samplesheet.csv \
    --results_dir ${results_dir} \
    --temp_dir ${results_dir}/tmp \
    > ${results_dir}/config_atac.yaml

snakemake -p --snakefile ${repo_dir}/ATACseq.smk --configfile ${results_dir}/config_atac.yaml \
    --latency-wait 60 --jobs 144 --use-conda \
    --cluster-config ${repo_dir}/config/cluster_config.yaml --cluster-status ${repo_dir}/scripts/slurm_status.py \
    --cluster '\
        sbatch --job-name={cluster.name} --account={cluster.account} --partition={cluster.partition} --nodes={cluster.nodes} \
        --ntasks-per-node={cluster.ntask} --mem={cluster.memory} --time={cluster.time} --parsable --output=logs/%x-%j.out'
