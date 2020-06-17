conda_env_dir=/nfs/turbo/<lab-turbo>/envs/atac_chip_pipeline

repo_dir=/nfs/med-bfx-activeprojects/Ryan_rjhryan_CU1/ChIPseq_ATACseq_pipelines
results_dir=${repo_dir}/test_sra_chip_tf_pe

mkdir -p ${results_dir}
cd ${results_dir}

conda activate ${conda_env_dir}

${repo_dir}/scripts/config_creator.py \
    --general_input ${repo_dir}/config/ChIP_TF_general_pe.yaml \
    --per_lib_input ${repo_dir}/data/sra_chip_test_data/sra_chip_samplesheet.csv \
    --results_dir ${results_dir} \
    --temp_dir ${results_dir}/tmp \
    > ${results_dir}/config_tf_pe.yaml

snakemake -p --snakefile ${repo_dir}/ChIPseq_pe.smk --configfile config_tf_pe.yaml \
    --latency-wait 60 --jobs 144 --use-conda \
    --cluster-config ${repo_dir}/config/cluster_config.yaml --cluster-status ${repo_dir}/scripts/slurm_status.py \
    --cluster '\
        sbatch --job-name={cluster.name} --account={cluster.account} --partition={cluster.partition} --nodes={cluster.nodes} \
        --ntasks-per-node={cluster.ntask} --mem={cluster.memory} --time={cluster.time} --parsable --output=logs/%x-%j.out'
