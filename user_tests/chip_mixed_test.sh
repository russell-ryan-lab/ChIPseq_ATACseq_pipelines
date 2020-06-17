conda_env_dir=/nfs/turbo/<lab-turbo>/envs/atac_chip_pipeline

repo_dir=/nfs/med-bfx-activeprojects/Ryan_rjhryan_CU1/ChIPseq_ATACseq_pipelines
result_dir=${repo_dir}/chip_test_mixed2

mkdir -p ${result_dir}
cd ${result_dir}

conda activate ${conda_env_dir}

# Create the configs for the input and IP, run paired-end and single-end, respectively
${repo_dir}/scripts/config_creator.py \
    --general_input ${repo_dir}/config/ChIP_histone_general_se.yaml \
    --per_lib_input ${repo_dir}/data/sra_mixed_chip_test_data/sra_chip_histone_input.csv \
    --results_dir ${result_dir} \
    --temp_dir ${result_dir}/tmp \
    > ${result_dir}/config_histone_input.yaml

${repo_dir}/scripts/config_creator.py \
    --general_input ${repo_dir}/config/ChIP_histone_general_pe.yaml \
    --per_lib_input ${repo_dir}/data/sra_mixed_chip_test_data/sra_chip_histone_IP.csv \
    --results_dir ${result_dir} \
    --temp_dir ${result_dir}/tmp \
    > ${result_dir}/config_histone_IP.yaml

# Run input through SE pipeline with SE params, and IP through PE pipeline with PE params. These can run simultaneously.
# On GL (screen 0)
snakemake -p --snakefile ${repo_dir}/ChIPseq_se.smk --configfile ${result_dir}/config_histone_input.yaml \
    --latency-wait 60 --jobs 144 --cluster-status ${repo_dir}/scripts/slurm_status.py --cluster-config ${repo_dir}/config/cluster_config.yaml --use-conda \
    --cluster 'sbatch --job-name={cluster.name} --account={cluster.account} --partition={cluster.partition} --nodes={cluster.nodes} --ntasks-per-node={cluster.ntask} --mem={cluster.memory} --time={cluster.time} --output=logs/%x-%j.out --parsable'

# On GL (screen 1) (have to re-set the variables)
snakemake -p --snakefile ${repo_dir}/ChIPseq_pe.smk --configfile ${result_dir}/config_histone_IP.yaml \
    --latency-wait 60 --jobs 144 --cluster-status ${repo_dir}/scripts/slurm_status.py --cluster-config ${repo_dir}/config/cluster_config.yaml --use-conda \
    --cluster 'sbatch --job-name={cluster.name} --account={cluster.account} --partition={cluster.partition} --nodes={cluster.nodes} --ntasks-per-node={cluster.ntask} --mem={cluster.memory} --time={cluster.time} --output=logs/%x-%j.out --parsable'


# Create the config that will be used for homer_only workflow
${repo_dir}/scripts/config_creator.py \
    --general_input ${repo_dir}/config/ChIP_histone_general_se.yaml \
    --per_lib_input ${repo_dir}/data/sra_mixed_chip_test_data/sra_chip_histone_all.csv \
    --results_dir ${result_dir} \
    --temp_dir ${result_dir}/tmp \
    --homer_only \
    > ${result_dir}/config_histone_homer.yaml

snakemake -p --snakefile ${repo_dir}/homer_only.smk --configfile ${result_dir}/config_histone_homer.yaml \
    --latency-wait 60 --jobs 144 --cluster-status ${repo_dir}/scripts/slurm_status.py --cluster-config ${repo_dir}/config/cluster_config.yaml --use-conda \
    --cluster 'sbatch --job-name={cluster.name} --account={cluster.account} --partition={cluster.partition} --nodes={cluster.nodes} --ntasks-per-node={cluster.ntask} --mem={cluster.memory} --time={cluster.time} --output=logs/%x-%j.out --parsable'
