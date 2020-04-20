base_dir=/nfs/med-bfx-activeprojects/Ryan_rjhryan_CU1
test_dir=${base_dir}/atac_test_run
repo_dir=${base_dir}/ChIPseq_ATACseq_pipelines

mkdir -p ${test_dir}/tmp

cd ${test_dir}

${repo_dir}/ngs_rawdata_config_creator.py \
    --general_input ${repo_dir}/example/ATAC_general.json \
    --per_lib_input ${repo_dir}/tests/atac_test_samplesheet.csv \
    --results_dir ${test_dir} \
    --temp_dir ${test_dir}/tmp \
    > ${test_dir}/config.json

# On comps
snakemake --snakefile ${repo_dir}/Snakefile_ATACseq --configfile config.json --use-conda

# On GL
snakemake -p --snakefile ${repo_dir}/Snakefile_ATACseq --configfile ${test_dir}/config.json \
    --latency-wait 60 --jobs 144 --cluster-config ${repo_dir}/cluster_config.json --use-conda \
    --cluster 'sbatch --job-name={cluster.name} --account={cluster.account} --partition={cluster.partition} --nodes={cluster.nodes} --ntasks-per-node={cluster.ntask} --mem={cluster.memory} --time={cluster.time} --output=%x-%j.out'
