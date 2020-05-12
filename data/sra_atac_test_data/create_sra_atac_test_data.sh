# module load Bioinformatics
# module load sratoolkit
#
# project_dir=/nfs/med-bfx-activeprojects/Ryan_rjhryan_CU1
# data_dir=${project_dir}/ChIPseq_ATACseq_pipelines/data/sra_atac_test_data
#
# mkdir -p ${data_dir}
#
# cd ${data_dir}
#
# cat SRR_Acc_List.txt | xargs -L 1 -P 3 -I {} sh -c "fastq-dump --outdir ./ --split-files {}"

# The purpose of this test data is to mimic an SRA download (with empty data) where there may be
# multiple runs per sample (*401, *402) and to rearrange the files into the correct
# tree strcture based on SraRunTable.txt
# NOTE: This data is not run through the pipeline, only used to test correct tree structure

project_dir=/nfs/med-bfx-activeprojects/Ryan_rjhryan_CU1
data_dir=${project_dir}/ChIPseq_ATACseq_pipelines/data/sra_atac_test_data

mkdir -p ${data_dir}

touch ${data_dir}/SRR5799398_1.fastq.gz
touch ${data_dir}/SRR5799398_2.fastq.gz
touch ${data_dir}/SRR5799401_1.fastq.gz
touch ${data_dir}/SRR5799401_2.fastq.gz
touch ${data_dir}/SRR5799402_1.fastq.gz
touch ${data_dir}/SRR5799402_2.fastq.gz
touch ${data_dir}/SRR5799447_1.fastq.gz
touch ${data_dir}/SRR5799447_2.fastq.gz
