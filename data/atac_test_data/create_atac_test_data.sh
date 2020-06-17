module purge
module load singularity/3.5.2

singularity shell -B /nfs/med-bfx-activeprojects /nfs/turbo/epicore-active/common/singularity/seq_utilities.simg

project_dir=/nfs/med-bfx-activeprojects/Ryan_rjhryan_CU1

full_dir=${project_dir}/atac_test_run
data_dir=${project_dir}/ChIPseq_ATACseq_pipelines/data/atac_test_data
seqcore_dir=/nfs/med-bfx-dnaseqcore/Run_testatacseq/epicore

mkdir -p ${data_dir}
cd ${data_dir}

###
###

# Find reads
# NOTES
# column 7 of the result of bedtools intersect has the readnames
# a /1 or /2 is introduced into the readname and has to be cut out
find ${full_dir}/aligned -name "*merged.bam" | sort | xargs -I{} basename {} ".merged.bam" | \
    xargs -n 1 -P 18 -I{} sh -c '\
        bedtools intersect \
            -a ${3}/atac_test_data_regions.bed \
            -b "${2}/aligned/${1}.merged.bam" \
            -sorted \
            -g /nfs/turbo/epicore-active/genomes/Mus_musculus/UCSC/mm10/mm10.chrom.sizes.fa.sort \
            -wb \
        | cut -f 7 | cut -d '/' -f 1 | sort -T . | uniq > "$1.read_names"' -- {} ${full_dir} ${data_dir}

# Create new fastqs
find ${full_dir}/aligned -name "*merged.bam" | sort | xargs -I{} basename {} ".merged.bam" | \
    xargs -n 1 -P 18 -I{} sh -c 'seqtk subseq "${3}/Sample_${1}/${1}_L001_R1.fastq.gz" "${2}/${1}.read_names" > "${2}/${1}_L001_R1.fastq"' -- {} ${data_dir}
find ${full_dir}/aligned -name "*merged.bam" | sort | xargs -I{} basename {} ".merged.bam" | \
    xargs -n 1 -P 18 -I{} sh -c 'seqtk subseq "${3}/Sample_${1}/${1}_L001_R2.fastq.gz" "${2}/${1}.read_names" > "${2}/${1}_L001_R2.fastq"' -- {} ${data_dir} ${seqcore_dir}

# Zip them
gzip *fastq

# Delete .read_name
rm *read_names

# Create the test directory and move files into a MAGC-like structure
for i in `seq 1 4`
do
    mkdir -p ${data_dir}/Sample_${i}
    mv ${i}_L001_R*.fastq.gz ${data_dir}/Sample_${i}/
done
