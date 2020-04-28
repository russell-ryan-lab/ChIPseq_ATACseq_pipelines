module purge
module load singularity/3.5.2

singularity shell -B /nfs/med-bfx-activeprojects /nfs/turbo/epicore-active/common/singularity/seq_utilities.simg

################################################################################


ORIG_DIR=/nfs/med-bfx-activeprojects/Ryan_rjhryan_CU1/atac_test_run
CORE_DIR=/nfs/med-bfx-activeprojects/Ryan_rjhryan_CU1/ChIPseq_ATACseq_pipelines/data/atac_test_data

mkdir -p ${CORE_DIR}
cd ${CORE_DIR}

########################################

# Collect the read names from these bams.
# NOTE: 1. Somewhere a /1 or /2 is .
# NOTE: 2. We need to sort and take uniq of the read names because seqtk subseq
# will find read names in fastqs with multiplicity.
find ${ORIG_DIR}/aligned -name "*merged.bam" | sort | xargs -I{} basename {} ".merged.bam" | \
    xargs -n 1 -P 18 -I{} sh -c '\
        bedtools intersect \
            -a /nfs/med-bfx-activeprojects/Ryan_rjhryan_CU1/ChIPseq_ATACseq_pipelines/data/atac_test_data_regions.bed \
            -b "${2}/aligned/${1}.merged.bam" \
            -sorted \
            -g /nfs/turbo/epicore-active/genomes/Mus_musculus/UCSC/mm10/mm10.chrom.sizes.fa.sort \
            -wb \
        | cut -f 7 | cut -d '/' -f 1 | sort -T . | uniq > "$1.read_names"' -- {} ${ORIG_DIR}

# Convert to fastq
find ${ORIG_DIR}/aligned -name "*merged.bam" | sort | xargs -I{} basename {} ".merged.bam" | \
    xargs -n 1 -P 18 -I{} sh -c 'seqtk subseq "/nfs/med-bfx-dnaseqcore/Run_testatacseq/epicore/Sample_${1}/${1}_L001_R1.fastq.gz" "${2}/${1}.read_names" > "${2}/${1}_L001_R1.fastq"' -- {} ${CORE_DIR}
find ${ORIG_DIR}/aligned -name "*merged.bam" | sort | xargs -I{} basename {} ".merged.bam" | \
    xargs -n 1 -P 18 -I{} sh -c 'seqtk subseq "/nfs/med-bfx-dnaseqcore/Run_testatacseq/epicore/Sample_${1}/${1}_L001_R2.fastq.gz" "${2}/${1}.read_names" > "${2}/${1}_L001_R2.fastq"' -- {} ${CORE_DIR}

# Zip the fastq
gzip *fastq

# Delete .read_name
rm *read_names

# Create the test directory and move files into a MAGC-like structure
for i in `seq 1 4`
do
    mkdir -p ${CORE_DIR}/Sample_${i}
    mv ${i}_L001_R*.fastq.gz ${CORE_DIR}/Sample_${i}/
done
