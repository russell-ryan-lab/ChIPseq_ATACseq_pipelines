import os

#Note: ALIGN_DIR and CONCAT_READS_DIR variables must be set in parent snakefile


# Rules

rule bwa_mem:
    input:
        expand(os.path.join(CONCAT_READS_DIR, "{{sample}}_R{read}.fastq.gz"), read=[1,2], allow_missing=True)
    output:
        temp(os.path.join(ALIGN_DIR, "{sample}.sorted.bam"))
    params:
        index = lambda wildcards: config['bwa_index'][config['sample_genome'][wildcards.sample]],  #nested config call has identical functionality to config['bwa_index'][get_genome(wildcards.sample)]
        sort_tmp = os.path.join(config['tmpdir'], '{sample}.sort.tmp')
    threads: 8
    conda: "../envs/bwa.yaml"
    shell:
        "HALF_THREADS=$(( {threads} / 2 )); "
        "bwa mem -M -t $HALF_THREADS {params.index} {input} | samtools sort -@ $HALF_THREADS -O bam -T {params.sort_tmp} -o {output} -"
