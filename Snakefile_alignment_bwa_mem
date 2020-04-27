import os

#Note: ALIGN_DIR and CONCAT_READS_DIR variables must be set in parent snakefile


# Rules

rule bwa_mem:
    input:
        os.path.join(CONCAT_READS_DIR, "{library}_R{read}.fastq.gz")
    output:
        temp(os.path.join(ALIGN_DIR, "{library}.sorted.bam"))
    params:
        index = lambda wildcards: config['bwa_index'][config['lib_genome'][wildcards.library]],  #nested config call has identical functionality to config['bwa_index'][get_genome(wildcards.library)]
        sort_tmp = os.path.join(config['tmpdir'], '{library}.sort.tmp')
    threads: 8
    conda: "envs/bwa.yaml"
    shell:
        "bwa mem -M -t {threads} {params.index} {input} | samtools sort -m 1g -@ {threads} -O bam -T {params.sort_tmp} -o {output} -"
