import os

#Note: ALIGN_DIR and CONCAT_READS_DIR variables must be set in parent snakefile


# Rules

rule bwa_aln:
    input:
        os.path.join(CONCAT_READS_DIR, "{sample}_R{read}.fastq.gz")
    output:
        temp(os.path.join(ALIGN_DIR, "{sample}_R{read}.sai"))
    params:
        index = lambda wildcards: config['bwa_index'][config['sample_genome'][wildcards.sample]] #nested config call identical functionality to config['bwa_index'][get_genome(wildcards.sample)]
    threads: 8
    conda: "../envs/bwa.yaml"
    shell:
        "bwa aln -t {threads} {params.index} {input} > {output}"

rule bwa_sampe:
    input:
        fq_left = lambda wildcards: config['sample_paths'][wildcards.sample]['1'],
        fq_right = lambda wildcards: config['sample_paths'][wildcards.sample]['2'],
        sai_left = os.path.join(ALIGN_DIR, "{sample}_R1.sai"),
        sai_right = os.path.join(ALIGN_DIR, "{sample}_R2.sai")
    output:
        temp(os.path.join(ALIGN_DIR, "{sample}.aligned.sam"))
    params:
        index = lambda wildcards: config['bwa_index'][config['sample_genome'][wildcards.sample]]
    conda: "../envs/bwa.yaml"
    shell:
        "bwa sampe {params.index} {input.sai_left} {input.sai_right} {input.fq_left} {input.fq_right} > {output}"

rule picard_sort:
    input:
        os.path.join(ALIGN_DIR, "{sample}.aligned.sam")
    output:
        temp(os.path.join(ALIGN_DIR, "{sample}.sorted.bam"))
    conda: "../envs/picard.yaml"
    shell:
        "picard SortSam I={input} O={output} SO=coordinate VALIDATION_STRINGENCY=LENIENT"
