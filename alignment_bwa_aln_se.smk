import os

#Note: ALIGN_DIR and CONCAT_READS_DIR variables must be set in parent snakefile


# Rules

rule bwa_aln:
    input:
        os.path.join(CONCAT_READS_DIR, "{library}_R{read}.fastq.gz")
    output:
        temp(os.path.join(ALIGN_DIR, "{library}_R{read}.sai"))
    params:
        index = lambda wildcards: config['bwa_index'][config['lib_genome'][wildcards.library]] #nested config call identical functionality to config['bwa_index'][get_genome(wildcards.library)]
    threads: 8
    shell:
        "bwa aln -t {threads} {params.index} {input} > {output}"

rule bwa_samse:
    input:
        fq_left = lambda wildcards: config['lib_paths'][wildcards.library]['1'],
        sai_left = os.path.join(ALIGN_DIR, "{library}_R1.sai"),
    output:
        temp(os.path.join(ALIGN_DIR, "{library}.aligned.sam"))
    params:
        index = lambda wildcards: config['bwa_index'][config['lib_genome'][wildcards.library]]
    shell:
        "bwa samse {params.index} {input.sai_left} {input.fq_left} > {output}"

rule picard_sort:
    input:
        os.path.join(ALIGN_DIR, "{library}.aligned.sam")
    output:
        temp(os.path.join(ALIGN_DIR, "{library}.sorted.bam"))
    shell:
        "picard SortSam I={input} O={output} SO=coordinate VALIDATION_STRINGENCY=LENIENT"
