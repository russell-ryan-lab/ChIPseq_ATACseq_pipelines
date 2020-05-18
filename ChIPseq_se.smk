import os
import sys
import functools


# GENERIC DATA

INCLUDE_CHRS = {
    'hg19': ['chr{}'.format(i) for i in range(1, 23)] + ["chrX", "chrY"],
    'hg38': ['chr{}'.format(i) for i in range(1, 23)] + ["chrX", "chrY"],
    'mm9': ['chr{}'.format(i) for i in range(1, 20)] + ["chrX", "chrY"],
    'mm10': ['chr{}'.format(i) for i in range(1, 20)] + ["chrX", "chrY"],
    'rn4': ['chr{}'.format(i) for i in range(1, 21)] + ["chrX", "chrY"],
    'rn5': ['chr{}'.format(i) for i in range(1, 21)] + ["chrX", "chrY"]
}


# Helper functions

def get_genome(sample):
    return(config['sample_genome'][sample])


# RESULT PATHS

prefix_results = functools.partial(os.path.join, config['results_dir'])
CONCAT_READS_DIR = prefix_results('concat_reads')
ALIGN_DIR = prefix_results('aligned')
PRUNE_DIR = prefix_results('pruned')
DISP_DIR = prefix_results('display_tracks')
HOMERTAG_DIR = prefix_results('tag')
HOMERPEAK_DIR = prefix_results('peaks')
HOMERMOTIF_DIR = prefix_results('homer_motifs')


SCRIPTS_DIR = os.path.join(workflow.basedir, 'scripts')


# Set workdir - Snakemake will be run from this location.
workdir:
    config['results_dir']


# Rules

rule all:
    input:
        expand(os.path.join(DISP_DIR, "{sample}.tdf"), sample=config['sample_paths'].keys()), #Create tdfs for all samples
        expand(os.path.join(DISP_DIR, "{sample}.1m.bw"), sample=config['sample_paths'].keys()), #Create bigwigs for all samples
        expand(os.path.join(HOMERPEAK_DIR, "{paramset}", "{sample}_BLfiltered.hpeaks"), paramset=config['homer_findPeaks_params'].keys(), sample=config['sample_input'].keys()), #Call peaks for all samples with matched inputs
        expand(os.path.join(HOMERMOTIF_DIR, "{sample}"), sample=config['sample_homer_fmg_genome'].keys()), #Homermotifs for all samples with a specified genome for homer findMotifsGenome


include:
    "rules/alignment_bwa_aln_se.smk"

rule concatenate_reads:
    input:
        lambda wildcards: config['sample_paths'][wildcards.sample][wildcards.read]
    output:
        os.path.join(CONCAT_READS_DIR, "{sample}_R{read}.fastq.gz")
    shell:
        "cat {input} > {output}"

rule mark_duplicates:
    input:
        os.path.join(ALIGN_DIR, "{sample}.sorted.bam")
    output:
        bam = os.path.join(ALIGN_DIR, "{sample}.mrkdup.bam"),
        metric = os.path.join(ALIGN_DIR, "{sample}.mrkdup.metric")
    params:
        tmpdir = config['tmpdir']
    conda: "envs/picard.yaml"
    shell:
        "export JAVA_OPTIONS=-Xmx12g ; "
        "picard MarkDuplicates I={input} O={output.bam} "
        "METRICS_FILE={output.metric} "
        "ASSUME_SORTED=True "
        "VALIDATION_STRINGENCY=LENIENT "
        "TMP_DIR={params.tmpdir}"

rule index_dupmarked_bams:
    input:
        os.path.join(ALIGN_DIR, "{sample}.mrkdup.bam")
    output:
        os.path.join(ALIGN_DIR, "{sample}.mrkdup.bai")
    conda: "envs/samtools.yaml"
    shell:
        "samtools index {input} {output}"

rule samtools_prune:
    input:
        bam = os.path.join(ALIGN_DIR, "{sample}.mrkdup.bam"),
        bai = os.path.join(ALIGN_DIR, "{sample}.mrkdup.bai")
    output:
        bam = temp(os.path.join(PRUNE_DIR, "{sample}.pruned.bam"))
    params:
        incl_chr = lambda wildcards: INCLUDE_CHRS[get_genome(wildcards.sample)],
        flags = config['samtools_prune_flags']
    conda: "envs/samtools.yaml"
    shell:
        "samtools view -b {params.flags} {input.bam} {params.incl_chr} > {output.bam}"

rule index_pruned:
    input:
        bam = os.path.join(PRUNE_DIR, "{sample}.pruned.bam")
    output:
        bai = os.path.join(PRUNE_DIR, "{sample}.pruned.bai")
    conda: "envs/samtools.yaml"
    shell:
        "samtools index {input.bam} {output.bai}"

rule igvtools_count_tdf:
    input:
        os.path.join(PRUNE_DIR, "{sample}.pruned.bam")
    output:
        os.path.join(DISP_DIR, "{sample}.tdf")
    params:
        genome = lambda wildcards: get_genome(wildcards.sample),
        args = config['igvtools_count_params']
    conda: "envs/igvtools.yaml"
    shell:
        "igvtools count {params.args} {input} {output} {params.genome}"

rule deeptools_bamcoverage_bw:
    input:
        bam = os.path.join(PRUNE_DIR, "{sample}.pruned.bam"),
        bai = os.path.join(PRUNE_DIR, "{sample}.pruned.bai")
    output:
        os.path.join(DISP_DIR, "{sample}.1m.bw")
    params:
        blacklist = lambda wildcards: config['blacklist'][get_genome(wildcards.sample)],
        args = config['deeptools_bamcoverage_params']
    conda: "envs/deeptools.yaml"
    shell:
        "bamCoverage --bam {input.bam} -o {output} -bl {params.blacklist} {params.args}"

rule makeTagDirectory:
    input:
        os.path.join(PRUNE_DIR, "{sample}.pruned.bam")
    output:
        directory(os.path.join(HOMERTAG_DIR, "{sample}"))
    params:
        genome = lambda wildcards: get_genome(wildcards.sample),
        params = config['makeTagDir_params']
    shell:
        "makeTagDirectory {output} {params.params} -genome {params.genome} {input}"

rule findPeaks:
    input:
        sample = os.path.join(HOMERTAG_DIR, "{sample}"),
        input = lambda wildcards: os.path.join(HOMERTAG_DIR, config['sample_input'][wildcards.sample])
    output:
        os.path.join(HOMERPEAK_DIR, "{paramset}", "{sample}.all.hpeaks")
    params:
        config['homer_findPeaks_params']['{paramset}']
    shell:
        "findPeaks {input.sample} -i {input.input} {params} -o {output}"

rule pos2bed:
    input:
        os.path.join(HOMERPEAK_DIR, "{paramset}", "{sample}.all.hpeaks")
    output:
        os.path.join(HOMERPEAK_DIR, "{paramset}", "{sample}.all.bed")
    shell:
        "pos2bed.pl {input} > {output}"

rule blacklist_filter_bed:
    input:
        os.path.join(HOMERPEAK_DIR, "{paramset}", "{sample}.all.bed"),
    output:
        os.path.join(HOMERPEAK_DIR, "{paramset}", "{sample}_BLfiltered.bed"),
    params:
        blacklist = lambda wildcards: config['blacklist'][get_genome(wildcards.sample)]
    conda: "envs/bedtools.yaml"
    shell:
        "bedtools intersect -a {input} -b {params.blacklist} -v > {output}"

rule keepBedEntriesInHpeaks:
    input:
        filtbed = os.path.join(HOMERPEAK_DIR, "{paramset}", "{sample}_BLfiltered.bed"),
        allhpeaks = os.path.join(HOMERPEAK_DIR, "{paramset}", "{sample}.all.hpeaks")
    output:
        os.path.join(HOMERPEAK_DIR, "{paramset}", "{sample}_BLfiltered.hpeaks")
    conda: "envs/pysam.yaml"
    shell:
        "python {SCRIPTS_DIR}/keepBedEntriesInHpeaks.py -i {input.allhpeaks} -b {input.filtbed} -o {output}"

rule findMotifsGenome:
    input:
        os.path.join(HOMERPEAK_DIR, "{paramset}", "{sample}_BLfiltered.hpeaks")
    output:
        directory(os.path.join(HOMERMOTIF_DIR, "{paramset}", "{sample}"))
    params:
        genome = lambda wildcards: config['sample_homer_fmg_genome'][wildcards.sample],
        params = config['homer_fmg_params'],
        tmpdir = config['tmpdir']
    shell:
        "findMotifsGenome.pl {input} {params.genome} {output} {params.params} -preparsedDir {params.tmpdir}"
