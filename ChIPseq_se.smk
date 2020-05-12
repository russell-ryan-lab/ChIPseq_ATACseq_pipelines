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

def get_genome(library):
    return(config['lib_genome'][library])


# RESULT PATHS

prefix_results = functools.partial(os.path.join, config['results_dir'])
ALIGN_DIR = prefix_results('aligned')
PRUNE_DIR = prefix_results('pruned')
DISP_DIR = prefix_results('display_tracks')
HOMERTAG_DIR = prefix_results('tag')
HOMERPEAK_DIR = prefix_results('peaks')
HOMERMOTIF_DIR = prefix_results('homer_motifs')


SCRIPTS_DIR = os.path.join(os.getcwd(), 'scripts')


# Set workdir - If running on Flux cluster, logs will be placed in this location
workdir:
    config['flux_log_dir']


# Rules

rule all:
    input:
        expand(os.path.join(DISP_DIR, "{library}.tdf"), library=config['lib_paths'].keys()), #Create tdfs for all samples
        expand(os.path.join(DISP_DIR, "{library}.1m.bw"), library=config['lib_paths'].keys()), #Create bigwigs for all samples
        expand(os.path.join(HOMERPEAK_DIR, "{library}_BLfiltered.hpeaks"), library=config['lib_input'].keys()), #Call peaks for all samples with matched inputs
        expand(os.path.join(HOMERMOTIF_DIR, "{library}"), library=config['lib_homer_fmg_genome'].keys()), #Homermotifs for all samples with a specified genome for homer findMotifsGenome


include:
    "alignment_bwa_aln_se.smk"

rule mark_duplicates:
    input:
        os.path.join(ALIGN_DIR, "{library}.merged.bam")
    output:
        bam = os.path.join(ALIGN_DIR, "{library}.mrkdup.bam"),
        metric = os.path.join(ALIGN_DIR, "{library}.mrkdup.metric")
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
        os.path.join(ALIGN_DIR, "{library}.mrkdup.bam")
    output:
        os.path.join(ALIGN_DIR, "{library}.mrkdup.bai")
    conda: "envs/samtools.yaml"
    shell:
        "samtools index {input} {output}"

rule samtools_prune:
    input:
        bam = os.path.join(ALIGN_DIR, "{library}.mrkdup.bam"),
        bai = os.path.join(ALIGN_DIR, "{library}.mrkdup.bai")
    output:
        bam = temp(os.path.join(PRUNE_DIR, "{library}.pruned.bam"))
    params:
        incl_chr = lambda wildcards: INCLUDE_CHRS[get_genome(wildcards.library)],
        flags = config['samtools_prune_flags']
    conda: "envs/samtools.yaml"
    shell:
        "samtools view -b {params.flags} {input.bam} {params.incl_chr} > {output.bam}"

rule index_pruned:
    input:
        bam = os.path.join(PRUNE_DIR, "{library}.pruned.bam")
    output:
        bai = os.path.join(PRUNE_DIR, "{library}.pruned.bai")
    conda: "envs/samtools.yaml"
    shell:
        "samtools index {input.bam} {output.bai}"

rule igvtools_count_tdf:
    input:
        os.path.join(PRUNE_DIR, "{library}.pruned.bam")
    output:
        os.path.join(DISP_DIR, "{library}.tdf")
    params:
        genome = lambda wildcards: get_genome(wildcards.library),
        args = config['igvtools_count_params']
    conda: "envs/igvtools.yaml"
    shell:
        "igvtools count {params.args} {input} {output} {params.genome}"

rule deeptools_bamcoverage_bw:
    input:
        bam = os.path.join(PRUNE_DIR, "{library}.pruned.bam"),
        bai = os.path.join(PRUNE_DIR, "{library}.pruned.bai")
    output:
        os.path.join(DISP_DIR, "{library}.1m.bw")
    params:
        blacklist = lambda wildcards: config['blacklist'][get_genome(wildcards.library)],
        args = config['deeptools_bamcoverage_params']
    conda: "envs/deeptools.yaml"
    shell:
        "bamCoverage --bam {input.bam} -o {output} -bl {params.blacklist} {params.args}"

rule makeTagDirectory:
    input:
        os.path.join(PRUNE_DIR, "{library}.pruned.bam")
    output:
        directory(os.path.join(HOMERTAG_DIR, "{library}"))
    params:
        genome = lambda wildcards: get_genome(wildcards.library),
        params = config['makeTagDir_params']
    conda: "envs/homer.yaml"
    shell:
        "makeTagDirectory {output} {params.params} -genome {params.genome} {input}"

rule findPeaks:
    input:
        sample = os.path.join(HOMERTAG_DIR, "{library}"),
        input = lambda wildcards: os.path.join(HOMERTAG_DIR, config['lib_input'][wildcards.library])
    output:
        os.path.join(HOMERPEAK_DIR, "{library}.all.hpeaks")
    params:
        config['homer_findPeaks_params']
    conda: "envs/homer.yaml"
    shell:
        "findPeaks {input.sample} -i {input.input} {params} -o {output}"

rule pos2bed:
    input:
        os.path.join(HOMERPEAK_DIR, "{library}.all.hpeaks")
    output:
        os.path.join(HOMERPEAK_DIR, "{library}.all.bed")
    conda: "envs/homer.yaml"
    shell:
        "pos2bed.pl {input} > {output}"

rule blacklist_filter_bed:
    input:
        os.path.join(HOMERPEAK_DIR, "{library}.all.bed"),
    output:
        os.path.join(HOMERPEAK_DIR, "{library}_BLfiltered.bed"),
    params:
        blacklist = lambda wildcards: config['blacklist'][get_genome(wildcards.library)]
    conda: "envs/bedtools.yaml"
    shell:
        "bedtools intersect -a {input} -b {params.blacklist} -v > {output}"

rule keepBedEntriesInHpeaks:
    input:
        filtbed = os.path.join(HOMERPEAK_DIR, "{library}_BLfiltered.bed"),
        allhpeaks = os.path.join(HOMERPEAK_DIR, "{library}.all.hpeaks")
    output:
        os.path.join(HOMERPEAK_DIR, "{library}_BLfiltered.hpeaks")
    conda: "envs/pysam.yaml"
    shell:
        "python {SCRIPTS_DIR}/keepBedEntriesInHpeaks.py -i {input.allhpeaks} -b {input.filtbed} -o {output}"

rule findMotifsGenome:
    input:
        os.path.join(HOMERPEAK_DIR, "{library}_BLfiltered.hpeaks")
    output:
        directory(os.path.join(HOMERMOTIF_DIR, "{library}"))
    params:
        genome = lambda wildcards: config['lib_homer_fmg_genome'][wildcards.library],
        params = config['homer_fmg_params']
    conda: "envs/homer.yaml"
    shell:
        "findMotifsGenome.pl {input} {params.genome} {output} {params.params}"
