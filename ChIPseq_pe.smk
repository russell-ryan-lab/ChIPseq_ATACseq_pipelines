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

# Conditional inputs to all rule
if config.get('sample_input'):
    if config.get('sample_homer_fmg_genome'):
        all_input = [
            expand(os.path.join(DISP_DIR, "{sample}.tdf"), sample=config['sample_paths'].keys()), #Create tdfs for all samples
            expand(os.path.join(DISP_DIR, "{sample}.1m.bw"), sample=config['sample_paths'].keys()), #Create bigwigs for all samples
            expand(os.path.join(HOMERPEAK_DIR, "{paramset}", "{sample}.BLfiltered.hpeaks"), paramset=config['homer_findPeaks_params'].keys(), sample=config['sample_input'].keys()), #Call peaks for all samples with matched inputs
            expand(os.path.join(HOMERMOTIF_DIR, "{paramset}", "{sample}"), paramset=config['homer_findPeaks_params'].keys(), sample=config['sample_homer_fmg_genome'].keys()), #Homermotifs for all samples with a specified genome for homer findMotifsGenome
        ]
    else:
        all_input = [
            expand(os.path.join(DISP_DIR, "{sample}.tdf"), sample=config['sample_paths'].keys()), #Create tdfs for all samples
            expand(os.path.join(DISP_DIR, "{sample}.1m.bw"), sample=config['sample_paths'].keys()), #Create bigwigs for all samples
            expand(os.path.join(HOMERPEAK_DIR, "{paramset}", "{sample}.BLfiltered.hpeaks"), paramset=config['homer_findPeaks_params'].keys(), sample=config['sample_input'].keys()), #Call peaks for all samples with matched inputs
        ]
else:
    all_input = [
        expand(os.path.join(DISP_DIR, "{sample}.tdf"), sample=config['sample_paths'].keys()), #Create tdfs for all samples
        expand(os.path.join(DISP_DIR, "{sample}.1m.bw"), sample=config['sample_paths'].keys()), #Create bigwigs for all samples
        expand(os.path.join(HOMERPEAK_DIR, "{paramset}", "{sample}.BLfiltered.hpeaks"), paramset=config['homer_findPeaks_params'].keys(), sample=config['sample_input'].keys()), #Call peaks for all samples with matched inputs
        expand(os.path.join(HOMERMOTIF_DIR, "{paramset}", "{sample}"), paramset=config['homer_findPeaks_params'].keys(), sample=config['sample_homer_fmg_genome'].keys()), #Homermotifs for all samples with a specified genome for homer findMotifsGenome
    ]

# Print pipeline version number to log
include: "version.smk"
version_string = "Pipeline version: {}\n".format(version)
logger.logger.info(version_string)

# Rules

rule all:
    input:
        all_input

include: "rules/alignment_bwa_aln_pe.smk"
include: "rules/homer.smk"

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

# samtools is available in the parent environment atac_chip_pipeline
rule index_dupmarked_bams:
    input:
        os.path.join(ALIGN_DIR, "{sample}.mrkdup.bam")
    output:
        os.path.join(ALIGN_DIR, "{sample}.mrkdup.bai")
    shell:
        "samtools index {input} {output}"

# samtools is available in the parent environment atac_chip_pipeline
rule samtools_prune:
    input:
        bam = os.path.join(ALIGN_DIR, "{sample}.mrkdup.bam"),
        bai = os.path.join(ALIGN_DIR, "{sample}.mrkdup.bai")
    output:
        bam = temp(os.path.join(PRUNE_DIR, "{sample}.stpruned.bam"))
    params:
        incl_chr = lambda wildcards: INCLUDE_CHRS[get_genome(wildcards.sample)],
        flags = config['samtools_prune_flags']
    shell:
        "samtools view -b {params.flags} {input.bam} {params.incl_chr} > {output.bam}"

# samtools is available in the parent environment atac_chip_pipeline
rule namesort_st_pruned:
    input:
        os.path.join(PRUNE_DIR, "{sample}.stpruned.bam")
    output:
        temp(os.path.join(PRUNE_DIR, "{sample}.ns.bam"))
    shell:
        "samtools sort -n -o {output} {input}"

rule X0_pair_filter:
    input:
        os.path.join(PRUNE_DIR, "{sample}.ns.bam")
    output:
        temp(os.path.join(PRUNE_DIR, "{sample}.x0_filtered.bam"))
    params:
        config['X0_pair_filter_params']
    conda: "envs/pysam.yaml"
    shell:
        "python {SCRIPTS_DIR}/X0_pair_filter.py {params} -b {input} -o {output}"

# samtools is available in the parent environment atac_chip_pipeline
rule coordsort_index_final_pruned:
    input:
        os.path.join(PRUNE_DIR, "{sample}.x0_filtered.bam")
    output:
        bam = os.path.join(PRUNE_DIR, "{sample}.pruned.bam"),
        bai = os.path.join(PRUNE_DIR, "{sample}.pruned.bai")
    shell:
        "samtools sort -o {output.bam} {input} ;"
        "samtools index {output.bam} {output.bai}"

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
        blacklist = lambda wildcards: "-bl {}".format(config['blacklist'][get_genome(wildcards.sample)]) if config.get('deeptools_bamcoverage_use_blacklist') == True else '', # deeptools_bamcoverage_use_blacklist should be True or False. If True, add `-bl /path/to/blacklist` to command. Otherwise empty string
        args = config['deeptools_bamcoverage_params']
    conda: "envs/deeptools.yaml"
    shell:
        "bamCoverage --bam {input.bam} -o {output} {params.blacklist} {params.args}"
