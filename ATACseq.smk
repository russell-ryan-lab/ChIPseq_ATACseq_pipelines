import os
import sys
import functools


# GENERIC DATA

ORGANISMS = {
    'rn4': 'rat',
    'rn5': 'rat',
    'mm9': 'mouse',
    'mm10': 'mouse',
    'hg19': 'human',
    'hg38': 'human'
}

INCLUDE_CHRS = {
    'hg19': ['chr{}'.format(i) for i in range(1, 23)] + ["chrX", "chrY"],
    'hg38': ['chr{}'.format(i) for i in range(1, 23)] + ["chrX", "chrY"],
    'mm9': ['chr{}'.format(i) for i in range(1, 20)] + ["chrX", "chrY"],
    'mm10': ['chr{}'.format(i) for i in range(1, 20)] + ["chrX", "chrY"],
    'rn4': ['chr{}'.format(i) for i in range(1, 21)] + ["chrX", "chrY"],
    'rn5': ['chr{}'.format(i) for i in range(1, 21)] + ["chrX", "chrY"]
}

MACS2_GENOME_SIZE = {
    'rn4': 'mm',
    'rn5': 'mm',
    'mm9': 'mm',
    'mm10': 'mm',
    'hg19': 'hs',
    'hg38': 'hs'
}

# Helper functions

def get_genome(library):
    return(config['lib_genome'][library])


# RESULT PATHS

prefix_results = functools.partial(os.path.join, config['results_dir'])
CONCAT_READS_DIR = prefix_results('concat_reads')
ALIGN_DIR = prefix_results('aligned')
PRUNE_DIR = prefix_results('pruned')
MACS2_DIR = prefix_results('macs2')
ATAQV_DIR = prefix_results('ataqv')
DISP_DIR = prefix_results('bigwig')


# Set workdir - Snakemake will be run from this location.
workdir:
    config['results_dir']


# Rules

rule all:
    input:
        expand(os.path.join(MACS2_DIR, "{library}_BLfiltered.narrowPeak"), library=config['lib_paths'].keys()),
        expand(os.path.join(MACS2_DIR, "{library}_summits_BLfiltered.bed"), library=config['lib_paths'].keys()),
        expand(os.path.join(ATAQV_DIR, "{library}.ataqv.out"), library=config['lib_paths'].keys()),
        expand(os.path.join(DISP_DIR, "{library}.1m.bw"), library=config['lib_paths'].keys()),


include:
    "alignment_bwa_aln_pe.smk"


rule concatenate_reads:
    input:
        lambda wildcards: config['lib_paths'][wildcards.library][wildcards.read]
    output:
        os.path.join(CONCAT_READS_DIR, "{library}_R{read}.fastq.gz")
    shell:
        "cat {input} > {output}"

rule mark_duplicates:
    input:
        # This input comes from rule included alignment module e.g. Snakefile_alignment_bwa_aln_pe
        os.path.join(ALIGN_DIR, "{library}.sorted.bam")
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

rule prune:
    input:
        bam = os.path.join(ALIGN_DIR, "{library}.mrkdup.bam"),
        bai = os.path.join(ALIGN_DIR, "{library}.mrkdup.bai")
    output:
        bam = os.path.join(PRUNE_DIR, "{library}.pruned.bam"),
        bai = os.path.join(PRUNE_DIR, "{library}.pruned.bai")
    params:
        incl_chr = lambda wildcards: INCLUDE_CHRS[get_genome(wildcards.library)],
        flags = config['samtools_prune_flags']
    conda: "envs/samtools.yaml"
    shell:
        "samtools view -b {params.flags} {input.bam} {params.incl_chr} > {output.bam} ; "
        "samtools index {output.bam} {output.bai}"

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

rule peaks:
    input:
        os.path.join(PRUNE_DIR, "{library}.pruned.bam")
    output:
        temp(os.path.join(MACS2_DIR, "{library}_peaks.narrowPeak")),
        os.path.join(MACS2_DIR, "{library}_peaks.xls"),
        temp(os.path.join(MACS2_DIR, "{library}_summits.bed"))
    params:
        name = "{library}",
        genome = lambda wildcards: MACS2_GENOME_SIZE[get_genome(wildcards.library)],
        outdir = MACS2_DIR
    conda: "envs/macs2.yaml"
    shell:
        "macs2 callpeak -t {input} --outdir {params.outdir} -n {params.name} -f BAMPE -g {params.genome} --keep-dup all"

rule blacklist_filter:
    input:
        narrowpeak = os.path.join(MACS2_DIR, "{library}_peaks.narrowPeak"),
        summits = os.path.join(MACS2_DIR, "{library}_summits.bed")
    output:
        narrowpeak = os.path.join(MACS2_DIR, "{library}_BLfiltered.narrowPeak"),
        summits = os.path.join(MACS2_DIR, "{library}_summits_BLfiltered.bed")
    params:
        blacklist = lambda wildcards: config['blacklist'][get_genome(wildcards.library)]
    conda: "envs/bedtools.yaml"
    shell:
        "bedtools intersect -a {input.narrowpeak} -b {params.blacklist} -v > {output.narrowpeak} ; "
        "bedtools intersect -a {input.summits} -b {params.blacklist} -v > {output.summits} ; "

rule ataqv:
    input:
        md_bam = os.path.join(ALIGN_DIR, "{library}.mrkdup.bam"),
        peaks = os.path.join(MACS2_DIR, '{library}_peaks.narrowPeak')
    output:
        metrics = os.path.join(ATAQV_DIR, '{library}.ataqv.json.gz'),
        stdout_destination = os.path.join(ATAQV_DIR, '{library}.ataqv.out')
    params:
        samplename = lambda wildcards: config['lib_samplename'][wildcards.library],
        organism = lambda wildcards: ORGANISMS[get_genome(wildcards.library)],
        tss_file = lambda wildcards: config['tss'][get_genome(wildcards.library)],
        blacklist = lambda wildcards: config['blacklist'][get_genome(wildcards.library)]
    conda: "envs/ataqv.yaml"
    shell:
        "ataqv --peak-file {input.peaks} --name {params.samplename} --metrics-file {output.metrics} --excluded-region-file {params.blacklist} --tss-file {params.tss_file} --ignore-read-groups {params.organism} {input.md_bam} > {output.stdout_destination}"
