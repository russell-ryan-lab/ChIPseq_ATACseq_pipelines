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

def get_genome(sample):
    return(config['sample_genome'][sample])


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
        expand(os.path.join(MACS2_DIR, "{sample}.BLfiltered.narrowPeak"), sample=config['sample_paths'].keys()),
        expand(os.path.join(MACS2_DIR, "{sample}.summits.BLfiltered.bed"), sample=config['sample_paths'].keys()),
        expand(os.path.join(ATAQV_DIR, 'viewer', 'index.html')),
        expand(os.path.join(DISP_DIR, "{sample}.1m.bw"), sample=config['sample_paths'].keys()),


include:
    "alignment_bwa_aln_pe.smk"


rule concatenate_reads:
    input:
        lambda wildcards: config['sample_paths'][wildcards.sample][wildcards.read]
    output:
        os.path.join(CONCAT_READS_DIR, "{sample}_R{read}.fastq.gz")
    shell:
        "cat {input} > {output}"

rule mark_duplicates:
    input:
        # This input comes from rule included alignment module e.g. Snakefile_alignment_bwa_aln_pe
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

rule prune:
    input:
        bam = os.path.join(ALIGN_DIR, "{sample}.mrkdup.bam"),
        bai = os.path.join(ALIGN_DIR, "{sample}.mrkdup.bai")
    output:
        bam = os.path.join(PRUNE_DIR, "{sample}.pruned.bam"),
        bai = os.path.join(PRUNE_DIR, "{sample}.pruned.bai")
    params:
        incl_chr = lambda wildcards: INCLUDE_CHRS[get_genome(wildcards.sample)],
        flags = config['samtools_prune_flags']
    conda: "envs/samtools.yaml"
    shell:
        "samtools view -b {params.flags} {input.bam} {params.incl_chr} > {output.bam} ; "
        "samtools index {output.bam} {output.bai}"

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

rule peaks:
    input:
        os.path.join(PRUNE_DIR, "{sample}.pruned.bam")
    output:
        temp(os.path.join(MACS2_DIR, "{sample}.peaks.narrowPeak")),
        os.path.join(MACS2_DIR, "{sample}.peaks.xls"),
        temp(os.path.join(MACS2_DIR, "{sample}.summits.bed"))
    params:
        name = "{sample}",
        genome = lambda wildcards: MACS2_GENOME_SIZE[get_genome(wildcards.sample)],
        outdir = MACS2_DIR
    conda: "envs/macs2.yaml"
    shell:
        "macs2 callpeak -t {input} --outdir {params.outdir} -n {params.name} -f BAMPE -g {params.genome} --keep-dup all ; "
        # aligning filename formats - separated by '.'s preferred
        "mv {params.outdir}/{wildcards.sample}_peaks.xls {params.outdir}/{wildcards.sample}.peaks.xls ; "
        "mv {params.outdir}/{wildcards.sample}_peaks.narrowPeak {params.outdir}/{wildcards.sample}.peaks.narrowPeak ; "
        "mv {params.outdir}/{wildcards.sample}_summits.bed {params.outdir}/{wildcards.sample}.summits.bed"

rule blacklist_filter:
    input:
        narrowpeak = os.path.join(MACS2_DIR, "{sample}.peaks.narrowPeak"),
        summits = os.path.join(MACS2_DIR, "{sample}.summits.bed")
    output:
        narrowpeak = os.path.join(MACS2_DIR, "{sample}.BLfiltered.narrowPeak"),
        summits = os.path.join(MACS2_DIR, "{sample}.summits.BLfiltered.bed")
    params:
        blacklist = lambda wildcards: config['blacklist'][get_genome(wildcards.sample)]
    conda: "envs/bedtools.yaml"
    shell:
        "bedtools intersect -a {input.narrowpeak} -b {params.blacklist} -v > {output.narrowpeak} ; "
        "bedtools intersect -a {input.summits} -b {params.blacklist} -v > {output.summits} ; "

rule ataqv_get_stats:
    input:
        md_bam = os.path.join(ALIGN_DIR, "{sample}.mrkdup.bam"),
        peaks = os.path.join(MACS2_DIR, '{sample}.peaks.narrowPeak')
    output:
        metrics = os.path.join(ATAQV_DIR, '{sample}.ataqv.json.gz'),
        stdout_destination = os.path.join(ATAQV_DIR, '{sample}.ataqv.out')
    params:
        samplename = '{wildcards.sample}',
        organism = lambda wildcards: ORGANISMS[get_genome(wildcards.sample)],
        tss_file = lambda wildcards: config['tss'][get_genome(wildcards.sample)],
        blacklist = lambda wildcards: config['blacklist'][get_genome(wildcards.sample)]
    conda: "envs/ataqv.yaml"
    shell:
        "ataqv --peak-file {input.peaks} --name {params.samplename} --metrics-file {output.metrics} --excluded-region-file {params.blacklist} --tss-file {params.tss_file} --ignore-read-groups {params.organism} {input.md_bam} > {output.stdout_destination}"

rule ataqv_make_viewer:
    input:
        metrics = expand(os.path.join(ATAQV_DIR, '{sample}.ataqv.json.gz'), sample=config['sample_paths'].keys())
    output:
        index = os.path.join(ATAQV_DIR, 'viewer', 'index.html')
    params:
        output_dir = os.path.join(ATAQV_DIR, 'viewer')
    conda: "envs/ataqv.yaml"
    shell:
        "mkarv --force {params.output_dir} {input.metrics}" # Use --force to overwrite the viewer/ directory or snakemake complains
