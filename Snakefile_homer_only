import os
import sys
import functools


# Helper functions

def get_genome(library):
    return(config['lib_genome'][library])


# RESULT PATHS

prefix_results = functools.partial(os.path.join, config['results_dir'])
HOMERPEAK_DIR = prefix_results('peaks')
HOMERMOTIF_DIR = prefix_results('homer_motifs')


# Load modules

#The loaded software versions are explicitly stated here.
#If additional functions are added above without explicit versions in the string, add inline comments to indicate
software_strings = [
    "module load Bioinformatics ; "
    "module load samtools/1.9 ;",
    "module load bedtools2/2.25.0 ;",
    "module load homer/4.8 ;"
]
shell.prefix("".join(software_strings))


# Set workdir - If running on Flux cluster, logs will be placed in this location
workdir:
    config['flux_log_dir']


# Rules

rule all:
    input:
        expand(os.path.join(HOMERPEAK_DIR, "{library}_BLfiltered.hpeaks"), library=config['lib_sample_tagdir'].keys()),
        expand(os.path.join(HOMERMOTIF_DIR, "{library}"), library=config['lib_sample_tagdir'].keys())



rule findPeaks:
    input:
        sample = lambda wildcards: config['lib_sample_tagdir'][wildcards.library],
        input = lambda wildcards: config['lib_input_tagdir'][wildcards.library]
    output:
        os.path.join(HOMERPEAK_DIR, "{library}.all.hpeaks")
    params:
        config['homer_findPeaks_params']
    shell:
        "findPeaks {input.sample} -i {input.input} {params} -o {output}"

rule pos2bed:
    input:
        os.path.join(HOMERPEAK_DIR, "{library}.all.hpeaks")
    output:
        os.path.join(HOMERPEAK_DIR, "{library}.all.bed")
    shell:
        "pos2bed.pl {input} > {output}"

rule blacklist_filter_bed:
    input:
        os.path.join(HOMERPEAK_DIR, "{library}.all.bed"),
    output:
        os.path.join(HOMERPEAK_DIR, "{library}_BLfiltered.bed"),
    params:
        blacklist = lambda wildcards: config['blacklist'][get_genome(wildcards.library)]
    shell:
        "bedtools intersect -a {input} -b {params.blacklist} -v > {output}"

rule keepBedEntriesInHpeaks:
    input:
        filtbed = os.path.join(HOMERPEAK_DIR, "{library}_BLfiltered.bed"),
        allhpeaks = os.path.join(HOMERPEAK_DIR, "{library}.all.hpeaks")
    output:
        os.path.join(HOMERPEAK_DIR, "{library}_BLfiltered.hpeaks")
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
    shell:
        "findMotifsGenome.pl {input} {params.genome} {output} {params.params}"
