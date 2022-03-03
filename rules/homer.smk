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

# rule findPeaks_noinput:
#     input:
#         sample = os.path.join(HOMERTAG_DIR, "{sample}"),
#     output:
#         os.path.join(HOMERPEAK_DIR, "{paramset}", "{sample}.all.hpeaks")
#     params:
#         lambda wildcards: config['homer_findPeaks_params'][wildcards.paramset]
#     shell:
#         "findPeaks {input.sample} {params} -o {output}"

# rule findPeaks:
#     input:
#         sample = os.path.join(HOMERTAG_DIR, "{sample}"),
#         input = lambda wildcards: os.path.join(HOMERTAG_DIR, config['sample_input'][wildcards.sample])
#     output:
#         os.path.join(HOMERPEAK_DIR, "{paramset}", "{sample}.all.hpeaks")
#     params:
#         lambda wildcards: config['homer_findPeaks_params'][wildcards.paramset]
#     shell:
#         "findPeaks {input.sample} -i {input.input} {params} -o {output}"

rule findPeaks:
    input:
        sample = os.path.join(HOMERTAG_DIR, "{sample}"),
        input = lambda wildcards: os.path.join(HOMERTAG_DIR, config['sample_input'][wildcards.sample])
    output:
        os.path.join(HOMERPEAK_DIR, "{paramset}", "{sample}.all.hpeaks")
    params:
        lambda wildcards: config['homer_findPeaks_params'][wildcards.paramset]
    shell:
        "findPeaks {input.sample} {params} -o {output}" if input.sample == input.input else "findPeaks {input.sample} -i {input.input} {params} -o {output}"

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
        os.path.join(HOMERPEAK_DIR, "{paramset}", "{sample}.BLfiltered.bed"),
    params:
        blacklist = lambda wildcards: config['blacklist'][get_genome(wildcards.sample)]
    conda: "../envs/bedtools.yaml"
    shell:
        "bedtools intersect -a {input} -b {params.blacklist} -v > {output}"

rule keepBedEntriesInHpeaks:
    input:
        filtbed = os.path.join(HOMERPEAK_DIR, "{paramset}", "{sample}.BLfiltered.bed"),
        allhpeaks = os.path.join(HOMERPEAK_DIR, "{paramset}", "{sample}.all.hpeaks")
    output:
        os.path.join(HOMERPEAK_DIR, "{paramset}", "{sample}.BLfiltered.hpeaks")
    conda: "../envs/pysam.yaml"
    shell:
        "python {SCRIPTS_DIR}/keepBedEntriesInHpeaks.py -i {input.allhpeaks} -b {input.filtbed} -o {output}"

rule findMotifsGenome:
    input:
        os.path.join(HOMERPEAK_DIR, "{paramset}", "{sample}.BLfiltered.hpeaks")
    output:
        directory(os.path.join(HOMERMOTIF_DIR, "{paramset}", "{sample}"))
    params:
        genome = lambda wildcards: config['sample_homer_fmg_genome'][wildcards.sample],
        params = config['homer_fmg_params']
    shell:
        "findMotifsGenome.pl {input} {params.genome} {output} {params.params}"
