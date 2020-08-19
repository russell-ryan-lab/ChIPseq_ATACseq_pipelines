import os
import sys
import functools


# Helper functions

def get_genome(sample):
    return(config['sample_genome'][sample])


# RESULT PATHS

prefix_results = functools.partial(os.path.join, config['results_dir'])
PRUNE_DIR = prefix_results('pruned')
HOMERTAG_DIR = prefix_results('tag')
HOMERPEAK_DIR = prefix_results('peaks')
HOMERMOTIF_DIR = prefix_results('homer_motifs')

SCRIPTS_DIR = os.path.join(workflow.basedir, 'scripts')

# Set workdir - Snakemake will be run from this location.
workdir:
    config['results_dir']

# Error-checking here produces informative error messages for missing config keys
for key in ['sample_input', 'sample_homer_fmg_genome']:
    if not config.get(key):
        msg = "\nMissing required config value: {}\n".format(key)
        raise(RuntimeError(msg))

# Print pipeline version number to log
include: "version.smk"
version_string = "Pipeline version: {}\n".format(version)
logger.logger.info(version_string)

# Rules

rule all:
    input:
        expand(os.path.join(HOMERPEAK_DIR, "{paramset}", "{sample}.BLfiltered.hpeaks"), paramset=config['homer_findPeaks_params'].keys(), sample=config['sample_input'].keys()),
        expand(os.path.join(HOMERMOTIF_DIR, "{paramset}", "{sample}"), paramset=config['homer_findPeaks_params'].keys(), sample=config['sample_input'].keys())


include: 'rules/homer.smk'
