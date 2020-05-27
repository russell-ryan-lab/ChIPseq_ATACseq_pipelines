#!/usr/bin/env python3

import argparse
import os
import re
import shutil
import warnings


THIS_SCRIPT_DIR = os.path.realpath(os.path.dirname(__file__))
BASE_DIR = os.path.abspath(os.path.join(THIS_SCRIPT_DIR, '..'))

files_to_modify_basepaths = [
    os.path.join(BASE_DIR, 'user_tests', 'atac_test.sh'),
    os.path.join(BASE_DIR, 'user_tests', 'chip_test.sh'),
    os.path.join(BASE_DIR, 'user_tests', 'chip_mixed_test.sh'),
    os.path.join(BASE_DIR, 'config', 'ATAC_general.yaml'),
    os.path.join(BASE_DIR, 'config', 'ChIP_histone_general_pe.yaml'),
    os.path.join(BASE_DIR, 'config', 'ChIP_histone_general_se.yaml'),
    os.path.join(BASE_DIR, 'config', 'ChIP_TF_general_pe.yaml'),
    os.path.join(BASE_DIR, 'config', 'ChIP_TF_general_se.yaml'),
    os.path.join(BASE_DIR, 'data', 'atac_test_data', 'atac_test_samplesheet.csv'),
    os.path.join(BASE_DIR, 'data', 'sra_chip_test_data', 'sra_chip_samplesheet.csv'),
    os.path.join(BASE_DIR, 'data', 'sra_mixed_chip_test_data', 'sra_chip_histone_input.csv'),
    os.path.join(BASE_DIR, 'data', 'sra_mixed_chip_test_data', 'sra_chip_histone_IP.csv')
]

files_to_modify_envpaths = [
    os.path.join(BASE_DIR, 'user_tests', 'atac_test.sh'),
    os.path.join(BASE_DIR, 'user_tests', 'chip_test.sh'),
    os.path.join(BASE_DIR, 'user_tests', 'chip_mixed_test.sh')
]

cluster_config = os.path.join(BASE_DIR, 'config', 'cluster_config.yaml')


def replace_lines(file_lines, match_string, replacement):
    '''Takes array of lines of text read in from input file, regex string, and replacement,
    Returns array of lines with modified basepath.'''
    file_lines = [re.sub(match_string, replacement, line) for line in file_lines]
    return(file_lines)

if __name__ == '__main__':

    parser = argparse.ArgumentParser(prog='python setup.py', description = "Changes values in CSV and .yaml files to simplify setup process for end users. Modifies base paths to match the directory of the repository, other values are filled in based on command-line arguments.")
    parser.add_argument('-a', '--hpc_account', required=True, help="Value to place in cluster_config.yaml for High Performance Compute (HPC) account.")
    parser.add_argument('-e', '--env_path', required=True, help="Path to atac_chip_pipeline conda environment.")
    parser.add_argument('-s', '--save_original', default = False, action = 'store_true', help="Keep the unmodified files with an '.orig' extension added on.")

    args = parser.parse_args()

    # Error checking - does env_path exist?
    if not os.path.isdir(args.env_path):
        msg = 'Error: Environment path {} is not a valid directory.'.format(args.env_path)
        raise RuntimeError(msg)
    # Does it look like a conda environment (with snakemake?)
    if not os.path.exists(os.path.join(args.env_path, 'bin', 'snakemake')):
        msg = 'Warning: Environment path {} does not look like the desired conda environment - Cannot find snakemake.'.format(args.env_path)
        warnings.warn(msg)

    all_files = set(files_to_modify_basepaths + files_to_modify_envpaths + [cluster_config])

    for filename in all_files:
        lines = []
        with open(filename, 'r') as file_handle:
            lines = list(file_handle) # Read file into memory
            if filename in files_to_modify_basepaths:
                lines = replace_lines(lines, '/nfs/med-bfx-activeprojects/Ryan_rjhryan_CU1/ChIPseq_ATACseq_pipelines', BASE_DIR)
            if filename in files_to_modify_envpaths:
                lines = replace_lines(lines, '/nfs/turbo/<lab-turbo>/envs/atac_chip_pipeline', args.env_path)
            if filename == cluster_config:
                lines = replace_lines(lines, 'cgates1', args.hpc_account)

        # Optionally copy original file first
        if args.save_original:
            copy_filename = filename + '.orig'
            shutil.copyfile(filename, copy_filename)

        with open(filename, 'w') as file_handle_w:
            file_handle_w.writelines(lines)
