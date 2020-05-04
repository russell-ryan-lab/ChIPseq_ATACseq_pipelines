#! /usr/bin/env python3

import argparse
from collections import defaultdict
import glob
import os
import re
import pandas as pd

def get_fastqs(fastq_dir):
    fastq_glob = glob.glob(os.path.join(fastq_dir, "*.fastq*"))
    fq_glob = glob.glob(os.path.join(fastq_dir, "*.fq*"))

    all_fastqs = fastq_glob + fq_glob

    if not all_fastqs:
        msg = "Error: {} does not contain any fastqs.".format(fastq_dir)
        raise RuntimeError(msg)

    return(all_fastqs)

def group_fastqs_by_metadata(fastqs, metadata, indiv_col, group_col, strip_regex):
    # Strip filenames to get remaining identifiers
    fq_bnames_dict = defaultdict(list)
    fq_filenames = [os.path.basename(x) for x in fastqs]
    for fn in fq_filenames:
        basename = re.sub(strip_regex, "", fn)
        fq_bnames_dict[basename].append(fn)

    grouping_dict = defaultdict(list)
    groups = metadata[group_col].unique().tolist()
    # Change index to indiv_col so that it will return desired index vals (from indiv_col) from group lookup
    metadata.set_index(indiv_col, inplace=True)
    for group in groups:
        for runs in metadata.groupby(group_col).groups[group].to_list():
            grouping_dict[group].extend(fq_bnames_dict[runs])

    return(grouping_dict)

def validate_input_arguments(metadata, args):
    required_cols = set([args.group_col, args.indiv_col])
    actual_cols = set(metadata.columns)
    missing_cols = required_cols - actual_cols

    if missing_cols:
        msg = "Metadata file {} is missing required columns: {}".format(args.metadata_csv, missing_cols)
        raise RuntimeError(msg)

    # Check directory exists and has fastqs
    if not os.path.isdir(args.fastq_dir):
        msg = "Error: Is {} a valid directory?".format(args.fastq_dir)
        raise RuntimeError(msg)

if __name__ == '__main__':

    parser = argparse.ArgumentParser(prog='python prepare_fastq_inputs.py', description = "Rearranges fastq files into a nested folder structure where each sample has a folder with its constituent fastq files.")
    parser.add_argument('-d', '--fastq_dir', required=True, help="Directory of fastq files (all samples)")
    parser.add_argument('-m', '--metadata_csv', required=True, help="CSV file with sample grouping information. Two columns will be used for grouping - for instance 'sample', and 'lane', where a sample may have multiple fastqs from several lanes. In the case of SRA metadata, they could be grouped by 'Run' and 'BioSample', where there may be multiple runs per biosample.")
    parser.add_argument('-g', '--group_col', default = 'BioSample', help="Column to use for defining sample names. Default 'BioSample'")
    parser.add_argument('-i', '--indiv_col', default = 'Run', help="Column to use for defining constituent parts of samples. Default 'Run'")
    parser.add_argument('-r', '--strip_regex', default = '_R*[12](\.fastq|\.fq)(\.gz)*$', help="Regular expression to strip from all input filenames. Default '_R*[12](\.fastq|\.fq)(\.gz)*$'")

    args = parser.parse_args()

    # Read in the file and verify needed columns exist
    metadata = pd.read_csv(args.metadata_csv)

    validate_input_arguments(metadata, args)

    fastqs = get_fastqs(args.fastq_dir)

    # Make the groupings and move the fastqs based on groupings
    groups = group_fastqs_by_metadata(fastqs, metadata, args.indiv_col, args.group_col, args.strip_regex)
    # Move the fastqs based on groupings separately?

    for key in groups:
        # Create the tree structure
        new_dir = os.path.join(args.fastq_dir, key)
        os.mkdir(new_dir)

        for file in groups[key]:
            # Create the hardlinks in the tree structure
            src_file = os.path.join(args.fastq_dir, file)
            dest_file = os.path.join(args.fastq_dir, key, file)

            os.link(src_file, dest_file)
