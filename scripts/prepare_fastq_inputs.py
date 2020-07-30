#! /usr/bin/env python3

import argparse
import binascii
from collections import defaultdict
import glob
import gzip
import os
import re
import pandas as pd
import warnings

def get_fastqs(fastq_dir):
    """
    Takes a path containing fastqs (ending in '.fastq.gz') and returns a list of the paths to them.
    """
    all_fastqs = glob.glob(os.path.join(fastq_dir, "*.fastq.gz"))

    if not all_fastqs:
        msg = "Error: Cannnot find any files ending in '.fastq.gz' in directory {}.".format(fastq_dir)
        raise RuntimeError(msg)

    return(all_fastqs)

def group_fastqs_by_metadata(fastqs, metadata, indiv_col, group_col, strip_regex):
    """
    Organize a folder of fastqs in a flat structure into a tree structure where the fastqs
    are grouped in folders according to the samples to which they belong.

    fastqs: The result of get_fastqs(), a list of paths to fastqs.
    metadata: A pandas DataFrame relating samples labels to the consituent part labels.
    indiv_col: Column to use for defining constituent parts of samples.
    group_col: Column to use for defining sample names.
    strip_regex: Regular expression to strip from all input filenames.
    """
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
        grouping_dict[group].sort()

    return(grouping_dict)

def validate_input_arguments(metadata, args):
    """
    Checks for existence of the fastq_dir passed to __main__, and that the group_col and indiv_col columns
    are in the metadata DataFrame.
    """
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

def add_readnum_to_filename(filename):
    already_r_match = re.match(r'.*_R([12])(?=[_\.]).*\.fastq\.gz', filename)
    if already_r_match:
        msg = "Warning: filename appears to already have _R preceding read number: {}. Using --add_R flag may be unnecessary.".format(filename)
        warnings.warn(msg)

    sra_type_match = re.match('.*_[12]\.fastq\.gz', filename)
    if sra_type_match:
        modified_filename = re.sub(r"_([12])\.fastq\.gz", "_R\\1.fastq.gz", filename)
    else:
        modified_filename = re.sub(r"(.*)\.fastq\.gz", "\\1_R1.fastq.gz", filename)

    return(modified_filename)

def check_problematic_readnames(filename):
    ''' Very basic check within first line of gzipped fastq files for problematic read names -
.1 & .2 suffixes on read names cause problems for bwa sampe.'''

    with gzip.open(filename, 'rt') as fh:
        for line in fh:
            first_line = line.strip()
            first_readname = first_line.split(" ")[0]
            break

    if re.match('\S+\.[12] ', first_readname):
        msg = "Read name will cause issues with sampe due to .1 or .2 suffix: {}\nIf using fastq-dump, using the -F flag of fastq-dump will give original headers, preventing this issue."
        warnings.warn(msg)



if __name__ == '__main__':

    parser = argparse.ArgumentParser(prog='python prepare_fastq_inputs.py', description = "Rearranges fastq files into a nested folder structure where each sample has a folder with its constituent fastq files.")
    parser.add_argument('-d', '--fastq_dir', required=True, help="Directory of fastq files (all samples)")
    parser.add_argument('-m', '--metadata_csv', required=True, help="CSV file with sample grouping information. Two columns will be used for grouping - for instance 'sample', and 'lane', where a sample may have multiple fastqs from several lanes. In the case of SRA metadata, they could be grouped by 'Run' and 'BioSample', where there may be multiple runs per biosample.")
    parser.add_argument('-g', '--group_col', default = 'BioSample', help="Column to use for defining sample names. Default 'BioSample'")
    parser.add_argument('-i', '--indiv_col', default = 'Run', help="Column to use for defining constituent parts of samples. Default 'Run'")
    parser.add_argument('-r', '--strip_regex', default = '_R?[12]_?[0-9]{0,3}\.fastq\.gz$', help="Regular expression to strip from all input filenames. The remaining filename text after stripping can then be matched to the values in the indiv_col of the metadata_csv. Default '_R?[12]_?[0-9]{0,3}\.fastq\.gz$'")
    parser.add_argument('--add_R', default = False, action = 'store_true', help = "Use this option to transform _1.fastq.gz to _R1.fastq.gz in preparation for config_creator.py")
    parser.add_argument('-n', '--dryrun', default = False, action = 'store_true', help="Perform a dry-run. Do not actually move any files.")

    args = parser.parse_args()

    # Read in the file and verify needed columns exist
    metadata = pd.read_csv(args.metadata_csv)

    validate_input_arguments(metadata, args)

    fastqs = get_fastqs(args.fastq_dir)

    # If fastqs are SRR by filename, error-checking and warnings
    if any([re.match("^SRR\d+_\d\.(fastq|fq)", fq) for fq in fastqs]):
        # Warn if --add_R argument not used for SRA data
        if not args.add_R:
            msg = "Fastq files appear to appear to be from SRA. Make sure to use the --add_R argument as needed"
            warnings.warn(msg)

        # # Check that readnames within SRA fastq files are not problematic (.1 & .2 in read names cause bwa sampe to throw errors)
        # # randomly check 2 fastq files
        # if any(check_problematic_readnames(random.sample(fastqs, 2))):


    # Make the groupings and move the fastqs based on groupings
    groups = group_fastqs_by_metadata(fastqs, metadata, args.indiv_col, args.group_col, args.strip_regex)
    # Move the fastqs based on groupings separately?

    # # Check if fastqs are SRR by filename. If so, :
    # for group in groups.values():
    #     if any([re.match("^SRR.*\.fastq.gz", fq) for fq in group]):
    #         msg_fmt = "# Warning: Combining SRA fastqs can cause bwa sampe to throw errors. To remedy this, run the following: \n" + \
    #         "for f in {} ; do BNAME=$(basename $f '.fastq.gz') ; echo $f ; gunzip -c $f | sed -E 's/(^[@+]SRR[0-9]+\.[0-9]+)\.[12]/\1/' | gzip -c > ${{BNAME.fastq.gz\n"
    #         msg_fmt.format(" ".join(group))
    #         import pdb; pdb.set_trace()

    for key in groups:
        # Create the tree structure
        new_dir = os.path.join(args.fastq_dir, key)
        if not args.dryrun:
            os.makedirs(new_dir, exist_ok=True)

        for src_file in groups[key]:
            if args.add_R:
                dest_file = add_readnum_to_filename(src_file)
            else:
                dest_file = src_file
            # Create the hardlinks in the tree structure
            src_path = os.path.join(args.fastq_dir, src_file)
            dest_path = os.path.join(args.fastq_dir, key, dest_file)
            if args.dryrun:
                print('{} -> {}'.format(src_path, dest_path))
            else:
                os.link(src_path, dest_path)
