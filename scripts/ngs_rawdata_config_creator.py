#!/usr/bin/env python3

import argparse
from collections import defaultdict
import glob
import json
import os
import pandas as pd
import re
import sys
import yaml


#https://stackoverflow.com/a/5369814
class Tree(defaultdict):
    def __init__(self, value=None):
        super(Tree, self).__init__(Tree)
        self.value = value


def parse_per_lib(lib_table):
    """
    Takes a pandas df as input, and returns a dict of dicts reflecting the per-library data
    """
    if not 'lib' in lib_table.columns:
        raise RuntimeError("Missing required 'lib' column specifying library/sample IDs")
    if not 'basepath' in lib_table.columns:
        raise RuntimeError("Missing required 'basepath' column")

    per_lib_dict = dict()

    lib_basepaths = lib_table[['lib', 'basepath']].values

    per_lib_dict['lib_paths'] = assign_libpaths(lib_basepaths)

    othercols = lib_table.columns.drop(['lib', 'basepath'])
    for c in othercols:
        two_cols = ['lib', c]
        combined_name = "_".join(two_cols)
        col_dict = dict(lib_table[two_cols].dropna().values)
        per_lib_dict[combined_name] = col_dict

    return(per_lib_dict)


def assign_libpaths(lib_basepaths):
    """
    Takes an array of arrays as input, each inner array having two items,
    library name and basepath. Calls basepath_to_filepathsdict for each library
    and adds each libraries' filepathsdict to the libpaths dict
    """
    libpaths_dict = dict()
    for row in lib_basepaths:
        lib, path = row
        if args.simulate_single_lane:
            #If simulating single lane and capture regex is default, replace it (disregard lane)
            if args.capture_regex == ".*_L(\d+)_R(\d+).*\.fastq\.gz":
                capture_regex = ".*_R(\d+).*\.fastq\.gz"
            #Otherwise, use user-specified capture_regex
            else:
                capture_regex = args.capture_regex
            #If samples do not contain lane information, use different capture regex which only captures read number
            libpaths_dict[lib] = basepath_to_filepathsdict(path, args.file_glob, capture_regex)
        else:
            libpaths_dict[lib] = basepath_to_filepathsdict(path, args.file_glob, args.capture_regex)
    return(libpaths_dict)


def basepath_to_filepathsdict(basepath, glob_regex, capture_regex):
    """
    Takes a basepath, and two regex strings, uses glob to find corresponding fastq files
    Then subclassifies these files by lane (readgroup) and read orientation
    and places these into a dict
    """
    all_fastqs = glob.glob(os.path.join(basepath, glob_regex))

    if len(all_fastqs) == 0:
        raise RuntimeError("Input files not found in the directory " + basepath + "\nNote that inputs are found using the following shell glob: " + glob_regex)

    nested_dict = lambda: defaultdict(nested_dict)
    readgroups = nested_dict()

    for fq in all_fastqs:
        basename = os.path.basename(fq)
        rmatch = re.match(capture_regex, basename)
        if not rmatch:
            msg_fmt = ("\nFile {} did not match regular expression {}. "
                "Could not capture desired group(s).\n"
                "If filenames resemble sample_R1.fastq.gz, use --simulate_single_lane flag. "
                "Note that all input files should be formatted consistently. "
                "File glob and capture regex can be controlled with --file_glob and --capture_regex if desired.")
            raise RuntimeError(msg_fmt.format(basename, capture_regex))
        if rmatch.group(0) == basename:
            #If samples do not contain lane information, treat them all as lane 001
            if args.simulate_single_lane:
                lane = "001"
                read = rmatch.group(1)
                if not len(all_fastqs) in [1,2]:
                    raise ValueError("--simulate_single_lane flag can only be used if there are one or two fastq's per sample. There are " + str(len(all_fastqs)) + " in " + basepath)
            else:
                lane = rmatch.group(1)
                read = rmatch.group(2)
            #Add fastq to dict
            readgroups[lane][read] = fq

    return(readgroups)

def read_input(input_filename):
    if input_filename.endswith('.yaml') or input_filename.endswith('.yml'):
        try:
            with open(input_filename) as infile:
                config_dict = yaml.load(infile, Loader=yaml.SafeLoader)
        except:
            print("Could not load input file {}. Assuming YAML input.".format(input_filename))

    elif input_filename.endswith('.json'):
        try:
            with open(input_filename) as infile:
                config_dict = json.load(infile)
        except:
            print("Could not load input file {}. Assuming JSON input.".format(input_filename))

    return config_dict


if __name__ == '__main__':

    parser = argparse.ArgumentParser(prog='python ngs_rawdata_config_creator.py', description = "")
    parser.add_argument('-g', '--general_input', required=True, help="json or yaml file with general config information (results location, reference paths, etc)")
    parser.add_argument('-p', '--per_lib_input', required=True, help="CSV file with per-lib information")
    parser.add_argument('-r', '--results_dir', required=True, help="Results basepath to use in the config")
    parser.add_argument('-l', '--log_dir', help="Log directory to use in the config. Defaults to results_dir/logs")
    parser.add_argument('-t', '--temp_dir', required=True, help="Temporary directory basepath to use in the config")
    parser.add_argument('-s', '--simulate_single_lane', action='store_true', help="If sample fastq's don't contain lane information, treat them all as a single lane")
    parser.add_argument('--file_glob', help="Override default file glob of '*.fastq.gz'", default='*.fastq.gz')
    parser.add_argument('--capture_regex', help="Override default capture regex of '.*_L(\d+)_R(\d+).*\.fastq\.gz'", default='.*_L(\d+)_R(\d+).*\.fastq\.gz')

    args = parser.parse_args()

    config_dict = read_input(args.general_input)

    if not args.log_dir:
        log_dir = os.path.join(args.results_dir, "logs")
    else:
        log_dir = args.log_dir

    config_dict.update({'results' : args.results_dir, 'flux_log_dir' : log_dir, 'tmpdir' : args.temp_dir})

    per_lib = parse_per_lib(pd.read_csv(args.per_lib_input, dtype=str))

    config_dict.update(per_lib)

    config_dict = json.loads(json.dumps(config_dict)) #Standardize dict type throughout object by using json as intermediate

    yaml.dump(config_dict, sys.stdout, default_flow_style=False)
