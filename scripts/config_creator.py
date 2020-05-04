#!/usr/bin/env python3

import argparse
from collections import defaultdict
import glob
import json
import jsonschema
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

class _InputValidationError(Exception):
    """Raised for problematic input data."""

    def __init__(self, msg, *args):
        super(_InputValidationError, self).__init__(msg, *args)


def parse_per_lib(sample_table):
    """
    Takes a pandas df as input, and returns a dict of dicts reflecting the per-library data
    """

    validate_sample_table(sample_table)

    per_lib_dict = dict()

    lib_basepaths = sample_table[['sample', 'basepath']].values

    per_lib_dict['sample_paths'] = assign_libpaths(lib_basepaths)

    othercols = sample_table.columns.drop(['lib', 'basepath', 'sample'])
    for c in othercols:
        two_cols = ['sample', c]
        combined_name = "_".join(two_cols)
        col_dict = dict(sample_table[two_cols].dropna().values)
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
        sample, path = row
        libpaths_dict[sample] = basepath_to_filepathsdict(path, args.file_glob, args.capture_regex)
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

    readgroups = defaultdict(list)

    for fq in all_fastqs:
        basename = os.path.basename(fq)
        rmatch = re.match(capture_regex, basename)
        if not rmatch:
            msg_fmt = ("\nFile {} did not match regular expression {}. "
                "Could not capture desired group(s).\n"
                "Note that all input files should be formatted consistently. "
                "File glob and capture regex can be controlled with --file_glob and --capture_regex if desired.")
            raise RuntimeError(msg_fmt.format(basename, capture_regex))
        if rmatch.group(0) == basename:
            readnum = rmatch.group(1)
            #Add fastq to dict
            readgroups[readnum].append(fq)

    for key in readgroups:
        readgroups[key].sort()

    return(readgroups)

def read_input(input_filename):
    """Takes input filename as argument, chooses to load json or yaml based on file extension,
    and handles any errors that occur as the input file is loaded. Returns the loaded config dictionary.
    """
    if input_filename.endswith('.yaml') or input_filename.endswith('.yml'):
        try:
            with open(input_filename) as infile:
                config_dict = yaml.load(infile, Loader=yaml.SafeLoader)
        except (yaml.parser.ParserError, yaml.scanner.ScannerError) as e:
            print(e)
            msg = "Error: Error loading YAML. Assuming YAML input based on file extension."
            sys.exit(str(msg))

    elif input_filename.endswith('.json'):
        try:
            with open(input_filename) as infile:
                config_dict = json.load(infile)
        except json.decoder.JSONDecodeError as e:
            print(e)
            msg = "Error: Error loading JSON. Assuming JSON input based on file extension."
            sys.exit(str(msg))

    else:
        msg = "Error: Cannot assume YAML or JSON based on file extension. Filename is {}.".format(input_filename)
        sys.exit(str(msg))

    return config_dict

def validate_config_with_schema(config_dict, schema_filename):
    """Takes config dictionary and schema filename as arguments. Loads the schema file (assumes YAML format),
    and uses jsonschema.validate to compare the config dictionary to the schema. Handles any
    errors thrown by the validator.
    """
    try:
        pipeline_basedir = os.path.realpath(os.path.dirname(os.path.dirname(__file__)))
        schema_filename = os.path.join(pipeline_basedir, 'tests', 'pipeline_config_schema.yaml')
        with open(schema_filename) as infile:
            schema_dict = yaml.load(infile, Loader=yaml.SafeLoader)
    except:
        msg = "Error loading schema file {}. Expecting YAML.".format(schema_filename)

    try:
        jsonschema.validate(config_dict, schema_dict)
    except Exception as e:
        #print(e)
        msg = "Error validating config against schema:\nSchema file: {}\nReason: {}".format(schema_filename, e.message)
        sys.exit(str(msg))

def validate_sample_table(sample_table):
    """Takes pandas dataframe as input (loaded from per_lib_input), verifies that required colums are present,
    and that samplenames don't contain invalid characters. Only alphanumeric and _ are allowed.
    """
    required_cols = set(['lib', 'sample', 'basepath'])
    actual_cols = set(sample_table.columns)
    missing_cols = required_cols - actual_cols
    if missing_cols:
        raise RuntimeError("Missing required {} column(s)".format(missing_cols))
    # Error checking for samplenames - can only contain alphanumeric and _
    samplenames = sample_table['sample'].to_list()
    goodname_matches = [re.match(r'^[a-zA-Z0-9_]+$', x) for x in samplenames]
    problematic_samples = []
    for i, name in enumerate(goodname_matches):
        if name == None:
            problematic_samples.append(samplenames[i])
    if problematic_samples:
        raise RuntimeError("Samplenames contain invalid characters: {}".format(problematic_samples))


if __name__ == '__main__':

    parser = argparse.ArgumentParser(prog='python ngs_rawdata_config_creator.py', description = "")
    parser.add_argument('-g', '--general_input', required=True, help="json or yaml file with general config information (results location, reference paths, etc)")
    parser.add_argument('-p', '--per_lib_input', required=True, help="CSV file with per-lib information")
    parser.add_argument('-r', '--results_dir', required=True, help="Results basepath to use in the config")
    parser.add_argument('-t', '--temp_dir', required=True, help="Temporary directory basepath to use in the config")
    parser.add_argument('--file_glob', help="Override default file glob of '*.fastq.gz'", default='*.fastq.gz')
    parser.add_argument('--capture_regex', help="Override default capture regex of '.*_R(\d+).*\.fastq\.gz'", default='.*_R(\d+).*\.fastq\.gz')

    args = parser.parse_args()

    config_dict = read_input(args.general_input)

    config_dict.update({'results_dir' : args.results_dir, 'tmpdir' : args.temp_dir})

    per_lib = parse_per_lib(pd.read_csv(args.per_lib_input, dtype=str))

    config_dict.update(per_lib)

    config_dict = json.loads(json.dumps(config_dict)) #Standardize dict type throughout object by using json as intermediate

    pipeline_basedir = os.path.realpath(os.path.dirname(os.path.dirname(__file__)))
    schema_filename = os.path.join(pipeline_basedir, 'tests', 'pipeline_config_schema.yaml')
    validate_config_with_schema(config_dict, schema_filename)

    yaml.dump(config_dict, sys.stdout, default_flow_style=False)
