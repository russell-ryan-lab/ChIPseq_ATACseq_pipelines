#!/usr/bin/env python3

import argparse
from collections import defaultdict
import glob
import json
import os
import pandas as pd
import re


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
			#If samples do not contain lane information, use different capture regex which only captures read number
			libpaths_dict[lib] = basepath_to_filepathsdict(path, "*.fastq.gz", ".*_R(\d+).*\.fastq\.gz")
		else:
			libpaths_dict[lib] = basepath_to_filepathsdict(path, "*.fastq.gz", ".*_L(\d+)_R(\d+).*\.fastq\.gz")
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

	readgroups = Tree()

	for fq in all_fastqs:
		basename = os.path.basename(fq)
		rmatch = re.match(capture_regex, basename)
		if rmatch.group(0) == basename:
			#If samples do not contain lane information, treat them all as lane 001
			if args.simulate_single_lane:
				lane = "001"
				read = rmatch.group(1)
				if len(all_fastqs) != 2:
					raise ValueError("--simulate_single_lane flag can only be used if there are two fastq's per sample. There are " + str(len(all_fastqs)) + " in " + basepath)
			else:
				lane = rmatch.group(1)
				read = rmatch.group(2)
			#Add fastq to dict
			readgroups[lane][read] = fq

	return(readgroups)



if __name__ == '__main__':

	parser = argparse.ArgumentParser(prog='python ngs_rawdata_config_creator.py', description = "")
	parser.add_argument('-g', '--general_input', required=True, help="json file with general config information (results location, reference paths, etc)")
	parser.add_argument('-p', '--per_lib_input', required=True, help="CSV file with per-lib information")
	parser.add_argument('-r', '--results_dir', required=True, help="Results basepath to use in the config")
	parser.add_argument('-l', '--log_dir', help="Log directory to use in the config. Defaults to results_dir/logs")
	parser.add_argument('-t', '--temp_dir', required=True, help="Temporary directory basepath to use in the config")
	parser.add_argument('-s', '--simulate_single_lane', action='store_true', help="If sample fastq's don't contain lane information, treat them all as a single lane")


	args = parser.parse_args()

	with open(args.general_input) as general:
		config_dict = json.load(general)

	if not args.log_dir:
		log_dir = os.path.join(args.results_dir, "logs")
	else:
		log_dir = args.log_dir

	config_dict.update({'results' : args.results_dir, 'flux_log_dir' : log_dir, 'tmpdir' : args.temp_dir})

	per_lib = parse_per_lib(pd.read_csv(args.per_lib_input, dtype=str))

	config_dict.update(per_lib)

	print(json.dumps(config_dict, sort_keys=True, indent=4, separators=(',', " : ")))
