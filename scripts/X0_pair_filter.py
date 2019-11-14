#!/usr/bin/env python3

import argparse
import pysam
import sys


def count_mark_pair(thres, read1, read2):
	'''
	Counts (0, 1, or 2) failed reads, and adds custom tags to reads for type of failure
	inputs: thres - X0 threshold (int). r1 & r2 are read1 and read2 from a pair of alignments
	output: failed_count is count of reads which didn't meet threshold, r1_arr & r2_arr are as above, with addt'l tags as needed
	'''
	failed_count = 0 #Count of reads which didn't meet threshold

	if not read1.has_tag('X0'):
		read1.set_tag("X?", "NO_X0")
		failed_count +=1
	elif read1.get_tag('X0') > thres:
		read1.set_tag("X?", "HI_X0")
		failed_count +=1

	if not read2.has_tag('X0'):
		read2.set_tag("X?", "NO_X0")
		failed_count +=1
	elif read2.get_tag('X0') > thres:
		read2.set_tag("X?", "HI_X0")
		failed_count +=1

	return(failed_count, read1, read2)


def decide_fate(fail_count, mode):
	'''
	Decides whether or not to keep the pair
	inputs: fail_count - count of reads which didn't meet threshold. mode - either "one" or "both", describing the rule for which pairs to keep
	output: Boolean. True - Keep the pair. False - Discard the pair
	'''
	if mode == "both":
		if fail_count == 2:
			return(False)
		else: #0 or 1 failed
			return(True)
	elif mode == "one":
		if fail_count == 0:
			return(True)
		else: #1 or 2 failed
			return(False)


if __name__ == "__main__":

	parser = argparse.ArgumentParser(description="X0 threshold paired filter. Uses mapped, paired, name-sorted alignment file to process read pairs, and removes the pair if one or both (configurable) reads fall below the X0 threshold.")
	parser.add_argument("-b", "--bam_input", required=True, help="Name-sorted alignment file (SAM/BAM) input.")
	parser.add_argument("-o", "--bam_output", required=True, help="Output filename")
	parser.add_argument("-x", "--X0_threshold", type=int, default=1, required=True, help="X0 threshold, it will filter out all reads with X0 tag  greater than x0_threshold. x0_threshold of 1 will keep only uniquely mapped reads. Default = 1")
	parser.add_argument("-m", "--pair_mode", choices=["one", "both"], default="both", required=True, help='Paired filtering mode. Remove the pair if one or both reads fall below threshold. Default = both')

	args = parser.parse_args()

	#Initialize variables
	read1 = None
	name1 = None

	#Open files and perform error checking
	alignments = pysam.AlignmentFile(args.bam_input)

	if alignments.is_read == False:
		raise RuntimeError("Input is not readable")
	if alignments.header['HD']['SO'] != 'queryname':
		raise RuntimeError("Input is not name sorted")

	output = pysam.AlignmentFile(args.bam_output, mode="wb", template=alignments)

	#Iterate & process reads
	for entry in alignments:
		#First observed in putative pair
		if (name1 == None):
			read1 = entry
			name1 = entry.query_name
		#Second observed in putative pair
		else:
			read2 = entry
			name2 = entry.query_name

			#This is a pair
			if (name1 == name2):
				#counting and tagging failed reads as pairs
				failct, mod_read1, mod_read2 = count_mark_pair(args.X0_threshold, read1, read2)

				#Decide whether or not to keep the pair
				fate = decide_fate(failct, args.pair_mode)

				#As fate would have it, the pair is kept
				if fate:
					output.write(read1)
					output.write(read2)

				#Reset variables after the current pair has been processed
				read1 = None
				name1 = None

			#Putative pair was not a pair
			else:
				print("Names do not match. Skipping " + name1 + ". If this error occurs frequently, verify that the input is paired and name-sorted", file=sys.stderr)
				#Reset variables to set current read as read1
				read1 = entry
				name1 = entry.query_name
