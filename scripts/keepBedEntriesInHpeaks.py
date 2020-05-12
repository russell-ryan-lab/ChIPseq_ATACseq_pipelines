#! /usr/bin/env python3

import argparse
import os
import pysam
import sys


def keepBedEntriesInHpeaks(hpeaks_in_fn, hpeaks_out_fn, bed_fn):
  ''' this function takes as an input a bed file and a hpeaks file,
  both with the same col4 (id) and outputs a hpeaks file WITH matching
  entries in the bed file '''
  hpeaks_in = open(hpeaks_in_fn, 'r')
  hpeaks_out = open(hpeaks_out_fn, 'w')
  bed = open(bed_fn, 'r')
  bed_dict = {}
  nlines_bed = 0
  for line in bed:
    cols = line.strip().split('\t')
    uniq_id = cols[3]
    bed_dict[uniq_id] = line
    nlines_bed+=1
  bed.close()
  line_count = 1
  for line in hpeaks_in:
    if line[0] == '#':
      if line_count == 4:
        hpeaks_out.write('# peaks remaining after blacklist filtering %s \n' % (nlines_bed) )
      else:
        hpeaks_out.write(line)
    cols = line.strip().split('\t')
    uniq_id = cols[0]
    if uniq_id in bed_dict:
      hpeaks_out.write(line)
    line_count+=1
  hpeaks_in.close()
  return 0


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--hpeaks_in_fn', help='hpeaks_in_fn')
    parser.add_argument('-o', '--hpeaks_out_fn', help='hpeaks_out_fn')
    parser.add_argument('-b', '--bed_fn', help='bed_fn')
    args = parser.parse_args()
    keepBedEntriesInHpeaks(args.hpeaks_in_fn, args.hpeaks_out_fn, args.bed_fn)
