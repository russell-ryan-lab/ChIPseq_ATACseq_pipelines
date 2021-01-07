#!/bin/bash

# divides reads from a single fastq file into multiple fastq files
# either based on desired number of files (method "files") or desired number of reads per file (method "reads")
# input files are expected to be named like sample2_L001_R1.fastq.gz
# output file names increment in the lane (e.g. L001) part of the file name
# use param 2 (OUT_ID) to change sample ID of output files

INFILE=$1     # e.g. sample2_L001_R1.fastq.gz
OUT_ID=$2     # e.g. sample2 -- the sample id/name for the outfile
OUT_LANE=$3   # e.g. 1       -- the starting lane number for the outfile, later formatted as L001
OUT_MATE=$4   # e.g. 1        -- 1 or 2, later formatted as R1/R2
METHOD=$5     # one of: files, reads
              # files: value is the number of files to create (divides reads evenly)
              # reads: value is the number of reads per file (creates unknown number of files)
VALUE=$6      # the value attached to specified method

if [[ "$1" == "" || "$2" == "" || "$3" == "" || "$4" == "" || "$5" == "" || "$6" == "" ]]
then
	echo "incorrect number of arguments."
	echo "usage: . split_fastq.sh [infile.fastq.gz] [sample] [lane] [mate] [method] [value]"
	echo "sample - a sample id, e.g. sample1"
	echo "lane   - starting lane number, e.g. 1"
	echo "mate   - mate number, either 1 or 2"
	echo "method - files or reads. files if specifying number of output files, reads if specifying number of reads per file"
	echo "value  - the integer value attached to the method"
	echo ""
	echo "  NOTE: leftover reads are added to last file."
else

	READ_COUNT=$(zcat $INFILE | gawk '(FNR % 4 ==1)' | wc -l)

	if [[ $METHOD == "files" ]]
	then
		DIVIDE_EACH=$(echo "$READ_COUNT" | gawk -v val="$VALUE" '{print int($0/val)}')
	else
		DIVIDE_EACH=$VALUE
	fi

	zcat $INFILE \
	| gawk \
		-v div=$DIVIDE_EACH \
		-v id=$OUT_ID \
		-v lane=$OUT_LANE \
		-v mate=$OUT_MATE \
		-v reads=$READ_COUNT \
		'BEGIN{
			readct=0
			reads_remaining=reads
			filect=lane-1
		}{
			if(FNR % 4 == 1){
				readct++
				reads_remaining--
			}
			if(FNR % 4 == 1 && readct % div == 1 && reads_remaining >= div){
				filect++
			}
                	filect_leading_zeros = sprintf("%03d", filect)
	                outfile = id "_L" filect_leading_zeros "_R" mate ".fastq"
			print $0 > outfile
		}'

fi
