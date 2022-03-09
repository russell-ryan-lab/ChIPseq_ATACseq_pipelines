# v1.1.0 (03/08/2022)

- Tweaked SE and PE ChIP-Seq pipelines to allow for peak calling and motif analysis on samples without an Input control (e.g. public data without controls). No input peak calling can be accessed by supplying the sample's name as its own input in the sample sheet.
- Updated README to reflect new features.
- Fixed outdated/missing information in README.
- Iterated version to 1.1.0

# v1.0.1 (01/07/2021)

- updated igvtools version from 2.3.93 to 2.5.3 (igvtools.yaml)
- bugfix: file naming issue in sampe rule that caused incorrect output from ATACV (alignment_bwa_aln_pe.smk)
- README reflects Python 3.7 requirement for conda environment. This fix is only relevant to developers.  Python 3.8 caused errors in nosetests.

# v1.0.0 (09/04/2020)

Initial release of combined ChIP-seq and ATAC-seq pipeline.

The ChIPseq pipeline has two portions. The first portion is similar to the ATAC-seq pipeline, in that it performs alignment, filtering, and produces bigwig files (though the tools and procedures sometimes differ from ATACseq). The second portion takes filtered bam files, and performs peak-calling and motif enrichment steps. These two portions can be run separately. See below for details.

The ATACseq pipeline takes raw fastq files as input, performs alignment, filtering, and QC, calls peaks using MACS2, and produces bigwig files for viewing.
