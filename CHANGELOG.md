# v1.0.0 (09/04/2020)

Initial release of combined ChIP-seq and ATAC-seq pipeline.

The ChIPseq pipeline has two portions. The first portion is similar to the ATAC-seq pipeline, in that it performs alignment, filtering, and produces bigwig files (though the tools and procedures sometimes differ from ATACseq). The second portion takes filtered bam files, and performs peak-calling and motif enrichment steps. These two portions can be run separately. See below for details.

The ATACseq pipeline takes raw fastq files as input, performs alignment, filtering, and QC, calls peaks using MACS2, and produces bigwig files for viewing.
