## Fastq inputs

In the script config_creator.py, one of the inputs is a CSV file which has a row for each sample, and among other things, a column containing a directory path which contains the FASTQs for a given sample. This configuration offers great flexibility and control. Various scenarios are exemplified here where we can take various actions to prepare inputs so that they are pipeline-ready.

### A note on filename restrictions

One requirement for the input fastq files is that the read number must be discernable based on the filename (.e.g `_R1.fastq.gz`, `_1.fastq`, etc.), so that they can be appropriately placed as read 1 or read 2 in the configuration file. Since all three scenarios outlined here ultimately rely on the same config_creator script, they will all have this same restriction. [FIXME: recommendation for single-end reads? config_creator.py is already flexible enough to handle this with some creativity. Should we just define this and recommend it? Or make an argument/handler for it?]


### Fastq files already in separate directories for each sample

Example: Sequencing results from UMich Advanced Genomics Core

This is the simplest case, where no extra steps need to be performed. The directories of fastq files can be fed to config_creator.py via the sample info CSV as described [FIXME: is this still in the main readme?]. The config creator script will gather the paths to the input fastqs to all of a sample's read1 or read2 inputs. The resulting information in the configuration file will eventually be used in the pipeline to regard these inputs as a sample.

### Fastq files for all samples in a single directory, (somewhat) similarly named

Example: Sequencing runs downloaded from NCBI Sequence Read Archive

This is another common case, where the fastqs from multiple samples reside within a single directory. For this scenario, an additional step is taken to place the fastq files into separate subdirectories on a per-sample basis; we've created the prepare_fastq_inputs.py script to simplify this process. After the fastqs are placed in subdirectories, they can be fed to config_creator.py via the sample info CSV as described [FIXME: is this still in the main README?]
In this case, an additional requirement is that filenames must be similar enough that they can be parsed uniformly to identify files to be grouped together as a sample. With this script, a table is supplied which maps sample names to the unique portion of the fastq filenames (the portion that remains after uniformly parsing them). In case of SRA files, we can group BioSamples which are comprised of multiple sequencing runs (multiple SRR accessions). In a more general use-case, we could group samples that have fastqs performed over multiple lanes of sequencing.

Note: prepare_fastq_inputs.py hardlinks the files into the proper locations, instead of moving them. This preserves the original flat structure in the given fastq directory in case that is desired, but also introduces in that location the nested structure needed for config_creator.py. These links do not consume any extra storage space.

### Fastq files in heterogeneous directory structure or with filenames too disparate to parse uniformly

Example: Some samples from a collaborator, some samples internal, etc.

This last scenario should be uncommon, but with the inherent flexibility in this setup it is still not too difficult to support it. Once again, the read number information must be discernible from the filename [FIXME: single-end read caveat]. However, that's the only restriction. The key objective is to get each sample's fastq files into their own directory. If you are able to manually create directories and move the appropriate fastqs into those locations, then that is the extent of the work to be done. After that, config_creator.py can be run exactly as described above, by supplying the appropriate fastq directory for each sample in the CSV.
