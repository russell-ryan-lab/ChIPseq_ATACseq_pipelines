# Reference Data Sources

## blacklist

Data for versions 2 and 1 of ENCODE blacklists were downloaded from the Boyle Lab at the University of Michigan (https://github.com/Boyle-Lab/Blacklist/tree/master/lists).

## chromsize

Chromosome sizes were downloaded using the `fetchChromSizes` utility from the University of California Santa Cruz.

## tss

TSSs are computed using code from the [ataqv](https://github.com/ParkerLab/ataqv/tree/master/data/tss) repository.

### hg19
```
mysql --host=genome-mysql.cse.ucsc.edu --user=genome -D hg19 -e "SELECT * FROM refGene" > hg19.tss.refseq.tsv
Rscript ataqv/data/tss/make_tss.R --protein-coding --ucsc-refgene hg19.tss.refseq.tsv --out hg19.tss.refseq.bed
gzip hg19.tss.refseq.bed
```

### hg38
```
mysql --host=genome-mysql.cse.ucsc.edu --user=genome -D hg38 -e "SELECT * FROM refGene" > hg38.tss.refseq.tsv
Rscript ataqv/data/tss/make_tss.R --protein-coding --ucsc-refgene hg38.tss.refseq.tsv --out hg38.tss.refseq.bed
gzip hg38.tss.refseq.bed
```

# mm10
```
mysql --host=genome-mysql.cse.ucsc.edu --user=genome -D mm10 -e "SELECT * FROM refGene" > mm10.tss.refseq.tsv
Rscript ataqv/data/tss/make_tss.R --protein-coding --ucsc-refgene mm10.tss.refseq.tsv --out mm10.tss.refseq.bed
gzip mm10.tss.refseq.bed
```
