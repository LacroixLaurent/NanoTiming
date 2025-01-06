# Generating the refgen and the annotations  
***
### script to generate the final gff3 and the seqinf  

The genome assembly was generated as explained in [Theulot et al., 2024](https://doi.org/10.1038/s41467-024-55520-3). The original GFF file was produced by LRSDAY (v1.7).

Gene to name correspondence table was created from S288C_reference_genome_R64-3-1_20210421.gff file downloaded from [SGD](http://sgd-archive.yeastgenome.org/sequence/S288C_reference/genome_releases/S288C_reference_genome_R64-3-1_20210421.tgz). 

One of the X element mapped by LRSDAY on chrV was annotated as X_element_partial and we manually change this annotation to X_element in order to have at least one X element by telomere. 

Telomeric repeat were annotated by [telofinder](https://github.com/GillesFischerSorbonne/telofinder). 

***

ARS annotation was performed as the following:  
- ARS positions were downloaded from [OriDB](http://cerevisiae.oridb.org/) on December 2023 and renamed so that each ARS had a unique name.  
- corresponding DNA sequence were extracted from the [sacCer1 reference genome](https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE36045&format=file&file=GSE36045%5FsacCer1%5FOct2003%2Efa%2Egz) used for OriDB.
- the DNA sequences were mapped on our reference genome using bwa mem (v0.7.17-r1198-dirty) and the resulting bam file was converted to a bed file using bedtools bamtobed (v2.26.0).  
- ARS that map on a different chromosome between OriDB and our genome were discarded

rDNA sequences were extracted from the S288C reference genome (R64-3-1) and mapped on the BT1genome using bwa mem (v0.7.17-r1198-dirty) and the resulting bam file was converted to a bed file using bedtools bamtobed (v2.26.0). 

pBL-hsvTKco-hENT1co sequences were mapped on the BT1genome using bwa mem (v0.7.17-r1198-dirty) and the resulting bam file was converted to a bed file using bedtools bamtobed (v2.26.0)). 

All these features were combined in the BT1multiUra.gff3 file

The seqinfo file was generated using the following R code

```R
suppressMessages(library(GenomicRanges))  
library(Biostrings)  
BT1_genome <- readDNAStringSet("BT1_multiUra.fa.zip")  

seqinf <- Seqinfo(  
	seqnames=names(BT1_genome),  
	seqlengths=lengths(BT1_genome),  
	isCircular=c(rep(F,16),T),  
	genome="BT1multiUra")  
seqnames(seqinf)[17] <- "chrM"  

saveRDS(seqinf,file="seqinfBT1multiUra.rds")  
```
