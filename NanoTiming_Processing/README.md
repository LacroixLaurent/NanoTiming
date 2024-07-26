# Script to parse mod_mappings.bam file and compute NanoTiming  
***
### mod_mappings file are splitted into 50000 reads and processed in parallel on a our bioinformatic cluster.

Submission files provided are based on HTCondor submission architecture and should be adapted to other environment.  
The first script (*nanoT_prep*) split the mod_mappings.bam file into smaller bam files.  
The second script (*nanoT_parsing*) extract the information from the Ml and Mm fileds of the bam file. This script also compute the mean of the BrdU signal per 1kb bin. 
The third script (*nanoT_merging*) merge the per read 1kb binned signal into genomic 1kb bin.  
The *nanoT_bam2bed* script created a compressed bed file containing all the chrom, start and end of the reads in order to check the genomic coverage of the experiment.  
The *nanoT_cleaning* is used to clean data at the end of the procedure.  
The *nanoT_ParsingFunction.r* script contains the function used in the *nanoT_parsing* script.  

Examples of the output files are in the Zenodo repository under the DOI [10.5281/zenodo.12668295](https://doi.org/10.5281/zenodo.12668295)

To generate session_info file

```R
library(sessioninfo)  
session_info() %>% capture.output(file="nanotiming_figures_session_info.txt")  
``` 

***
