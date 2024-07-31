### R script to parse megalodon data on the cluster
### sub as list

source("nanoT_Parsingfunction.r")
`%+%` <- paste0
library(stringr)
library(dplyr)

args <- commandArgs(trailingOnly=T)
bamfile <- args[1]

NanoTimingDir <- "YourChoice"
Exp <- "WT_rep3"

path2work <- NanoTimingDir %+% "/" %+% Exp
path_in <- path2work %+% "/mega/"

path_out <- path2work %+% "/Parsed/"
system("mkdir -p " %+% path_out)

nc <- 12L

bamf <- bamfile
i <- basename(bamf) %>% str_remove(".bam") %>% str_remove("mod_splitted_")
outf <- paste0(path_out,Exp,"_meg_",i)
parsing_BrdUbin(bam.in=bamf,out.file=outf,ncores=nc,bin.size=1000,savefile=T)
