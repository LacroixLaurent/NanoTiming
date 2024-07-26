### script to compute binned nanoT genomic profile

suppressMessages(library(tidyverse))
`%+%` <- paste0

NanoTimingDir <- "YourChoice"
Exp <- "WT_rep3"

path2work <- NanoTimingDir %+% "/" %+% Exp
path_in <- path2work %+% "/Parsed/"

file_list <- dir(path_in,full.names=T)

alldata <- do.call(bind_rows,lapply(file_list, function(x) {
		readRDS(x)})) %>%filter(chrom!="chrM_flye" & chrom!="chrM")

nanoT <- alldata %>%
		unnest(cols = c(signalr))%>%
		filter(signal>0.02)%>%
		group_by(chrom,positions)%>%
		summarise(mean_br_bin=mean(signal))
saveRDS(nanoT,file=paste0(path2work,"/",Exp,"_",RefGen,"_nanoT.rds"))
saveRDS(alldata,file=paste0(path2work,"/",Exp,"_",RefGen,"_nanoT_alldata.rds"))

q("no")
