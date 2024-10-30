#### Supplementary Figure 23
### Noise analysis
### using i vs i+1 dif with and without intraread correlation broken

suppressMessages(library(rtracklayer))
suppressMessages(library(GenomicRanges))
suppressMessages(library(tidyverse))
library(patchwork)
theme_set(theme_bw())
mypal <- c(paletteer::paletteer_d("ggthemes::Classic_20"),"grey40","black","gold","greenyellow")

`%+%` <- paste0

setwd("/Users/ll/work/RStudioProjects/NanoTiming")
path2fig <- "Figures_Article/Figures_pdf/"
path2data <- "Figures_Article/Figures_data/"

source("Figures_Article/rescaling_function.r")


## load genomic informations
chrom_order <- c("chrI","chrII","chrIII","chrIV","chrV","chrVI","chrVII","chrVIII","chrIX","chrX","chrXI","chrXII","chrXIII","chrXIV","chrXV","chrXVI")
seqinf <- readRDS("Reference_Genome/seqinfBT1multiUra.rds")
seqlevels(seqinf)=seqlevels(seqinf)[1:16]

chrom_sizes <- as.data.frame(seqinf)
chrom_sizes$chrom <- rownames(chrom_sizes)
chrom_sizes <- as_tibble(chrom_sizes) %>%
	select(chrom,seqlengths) %>%
	mutate(chrom=factor(chrom,levels=chrom_order))
genome_size <- chrom_sizes %>% pull(seqlengths) %>% sum

### create a background to fill empty bins with NA
bs=1000
bingen0 <- tileGenome(seqinf,tilewidth=bs, cut.last.tile.in.chrom=T)

backgr <- as_tibble(bingen0) %>%
	#dplyr::rename(positions=start,chrom=seqnames)%>%
	mutate(y=NA) %>%
	select(seqnames,start,end,y)


# Import precomputed data
### MFA-seq
toplot_mfa <- readRDS(paste0(path2data,"data_SupFig23_mfa.rds")) %>% mutate(type="MFA-seq")

### sort-seq
toplot_sort <- readRDS(paste0(path2data,"data_SupFig23_sortseq.rds")) %>% mutate(type="sort-seq")

### nanoT
nanoT_noise <- readRDS(paste0(path2data,"data_SupFig23_nanoT.rds")) %>% mutate(type="Nanotiming")

### nanoT_uncor
nanoT_noise_uncor <- readRDS(paste0(path2data,"data_SupFig23_nanoT_uncor.rds")) %>% mutate(type="Nanotiming-uncor")

### nanoT_shuffled
nanoTrandom_noise <- readRDS(paste0(path2data,"data_SupFig23_nanoT_shuf.rds")) %>% mutate(type="Nanotiming-shuffled")


### plot
sort_med <- toplot_sort %>% group_by(mapGb,type) %>% summarize(med=median(noisevar,na.rm=T),nb=median(nbin,na.rm=T)) %>% ungroup
mfa_med <- toplot_mfa %>% group_by(mapGb,type) %>% summarize(med=median(noisevar,na.rm=T),nb=median(nbin,na.rm=T)) %>% ungroup
nanoT_med <- nanoT_noise %>% group_by(mapGb,type) %>% summarize(med=median(noisevar,na.rm=T),nb=median(nbin,na.rm=T)) %>% ungroup
nanoT_med_uncor <- nanoT_noise_uncor %>% group_by(mapGb,type) %>% summarize(med=median(noisevar,na.rm=T),nb=median(nbin,na.rm=T)) %>% ungroup
nanoTrand_med <- nanoTrandom_noise %>% group_by(mapGb,type) %>% summarize(med=median(noisevar,na.rm=T),nb=median(nbin,na.rm=T)) %>% ungroup

boxline <- 0.2
p1<- ggplot(nanoT_noise)+
	geom_line(data=nanoT_med,aes(x=mapGb,y=med,col=type))+
	geom_boxplot(aes(x=mapGb,y=noisevar,group=mapGb,col=type),outlier.shape = NA,show.legend=F,size=boxline)+
	geom_line(data=nanoT_med_uncor,aes(x=mapGb,y=med,col=type))+
	geom_boxplot(data=nanoT_noise_uncor,aes(x=mapGb,y=noisevar,group=mapGb,col=type),outlier.shape = NA,show.legend=F,size=boxline)+
	geom_line(data=sort_med,aes(x=mapGb,y=med,col=type))+
	geom_boxplot(data=toplot_sort,aes(x=mapGb,y=noisevar,group=mapGb,col=type),outlier.shape = NA,show.legend=F,size=boxline)+
	geom_line(data=mfa_med,aes(x=mapGb,y=med,col=type))+
	geom_boxplot(data=toplot_mfa,aes(x=mapGb,y=noisevar,group=mapGb,col=type),outlier.shape = NA,show.legend=F,size=boxline)+
	geom_line(data=nanoTrand_med,aes(x=mapGb,y=med,col=type))+
	geom_boxplot(data=nanoTrandom_noise,aes(x=mapGb,y=noisevar,group=mapGb,col=type),outlier.shape = NA,show.legend=F,size=boxline)+
	scale_x_log10()+
	scale_y_log10()+
	xlab("Number of mapped bases")+
	ylab("Noise estimator")
p2<- ggplot(nanoT_noise)+
	geom_line(data=nanoT_med,aes(x=mapGb,y=nb/123.60,col=type))+
	geom_boxplot(aes(x=mapGb,y=nbin/123.60,group=mapGb,col=type),outlier.shape = NA,show.legend=F,size=boxline)+
	geom_line(data=nanoT_med_uncor,aes(x=mapGb,y=nb/123.60,col=type))+
	geom_boxplot(data=nanoT_noise_uncor,aes(x=mapGb,y=nbin/123.60,group=mapGb,col=type),outlier.shape = NA,show.legend=F,size=boxline)+
	geom_line(data=sort_med,aes(x=mapGb,y=nb/123.60,col=type))+
	geom_boxplot(data=toplot_sort,aes(x=mapGb,y=nbin/123.60,group=mapGb,col=type),outlier.shape = NA,show.legend=F,size=boxline)+
	geom_line(data=mfa_med,aes(x=mapGb,y=nb/123.60,col=type))+
	geom_boxplot(data=toplot_mfa,aes(x=mapGb,y=nbin/123.60,group=mapGb,col=type),outlier.shape = NA,show.legend=F,size=boxline)+
	geom_line(data=nanoTrand_med,aes(x=mapGb,y=nb/123.60,col=type))+
	geom_boxplot(data=nanoTrandom_noise,aes(x=mapGb,y=nbin/123.60,group=mapGb,col=type),outlier.shape = NA,show.legend=F,size=boxline)+
	scale_x_log10()+
	coord_cartesian(ylim=c(0,100))+
	xlab("Number of mapped bases")+
	ylab("Genome covered (%)")

p1/p2 + plot_layout(guides="collect")& scale_color_manual("",breaks=c("MFA-seq","sort-seq","Nanotiming","Nanotiming-uncor","Nanotiming-shuffled"),values=mypal[c(5,1,7,3,13)])& plot_annotation(title="Supplementary Figure 23",tag_levels="a")
ggsave(paste0(path2fig,"SupFig23.pdf"),h=7,w=6)




toreport <- bind_rows(sort_med,mfa_med,nanoT_med,nanoT_med_uncor,nanoTrand_med)
toreport %>% filter(mapGb==1e9)
#       mapGb type             med     nb
# 1 1000000000 sort-seq            0.0111  11308.
# 2 1000000000 MFA-seq             0.0377  11292
# 3 1000000000 Nanotiming          0.00293 12358
# 4 1000000000 Nanotiming-uncor    0.0126  12358
# 5 1000000000 Nanotiming-shuffled 0.129   12358

