### Supplementary Figure 26
### comparison nanotrace at 5,10 and 20µM

suppressMessages(library(GenomicRanges))
suppressMessages(library(rtracklayer))
suppressMessages(library(tidyverse))
library(patchwork)
library(ggprism)
library(ggcorrplot)
theme_set(theme_bw())
mypal <- c(paletteer::paletteer_d("ggthemes::Classic_20"),"grey40","black","gold","greenyellow")
mypal1 <- mypal[c(1,7)]
mypal2 <- mypal[c(2,19,1,7,8,3,4,13,14)]
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

fullGFF <- import("Reference_Genome/BT1multiUra.gff3") %>% as_tibble()
ARS <- fullGFF %>% filter(Conf=="Confirmed",name!="ARS1216.5",type %in% c("ORI"))
CEN <- fullGFF %>% filter(type=="centromere")
rDNA <- fullGFF %>% filter(type=="rRNA")

### loading inputs

### import nanot data
nanot_wt <- bind_rows(
	readRDS(path2data %+% "nanoT_WT_rep1.rds") %>% mutate(filename="wt_rep1_5")%>%ungroup() %>%
		mutate(mod=1+myscaling0(mean_br_bin,infq=0.005,supq=0.995)),
	readRDS(path2data %+% "nanoT_WT_10.rds") %>% mutate(filename="wt_rep1_10")%>%ungroup() %>%
		mutate(mod=1+myscaling0(mean_br_bin,infq=0.005,supq=0.995)),
	readRDS(path2data %+% "nanoT_WT_20.rds") %>% mutate(filename="wt_rep1_20")%>%ungroup() %>%
		mutate(mod=1+myscaling0(mean_br_bin,infq=0.005,supq=0.995)))

nanot2 <- do.call(bind_rows,nanot_wt) %>% mutate(chrom=factor(chrom,levels=c("chrI","chrII","chrIII","chrIV","chrV","chrVI","chrVII","chrVIII","chrIX","chrX","chrXI","chrXII","chrXIII","chrXIV","chrXV","chrXVI"))) %>% ungroup %>% mutate(positions=positions+500)
#### add a capping to chrom size

nanot3 <- left_join(nanot2,chrom_sizes, by="chrom") %>% mutate(positions=map2_dbl(positions,seqlengths, function(x,y) {res=ifelse(x<y,x,y); return(res)}))

### reorder filename

input <- nanot3%>%
	mutate(filename=factor(filename,levels=c("wt_rep1_5","wt_rep1_10","wt_rep1_20"))) %>% select(-c(seqlengths,mean_br_bin,BrdU))

SortSeqWT <- readRDS(path2data %+% "sortseq_WT.rds") %>%
	mutate(mod=1+myscaling0(timing,infq=0.005,supq=0.995))%>%
	mutate(filename="sort-seq") %>%
	mutate(mut="WT") %>%
	mutate(positions=positions+500) %>%
	left_join(.,chrom_sizes, by="chrom") %>%
	mutate(positions=map2_dbl(positions,seqlengths, function(x,y) {res=ifelse(x<y,x,y); return(res)})) %>%
	mutate(chrom=factor(chrom,levels=chrom_order)) %>% select(-seqlengths)

ars2plot <- ARS %>% rename(chrom=seqnames)
cen2plot <- CEN %>% rename(chrom=seqnames)
rdna2plot <- rDNA %>% rename(chrom=seqnames) %>% group_by(chrom,type) %>% summarise(start=min(start),end=max(end))

mypal3 <- mypal[c(1,2,3,22)]

input2 <- bind_rows(input,SortSeqWT) %>% mutate(filename=factor(filename,levels=c("wt_rep1_5","wt_rep1_10","wt_rep1_20","sort-seq")))
miny=0.8
maxy=2.2

pl <- ggplot(input2)+
	geom_vline(data=ars2plot,aes(xintercept=(start+end)/2),col=mypal[9],show.legend=F,linewidth=0.2)+
	geom_vline(data=cen2plot,aes(xintercept=(start+end)/2),col=mypal[5],show.legend=F)+
	geom_rect(data=rdna2plot,aes(xmin=start,xmax=end,ymin=0.9,ymax=1),col=NA,fill=mypal[21],show.legend=F,alpha=0.6)+
	geom_point(aes(x=positions,y=mod,col=filename,group=filename),shape=16,size=0.2)+
	scale_color_manual("",values=mypal3)+
	ylab("RT")+
	xlab("Genomic position (kb)")+
	coord_cartesian(ylim=c(miny,maxy),expand=F)+
	scale_x_continuous(
		guide = "prism_minor",
		labels=scales::unit_format(big.mark ="",suffix="",scale=1e-3,sep=""),
		limits=c(1,seqlengths(seqinf)[4]),
		breaks=seq(0,seqlengths(seqinf)[4],100000),
		minor_breaks=seq(0,seqlengths(seqinf)[4],50000),
		expand=c(0,0))+
	facet_wrap(~chrom,ncol=1,scales="fixed",strip.position = "right")+
	theme(legend.position="top")+
	guides(color = guide_legend(override.aes = list(size = 2)))+
	ggtitle("Supplementary Figure 26")

quartz(file=paste0(path2fig,"SupFig26.pdf"),height=14,width=12,type="pdf")
pl
dev.off()

