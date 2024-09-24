#### Figure 2
### comparison nanotrace vs sortseq on chrXII

suppressMessages(library(GenomicRanges))
suppressMessages(library(rtracklayer))
suppressMessages(library(tidyverse))
library(patchwork)
library(ggprism)
theme_set(theme_bw())
mypal <- c(paletteer::paletteer_d("ggthemes::Classic_20"),"grey40","black","gold","greenyellow")
mypal1 <- mypal[c(1,7,15)]
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
file2load <- c(
		"WT_rep1",
		"WT_rep2",
		"WT_rep3",
		"WT_rep4",
		"WT_rep5",
		"WT_rep6",
		"FKH1_rep1",
		"FKH1_rep2",
		"FKH1_rep3",
		"CTF19_rep1",
		"CTF19_rep2",
		"CTF19_rep3",
		"RIF1_rep1",
		"RIF1_rep2",
		"RIF1_rep3",
		"KU70_rep1",
		"KU70_rep2",
		"KU70_rep3"
		)

file_order <- c(
		"FKH1_rep1",
		"FKH1_rep2",
		"FKH1_rep3",
		"CTF19_rep1",
		"CTF19_rep2",
		"CTF19_rep3",
		"RIF1_rep1",
		"RIF1_rep2",
		"RIF1_rep3",
		"KU70_rep1",
		"KU70_rep2",
		"KU70_rep3",
		"WT_rep1",
		"WT_rep2",
		"WT_rep3",
		"WT_rep4",
		"WT_rep5",
		"WT_rep6"
	)

nanot <- lapply(file2load, function(x)
	{
	res <- readRDS(path2data %+% "nanoT_" %+% x %+% ".rds") %>%
		ungroup %>%
		mutate(mod=1+myscaling0(mean_br_bin,infq=0.005,supq=0.995))
	return(res)
	}
	)
nanot[[1]] <- nanot[[1]] %>% select(-BrdU)
nanot2 <- do.call(bind_rows,nanot) %>% mutate(chrom=factor(chrom,levels=chrom_order)) %>% ungroup %>% mutate(positions=positions+500)
#### add a capping to chrom size
nanot3 <- left_join(nanot2,chrom_sizes, by="chrom") %>% mutate(positions=map2_dbl(positions,seqlengths, function(x,y) {res=ifelse(x<y,x,y); return(res)})) %>% select(-seqlengths)

### Import sortseq data with scaling (and MFAseq also)
SortSeqWT <- readRDS(path2data %+% "sortseq_WT.rds") %>%
	mutate(mod=1+myscaling0(timing,infq=0.005,supq=0.995))%>%
	mutate(filename="sort-seq") %>%
	mutate(mut="WT") %>%
	mutate(positions=positions+500) %>%
	left_join(.,chrom_sizes, by="chrom") %>%
	mutate(positions=map2_dbl(positions,seqlengths, function(x,y) {res=ifelse(x<y,x,y); return(res)})) %>%
	mutate(chrom=factor(chrom,levels=chrom_order)) %>% select(-seqlengths)


MFAseq <- readRDS(path2data %+% "MFAseq_WT.rds") %>%
	mutate(mod=1+myscaling0(ratio,infq=0.005,supq=0.995)) %>%
	mutate(filename="MFA-seq") %>%
	mutate(mut="WT")  %>%
	mutate(positions=positions+500) %>%
	left_join(.,chrom_sizes, by="chrom") %>%
	mutate(positions=map2_dbl(positions,seqlengths, function(x,y) {res=ifelse(x<y,x,y); return(res)})) %>%
	mutate(chrom=factor(chrom,levels=chrom_order)) %>% select(-seqlengths)

### recode filename
input <- nanot3 %>%
	mutate(filename=factor(filename,levels=file_order)) %>% mutate(filename=fct_recode(filename,"fkh1∆_rep1"="FKH1_rep1","fkh1∆_rep2"="FKH1_rep2","fkh1∆_rep3"="FKH1_rep3","ctf19∆_rep1"="CTF19_rep1","ctf19∆_rep2"="CTF19_rep2","ctf19∆_rep3"="CTF19_rep3","rif1∆_rep1"="RIF1_rep1","rif1∆_rep3"="RIF1_rep3","rif1∆_rep2"="RIF1_rep2","yku70∆_rep1"="KU70_rep1","yku70∆_rep2"="KU70_rep2","yku70∆_rep3"="KU70_rep3",wt_rep1="WT_rep1",wt_rep2="WT_rep2",wt_rep3="WT_rep3",wt_rep4="WT_rep4",wt_rep5="WT_rep5",wt_rep6="WT_rep6"))

### ChrXII nanotraces
chr2plot <-"chrXII"

ars2plot <- ARS %>% rename(chrom=seqnames) %>% filter(chrom==chr2plot)
cen2plot <- CEN %>% rename(chrom=seqnames) %>% filter(chrom==chr2plot)
rdna2plot <- rDNA %>% rename(chrom=seqnames) %>% filter(chrom==chr2plot) %>% group_by(chrom,type) %>% summarise(start=min(start),end=max(end))

input2 <- bind_rows(input %>% filter(filename=="wt_rep1"),SortSeqWT)
miny=0.8
maxy=2.2

pl1 <- ggplot(input2 %>% filter(chrom==chr2plot))+
	geom_vline(data=ars2plot,aes(xintercept=(start+end)/2),col=mypal[9],show.legend=F,linewidth=0.2)+
	geom_vline(data=cen2plot,aes(xintercept=(start+end)/2),col=mypal[5],show.legend=F)+
	geom_rect(data=rdna2plot,aes(xmin=start,xmax=end,ymin=0.9,ymax=1),col=NA,fill=mypal[21],show.legend=F,alpha=0.6)+
	geom_point(aes(x=positions,y=mod,col=filename,group=filename),shape=16,size=0.2)+
	scale_color_manual("",values=mypal1,labels = c("Relative copy number by sort-seq", "Mean BrdU content"))+
	ylab("RT")+
	xlab("Genomic position on chrXII (kb)")+
	coord_cartesian(ylim=c(miny,maxy),expand=F)+
	scale_x_continuous(
		guide = "prism_minor",
		labels=scales::unit_format(big.mark ="",suffix="",scale=1e-3,sep=""),
		limits=c(1,seqlengths(seqinf)[12]),
		breaks=seq(0,seqlengths(seqinf)[12],100000),
		minor_breaks=seq(0,seqlengths(seqinf)[12],50000),
		expand=c(0,0))+
	theme(legend.position="top",axis.title.x=element_blank())+
	guides(color = guide_legend(override.aes = list(size = 2)))+
	annotate("text",x=10000,y=2.1,label="Early",hjust=0,fontface="italic",size=3)+
	annotate("text",x=10000,y=0.9,label="Late",hjust=0,fontface="italic",size=3)



input2 <- bind_rows(input %>% filter(filename=="wt_rep1"),MFAseq)
miny=0.8
maxy=2.2

pl2 <- ggplot(input2 %>% filter(chrom==chr2plot))+
	geom_vline(data=ars2plot,aes(xintercept=(start+end)/2),col=mypal[9],show.legend=F,linewidth=0.2)+
	geom_vline(data=cen2plot,aes(xintercept=(start+end)/2),col=mypal[5],show.legend=F)+
	geom_rect(data=rdna2plot,aes(xmin=start,xmax=end,ymin=0.9,ymax=1),col=NA,fill=mypal[21],show.legend=F,alpha=0.6)+
	geom_point(aes(x=positions,y=mod,col=filename,group=filename),shape=16,size=0.2)+
	scale_color_manual("",values=mypal1,labels = c("Relative copy number by MFA-seq", "Mean BrdU content"))+
	ylab("RT")+
	xlab("Genomic position on chrXII (kb)")+
	coord_cartesian(ylim=c(miny,maxy),expand=F)+
	scale_x_continuous(
		guide = "prism_minor",
		labels=scales::unit_format(big.mark ="",suffix="",scale=1e-3,sep=""),
		limits=c(1,seqlengths(seqinf)[12]),
		breaks=seq(0,seqlengths(seqinf)[12],100000),
		minor_breaks=seq(0,seqlengths(seqinf)[12],50000),
		expand=c(0,0))+
	theme(legend.position="top",axis.title.x=element_blank())+
	guides(color = guide_legend(override.aes = list(size = 2)))

#### MUT vs WT

mut2plot <- c("WT","CTF19")

input2 <- bind_rows(input %>% filter(mut %in% mut2plot))
miny=0.8
maxy=2.2

pl3 <- ggplot(input2 %>% filter(chrom==chr2plot))+
	geom_vline(data=ars2plot,aes(xintercept=(start+end)/2),col=mypal[9],show.legend=F,linewidth=0.2)+
	geom_vline(data=cen2plot,aes(xintercept=(start+end)/2),col=mypal[5],show.legend=F)+
	geom_rect(data=rdna2plot,aes(xmin=start,xmax=end,ymin=0.9,ymax=1),col=NA,fill=mypal[21],show.legend=F,alpha=0.6)+
	geom_point(aes(x=positions,y=mod,col=filename,group=filename),shape=16,size=0.2)+
	scale_color_manual("",values=mypal2)+
	ylab("RT")+
	xlab("Genomic position on chrXII (kb)")+
	coord_cartesian(ylim=c(miny,maxy),expand=F)+
	scale_x_continuous(
		guide = "prism_minor",
		labels=scales::unit_format(big.mark ="",suffix="",scale=1e-3,sep=""),
		limits=c(1,seqlengths(seqinf)[12]),
		breaks=seq(0,seqlengths(seqinf)[12],100000),
		minor_breaks=seq(0,seqlengths(seqinf)[12],50000),
		expand=c(0,0))+
	theme(legend.position="top",axis.title.x=element_blank())+
	guides(color = guide_legend(nrow=1,override.aes = list(size = 2)))

mut2plot <- c("WT","RIF1")

input2 <- bind_rows(input %>% filter(mut %in% mut2plot))
miny=0.8
maxy=2.2

pl4 <- ggplot(input2 %>% filter(chrom==chr2plot))+
	geom_vline(data=ars2plot,aes(xintercept=(start+end)/2),col=mypal[9],show.legend=F,linewidth=0.2)+
	geom_vline(data=cen2plot,aes(xintercept=(start+end)/2),col=mypal[5],show.legend=F)+
	geom_rect(data=rdna2plot,aes(xmin=start,xmax=end,ymin=0.9,ymax=1),col=NA,fill=mypal[21],show.legend=F,alpha=0.6)+
	geom_point(aes(x=positions,y=mod,col=filename,group=filename),shape=16,size=0.2)+
	scale_color_manual("",values=mypal2)+
	ylab("RT")+
	xlab("Genomic position on chrXII (kb)")+
	coord_cartesian(ylim=c(miny,maxy),expand=F)+
	scale_x_continuous(
		guide = "prism_minor",
		labels=scales::unit_format(big.mark ="",suffix="",scale=1e-3,sep=""),
		limits=c(1,seqlengths(seqinf)[12]),
		breaks=seq(0,seqlengths(seqinf)[12],100000),
		minor_breaks=seq(0,seqlengths(seqinf)[12],50000),
		expand=c(0,0))+
	theme(legend.position="top",axis.title.x=element_blank())+
	guides(color = guide_legend(nrow=1,override.aes = list(size = 2)))

mut2plot <- c("WT","KU70")

input2 <- bind_rows(input %>% filter(mut %in% mut2plot))
miny=0.8
maxy=2.2

pl5 <- ggplot(input2 %>% filter(chrom==chr2plot))+
	geom_vline(data=ars2plot,aes(xintercept=(start+end)/2),col=mypal[9],show.legend=F,linewidth=0.2)+
	geom_vline(data=cen2plot,aes(xintercept=(start+end)/2),col=mypal[5],show.legend=F)+
	geom_rect(data=rdna2plot,aes(xmin=start,xmax=end,ymin=0.9,ymax=1),col=NA,fill=mypal[21],show.legend=F,alpha=0.6)+
	geom_point(aes(x=positions,y=mod,col=filename,group=filename),shape=16,size=0.2)+
	scale_color_manual("",values=mypal2)+
	ylab("RT")+
	xlab("Genomic position on chrXII (kb)")+
	coord_cartesian(ylim=c(miny,maxy),expand=F)+
	scale_x_continuous(
		guide = "prism_minor",
		labels=scales::unit_format(big.mark ="",suffix="",scale=1e-3,sep=""),
		limits=c(1,seqlengths(seqinf)[12]),
		breaks=seq(0,seqlengths(seqinf)[12],100000),
		minor_breaks=seq(0,seqlengths(seqinf)[12],50000),
		expand=c(0,0))+
	theme(legend.position="top",axis.title.x=element_blank())+
	guides(color = guide_legend(nrow=1,override.aes = list(size = 2)))

mut2plot <- c("WT","FKH1")

input2 <- bind_rows(input %>% filter(mut %in% mut2plot))
miny=0.8
maxy=2.2

pl6 <- ggplot(input2 %>% filter(chrom==chr2plot))+
	geom_vline(data=ars2plot,aes(xintercept=(start+end)/2),col=mypal[9],show.legend=F,linewidth=0.2)+
	geom_vline(data=cen2plot,aes(xintercept=(start+end)/2),col=mypal[5],show.legend=F)+
	geom_rect(data=rdna2plot,aes(xmin=start,xmax=end,ymin=0.9,ymax=1),col=NA,fill=mypal[21],show.legend=F,alpha=0.6)+
	geom_point(aes(x=positions,y=mod,col=filename,group=filename),shape=16,size=0.2)+
	scale_color_manual("",values=mypal2)+
	ylab("RT")+
	xlab("Genomic position on chrXII (kb)")+
	coord_cartesian(ylim=c(miny,maxy),expand=F)+
	scale_x_continuous(
		guide = "prism_minor",
		labels=scales::unit_format(big.mark ="",suffix="",scale=1e-3,sep=""),
		limits=c(1,seqlengths(seqinf)[12]),
		breaks=seq(0,seqlengths(seqinf)[12],100000),
		minor_breaks=seq(0,seqlengths(seqinf)[12],50000),
		expand=c(0,0))+
	theme(legend.position="top",axis.title.x=element_blank())+
	guides(color = guide_legend(nrow=1,override.aes = list(size = 2)))

### add WT_rep24 replicat
nanot24_2 <- readRDS(path2data %+% "nanoT_WT_24rep.rds") %>%
	group_by(Rep) %>%
	mutate(mod=1+myscaling0(mean_br_bin,infq=0.005,supq=0.995)) %>%
	ungroup %>%
	mutate(chrom=factor(chrom,levels=chrom_order)) %>%
	ungroup %>%
	mutate(positions=positions+500)
#### add a capping to chrom size
nanot24_3 <- left_join(nanot24_2,chrom_sizes, by="chrom") %>%
	mutate(positions=map2_dbl(positions,seqlengths, function(x,y) {res=ifelse(x<y,x,y); return(res)})) %>%
	select(-seqlengths)

mypal3 <- c(paletteer::paletteer_d("ggthemes::Classic_20"),"grey40","gold","darkblue","greenyellow","black")


input2 <- bind_rows(nanot24_3,SortSeqWT %>% rename(Rep=filename)) %>%
	mutate(Rep=factor(Rep,levels=c("rep" %+% 1:24,"sort-seq")))
miny=0.8
maxy=2.2

pl7 <- ggplot(input2 %>% filter(chrom==chr2plot))+
	geom_vline(data=ars2plot,aes(xintercept=(start+end)/2),col=mypal[9],show.legend=F,linewidth=0.2)+
	geom_vline(data=cen2plot,aes(xintercept=(start+end)/2),col=mypal[5],show.legend=F)+
	geom_rect(data=rdna2plot,aes(xmin=start,xmax=end,ymin=0.9,ymax=1),col=NA,fill=mypal[21],show.legend=F,alpha=0.6)+
	geom_point(aes(x=positions,y=mod,col=Rep,group=Rep),shape=16,size=0.2)+
	scale_color_manual("",values=mypal3)+
	ylab("RT")+
	xlab("Genomic position on chrXII (kb)")+
	coord_cartesian(ylim=c(miny,maxy),expand=F)+
	scale_x_continuous(
		guide = "prism_minor",
		labels=scales::unit_format(big.mark ="",suffix="",scale=1e-3,sep=""),
		limits=c(1,seqlengths(seqinf)[12]),
		breaks=seq(0,seqlengths(seqinf)[12],100000),
		minor_breaks=seq(0,seqlengths(seqinf)[12],50000),
		expand=c(0,0))+
	theme(legend.position="top")+
	guides(color = guide_legend(nrow=2,override.aes = list(size = 2)))

pl <- pl1/pl2/pl3/pl4/pl5/pl6/pl7 + plot_annotation(tag_levels = 'a',title="Figure 2")

quartz(file=paste0(path2fig,"Fig2.pdf"),height=14,width=10,type="pdf")
pl
dev.off()
