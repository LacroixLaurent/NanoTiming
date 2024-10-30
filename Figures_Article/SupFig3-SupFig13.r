### Supplementary Figures 3 to 13
### comparison nanotrace vs others_ All chromosomes from Figure 2

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
nanot3 <- left_join(nanot2,chrom_sizes, by="chrom") %>% mutate(positions=map2_dbl(positions,seqlengths, function(x,y) {res=ifelse(x<y,x,y); return(res)}))

### Import sortseq data with scaling (and MFAseq also)
SortSeqWT <- readRDS(path2data %+% "sortseq_WT.rds") %>%
	mutate(mod=1+myscaling0(timing,infq=0.005,supq=0.995))%>%
	mutate(filename="sort-seq") %>%
	mutate(mut="WT") %>%
	mutate(positions=positions+500) %>%
	left_join(.,chrom_sizes, by="chrom") %>%
	mutate(positions=map2_dbl(positions,seqlengths, function(x,y) {res=ifelse(x<y,x,y); return(res)})) %>%
	mutate(chrom=factor(chrom,levels=chrom_order)) %>% select(-seqlengths)

SortSeqCTF19 <- readRDS(path2data %+% "sortseq_CTF19.rds") %>%
	mutate(mod=1+myscaling0(timing,infq=0.005,supq=0.995))%>%
	mutate(filename="sort-seq_ctf19") %>%
	mutate(mut="CTF19") %>%
	mutate(positions=positions+500) %>%
	left_join(.,chrom_sizes, by="chrom") %>%
	mutate(positions=map2_dbl(positions,seqlengths, function(x,y) {res=ifelse(x<y,x,y); return(res)})) %>%
	mutate(chrom=factor(chrom,levels=chrom_order)) %>% select(-seqlengths)

SortSeqRIF1 <- readRDS(path2data %+% "sortseq_RIF1.rds") %>%
	mutate(mod=1+myscaling0(timing,infq=0.005,supq=0.995))%>%
	mutate(filename="sort-seq_rif1") %>%
	mutate(mut="RIF1") %>%
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

### rename filename
input <- nanot3 %>%
	mutate(filename=factor(filename,levels=file_order)) %>% mutate(filename=fct_recode(filename,"fkh1∆_rep1"="FKH1_rep1","fkh1∆_rep2"="FKH1_rep2","fkh1∆_rep3"="FKH1_rep3","ctf19∆_rep1"="CTF19_rep1","ctf19∆_rep2"="CTF19_rep2","ctf19∆_rep3"="CTF19_rep3","rif1∆_rep1"="RIF1_rep1","rif1∆_rep3"="RIF1_rep3","rif1∆_rep2"="RIF1_rep2","yku70∆_rep1"="KU70_rep1","yku70∆_rep2"="KU70_rep2","yku70∆_rep3"="KU70_rep3",wt_rep1="WT_rep1",wt_rep2="WT_rep2",wt_rep3="WT_rep3",wt_rep4="WT_rep4",wt_rep5="WT_rep5",wt_rep6="WT_rep6"))

### Correlation plots
## SupFig6
# WT 5µM correlation with SortSeq and MFAseq
toplot1 <- input %>% filter(mut=="WT") %>% select(c(1,2,4,6)) %>% pivot_wider(values_from=mod, names_from=filename, names_prefix="nanoT_")
tocor <- left_join(toplot1,SortSeqWT %>% select(c(1,3,4)) %>% rename("sort-seq_wt"=mod)) %>% left_join(.,MFAseq %>% select(c(1,3,4)) %>% rename("MFA-seq_wt"=mod))
cormat1 <- sapply(3:10, function(i) sapply((3:10), function(j) cor(tocor[,i],tocor[,j],use="pairwise.complete.obs",method="s")))
colnames(cormat1) <- rownames(cormat1) <- names(tocor)[3:10]
plS6 <- ggcorrplot(cormat1,lab=T,lab_size=4,digits=3)+
	scale_fill_gradient2(limit=c(0,1),low="blue",high="red",mid="white",midpoint=0.5)+
	ggtitle("Supplementary Figure 6",subtitle="Spearman's rank correlation coefficients")

quartz(file=paste0(path2fig,"SupFig6.pdf"),height=14,width=10,type="pdf")
plS6
dev.off()

## SupFig13
toplot1 <- input %>% filter(mut=="CTF19") %>% select(c(1,2,4,6)) %>% pivot_wider(values_from=mod, names_from=filename, names_prefix="nanoT_")
tocor <- left_join(toplot1,SortSeqCTF19 %>% select(c(1,3,4)) %>% rename("sort-seq_ctf19∆"=mod))
cormat1 <- sapply(3:6, function(i) sapply(3:6, function(j) cor(tocor[,i],tocor[,j],use="pairwise.complete.obs",method="s")))
colnames(cormat1) <- rownames(cormat1) <- names(tocor)[3:6]
plS13a <- ggcorrplot(cormat1,lab=T,lab_size=4,digits=3)+
	scale_fill_gradient2(limit=c(0,1),low="blue",high="red",mid="white",midpoint=0.5)+
	ggtitle("Spearman's rank correlation coefficients")

toplot1 <- input %>% filter(mut=="RIF1") %>% select(c(1,2,4,6)) %>% pivot_wider(values_from=mod, names_from=filename, names_prefix="nanoT_")
tocor <- left_join(toplot1,SortSeqRIF1 %>% select(c(1,3,4)) %>% rename("sort-seq_rif1∆"=mod))
cormat1 <- sapply(3:6, function(i) sapply((3:6), function(j) cor(tocor[,i],tocor[,j],use="pairwise.complete.obs",method="s")))
colnames(cormat1) <- rownames(cormat1) <- names(tocor)[3:6]
plS13b <- ggcorrplot(cormat1,lab=T,lab_size=4,digits=3)+
	scale_fill_gradient2(limit=c(0,1),low="blue",high="red",mid="white",midpoint=0.5)+
	ggtitle("Spearman's rank correlation coefficients")

quartz(file=paste0(path2fig,"SupFig13.pdf"),height=6,width=12,type="pdf")
plS13a + plS13b + plot_annotation(tag_levels = 'a',title="Supplementary Figure 13") + plot_layout(guides="collect")
dev.off()

### All Chr nanotraces
ars2plot <- ARS %>% rename(chrom=seqnames)
cen2plot <- CEN %>% rename(chrom=seqnames)
rdna2plot <- rDNA %>% rename(chrom=seqnames) %>% group_by(chrom,type) %>% summarise(start=min(start),end=max(end))

## SupFig5
mut2plot <- c("WT")
input2 <- input %>% filter(mut %in% mut2plot)
miny=0.8
maxy=2.2

pl <- ggplot(input2)+
	geom_vline(data=ars2plot,aes(xintercept=(start+end)/2),col=mypal[9],show.legend=F,linewidth=0.2)+
	geom_vline(data=cen2plot,aes(xintercept=(start+end)/2),col=mypal[5],show.legend=F)+
	geom_rect(data=rdna2plot,aes(xmin=start,xmax=end,ymin=0.9,ymax=1),col=NA,fill=mypal[21],show.legend=F,alpha=0.6)+
	geom_point(aes(x=positions,y=mod,col=filename,group=filename),shape=16,size=0.2)+
	scale_color_manual("",values=mypal2[-c(1,2,3)])+
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
	guides(color = guide_legend(nrow=1,override.aes = list(size = 2)))+
	ggtitle("Supplementary Figure 5")

quartz(file=paste0(path2fig,"SupFig5.pdf"),height=14,width=12,type="pdf")
pl
dev.off()


## SupFig3
input2 <- bind_rows(input %>% filter(filename=="wt_rep1"),SortSeqWT)
miny=0.8
maxy=2.2

pl <- ggplot(input2)+
	geom_vline(data=ars2plot,aes(xintercept=(start+end)/2),col=mypal[9],show.legend=F,linewidth=0.2)+
	geom_vline(data=cen2plot,aes(xintercept=(start+end)/2),col=mypal[5],show.legend=F)+
	geom_rect(data=rdna2plot,aes(xmin=start,xmax=end,ymin=0.9,ymax=1),col=NA,fill=mypal[21],show.legend=F,alpha=0.6)+
	geom_point(aes(x=positions,y=mod,col=filename,group=filename),shape=16,size=0.2)+
	scale_color_manual("",values=mypal1,labels = c("Relative copy number by sort-seq", "Mean BrdU content"))+
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
	ggtitle("Supplementary Figure 3")

quartz(file=paste0(path2fig,"SupFig3.pdf"),height=14,width=12,type="pdf")
pl
dev.off()

## SupFig4
input2 <- bind_rows(input %>% filter(filename=="wt_rep1"),MFAseq)
miny=0.8
maxy=2.2

pl <- ggplot(input2)+
	geom_vline(data=ars2plot,aes(xintercept=(start+end)/2),col=mypal[9],show.legend=F,linewidth=0.2)+
	geom_vline(data=cen2plot,aes(xintercept=(start+end)/2),col=mypal[5],show.legend=F)+
	geom_rect(data=rdna2plot,aes(xmin=start,xmax=end,ymin=0.9,ymax=1),col=NA,fill=mypal[21],show.legend=F,alpha=0.6)+
	geom_point(aes(x=positions,y=mod,col=filename,group=filename),shape=16,size=0.2)+
	scale_color_manual("",values=mypal1,labels = c("Relative copy number by MFA-seq", "Mean BrdU content"))+
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
	ggtitle("Supplementary Figure 4")

quartz(file=paste0(path2fig,"SupFig4.pdf"),height=14,width=12,type="pdf")
pl
dev.off()

## SupFig11
input2 <- bind_rows(input %>% filter(filename=="ctf19∆_rep3"),SortSeqCTF19) %>%
mutate(filename=factor(filename,levels=rev(unique(filename))))
miny=0.8
maxy=2.2

pl <- ggplot(input2)+
	geom_vline(data=ars2plot,aes(xintercept=(start+end)/2),col=mypal[9],show.legend=F,linewidth=0.2)+
	geom_vline(data=cen2plot,aes(xintercept=(start+end)/2),col=mypal[5],show.legend=F)+
	geom_rect(data=rdna2plot,aes(xmin=start,xmax=end,ymin=0.9,ymax=1),col=NA,fill=mypal[21],show.legend=F,alpha=0.6)+
	geom_point(aes(x=positions,y=mod,col=filename,group=filename),shape=16,size=0.2)+
	scale_color_manual("",values=mypal1,labels =c("Relative copy number by sort-seq ctf19∆","Mean BrdU content ctf19∆"))+
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
	ggtitle("Supplementary Figure 11")

quartz(file=paste0(path2fig,"SupFig11.pdf"),height=14,width=12,type="pdf")
pl
dev.off()

## SupFig12
input2 <- bind_rows(input %>% filter(filename=="rif1∆_rep1"),SortSeqRIF1)%>%
mutate(filename=factor(filename,levels=rev(unique(filename))))
miny=0.8
maxy=2.2

pl <- ggplot(input2)+
	geom_vline(data=ars2plot,aes(xintercept=(start+end)/2),col=mypal[9],show.legend=F,linewidth=0.2)+
	geom_vline(data=cen2plot,aes(xintercept=(start+end)/2),col=mypal[5],show.legend=F)+
	geom_rect(data=rdna2plot,aes(xmin=start,xmax=end,ymin=0.9,ymax=1),col=NA,fill=mypal[21],show.legend=F,alpha=0.6)+
	geom_point(aes(x=positions,y=mod,col=filename,group=filename),shape=16,size=0.2)+
	scale_color_manual("",values=mypal1,labels =c("Relative copy number by sort-seq rif1∆","Mean BrdU content rif1∆"))+
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
	ggtitle("Supplementary Figure 12")

quartz(file=paste0(path2fig,"SupFig12.pdf"),height=14,width=12,type="pdf")
pl
dev.off()

#### MUT vs WT
## SupFig7
mut2plot <- c("WT","CTF19")

input2 <- bind_rows(input %>% filter(mut %in% mut2plot))
miny=0.8
maxy=2.2

pl <- ggplot(input2)+
	geom_vline(data=ars2plot,aes(xintercept=(start+end)/2),col=mypal[9],show.legend=F,linewidth=0.2)+
	geom_vline(data=cen2plot,aes(xintercept=(start+end)/2),col=mypal[5],show.legend=F)+
	geom_rect(data=rdna2plot,aes(xmin=start,xmax=end,ymin=0.9,ymax=1),col=NA,fill=mypal[21],show.legend=F,alpha=0.6)+
	geom_point(aes(x=positions,y=mod,col=filename,group=filename),shape=16,size=0.2)+
	scale_color_manual("",values=mypal2)+
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
	guides(color = guide_legend(nrow=1,override.aes = list(size = 2)))+
	ggtitle("Supplementary Figure 7")

quartz(file=paste0(path2fig,"SupFig7.pdf"),height=14,width=12,type="pdf")
pl
dev.off()

## SupFig8
mut2plot <- c("WT","RIF1")

input2 <- bind_rows(input %>% filter(mut %in% mut2plot))
miny=0.8
maxy=2.2

pl <- ggplot(input2)+
	geom_vline(data=ars2plot,aes(xintercept=(start+end)/2),col=mypal[9],show.legend=F,linewidth=0.2)+
	geom_vline(data=cen2plot,aes(xintercept=(start+end)/2),col=mypal[5],show.legend=F)+
	geom_rect(data=rdna2plot,aes(xmin=start,xmax=end,ymin=0.9,ymax=1),col=NA,fill=mypal[21],show.legend=F,alpha=0.6)+
	geom_point(aes(x=positions,y=mod,col=filename,group=filename),shape=16,size=0.2)+
	scale_color_manual("",values=mypal2)+
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
	guides(color = guide_legend(nrow=1,override.aes = list(size = 2)))+
	ggtitle("Supplementary Figure 8")

quartz(file=paste0(path2fig,"SupFig8.pdf"),height=14,width=12,type="pdf")
pl
dev.off()

## SupFig9
mut2plot <- c("WT","KU70")

input2 <- bind_rows(input %>% filter(mut %in% mut2plot))
miny=0.8
maxy=2.2

pl <- ggplot(input2)+
	geom_vline(data=ars2plot,aes(xintercept=(start+end)/2),col=mypal[9],show.legend=F,linewidth=0.2)+
	geom_vline(data=cen2plot,aes(xintercept=(start+end)/2),col=mypal[5],show.legend=F)+
	geom_rect(data=rdna2plot,aes(xmin=start,xmax=end,ymin=0.9,ymax=1),col=NA,fill=mypal[21],show.legend=F,alpha=0.6)+
	geom_point(aes(x=positions,y=mod,col=filename,group=filename),shape=16,size=0.2)+
	scale_color_manual("",values=mypal2)+
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
	guides(color = guide_legend(nrow=1,override.aes = list(size = 2)))+
	ggtitle("Supplementary Figure 9")

quartz(file=paste0(path2fig,"SupFig9.pdf"),height=14,width=12,type="pdf")
pl
dev.off()

## SupFig10
mut2plot <- c("WT","FKH1")

input2 <- bind_rows(input %>% filter(mut %in% mut2plot))
miny=0.8
maxy=2.2

pl <- ggplot(input2)+
	geom_vline(data=ars2plot,aes(xintercept=(start+end)/2),col=mypal[9],show.legend=F,linewidth=0.2)+
	geom_vline(data=cen2plot,aes(xintercept=(start+end)/2),col=mypal[5],show.legend=F)+
	geom_rect(data=rdna2plot,aes(xmin=start,xmax=end,ymin=0.9,ymax=1),col=NA,fill=mypal[21],show.legend=F,alpha=0.6)+
	geom_point(aes(x=positions,y=mod,col=filename,group=filename),shape=16,size=0.2)+
	scale_color_manual("",values=mypal2)+
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
	guides(color = guide_legend(nrow=1,override.aes = list(size = 2)))+
	ggtitle("Supplementary Figure 10")

quartz(file=paste0(path2fig,"SupFig10.pdf"),height=14,width=12,type="pdf")
pl
dev.off()

