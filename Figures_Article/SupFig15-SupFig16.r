#### Supplementary Figures15 and 16
#### Nanotiming analysis at X near Y' or not for KU70

suppressMessages(library(GenomicRanges))
suppressMessages(library(rtracklayer))
suppressMessages(library(tidyverse))
library(patchwork)
library(ggprism)
library(ggdist)
library(ggh4x)

theme_set(theme_bw(base_size=12))
mypal <- c(paletteer::paletteer_d("ggthemes::Classic_20"),"grey40","black","gold","greenyellow")
mypal1 <- mypal[c(1,3:20)]
mypal2 <- mypal[c(7,8,3,1,2,19,4,13,14,20,9,10)]
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

bs=1000
bingen <- tileGenome(seqinf,tilewidth=bs, cut.last.tile.in.chrom=T)

backgr <- as_tibble(bingen) %>% dplyr::rename(positions=start,chrom=seqnames)%>% mutate(y=NA) %>% select(chrom,positions,y)

## load TEL and TEL_ELMNT
GFF <- import("Reference_Genome/BT1multiUra.gff3")
seqinfo(GFF) <- seqinf
TEL_BT1 <- GFF[GFF$type=="Tel_repeat"]
XY_BT1 <- GFF[GFF$type %in% c("X_element","Y_prime_element","X_element_partial")]

### import nanot data
file2load <- c(
		"WT_rep1",
		"WT_rep2",
		"WT_rep3",
		"WT_rep4",
		"WT_rep5",
		"WT_rep6",
		"KU70_rep1",
		"KU70_rep2",
		"KU70_rep3"
		)

file_order <- c(
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
		mutate(mod=1+myscaling0(mean_br_bin,infq=0.005,supq=0.995)) %>%
		left_join(chrom_sizes, by="chrom") %>%
		mutate(end=map2_dbl(positions,seqlengths, function(x,y) {
			res=ifelse(x+999<y,x+999,y); return(res)}))
			return(res)
			}
	)
nanot[[1]] <- nanot[[1]] %>% select(-BrdU)
### sort X
Xelnt <- XY_BT1[XY_BT1$type=="X_element" | (XY_BT1$type=="X_element_partial" & seqnames(XY_BT1)=="chrV")]
Ypelnt <- XY_BT1[XY_BT1$type=="Y_prime_element"]
distXY <- as_tibble(distanceToNearest(Xelnt,Ypelnt))

Xelnt$distY <- full_join(tibble(queryHits=1:31),distXY)$distance

### SupFig15a
strip_col <- mypal[c(10,6,6,6,10,6,10,6,6,6,6,10,10,6,6,10,6,10,6,10,10,10,6,6,10,6,6,6,10,6,6,6)]
strip <- strip_themed(background_x = elem_list_rect(fill = strip_col),text_x=element_text(size=12))

Telo <- flank(TEL_BT1[TEL_BT1$name=="term"],1000)

Telo$XorY <- ifelse(Xelnt$distY>20000 | is.na(Xelnt$distY),"X","XY'")
Telo$LorR <- ifelse(start(Telo)<20000,"L","R")
Telo$chrnum <- as.numeric(as.roman(str_remove(as.character(seqnames(Telo)),"chr")))
Telo$TeloName <- ifelse(Telo$chrnum<10,paste("TEL\n0",Telo$chrnum,Telo$LorR,sep=""),paste("TEL\n",Telo$chrnum,Telo$LorR,sep=""))
mcols(Telo) <- mcols(Telo)[c(2,11,12,14)]

reslist3 <- lapply(seq_along(nanot), function(i)
	{
	nanotGR <- with(nanot[[i]], GRanges(seqnames=chrom,ranges=IRanges(start=positions,end=end),strand="*",raw=mean_br_bin,mod=mod,seqinfo=seqinf))
	muta <- nanot[[i]]$mut[1]
	Telo$nanoT <- sapply(seq_along(Telo), function(x)  mean(nanotGR[overlapsAny(nanotGR,Telo[x])]$mod,na.rm=T))
	res <- as_tibble(Telo) %>%
		mutate(filename=file2load[i]) %>%
		mutate(mut=muta)
	})

toplot2 <- do.call(bind_rows,reslist3) %>%
	mutate(filename=factor(filename,levels=c("yku70∆_rep1"="KU70_rep1","yku70∆_rep2"="KU70_rep2","yku70∆_rep3"="KU70_rep3",wt_rep1="WT_rep1",wt_rep2="WT_rep2",wt_rep3="WT_rep3",wt_rep4="WT_rep4",wt_rep5="WT_rep5",wt_rep6="WT_rep6"))) %>%
	mutate(TeloName=factor(TeloName,levels=reslist3[[1]]$TeloName)) %>%
	mutate(mut=fct_recode(mut,"yku70∆"="KU70","wt"="WT")) %>%
	mutate(mut=factor(mut,levels=c("wt","yku70∆"))) %>%
	mutate(XorY=factor(XorY,levels=c("X","XY'")))

# add dummy data for missing TEL
dummy <- tibble(TeloName="TEL\n13R",XorY="XY'",nanoT=0.2,mut="yku70∆")%>%
	mutate(TeloName=factor(TeloName,levels=reslist3[[1]]$TeloName)) %>%
	mutate(mut=factor(mut,levels=c("wt","yku70∆"))) %>%
	mutate(XorY=factor(XorY,levels=c("X","XY'")))
toplot2bis <- bind_rows(toplot2,dummy)

fS15a <- ggplot(toplot2bis,aes(x=TeloName,y=nanoT,col=mut))+
	geom_boxplot(outlier.shape = NA) +
	geom_point(position=position_dodge(width=0.75),alpha=0.5,shape=16)+
	coord_cartesian(ylim=c(0.9,2.15))+
	scale_color_manual("",values=mypal[c(7,3)])+
	facet_wrap2("TeloName",scales="free_x",nrow=1,strip = strip)+
theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank())+
	scale_y_continuous(position = "right")


### SupFig15b
strip_col <- mypal[c(10,6)]
strip <- strip_themed(background_x = elem_list_rect(fill = strip_col),text_x=element_text(size=16))

fS15b <- ggplot(toplot2,aes(x=mut,y=nanoT,fill=mut,group=mut))+
	stat_slab(col=NA,alpha=0.9,scale=0.8,normalize="groups")+
	stat_pointinterval(.width=c(.5,.95),col="black",show.legend=F)+
	scale_fill_manual("",values=mypal[c(7,3)])+
	coord_cartesian(ylim=c(0.9,2.15))+
	ylab("RT at telomeres")+
	theme(axis.text.x = element_blank(),axis.title.x=element_blank(),axis.ticks.x=element_blank())+
	facet_wrap2("XorY",scales="free_x",nrow=1,strip = strip)


layout <- "
AAAAAAAAB
"
plo <- fS15a+theme(legend.position="none")+fS15b+
	plot_layout(design = layout,heights=c(2))+
	plot_annotation(tag_levels = 'a',title="Supplementary Figure 15")

quartz(file=paste0(path2fig,"SupFig15.pdf"),height=4,width=22,type="pdf")
plo
dev.off()

### figS16
ARS_BT1 <- GFF[GFF$type=="ORI"] %>% as_tibble() %>% filter(Conf=="Confirmed",name!="ARS1216.5")
CEN_BT1 <- GFF[GFF$type=="centromere"] %>% as_tibble()

### load TEL and TEL_ELMNT
tel_BT1 <- TEL_BT1 %>% as_tibble()
xy_BT1 <- XY_BT1 %>% as_tibble() %>% mutate(type=case_when(type=="X_element_partial"~"X_element",T~type))
#
nanot <- lapply(file2load, function(x)
	{
	res <- readRDS(path2data %+% "nanoT_" %+% x %+% ".rds") %>%
		ungroup %>%
		mutate(mod=1+myscaling0(mean_br_bin,infq=0.005,supq=0.995)) %>%
		full_join(backgr,by = join_by(chrom, positions)) %>%
		arrange(chrom,positions) %>%
		select(-y)
	mut0 <- res$mut[1]
	filename0 <- res$filename[1]
	res <- res %>% mutate(mut=mut0) %>% mutate(filename=filename0)
	return(res)
	}
	)
nanot[[1]] <- nanot[[1]] %>% select(-BrdU)
nanot2 <- do.call(bind_rows,nanot) %>% mutate(chrom=factor(chrom,levels=chrom_order)) %>% ungroup %>% mutate(positions=positions+500)
#### add a capping to chrom size
nanot3 <- left_join(nanot2,chrom_sizes, by="chrom") %>% mutate(positions=map2_dbl(positions,seqlengths, function(x,y) {res=ifelse(x<y,x,y); return(res)}))
### rename
input <- nanot3 %>%
	mutate(filename=factor(filename,levels=file2load)) %>%
	mutate(mut=fct_recode(mut,"yku70∆"="KU70","wt"="WT")) %>%
	mutate(mut=factor(mut,levels=c("wt","yku70∆")))
### All Chr nanotraces
ars2plot <- ARS_BT1 %>% rename(chrom=seqnames)
cen2plot <- CEN_BT1 %>% rename(chrom=seqnames)
tel2plot <- tel_BT1 %>% rename(chrom=seqnames) %>% mutate(type="Tel_repeat") %>% filter(name=="term")
xy2plot <- xy_BT1  %>% rename(chrom=seqnames)

mypal2 <- mypal[c(7,3,9,11,13,17)]
names(mypal2) <- c("wt","yku70∆","ORI","Tel_repeat","X_element","Y_prime_element")
minxL <- tibble(chrom=seqnames(seqinf),minx=seqlengths(seqinf)-50000)
ARSteloR <- left_join(ars2plot,minxL) %>% filter(start>minx)  %>% mutate(chrom=factor(chrom,levels=chrom_order))
input2 <- input
miny=0.8
maxy=2.2
minx=0
maxx=50000
strip_col <- mypal[c(10,6,6,6,10,6,10,6,6,6,6,10,10,6,6,10,6,10,6,10,10,10,6,6,10,6,6,6,10,6,6,6)]
strip_colL <- strip_col[1+2*(0:15)]
strip_colR <- strip_col[2+2*(0:15)]

strip <- strip_themed(background_y = elem_list_rect(fill = strip_colL),text_y=element_text(size=8))

pl1 <- ggplot(input2)+
			geom_rect(data=tel2plot,aes(xmin=start,xmax=end,ymin=(miny+maxy)/2,ymax=maxy),fill=mypal[11],col=NA,show.legend=F)+
			geom_rect(data=xy2plot %>% filter(type=="X_element"),aes(xmin=start,xmax=end,ymin=miny,ymax=(miny+maxy)/2),fill=mypal[13],col=NA,show.legend=F)+
			geom_rect(data=xy2plot %>% filter(type=="Y_prime_element"),aes(xmin=start,xmax=end,ymin=miny,ymax=(miny+maxy)/2),fill=mypal[17],col=NA,show.legend=F)+
		geom_vline(data=ars2plot,aes(xintercept=(start+end)/2),col=mypal[9],show.legend=F,linewidth=0.2)+
		geom_vline(data=cen2plot,aes(xintercept=(start+end)/2),col=mypal[5],show.legend=F,alpha=0.6)+
	geom_line(aes(x=positions,y=mod,col=mut,group=filename),linewidth=0.1)+
	scale_color_manual("",values=mypal2)+
	scale_fill_manual("",values=mypal2)+
	ylab("RT")+
	xlab("Genomic position (kb)")+
	coord_cartesian(ylim=c(miny,maxy),expand=T)+
	scale_x_continuous(
		guide = "prism_minor",
		labels=scales::unit_format(big.mark ="",suffix="",scale=1e-3,sep=""),
		limits=c(minx,maxx),
		expand=expansion(add = c(1000,1000)))+
	facet_wrap2(~chrom,scales="free_x",ncol=1,strip.position = "left",strip = strip)+
	theme(legend.position="top",axis.text.y=element_blank(),axis.title.y=element_blank())+
	guides(color = guide_legend(nrow=1,override.aes = list(linewidth = 1)))+
	scale_y_continuous(position = "right")

strip <- strip_themed(background_y = elem_list_rect(fill = strip_colR),text_y=element_text(size=8))

pl2 <- ggplot(input2 %>% group_by(filename,chrom) %>% slice_tail(n=50) %>% ungroup)+
			geom_rect(data=tel2plot %>% filter(start>50000),aes(xmin=start,xmax=end,ymin=(miny+maxy)/2,ymax=maxy),fill=mypal[11],col=NA)+
			geom_rect(data=xy2plot %>% filter(type=="X_element" & start>50000),aes(xmin=start,xmax=end,ymin=miny,ymax=(miny+maxy)/2),fill=mypal[13],col=NA)+
			geom_rect(data=xy2plot %>% filter(type=="Y_prime_element" & start>50000),aes(xmin=start,xmax=end,ymin=miny,ymax=(miny+maxy)/2),fill=mypal[17],col=NA)+
	geom_vline(data=ARSteloR,aes(xintercept=(start+end)/2),col=mypal[9],show.legend=F,linewidth=0.2)+
	geom_line(aes(x=positions,y=mod,col=mut,group=filename),linewidth=0.1)+
	scale_color_manual("",values=mypal2)+
	scale_fill_manual("",values=mypal2)+
	ylab("RT")+
	xlab("Genomic position (kb)")+
	coord_cartesian(ylim=c(miny,maxy),expand=T)+
	scale_x_continuous(
		guide = "prism_minor",
		labels=scales::unit_format(big.mark ="",suffix="",scale=1e-3,sep=""),
		expand=expansion(add = c(1000,1000)))+
	facet_wrap2(~chrom,scales="free_x",ncol=1,strip.position = "right",strip = strip)+
	theme(legend.position="top")+
	guides(color = guide_legend(nrow=1,override.aes = list(linewidth = 1)))


pl <- pl1+pl2 +plot_layout(guides="collect") & theme(legend.position="top") & plot_annotation(title="Supplementary Figure 16")

quartz(file=paste0(path2fig,"SupFig16.pdf"),height=14,width=8,type="pdf")
pl
dev.off()

