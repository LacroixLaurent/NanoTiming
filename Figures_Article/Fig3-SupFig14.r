#### Figure 3 and SupFig14
#### Nanotiming analysis at X near Y' or not

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
		"RIF1_rep1",
		"RIF1_rep2",
		"RIF1_rep3"
		)

file_order <- c(
		"RIF1_rep1",
		"RIF1_rep2",
		"RIF1_rep3",
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
	})

nanot[[1]] <- nanot[[1]] %>% select(-BrdU)

### sort X
Xelnt <- XY_BT1[XY_BT1$type=="X_element" | (XY_BT1$type=="X_element_partial" & seqnames(XY_BT1)=="chrV")]
Ypelnt <- XY_BT1[XY_BT1$type=="Y_prime_element"]
distXY <- as_tibble(distanceToNearest(Xelnt,Ypelnt))

Xelnt$distY <- full_join(tibble(queryHits=1:31),distXY)$distance

### fig3a
strip_cola <- mypal[c(10,6,6,6,10,6,10,6,6,6,6,10,10,6,6,10,6,10,6,10,10,10,6,6,10,6,6,6,10,6,6,6)]
strip <- strip_themed(background_x = elem_list_rect(fill = strip_cola),text_x=element_text(size=4))

Telo <- flank(TEL_BT1[TEL_BT1$name=="term"],1000)

Telo$XorY <- ifelse(Xelnt$distY>20000 | is.na(Xelnt$distY),"X","XY'")
Telo$LorR <- ifelse(start(Telo)<20000,"L","R")
Telo$chrnum <- as.numeric(as.roman(str_remove(as.character(seqnames(Telo)),"chr")))
Telo$TeloName <- ifelse(Telo$chrnum<10,paste("TEL\n0",Telo$chrnum,Telo$LorR,sep=""),paste("TEL\n",Telo$chrnum,Telo$LorR,sep=""))
mcols(Telo) <- mcols(Telo)[c(2,11,12,14)]

reslist3 <- lapply(seq_along(nanot), function(i)
	{
	nanotGR <- with(nanot[[i]], GRanges(seqnames=chrom,ranges=IRanges(start=positions,end=end),strand="*",raw=mean_br_bin,mod=mod,seqinfo=seqinf))
	cov_nanoT <- coverage(nanotGR,weight=nanotGR$mod)
	muta <- nanot[[i]]$mut[1]
	Telo$nanoT <- sapply(seq_along(Telo), function(x)  mean(cov_nanoT[Telo[x]],na.rm=T))
	res <- as_tibble(Telo) %>%
		mutate(filename=file2load[i]) %>%
		mutate(mut=muta)
	})

toplot2 <- do.call(bind_rows,reslist3) %>%
	mutate(filename=factor(filename,levels=c("rif1∆_rep1"="RIF1_rep1","rif1∆_rep3"="RIF1_rep3","rif1∆_rep2"="RIF1_rep2",wt_rep1="WT_rep1",wt_rep2="WT_rep2",wt_rep3="WT_rep3",wt_rep4="WT_rep4",wt_rep5="WT_rep5",wt_rep6="WT_rep6"))) %>%
	mutate(TeloName=factor(TeloName,levels=reslist3[[1]]$TeloName)) %>%
	mutate(mut=fct_recode(mut,"rif1∆"="RIF1","wt"="WT")) %>%
	mutate(mut=factor(mut,levels=c("wt","rif1∆"))) %>%
	mutate(XorY=factor(XorY,levels=c("X","XY'")))

f3a <- ggplot(toplot2,aes(x=TeloName,y=nanoT,col=mut))+
	geom_boxplot(coef=NULL,outlier.shape = NA,size=0.2) +
	geom_point(position=position_dodge(width=0.75),alpha=0.5,shape=16,size=0.4)+
	coord_cartesian(ylim=c(0.9,2.15))+
	scale_color_manual("",values=mypal[c(7,1,3)])+
	facet_wrap2("TeloName",scales="free_x",nrow=1,strip = strip)+
	theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank())+
	scale_y_continuous(position = "right")

### fig3b
strip_colb <- mypal[c(10,6)]
strip <- strip_themed(background_x = elem_list_rect(fill = strip_colb),text_x=element_text(size=8))

f3b <- ggplot(toplot2,aes(x=mut,y=nanoT,fill=mut,group=mut))+
	stat_slab(col=NA,alpha=0.9,scale=0.8,normalize="groups")+
	stat_pointinterval(
		.width=c(.5,.95),
		col="black",
		show.legend=F,
		point_size = 0.1,
		interval_size_domain = c(0.5, 0.95),
		interval_size_range = c(0.17, 0.2))+
	scale_fill_manual("",values=mypal[c(7,1,3)])+
	coord_cartesian(ylim=c(0.9,2.15))+
	ylab("RT at telomeres")+
	theme(axis.text.x = element_blank(),
				axis.title.x=element_blank(),
				axis.ticks.x=element_blank())+
	facet_wrap2("XorY",scales="free_x",nrow=1,strip = strip)


### fig3c/S14
ARS_BT1 <- GFF[GFF$type=="ORI"] %>% as_tibble() %>% filter(Conf=="Confirmed",name!="ARS1216.5")
CEN_BT1 <- GFF[GFF$type=="centromere"] %>% as_tibble()

### load TEL and TEL_ELMNT
tel_BT1 <- TEL_BT1 %>% as_tibble()
xy_BT1 <- XY_BT1 %>% as_tibble() %>% mutate(type=case_when(type=="X_element_partial"~"X_element",T~type))

nanot <- lapply(file2load, function(x)
	{
	res <- readRDS(path2data %+% "nanoT_" %+% x %+% ".rds") %>%
		ungroup %>%
		mutate(mod=1+myscaling0(mean_br_bin,infq=0.005,supq=0.995))
	return(res)
	})

nanot[[1]] <- nanot[[1]] %>% select(-BrdU)
nanot2 <- do.call(bind_rows,nanot) %>% mutate(chrom=factor(chrom,levels=chrom_order)) %>% ungroup %>% mutate(positions=positions+500)
#### add a capping to chrom size
nanot3 <- left_join(nanot2,chrom_sizes, by="chrom") %>% mutate(positions=map2_dbl(positions,seqlengths, function(x,y) {res=ifelse(x<y,x,y); return(res)})) %>% select(-seqlengths)

### rename
input <- nanot3 %>%
	mutate(filename=factor(filename,levels=file_order)) %>%
	mutate(mut=fct_recode(mut,"rif1∆"="RIF1","wt"="WT")) %>%
	mutate(mut=factor(mut,levels=c("wt","rif1∆")))

### All Chr nanotraces
ars2plot <- ARS_BT1 %>% rename(chrom=seqnames)
cen2plot <- CEN_BT1 %>% rename(chrom=seqnames)
tel2plot <- tel_BT1 %>% rename(chrom=seqnames) %>% mutate(type="Tel_repeat") %>% filter(name=="term")
xy2plot <- xy_BT1  %>% rename(chrom=seqnames)

mypal2 <- mypal[c(7,1,3,9,11,13,17)]
names(mypal2) <- c("wt","rif1∆","ku70∆","ORI","Tel_repeat","X_element","Y_prime_element")
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

strip <- strip_themed(background_y = elem_list_rect(fill = strip_colL),text_y=element_text(size=6))


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
#	scale_x_continuous(limits=c(minx,maxx),expand = expansion(add = c(1000,1000))) +
	facet_wrap2(~chrom,scales="free_x",ncol=1,strip.position = "left",strip = strip)+
	theme(legend.position="top",axis.text.y=element_blank(),axis.title.y=element_blank())+
	guides(color = guide_legend(nrow=1,override.aes = list(linewidth = 1)))+
	scale_y_continuous(position = "right")

strip <- strip_themed(background_y = elem_list_rect(fill = strip_colR),text_y=element_text(size=6))

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


pl <- pl1+pl2 +plot_layout(guides="collect") & theme(legend.position="top") & plot_annotation(title="Supplementary Figure 14")

quartz(file=paste0(path2fig,"SupFig14.pdf"),height=14,width=8,type="pdf")
pl
dev.off()

#####
ars2plot <- ARS_BT1 %>% rename(chrom=seqnames)%>% filter(chrom %in% c("chrXI","chrXII","chrXIII","chrXIV","chrXV"))
cen2plot <- CEN_BT1 %>% rename(chrom=seqnames)%>% filter(chrom %in% c("chrXI","chrXII","chrXIII","chrXIV","chrXV"))
tel2plot <- tel_BT1 %>% rename(chrom=seqnames) %>% mutate(type="Tel_repeat") %>% filter(name=="term")%>% filter(chrom %in% c("chrXI","chrXII","chrXIII","chrXIV","chrXV"))
xy2plot <- xy_BT1  %>% rename(chrom=seqnames)%>% filter(chrom %in% c("chrXI","chrXII","chrXIII","chrXIV","chrXV"))
mypal2 <- mypal[c(7,1,9,11,13,17)]
names(mypal2) <- c("wt","rif1∆","ORI","Tel_repeat","X_element","Y_prime_element")
minxL <- tibble(chrom=seqnames(seqinf),minx=seqlengths(seqinf)-50000)
ARSteloR <- left_join(ars2plot,minxL) %>% filter(start>minx)  %>% mutate(chrom=factor(chrom,levels=c("chrI","chrII","chrIII","chrIV","chrV","chrVI","chrVII","chrVIII","chrIX","chrX","chrXI","chrXII","chrXIII","chrXIV","chrXV","chrXVI")))
input2 <- input %>% filter(chrom %in% c("chrXI","chrXII","chrXIII","chrXIV","chrXV"))
miny=0.8
maxy=2.2
minx=0
maxx=50000

strip_colL1 <- strip_colL[11:15]
strip_colR1 <- strip_colR[11:15]

strip <- strip_themed(background_y = elem_list_rect(fill = strip_colL1),text_y=element_text(size=5))

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
#	scale_x_continuous(limits=c(minx,maxx),expand = expansion(add = c(1000,1000))) +
	facet_wrap2(~chrom,scales="free_x",ncol=1,strip.position = "left",strip = strip)+
	theme(axis.text.y=element_blank(),axis.title.y=element_blank())+
	guides(color = guide_legend(nrow=2,override.aes = list(linewidth = 1)))+
	scale_y_continuous(position = "right")

strip <- strip_themed(background_y = elem_list_rect(fill = strip_colR1),text_y=element_text(size=5))

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
	guides(color = guide_legend(nrow=2,override.aes = list(linewidth = 1)))


pl <- pl1+pl2 & theme(legend.position="right",axis.text = element_text(size = 8))

layout <- "
AAAAAAAAB
CCCCCCCCC
"
plo <- f3a+theme(legend.position="none")+f3b+theme(legend.position="none")+pl+
	plot_layout(design = layout,heights=c(2,5),guides="collect")+
	plot_annotation(tag_levels = 'a',title="Figure 3")

#quartz(file=paste0(path2fig,"Fig3.pdf"),height=12,width=22,type="pdf")
#plo
#dev.off()


plA4 <- plo & theme(
	text = element_text(size = 8),  # Taille des polices
	legend.margin = margin(t = 0, r = 0, b = 0, l = 0),
	legend.text = element_text(size = 5),
	legend.key.size = unit(3, "mm"),
	axis.title = element_text(size = 7),
	axis.text = element_text(size = 4),
	axis.line=element_line(size=0.1),
	axis.ticks = element_line(size=0.2),
	plot.title = element_text(size = 9),
	plot.margin = unit(c(1, 1, 1, 1), "mm"),
	panel.spacing = unit(0.1, "lines"),
	panel.grid.minor = element_line(linewidth = 0.1),
	panel.grid.major = element_line(linewidth = 0.2),
	panel.border = element_rect(linewidth = 0.4)
)
quartz(file=paste0(path2fig,"Fig3.pdf"),width=7.05,height=3.85,type="pdf")
plA4
dev.off()
