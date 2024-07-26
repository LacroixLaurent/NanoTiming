### Script for Figures 4 and S15 to S19

suppressMessages(library(GenomicRanges))
suppressMessages(library(rtracklayer))
suppressMessages(library(tidyverse))
library(patchwork)
library(ggdist)
library(ggpubr)
library(ggpointdensity)
library(ggh4x)

`%+%` <- paste0
theme_set(theme_bw())
mypal <- c(paletteer::paletteer_d("ggthemes::Classic_20"),"red","black","gold","greenyellow")
mypal1 <- mypal[c(1,3:20)]
mypal2 <- mypal[c(7,1,3,5,9)]
mypal3 <- c(mypal,rev(RColorBrewer::brewer.pal(8,"Dark2")))

setwd("/Users/ll/work/RStudioProjects/NanoTiming")
path2fig <- "Figures_Article/Figures_pdf/"
path2data <- "Figures_Article/Figures_data/"
path2table <- "Figures_Article/Tables/"

source("Figures_Article/rescaling_function.r")

chrom_order <- c("chrI","chrII","chrIII","chrIV","chrV","chrVI","chrVII","chrVIII","chrIX","chrX","chrXI","chrXII","chrXIII","chrXIV","chrXV","chrXVI")

seqinf <- readRDS("Reference_Genome/seqinfBT1multiUra.rds")
seqlevels(seqinf)=seqlevels(seqinf)[1:16]

chrom_sizes <- as.data.frame(seqinf)
chrom_sizes$chrom <- rownames(chrom_sizes)
chrom_sizes <- as_tibble(chrom_sizes) %>%
	select(chrom,seqlengths) %>%
	mutate(chrom=factor(chrom,levels=chrom_order))

tel_order0 <- sapply(1:16, function(i) ifelse(i<10,paste0("TEL0",i),paste0("TEL",i)))
tel_order <- do.call(c,lapply(tel_order0, function(i) paste0(i,c("L","R"))))

mut_order <- c("WT","RIF1","KU70","CTF19","FKH1")
Exp_order <- c("WT_noBrdU","WT_rep1","WT_rep2","WT_rep3","WT_rep4","WT_rep5","WT_rep6","RIF1_rep1","RIF1_rep2","RIF1_rep3","KU70_rep1","KU70_rep2","KU70_rep3","CTF19_rep1","CTF19_rep2","CTF19_rep3","FKH1_rep1","FKH1_rep2","FKH1_rep3")

################################################################################
### load Raw telomeric data (can be skipped by using the intermediate result file TeloLengthDataNanoT.rds or TeloLengthDataNanoTseqFilter2.rds)

telo_data <- readRDS(path2data %+% "TeloDataRaw.rds")

### compute metrics on the telomeric sequences
TeloList <- split(telo_data,telo_data$filename)

## filtering function)
require(Biostrings)

H_inW <- function(sequ,w=50) {
	H <- double()
	for (i in 1:(nchar(sequ)-w+1)) {
		input <- sequ[i:(i+w-1)]
		freq <- dinucleotideFrequency(input,as.prob=T)
		H <- c(H,-sum(freq * log(freq), na.rm = TRUE))
		}
	return(H)
	}

freqdi_inW <- function(sequ,w=50) {
	freqdi_list <- lapply((1:(nchar(sequ)-w+1)), function(i) {
		input <- sequ[i:(i+w-1)]
		freq <- dinucleotideFrequency(input,as.prob=T)
		})
	return(freqdi_list)
	}

freqGG_inW <- function(freqdi) {
	freqGG <- sapply(freqdi, "[[",11)
	return(freqGG)
	}
freqTelo_inW <- function(freqdi) {
	freqGG <- sapply(freqdi, function(x) {
		res <- sum(x[c(11,12,15)])
		return(res)
		})
	return(freqGG)
	}

library(furrr)
options(future.rng.onMisuse = "ignore")
options(future.globals.maxSize = 4e9)
nc <- 8L

plan(multicore,workers=nc)

toto <- lapply(seq_along(TeloList), function(i)
	{
	input <- TeloList[[i]] %>%
	mutate(seqlen=map_dbl(Seq,function(x) if(!is.na(x)) {width(x)}else{0})) %>%
	mutate(seq_data=future_map2(Seq,Type, function(x,y) {
		res <- tibble(pos=integer(),H=double(),freqTelo=double(),freqGG=double())
		if (!is.na(x))
			{
			if (nchar(x)<=50 & nchar(x)>0)
				{
				if (y=="left")
					{sequ=complement(DNAString(x))}
					else
					{sequ=DNAString(x)}
				res1 <- dinucleotideFrequency(sequ,as.prob=T)
				H_res <- -sum(res1  * log(res1), na.rm = TRUE)
				freqGG_res <- res1[11]
				freqTelo_res <- sum(res1[c(11,12,15)])
				res <- tibble(
					pos=seq_along(H_res),H=H_res,freqTelo=freqTelo_res,freqGG=freqGG_res)
					}
			if (nchar(x)>50)
				{
				if (y=="left")
					{sequ=complement(DNAString(x))}
					else
					{sequ=DNAString(x)}
				w <- 50
				res1 <- freqdi_inW(sequ,w)
				H_res <- H_inW(sequ,w)
				freqGG_res <- freqGG_inW(res1)
				freqTelo_res <- freqTelo_inW(res1)
				res <- tibble(
					pos=seq_along(H_res),H=H_res,freqTelo=freqTelo_res,freqGG=freqGG_res)
				}
			}
		return(res)
		}))
	saveRDS(input,file="TeloL00_"%+% names(TeloList)[i] %+%"_"%+% today %+%".rds")
	return(input)
	})
telo_data2 <- do.call(bind_rows,toto)
#### this part requires a significant computing time.
telo_data2 <- do.call(bind_rows,lapply(Exp_order, function(x) readRDS("TeloL00_"%+% x %+%"_20240704.rds")))

### rescale the RT value using each experiment 99% quantile
rescaling_limit <- readRDS(path2data %+% "rescaling_limit.rds")

telo_data3 <- left_join(telo_data2,rescaling_limit) %>%
	mutate(meanBscaled=pmap_dbl(., function(meanB,infq,supq,...)
		{
		1+(meanB-infq)/(supq-infq)
		}))

#### an intermediate file is provided on Zenodo (REF)

telo_data3 <- readRDS(path2data %+% "TeloLengthDataNanoT.rds")

## Compute seqfilter2 data based on the telomeric sequence metrics

telo_data4 <- telo_data3 %>%
	mutate(seqfilter2=map_dbl(seq_data, function(x) {
		res <- 0
		if (nrow(x)>1)
			{
			maxlim <- ifelse(nrow(x)>25,nrow(x)-25,nrow(x))
			noGG <- sum(!between(x$freqGG[1:maxlim],0.05,0.4))
			noTelo <- sum(x$freqTelo[1:maxlim]<=0.79)
			noH <- sum(!between(x$H[1:maxlim],0.85,1.7))
			res <- max(noGG,noTelo,noH)/maxlim
			}
		if (nrow(x)==1)
			{
			res <- 1-as.numeric(
				min(x$freqGG)>0.05 &
				max(x$freqGG)<0.4 &
				min(x$freqTelo)>0.79 &
				min(x$H)>0.85 &
				max(x$H)<1.7
			)
			}
		return(res)
		})) %>%
	select(-c(Seq,seq_data))

################################################################################

### to fit on github, an intermediate data file without Seq and seq_data is provided
telo_data4 <- readRDS(path2data %+% "TeloLengthDataNanoTseqFilter2.rds")

#### Fig4a (no BrdU filter)
### I need to remove seqlen>10000 otherwise the stat_slab is ugly
toplot <- telo_data4 %>%
	filter(seqlen<10000) %>%
	filter(seqfilter2<0.1) %>%
	mutate(filename=factor(filename,levels=Exp_order)) %>%
	mutate(mut=factor(mut,levels=mut_order)) %>%
		mutate(mut=fct_recode(mut,"rif1∆"="RIF1","yku70∆"="KU70","cft19∆"="CTF19","fkh1∆"="FKH1","wt"="WT")) %>%
	mutate(filename=fct_recode(filename,"fkh1∆_rep1"="FKH1_rep1","fkh1∆_rep2"="FKH1_rep2","fkh1∆_rep3"="FKH1_rep3","ctf19∆_rep1"="CTF19_rep1","ctf19∆_rep2"="CTF19_rep2","ctf19∆_rep3"="CTF19_rep3","rif1∆_rep1"="RIF1_rep1","rif1∆_rep3"="RIF1_rep3","rif1∆_rep2"="RIF1_rep2","yku70∆_rep1"="KU70_rep1","yku70∆_rep2"="KU70_rep2","yku70∆_rep3"="KU70_rep3",wt_rep1="WT_rep1",wt_rep2="WT_rep2",wt_rep3="WT_rep3",wt_rep4="WT_rep4",wt_rep5="WT_rep5",wt_rep6="WT_rep6",wt_noBrdU="WT_noBrdU"))

totext <- toplot %>%
	group_by(filename,mut) %>%
	summarise(mea=round(mean(seqlen, na.rm=T)),med=round(median(seqlen, na.rm=T)),n=n(),.groups="drop") %>%
	ungroup
mypal2 <- mypal[c(7,1,3,17,19)]

f4a <- ggplot(toplot,aes(x=filename,y=seqlen,fill=mut))+
	coord_cartesian(ylim=c(0,1100))+
	theme(axis.text.x = element_text(angle = 45,hjust=1),axis.title.x=element_blank())+
	stat_slab(col=NA,alpha=0.9,scale=0.8,normalize="groups")+
	stat_pointinterval(.width=c(.5,.95),col="black",show.legend=F)+
	geom_text(data=totext,aes(x=filename,y=1100,label=med),col="red") +
	geom_text(data=totext,aes(x=filename,y=0,label=n),fontface="italic",col="black") +
	scale_fill_manual("",values=mypal2)+
	ylab("Telomeric repeat length (bp)")

### fig4b and S15-19

# keep replicated telomeres (meanB>0.02) with "good" telomeric repeats (seqfilter2<0.1)
toplot <- telo_data4 %>% filter(meanB>0.02 & seqfilter2<0.1 & Rep !="noBrdU") %>% mutate(teloname=factor(teloname,levels=tel_order))

### add a dummy point for KU70 and TEL13R
dum <- tibble(Chromsome="chrXIII",read_id="dummy",mut="KU70",meanBscaled=0.1,seqlen=266,teloname="TEL13R")
toplotb <- bind_rows(toplot,dum) %>% mutate(teloname=factor(teloname,levels=tel_order))

strip_col0 <- mypal[c(10,6,6,6,10,6,10,6,6,6,6,10,10,6,6,10,6,10,6,10,10,10,6,6,10,6,6,6,10,6,6,6)]
strip <- strip_themed(background_y = elem_list_rect(fill = strip_col0))

# FigS15
plWT_all <- ggplot(toplot %>% filter(mut=="WT"))+
	geom_pointdensity(aes(x=seqlen,y=meanBscaled),shape=16,size=1)+
	stat_cor(aes(x=seqlen,y=meanBscaled),show.legend=F,label.x=650,label.y=c(3),method="spearman",cor.coef.name="rho",size=2)+
	geom_text(data=toplot %>% filter(mut=="WT") %>% group_by(teloname) %>% summarise(n=n()),aes(label="n=" %+% n),x=50,y=3,size=2)+
	coord_cartesian(ylim=c(0.5,3.5))+
	scale_fill_continuous(type = "viridis",direction=-1,option="magma")+
	scale_color_continuous(type = "viridis",direction=-1,option="magma")+
	facet_wrap2("teloname",ncol=2,strip.position = "right",strip = strip)+
	ylab("RT")+
	xlab("Telomeric repeat length (bp)")+
	ggtitle("Figure S15",subtitle="wt")

# FigS16
plRIF1_all <- ggplot(toplot %>% filter(mut=="RIF1"))+
	geom_pointdensity(aes(x=seqlen,y=meanBscaled),shape=16,size=1)+
	geom_text(data=toplot %>% filter(mut=="RIF1") %>% group_by(teloname) %>% summarise(n=n()),aes(label="n=" %+% n),x=100,y=3,size=2)+
	stat_cor(aes(x=seqlen,y=meanBscaled),show.legend=F,label.x=1100,label.y=c(3),method="spearman",cor.coef.name="rho",size=2)+
	coord_cartesian(ylim=c(0.5,3.5))+
	scale_fill_continuous(type = "viridis",direction=-1,option="magma")+
	scale_color_continuous(type = "viridis",direction=-1,option="magma")+
	facet_wrap2("teloname",ncol=2,strip.position = "right",strip = strip)+
	ylab("RT")+
	xlab("Telomeric repeat length (bp)")+
	ggtitle("Figure S16",subtitle="rif1∆")

# FigS17
plKU70_all <- ggplot(toplotb %>% filter(mut=="KU70"))+
	geom_pointdensity(aes(x=seqlen,y=meanBscaled),shape=16,size=1)+
	geom_text(data=toplotb %>% filter(mut=="KU70",read_id!="dummy") %>% group_by(teloname) %>% summarise(n=n()),aes(label="n=" %+% n),x=20,y=3,size=2)+
	stat_cor(aes(x=seqlen,y=meanBscaled),show.legend=F,label.x=350,label.y=c(3),method="spearman",cor.coef.name="rho",size=2)+
	coord_cartesian(ylim=c(0.5,3.5))+
	scale_fill_continuous(type = "viridis",direction=-1,option="magma")+
	scale_color_continuous(type = "viridis",direction=-1,option="magma")+
	facet_wrap2("teloname",ncol=2,strip.position = "right",strip = strip)+
	ylab("RT")+
	xlab("Telomeric repeat length (bp)")+
	ggtitle("Figure S17",subtitle="yku70∆")

# FigS18
plCTF19_all <- ggplot(toplot %>% filter(mut=="CTF19"))+
	geom_pointdensity(aes(x=seqlen,y=meanBscaled),shape=16,size=1)+
	geom_text(data=toplot %>% filter(mut=="CTF19") %>% group_by(teloname) %>% summarise(n=n()),aes(label="n=" %+% n),x=50,y=3,size=2)+
	stat_cor(aes(x=seqlen,y=meanBscaled),show.legend=F,label.x=450,label.y=c(3),method="spearman",cor.coef.name="rho",size=2)+
	coord_cartesian(ylim=c(0.5,3.5))+
	scale_fill_continuous(type = "viridis",direction=-1,option="magma")+
	scale_color_continuous(type = "viridis",direction=-1,option="magma")+
	facet_wrap2("teloname",ncol=2,strip.position = "right",strip = strip)+
	ylab("RT")+
	xlab("Telomeric repeat length (bp)")+
	ggtitle("Figure S18",subtitle="cft19∆")


# FigS19
plFKH1_all <- ggplot(toplot %>% filter(mut=="FKH1"))+
	geom_pointdensity(aes(x=seqlen,y=meanBscaled),shape=16,size=1)+
	geom_text(data=toplot %>% filter(mut=="FKH1") %>% group_by(teloname) %>% summarise(n=n()),aes(label="n=" %+% n),x=50,y=3,size=2)+
	stat_cor(aes(x=seqlen,y=meanBscaled),show.legend=F,label.x=450,label.y=c(3),method="spearman",cor.coef.name="rho",size=2)+
	coord_cartesian(ylim=c(0.5,3.5))+
	scale_fill_continuous(type = "viridis",direction=-1,option="magma")+
	scale_color_continuous(type = "viridis",direction=-1,option="magma")+
	facet_wrap2("teloname",ncol=2,strip.position = "right",strip = strip)+
	ylab("RT")+
	xlab("Telomeric repeat length (bp)")+
	ggtitle("Figure S19",subtitle="fkh1∆")

quartz(file=path2fig %+% "FigS15-19.pdf",height=10,width=8,type="pdf")
plWT_all
plRIF1_all
plKU70_all
plCTF19_all
plFKH1_all
dev.off()


### fig4B
toplot_red <- telo_data4 %>% filter(meanB>0.02 & seqfilter2<0.1 & Rep !="noBrdU") %>%
	filter(Chromosome %in% c("chrVII","chrVIII","chrIX","chrX"))

strip_col <- strip_col0[13:20]
strip <- strip_themed(background_y = elem_list_rect(fill = strip_col))

plWT_red <- ggplot(toplot_red %>% filter(mut=="WT"))+
	geom_pointdensity(aes(x=seqlen,y=meanBscaled),shape=16,size=1)+
	geom_text(data=toplot_red %>% filter(mut=="WT") %>% group_by(teloname) %>% summarise(n=n()),aes(label="n=" %+% n),x=20,y=3,size=3)+
	stat_cor(aes(x=seqlen,y=meanBscaled),show.legend=F,label.x=450,label.y=c(3),method="spearman",cor.coef.name="rho",size=3)+
	coord_cartesian(ylim=c(0.5,3.5),xlim=c(0,600))+
	scale_fill_continuous(type = "viridis",direction=-1,option="magma")+
	scale_color_continuous(type = "viridis",direction=-1,option="magma")+
	facet_wrap2("teloname",ncol=2,strip.position = "right",strip = strip)+
	ylab("RT")+
	xlab("Telomeric repeat length (bp)")+
	ggtitle("wt")

plRIF1_red <- ggplot(toplot_red %>% filter(mut=="RIF1"))+
	geom_pointdensity(aes(x=seqlen,y=meanBscaled),shape=16,size=1)+
	geom_text(data=toplot_red %>% filter(mut=="RIF1") %>% group_by(teloname) %>% summarise(n=n()),aes(label="n=" %+% n),x=150,y=3,size=3)+
	stat_cor(aes(x=seqlen,y=meanBscaled),show.legend=F,label.x=1200,label.y=c(3),method="spearman",cor.coef.name="rho",size=3)+
	coord_cartesian(ylim=c(0.5,3.5))+
	scale_fill_continuous(type = "viridis",direction=-1,option="magma")+
	scale_color_continuous(type = "viridis",direction=-1,option="magma")+
	facet_wrap2("teloname",ncol=2,strip.position = "right",strip = strip)+
	ylab("RT")+
	xlab("Telomeric repeat length (bp)")+
	ggtitle("rif1∆")


plKU70_red <- ggplot(toplot_red %>% filter(mut=="KU70"))+
	geom_pointdensity(aes(x=seqlen,y=meanBscaled),shape=16,size=1)+
	geom_text(data=toplot_red %>% filter(mut=="KU70") %>% group_by(teloname) %>% summarise(n=n()),aes(label="n=" %+% n),x=25,y=3,size=3)+
	stat_cor(aes(x=seqlen,y=meanBscaled),show.legend=F,label.x=300,label.y=c(3),method="spearman",cor.coef.name="rho",size=3)+
	coord_cartesian(ylim=c(0.5,3.5))+
	scale_fill_continuous(type = "viridis",direction=-1,option="magma")+
	scale_color_continuous(type = "viridis",direction=-1,option="magma")+
	facet_wrap2("teloname",ncol=2,strip.position = "right",strip = strip)+
	ylab("RT")+
	xlab("Telomeric repeat length (bp)")+
	ggtitle("yku70∆")

plCTF19_red <- ggplot(toplot_red %>% filter(mut=="CTF19"))+
	geom_pointdensity(aes(x=seqlen,y=meanBscaled),shape=16,size=1)+
	geom_text(data=toplot_red %>% filter(mut=="CTF19") %>% group_by(teloname) %>% summarise(n=n()),aes(label="n=" %+% n),x=20,y=3,size=3)+
	stat_cor(aes(x=seqlen,y=meanBscaled),show.legend=F,label.x=450,label.y=c(3),method="spearman",cor.coef.name="rho",size=3)+
	coord_cartesian(ylim=c(0.5,3.5),xlim=c(0,600))+
	scale_fill_continuous(type = "viridis",direction=-1,option="magma")+
	scale_color_continuous(type = "viridis",direction=-1,option="magma")+
	facet_wrap2("teloname",ncol=2,strip.position = "right",strip = strip)+
	ylab("RT")+
	xlab("Telomeric repeat length (bp)")+
	ggtitle("ctf19∆")

plFKH1_red <- ggplot(toplot_red %>% filter(mut=="FKH1"))+
	geom_pointdensity(aes(x=seqlen,y=meanBscaled),shape=16,size=1)+
	geom_text(data=toplot_red %>% filter(mut=="FKH1") %>% group_by(teloname) %>% summarise(n=n()),aes(label="n=" %+% n),x=20,y=3,size=3)+
	stat_cor(aes(x=seqlen,y=meanBscaled),show.legend=F,label.x=450,label.y=c(3),method="spearman",cor.coef.name="rho",size=3)+
	coord_cartesian(ylim=c(0.5,3.5),xlim=c(0,600))+
	scale_fill_continuous(type = "viridis",direction=-1,option="magma")+
	scale_color_continuous(type = "viridis",direction=-1,option="magma")+
	facet_wrap2("teloname",ncol=2,strip.position = "right",strip = strip)+
	ylab("RT")+
	xlab("Telomeric repeat length (bp)")+
	ggtitle("fkh1∆")

pl4b <- plWT_red/plRIF1_red/plKU70_red/plCTF19_red/plFKH1_red


plo <- f4a/pl4b+
	plot_layout(heights=c(1,5))+
	plot_annotation(tag_levels = 'a',title="Figure 4")

quartz(file=paste0(path2fig,"Fig4.pdf"),height=20,width=12,type="pdf")
plo
dev.off()

### Generating Sup table
## add XY info to telo
telXY <- readRDS(path2data %+% "telXY.rds")

totable0 <- telo_data4 %>% left_join(telXY) %>%
	filter(meanB>0.02 & seqfilter2<0.1 & Rep!="noBrdU") %>%
	mutate(mut=factor(mut,levels=mut_order))%>%
	mutate(teloname=factor(teloname,levels=tel_order)) %>%
mutate(mut=fct_recode(mut,"rif1∆"="RIF1","yku70∆"="KU70","cft19∆"="CTF19","fkh1∆"="FKH1","wt"="WT"))

#group by mut
totable1 <- totable0 %>%
	group_by(mut) %>%
summarise(n=n(),rho=signif(cor(meanBscaled,seqlen,method="s",use="pairwise.complete.obs"),2),pval=signif(cor.test(meanBscaled,seqlen,method="s",use="pairwise.complete.obs")$p.value,2),TeloMedLen=signif(median(seqlen,na.rm=T),3),meanRT=signif(mean(meanBscaled,na.rm=T),3)) %>% dplyr::rename(strain=mut)

write_tsv(totable1,file=path2table %+% "tableS1.tsv")

#group by mut,xy
totable2 <- totable0 %>%
	group_by(mut,xy) %>%
summarise(n=n(),rho=signif(cor(meanBscaled,seqlen,method="s",use="pairwise.complete.obs"),2),pval=signif(cor.test(meanBscaled,seqlen,method="s",use="pairwise.complete.obs")$p.value,2),TeloMedLen=signif(median(seqlen,na.rm=T),3),meanRT=signif(mean(meanBscaled,na.rm=T),3)) %>% dplyr::rename(strain=mut)

write_tsv(totable2,file=path2table %+% "tableS2.tsv")

#group by mut,tel,xy
totable3 <- totable0 %>%
	group_by(mut,teloname,xy) %>%
summarise(n=n(),rho=signif(cor(meanBscaled,seqlen,method="s",use="pairwise.complete.obs"),2),pval=signif(cor.test(meanBscaled,seqlen,method="s",use="pairwise.complete.obs")$p.value,2),TeloMedLen=signif(median(seqlen,na.rm=T),3),meanRT=signif(mean(meanBscaled,na.rm=T),3)) %>% dplyr::rename(strain=mut)

write_tsv(totable3,file=path2table %+% "tableS3.tsv")
