#### Supplementary Figures 24 and 25
### comparison nanotrace 24 rep vs sortseq

suppressMessages(library(GenomicRanges))
suppressMessages(library(rtracklayer))
suppressMessages(library(tidyverse))
library(patchwork)
library(ggprism)
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

fullGFF <- import("Reference_Genome/BT1multiUra.gff3") %>% as_tibble()
ARS <- fullGFF %>% filter(Conf=="Confirmed",name!="ARS1216.5",type %in% c("ORI"))
CEN <- fullGFF %>% filter(type=="centromere")
rDNA <- fullGFF %>% filter(type=="rRNA")

### loading inputs

## Import sortseq data with scaling
SortSeqWT <- readRDS(path2data %+% "sortseq_WT.rds") %>%
	mutate(mod=1+myscaling0(timing,infq=0.005,supq=0.995))%>%
	mutate(filename="sort-seq") %>%
	mutate(mut="WT") %>%
	mutate(positions=positions+500) %>%
	left_join(.,chrom_sizes, by="chrom") %>%
	mutate(positions=map2_dbl(positions,seqlengths, function(x,y) {res=ifelse(x<y,x,y); return(res)})) %>%
	mutate(chrom=factor(chrom,levels=chrom_order)) %>% select(-seqlengths)


### load WT_24rep data
nanot24_2 <- readRDS(path2data %+% "nanoT_WT_24rep.rds") %>%
	group_by(Rep) %>%
	mutate(mod=1+myscaling0(mean_br_bin,infq=0.005,supq=0.995)) %>%
	ungroup %>%
	mutate(chrom=factor(chrom,levels=chrom_order)) %>%
	ungroup %>%
	mutate(positions=positions+500)

#### add a capping to chrom size
nanot24_3 <- left_join(nanot24_2,chrom_sizes, by="chrom") %>% mutate(positions=map2_dbl(positions,seqlengths, function(x,y) {res=ifelse(x<y,x,y); return(res)}))

### nanotraces
ars2plot <- ARS %>% rename(chrom=seqnames)
cen2plot <- CEN %>% rename(chrom=seqnames)
rdna2plot <- rDNA %>% rename(chrom=seqnames) %>% group_by(chrom,type) %>% summarise(start=min(start),end=max(end))

mypal3 <- c(paletteer::paletteer_d("ggthemes::Classic_20"),"grey40","gold","darkblue","greenyellow","black")

input2 <- bind_rows(nanot24_3,SortSeqWT %>% rename(Rep=filename)) %>%
	mutate(Rep=factor(Rep,levels=c("rep" %+% 1:24,"sort-seq")))
miny=0.8
maxy=2.2

pl <- ggplot(input2)+
	geom_vline(data=ars2plot,aes(xintercept=(start+end)/2),col=mypal[9],show.legend=F,linewidth=0.2)+
	geom_vline(data=cen2plot,aes(xintercept=(start+end)/2),col=mypal[5],show.legend=F)+
	geom_rect(data=rdna2plot,aes(xmin=start,xmax=end,ymin=0.9,ymax=1),col=NA,fill=mypal[21],show.legend=F,alpha=0.6)+
	geom_point(aes(x=positions,y=mod,col=Rep,group=Rep),shape=16,size=0.2)+
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
	guides(color = guide_legend(nrow=2,override.aes = list(size = 2)))+
	ggtitle("Supplementary Figure 24")

ggsave(plot=pl,file=paste0(path2fig,"SupFig24.pdf"),height=14,width=12)

### Corr
library(ggcorrplot)
toplot1 <- nanot24_3  %>% select(c(1,2,6,7)) %>% pivot_wider(values_from=mod, names_from=Rep, names_prefix="nanoT_",values_fn = mean)
tocor <- left_join(toplot1,SortSeqWT %>% select(c(1,3,4)) %>% rename("sort-seq"=mod))
cormat1 <- sapply(3:27, function(i) sapply((3:27), function(j) cor(tocor[,i],tocor[,j],use="pairwise.complete.obs",method="s")))
colnames(cormat1) <- rownames(cormat1) <- names(tocor)[3:27]
ggcorrplot(cormat1,lab=T,lab_size=4,digits=3)+
	scale_fill_gradient2(limit=c(0,1),low="blue",high="red",mid="white",midpoint=0.5)+
	ggtitle("Supplementary Figure 25",subtitle="Spearman's rank correlation coefficients")
ggsave(path2fig %+% "SupFig25.pdf",height=20,width=20)
