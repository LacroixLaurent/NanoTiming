#### Figure S21
### Noise analysis
### using i vs i+1 dif with and without intraread correlation broken


suppressMessages(library(GenomicRanges))
suppressMessages(library(tidyverse))
library(patchwork)
library(parallel)
RNGkind("L'Ecuyer-CMRG")
theme_set(theme_bw())
mypal <- c(paletteer::paletteer_d("ggthemes::Classic_20"),"grey40","black","gold","greenyellow")

`%+%` <- paste0

setwd("/Users/ll/work/RStudioProjects/NanoTiming")
path2fig <- "Figures_Article/Figures_pdf/"
path2data <- "Figures_Article/Figures_data/"

source("Figures_Article/rescaling_function.r")

ncores <- 8L

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
bingen <- tileGenome(seqinf,tilewidth=bs, cut.last.tile.in.chrom=T)

backgr <- as_tibble(bingen) %>% dplyr::rename(positions=start,chrom=seqnames)%>% mutate(y=NA) %>% select(chrom,positions,y)

# load data
alldata <- do.call(bind_rows,lapply(1:24, function(i) readRDS(path2data %+% "nanoT_WT_24rep"%+%i%+%"_alldata.rds")))
alldata2 <- alldata %>% mutate(readlength=end-start+1)
sum(alldata2$readlength)/genome_size
# [1] 3367.668
target_list <- c(1,3,10,20,30,50,70,100,200,300,600,1000,3000)

### generate data while breaking the intrareads correlation

set.seed(123)
bigres2 <- mclapply(1:10, function(i) {
	totest <- alldata2 %>%
		sample_frac() %>%
		mutate(sumcov=cumsum(readlength)/genome_size)
	totestb <- alldata2 %>%
		sample_frac() %>%
		mutate(sumcov=cumsum(readlength)/genome_size)
	res <- lapply(target_list, function(target_cov) {
		totest2 <- totest %>% filter(sumcov<=target_cov)
		totest2b <- totestb %>% filter(sumcov<=target_cov)
		nanot <- totest2 %>%
			unnest(cols = c(signalr))%>%
			filter(signal>0.02)%>%
			group_by(chrom,positions)%>%
			summarise(mean_br_bin=mean(signal),.groups="drop") %>%
			ungroup %>%
			mutate(mod=1+myscaling0(mean_br_bin,infq=0.005,supq=0.995)) %>%
			select(-mean_br_bin)
		nanotb <- totest2b %>%
			unnest(cols = c(signalr))%>%
			filter(signal>0.02)%>%
			group_by(chrom,positions)%>%
			summarise(mean_br_bin=mean(signal),.groups="drop") %>%
			ungroup %>%
			mutate(modb=1+myscaling0(mean_br_bin,infq=0.005,supq=0.995)) %>%
			select(-mean_br_bin)
		noise_totest <- full_join(nanot,nanotb,by = join_by(chrom, positions)) %>%
			full_join(backgr,by = join_by(chrom, positions)) %>%
			arrange(chrom,positions) %>%
			select(-y) %>%
			group_by(chrom) %>%
			mutate(dif=mod-lead(modb)) %>%
			ungroup %>%
			mutate(gen_cov=target_cov)
		noisevar <- noise_totest %>% pull(dif) %>% var(na.rm=T)
		nbin <- noise_totest %>% reframe(nb=sum(!is.na(dif)))
		resu <- tibble(gen_cov=target_cov,noisevar=noisevar,nbin=nbin,i=i)
		return(resu)
		})
	},mc.cores=ncores,mc.set.seed=T)

bigres <- do.call(bind_rows,bigres2)
saveRDS(bigres,file=path2data %+% "Data_FigS21.rds")
#bigres <- readRDS(path2data %+% "Data_FigS21.rds")

### add sortseq/MFAseq data
## load data
SortSeqWT <- readRDS(path2data %+% "sortseq_WT.rds") %>%
	mutate(mod=1+myscaling0(timing,infq=0.005,supq=0.995))%>%
	mutate(chrom=factor(chrom,levels=chrom_order))


MFAseq <- readRDS(path2data %+% "MFAseq_WT.rds") %>%
	mutate(mod=1+myscaling0(ratio,infq=0.005,supq=0.995)) %>%
	mutate(chrom=factor(chrom,levels=chrom_order))

## noise
noise_sortseqWT <- SortSeqWT %>% select(-timing) %>%
	full_join(backgr,by = join_by(chrom, positions)) %>%
	arrange(chrom,positions) %>%
	select(-y) %>%
	group_by(chrom) %>%
	mutate(dif=mod-lead(mod)) %>%
	ungroup

noise_MFAseq <- MFAseq %>% select(-ratio) %>%
	full_join(backgr,by = join_by(chrom, positions)) %>%
	arrange(chrom,positions) %>%
	select(-y) %>%
	group_by(chrom) %>%
	mutate(dif=mod-lead(mod)) %>%
	ungroup

noisevarWT <- noise_sortseqWT %>% summarise(var=var(dif,na.rm=T)) %>% pull(var)
noisevarMFA <- noise_MFAseq %>% summarise(var=var(dif,na.rm=T)) %>% pull(var)

text_tb <- tibble(labl=factor(c("sort-seq","MFA-seq")),ypos=c(noisevarWT,noisevarMFA),gen_cov=c(1678,8740))

pl <- ggplot(bigres,aes(x=gen_cov,y=noisevar,group=gen_cov,col="Nanotiming")) +
	stat_summary(fun=median,geom="line",aes(group=1),linewidth = 0.2)+
	geom_boxplot(outlier.shape = NA,show.legend=F,linewidth = 0.2)+
	geom_point(data=text_tb,aes(x=gen_cov,y=ypos,col=labl), shape=16,size=1,show.legend=F)+
	scale_x_log10()+
	xlab("Genome coverage (x)")+
	ylab("Noise estimator")+
	scale_color_manual("",values=mypal[c(3,1,7)])+
	scale_y_log10()+
	guides(color = guide_legend(override.aes = list(linewidth = 1)))+
	ggtitle("Figure S21")
ggsave(plot=pl,file=path2fig %+% "FigS21.pdf",height=3,width=5)
