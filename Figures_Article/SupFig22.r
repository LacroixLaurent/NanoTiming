#### Supplementary Figure 22
### comparison nanotrace vs sortseq vs nreads nanot or coverage

suppressMessages(library(GenomicRanges))
suppressMessages(library(rtracklayer))
suppressMessages(library(tidyverse))
library(patchwork)
library(ggprism)
library(parallel)
RNGkind("L'Ecuyer-CMRG")
theme_set(theme_bw())
mypal <- c(paletteer::paletteer_d("ggthemes::Classic_20"),"grey40","black","gold","greenyellow")
`%+%` <- paste0
setwd("/Users/ll/work/Ori/nanotiming/BioRxivPaper")
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


nanot_alldata <-  readRDS(path2data %+% "nanoT_WT_rep1_alldata.rds") %>% mutate(chrom=factor(chrom,levels=chrom_order))

SortSeq <- readRDS(path2data %+% "sortseq_WT.rds")

#### do the different subsampling
samp_size <- c(1000,5000,10000,30000,50000,80000,100000,150000,200000,300000)
nb_sampling <- 100
set.seed(37)
samp_data <- lapply(samp_size, function(samp) {
	mclapply(1:nb_sampling, function(i) {
		nano_data <- nanot_alldata %>%
		slice_sample(n=samp,replace=T) %>%
		unnest(cols = c(signalr)) %>%
		filter(signal>0.02) %>%
		group_by(chrom,positions) %>%
		summarise(mean_br_bin=mean(signal), .groups="drop") %>%
		mutate(iter=i) %>%
		mutate(samp_size=samp) %>% ungroup %>%
		mutate(mod=1+myscaling0(mean_br_bin,infq=0.005,supq=0.995)) %>%
		ungroup
		return(nano_data)
		},mc.cores=8L,mc.set.seed=T)
	})

nanot2 <- do.call(bind_rows,samp_data)

### do the grouping
input <- nanot2 %>%
	group_by(samp_size,iter)

### do spearman with sortseq
tocor <- full_join(input,SortSeq,join_by(chrom, positions)) %>%
	summarise(sp_cor=cor(mod,timing,use="pairwise.complete.obs",method="s"))
saveRDS(tocor,file=path2data %+% "Data_FigS22a.rds")
tocor <- readRDS(path2data %+% "Data_FigS22a.rds")

p1 <- ggplot(tocor,aes(x=samp_size,y=sp_cor,group=samp_size)) +
	stat_summary(fun=median,geom="line",aes(group=1),col="grey",linewidth = 0.2)+
	geom_boxplot(outlier.shape = NA,width=5000,linewidth = 0.2)+
	geom_hline(aes(yintercept=0.9), linetype=2)+
	scale_x_continuous(
			guide = "prism_minor",
			labels=scales::unit_format(big.mark ="",suffix="k",scale=1e-3,sep=" "))+
	xlab("Sample size")+
	ylab("Spearman's rank correlation coefficient")+
	coord_cartesian(ylim=c(0.35,1))


### try with coverage_x
alldata2 <- nanot_alldata %>% mutate(readlength=end-start+1)
sum(alldata2$readlength)/genome_size
#[1] 697.5959
target_list <- c(1,3,10,20,30,50,100,300,600)
set.seed(123)

bigres <- mclapply(1:100, function(i) {

	totest <- alldata2 %>% sample_frac() %>% mutate(sumcov=cumsum(readlength)/genome_size)
	res <- sapply(target_list, function(target_cov) {
		totest2 <- totest %>% filter(sumcov<=target_cov)
		nanot <- totest2 %>%
				unnest(cols = c(signalr))%>%
				filter(signal>0.02)%>%
				group_by(chrom,positions)%>%
				summarise(mean_br_bin=mean(signal),.groups="drop") %>%
				ungroup %>%
				mutate(mod=1+myscaling0(mean_br_bin,infq=0.005,supq=0.995))
		mycor <- full_join(nanot,SortSeq,join_by(chrom, positions)) %>%
	summarise(sp_cor=cor(mod,timing,use="pairwise.complete.obs",method="s"),.groups="drop")
		return(mycor)
		})
	tibble(target_list,res %>% unlist,i)
	},mc.cores=8L,mc.set.seed=T)

bigres2 <- do.call(bind_rows,bigres) %>% rename(spcor=2)
saveRDS(bigres2,file=path2data %+% "Data_FigS22b.rds")
bigres2 <- readRDS(path2data %+% "Data_FigS22b.rds")

p2 <- ggplot(bigres2,aes(x=target_list,y=spcor,group=target_list)) +
	geom_hline(aes(yintercept=0.9), linetype=2)+
	stat_summary(fun=median,geom="line",aes(group=1),col="grey",linewidth = 0.2)+
	geom_boxplot(outlier.shape = NA,linewidth = 0.2)+
	xlab("Genome coverage (x)")+
	ylab("Spearman's rank correlation coefficient")+
	scale_x_log10()+
	coord_cartesian(ylim=c(0.35,1))


pcomb <- p1+p2 + plot_annotation(tag_levels="a",title="Supplementary Figure 22") + plot_layout(axis_titles="collect")
ggsave(plot=pcomb,file=path2fig %+% "SupFig22.pdf",height=4,width=10)
