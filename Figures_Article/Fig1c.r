#### Figure 1c
#### comparing nanoT data to sortseq vs BrdU

suppressMessages(library(tidyverse))
library(ggpubr)
theme_set(theme_bw())
mypal <- c(paletteer::paletteer_d("ggthemes::Classic_20"),"grey40","black")
`%+%` <- paste0

setwd("/Users/ll/work/RStudioProjects/NanoTiming")
path2fig <- "Figures_Article/Figures_pdf/"
path2data <- "Figures_Article/Figures_data/"

### Import Sortseq data

SortSeq_WT <- readRDS(path2data %+% "sortseq_WT.rds")

#### Import nanoT data

nanoT05 <- readRDS(path2data %+% "nanoT_WT_rep1.rds")
nanoT10 <- readRDS(path2data %+% "nanoT_WT_10.rds")
nanoT20 <- readRDS(path2data %+% "nanoT_WT_20.rds")
nanoT100 <- readRDS(path2data %+% "nanoT_WT_100.rds")
nanoT1000 <- readRDS(path2data %+% "nanoT_WT_1000.rds")

### Plot
toplot_nano <- bind_rows(nanoT05,nanoT10,nanoT20,nanoT100,nanoT1000) %>% rename(raw=mean_br_bin)

toplot <- full_join(toplot_nano,SortSeq_WT) %>% mutate(BrdU=factor(BrdU,levels=c("5 µM","10 µM","20 µM","100 µM","1000 µM")))

f1c <- ggplot(toplot,aes(x=timing,y=raw,col=BrdU))+
	geom_point(size=0.2,alpha=0.2,shape=16)+
	scale_color_manual(values=mypal,guide = guide_legend(reverse = TRUE) )+
	geom_smooth(method="lm",se=F)+
	stat_cor(label.x=2.1,label.y = c(0,0.05,0.1,0.15,0.2),method="spearman",cor.coef.name="rho",show.legend=F,aes(label=after_stat(r.label)))+
	ylab("Mean BrdU content")+
	xlab("Relative copy number (from sort-seq)")+
	ggtitle("Figure 1c")

ggsave(plot=f1c,file=path2fig %+% "Fig1c.pdf",h=3,w=4.5)


f1c <- ggplot(toplot,aes(x=timing,y=raw,col=BrdU))+
#	geom_point(size=0.2,alpha=0.2,shape=16)+
	ggrastr::rasterise(geom_point(size=0.2,alpha=0.1,shape=16),dpi=600)+
	scale_color_manual(values=mypal,guide = guide_legend(reverse = TRUE) )+
	geom_smooth(method="lm",se=F)+
	stat_cor(label.x=2.1,label.y = c(0,0.05,0.1,0.15,0.2),method="spearman",cor.coef.name="rho",show.legend=F,aes(label=after_stat(r.label)))+
	ylab("Mean BrdU content")+
	xlab("Relative copy number (from sort-seq)")
	
#ggsave(plot=f1c,file=path2fig %+% "Fig1c.pdf",h=3,w=4.5)
ggsave(plot=f1c,file=path2fig %+% "Fig1c_rast.pdf",h=3,w=4.5)

