#### Figure S1
#### fitting nanoT data to sortseq vs BrdU

suppressMessages(library(tidyverse))
library(ggpubr)
library(patchwork)
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

## fitting

B0 <- 5
data1 <- left_join(nanoT05,SortSeq_WT) %>%
	mutate(T=B0*(1/mean_br_bin-1)) %>%
	mutate(trel=2-timing)
data2 <- data1 #%>% filter(between(trel,0,1))
x <- data2$trel
y <- data2$T
fit2 <- nls(y ~ T0*exp(k*x), start = list(T0=10, k = 0.3))

fitter2 <- function(T0,k, x) {
	res <- T0*exp(k*x)
	return(res)
	}

param <- coef(fit2)
data2$fitted2 <- fitter2(param[1],param[2],data2$trel)
data_wt1_05 <- data2


B0 <- 10
data1 <- left_join(nanoT10,SortSeq_WT) %>%
	mutate(T=B0*(1/mean_br_bin-1)) %>%
	mutate(trel=2-timing)
data2 <- data1 #%>% filter(between(trel,0,1))
x <- data2$trel
y <- data2$T
fit2 <- nls(y ~ T0*exp(k*x), start = list(T0=10, k = 0.3))

fitter2 <- function(T0,k, x) {
	res <- T0*exp(k*x)
	return(res)
	}

param <- coef(fit2)
data2$fitted2 <- fitter2(param[1],param[2],data2$trel)
data_wt1_10 <- data2

B0 <- 20
data1 <- left_join(nanoT20,SortSeq_WT) %>%
	mutate(T=B0*(1/mean_br_bin-1)) %>%
	mutate(trel=2-timing)#
data2 <- data1 #%>% filter(between(trel,0,1))
x <- data2$trel
y <- data2$T
fit2 <- nls(y ~ T0*exp(k*x), start = list(T0=10, k = 0.3))

fitter2 <- function(T0,k, x) {
	res <- T0*exp(k*x)
	return(res)
	}

param <- coef(fit2)
data2$fitted2 <- fitter2(param[1],param[2],data2$trel)
data_wt1_20 <- data2

toplot1 <- bind_rows(data_wt1_05 %>% mutate(B0=5),data_wt1_10 %>% mutate(B0=10),data_wt1_20 %>% mutate(B0=20))%>% mutate(BrdU=factor(BrdU,levels=c("5 µM","10 µM","20 µM")))

p1 <- ggplot(toplot1)+
	geom_point(aes(x=trel, y=T,col=BrdU), size=0.2,alpha=0.1)+
	geom_line(aes(x=trel,y=fitted2,col=BrdU),linewidth=1)+
	scale_colour_manual("BrdU",values=mypal)+
	xlab("Relative replication timing")+
	ylab("dTTP (µM)")+
	coord_cartesian(xlim=c(0,1),ylim=c(5,50))


toplot1R <- toplot1 %>%
	mutate(pB_fit=B0/(B0+fitted2))

p2 <- ggplot(toplot1R)+
	geom_point(aes(x=trel, y=mean_br_bin,col=BrdU), size=0.2,alpha=0.1)+
	geom_line(aes(x=trel,y=pB_fit,col=BrdU),linewidth=1)+
	scale_colour_manual("BrdU",values=mypal)+
	xlab("Relative replication timing")+
	ylab("Mean BrdU content")+
	coord_cartesian(xlim=c(0,1))

pl <- p1/p2+
	plot_layout(guides="collect")+
	plot_annotation(title="Supplementary Figure 1",tag_levels="a")

ggsave(plot=pl,file=path2fig %+% "SupFig1_v1.pdf",h=6,w=5)

pl <- p1/p2+
	plot_layout(guides="collect")+
	plot_annotation(title="Supplementary Figure 1",tag_levels="a")&
	facet_wrap("BrdU")

ggsave(plot=pl,file=path2fig %+% "SupFig1_v1bis.pdf",h=6,w=8)


p1 <- ggplot(toplot1 %>% filter(B0==5))+
	geom_point(aes(x=trel, y=T,col=BrdU), size=0.2,alpha=0.1)+
	geom_line(aes(x=trel,y=fitted2,col=BrdU),linewidth=1)+
	scale_colour_manual("BrdU",values=mypal)+
	xlab("Relative replication timing")+
	ylab("dTTP (µM)")+
	coord_cartesian(xlim=c(0,1),ylim=c(5,50))

p2 <- ggplot(toplot1R %>% filter(B0==5))+
	geom_point(aes(x=trel, y=mean_br_bin,col=BrdU), size=0.2,alpha=0.1)+
	geom_line(aes(x=trel,y=pB_fit,col=BrdU),linewidth=1)+
	scale_colour_manual("BrdU",values=mypal)+
	xlab("Relative replication timing")+
	ylab("Mean BrdU content")+
	coord_cartesian(xlim=c(0,1))

pl <- p1/p2+
	plot_layout(guides="collect")+
	plot_annotation(title="Supplementary Figure 1",tag_levels="a")

ggsave(plot=pl,file=path2fig %+% "SupFig1_v2.pdf",h=6,w=5)
