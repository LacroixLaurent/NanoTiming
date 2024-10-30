#### Figure S2
#### pseudo nanoT data to sortseq vs BrdU for MCM899

suppressMessages(library(tidyverse))
library(ggpubr)
library(patchwork)
theme_set(theme_bw())
mypal <- c(paletteer::paletteer_d("ggthemes::Classic_20"),"grey40","black")
`%+%` <- paste0

setwd("/Users/ll/work/RStudioProjects/NanoTiming")
path2fig <- "Figures_Article/Figures_pdf/"
path2data <- "Figures_Article/Figures_data/"
path2data0 <- "/Users/ll/work/Ori/nanotiming/nanoT_data_BOBY/"

### Import Sortseq data

SortSeq_WT <- readRDS(path2data %+% "sortseq_WT.rds")

#### Import pseudo_nanoT data

BOdata <- bind_rows(readRDS(path2data0 %+% "BO_run14_RefBT1multi_nanoT_alldata.rds"),readRDS(path2data0 %+%"BO_runpercent_RefBT1multi_nanoT_alldata.rds"))
BPdata <- bind_rows(readRDS(path2data0 %+%"BP_run14_RefBT1multi_nanoT_alldata.rds"),readRDS(path2data0 %+%"BP_runpercent_RefBT1multi_nanoT_alldata.rds"))
BQdata <- bind_rows(readRDS(path2data0 %+%"BQ_run14_RefBT1multi_nanoT_alldata.rds"),readRDS(path2data0 %+%"BQ_runpercent_RefBT1multi_nanoT_alldata.rds"))
BRdata <- bind_rows(readRDS(path2data0 %+%"BR_run14_RefBT1multi_nanoT_alldata.rds"),readRDS(path2data0 %+%"BR_runpercent_RefBT1multi_nanoT_alldata.rds"))
BSdata <- bind_rows(readRDS(path2data0 %+%"BS_run14_RefBT1multi_nanoT_alldata.rds"),readRDS(path2data0 %+%"BS_runpercent_RefBT1multi_nanoT_alldata.rds"))
BTdata <- bind_rows(readRDS(path2data0 %+%"BT_run14_RefBT1multi_nanoT_alldata.rds"),readRDS(path2data0 %+%"BT_runpercent_RefBT1multi_nanoT_alldata.rds"))
BUdata <- bind_rows(readRDS(path2data0 %+%"BU_run14_RefBT1multi_nanoT_alldata.rds"),readRDS(path2data0 %+%"BU_runpercent_RefBT1multi_nanoT_alldata.rds"))
BVdata <- bind_rows(readRDS(path2data0 %+%"BV_run14_RefBT1multi_nanoT_alldata.rds"),readRDS(path2data0 %+%"BV_runpercent_RefBT1multi_nanoT_alldata.rds"))
BWdata <- bind_rows(readRDS(path2data0 %+%"BW_run14_RefBT1multi_nanoT_alldata.rds"),readRDS(path2data0 %+%"BW_runpercent_RefBT1multi_nanoT_alldata.rds"))
BXdata <- bind_rows(readRDS(path2data0 %+%"BX_run14_RefBT1multi_nanoT_alldata.rds"),readRDS(path2data0 %+%"BX_runpercent_RefBT1multi_nanoT_alldata.rds"))
BYdata <- bind_rows(readRDS(path2data0 %+%"BY_run14_RefBT1multi_nanoT_alldata.rds"),readRDS(path2data0 %+%"BY_runpercent_RefBT1multi_nanoT_alldata.rds"))

data_list <- list(BOdata,BPdata,BQdata,BRdata,BSdata,BTdata,BUdata,BVdata,BWdata,BXdata,BYdata)

## compute nanotData
pB_list <- seq(0,100,10)
data_nanoT <- lapply(seq_along(data_list), function(i) {
	res1 <- data_list[[i]] %>%
		unnest(cols = c(signalr))%>%
		filter(signal>0.02)%>%
		group_by(chrom,positions)%>%
		summarise(mean_br_bin=mean(signal),nf=n()) %>%
		ungroup
	res2 <- data_list[[i]] %>%
		unnest(cols = c(signalr))%>%
		group_by(chrom,positions)%>%
		summarise(n0=n()) %>%
		ungroup
	full_join(res1,res2) %>% mutate(pB=pB_list[i] %>% as.factor)
})


toplot <- do.call(bind_rows,data_nanoT) %>% left_join(SortSeq_WT)


pl <- ggplot(toplot,aes(x=timing,y=mean_br_bin,col=pB))+
	geom_point(size=0.2,alpha=0.3,shape=16)+
	ylab("Mean BrdU content")+
	xlab("Relative copy number (from sort-seq)")+	scale_color_manual("BrdU (%)",values=mypal)+
	guides(color = guide_legend(override.aes = list(size = 2,alpha=1)))+
	coord_cartesian(xlim=c(0.8,2.2))+
	ggtitle("Supplementary Figure 2")

ggsave(plot=pl,file=path2fig %+% "SupFig2.pdf",h=5,w=6)
