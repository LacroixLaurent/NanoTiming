### Figure 1b
### plot BRdU content in individual reads

suppressMessages(library(tidyverse))
library(patchwork)
library(ggprism)
theme_set(theme_bw())
`%+%` <- paste0

setwd("/Users/ll/work/RStudioProjects/NanoTiming")
path2fig <- "Figures_Article/Figures_pdf/"
path2data <- "Figures_Article/Figures_data/"

nanot_data <- bind_rows(
	readRDS(path2data %+% "nanoT_WT_rep1_alldata.rds"),
	readRDS(path2data %+% "nanoT_WT_10_alldata.rds"),
	readRDS(path2data %+% "nanoT_WT_20_alldata.rds"),
	readRDS(path2data %+% "nanoT_WT_100_alldata.rds"),
	readRDS(path2data %+% "nanoT_WT_1000_alldata.rds")
	)

# Figure1b BrdU variation in reads
# filtering non substituted reads

set.seed(37)

toplot <- nanot_data %>%
	mutate(dose=as_factor(dose))%>%
	mutate(length=end-start+1)%>%
	filter(length>75000)%>%
	unnest(cols = c(signalr))%>%
	group_by(read_id)%>%
	mutate(median_b_read=median(signal,na.rm=T))%>%
	mutate(mean_b_read=mean(signal,na.rm=T))%>%
	filter(median_b_read> 0.02)%>%
	ungroup()%>%
	nest(cols = c(positions,signal))%>%
#subset of 50 reads per dose
	group_by(dose)%>%
	sample_n(50)%>%
	arrange(median_b_read)%>%
	mutate(track_i=row_number())%>%
	unnest(cols = c(cols))%>%
	mutate(x=positions-start) %>%
	ungroup()

f1b <- ggplot(toplot)+
	geom_line(aes(x=x,y=track_i,group=paste(read_id),color=signal),linewidth=1)+
	facet_wrap(~dose,ncol=5)+
	scale_color_viridis_c( option = "C",limits=c(0,1))+
	labs(x="Position along the read (kb)", y="", color="BrdU content")+
	scale_x_continuous(guide = "prism_minor",
		labels=scales::unit_format(big.mark ="",suffix="",scale=1e-3,sep=""),
		expand=expansion(mult=c(0,0.01)))+
	scale_y_continuous(expand=expansion(mult=c(0.01,0.015)))+
	theme_bw()+
	theme(axis.text.y=element_blank(),axis.ticks.y=element_blank())+
	ggtitle("Figure 1b")
ggsave(plot=f1b, file=path2fig %+% "Fig1b.pdf",w=5.5,h=3)


f1b <- ggplot(toplot)+
#	geom_line(aes(x=x,y=track_i,group=paste(read_id),color=signal),linewidth=1)+
	ggrastr::rasterise(geom_line(aes(x=x,y=track_i,group=paste(read_id),color=signal),linewidth=0.8),dpi=1200)+
	facet_wrap(~dose,ncol=5)+
	scale_color_viridis_c( option = "C",limits=c(0,1))+
	labs(x="Position along the read (kb)", y="", color="BrdU content")+
	scale_x_continuous(guide = "prism_minor",
		labels=scales::unit_format(big.mark ="",suffix="",scale=1e-3,sep=""),
		expand=expansion(mult=c(0,0.01)))+
	scale_y_continuous(expand=expansion(mult=c(0.01,0.015)))+
	theme_bw()+
	theme(axis.text.y=element_blank(),axis.ticks.y=element_blank())
	
#ggsave(plot=f1b, file=path2fig %+% "Fig1b.pdf",w=5.5,h=3)
ggsave(plot=f1b, file=path2fig %+% "Fig1b_rast.pdf",w=5.5,h=3)
