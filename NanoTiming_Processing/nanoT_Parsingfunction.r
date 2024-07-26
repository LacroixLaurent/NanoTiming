parsing_BrdUbin <- function(bam.in,out.file,ncores=1L,bin.size=1000,savefile=T,min_len=1,B_thr=0.02,n_thr=0)
{
require(dplyr)
require(stringr)
require(future.apply)
# import from bam: read_id, flag, chrom,start,end,seq,pos,prob
data0 <- try(system(paste0("samtools view ",bam.in," | cut -f 1,2,3,4,9,10,13,14"),intern=T))

plan(multicore, workers = ncores)

data1 <- future_lapply(data0, function(x) 
	{
	input1 <- strsplit(x,"\t")[[1]]
	res <- tibble(read_id=character(), chrom=character(), start=integer(), end=integer(), strand=character(), signalr=list())
	if ((as.numeric(input1[5])-1)>=min_len) {
		start <- as.numeric(input1[4])
		end <- start+as.numeric(input1[5])-1
		id <- input1[1]
		if (input1[2]=="0") {strand="+"} else {strand="-"}
		chrom <- input1[3]
		if (strand=="+") {
				seq <- input1[6]
			} else {
				seq <- as.character(Biostrings::reverseComplement(Biostrings::DNAString(input1[6])))
			}
				# the Biostrings library is not loaded but should be installed.
		pos <- input1[7] %>% str_sub(6L,-2L) %>% str_split(";") %>% unlist
		pos01 <- str_sub(pos[1],5L,-1L)
		pos1 <- as.numeric(c(str_split(pos01,",")[[1]]))
		prob <- input1[8] %>% str_remove("Ml:B:C,")
		prob0 <- as.numeric(c(str_split(prob,",")[[1]]))/255
		prob1 <- prob0
		Tpos <- str_locate_all(seq,"T")[[1]][,1]
		pos_out <- Tpos[cumsum(pos1+1)]
		prob_out <- prob1
		signal0 <- tibble(pos_out,prob_out)
		if(strand=="+") {
				positions= start - 1 + signal0 %>% pull(pos_out);
				Bprob=signal0 %>% pull(prob_out)
			}
		# translate read coordinates to genomic coordinates and reverse the signal if strand==-
		if(strand=="-") {
				positions= end + 1 - signal0 %>% pull(pos_out);
				Bprob=signal0 %>% pull(prob_out)
			}
		out <- tibble(positions,Bprob) %>% arrange(positions)
		## final bining
		bs <- bin.size
		out2 <- out %>%
			mutate(positions = floor((positions-1)/bs)*bs+1) %>%
			group_by(positions) %>%
			summarise(signal = mean(Bprob,na.rm=T), .groups = "drop")
		## compiling
		res <- tibble(read_id=id,chrom,start,end,strand,signalr=list(out2))
		}
	})
data_treated <- do.call(bind_rows,data1) %>% arrange(chrom,start)

	if (savefile) {
		saveRDS(data_treated, file=paste0(out.file,"_parsed_BrdU_binned_",bin.size,".rds"))
	}else{
		return(data_treated)
	}
}
