library(GenomicRanges)
#library(rtracklayer)
library(dplyr)
library(optparse)
library(tidyverse)


### COMMAND LINE OPTIONS ###
option_list = list(make_option(c("--output_dir"), type = "character", default = NULL, help = "Output file", metavar = "character"),
		   make_option(c("--LOLA_db"), type = "character", default = NULL, help = "", metavar = "character"),
		   make_option(c("--enrichment_tsv"), type = "character", default = NULL, help = "CRE_set", metavar = "character"),
		   make_option(c("--diff_cre"), type = "character", default = NULL, help = "", metavar = "character"),
		   make_option(c("--TF_list"), type = "character", default = NULL, help = "TF names, space separated ", metavar = "character")
)

opt = parse_args(OptionParser(option_list = option_list))
output_dir <- opt$output_dir
diff_cre_path <- opt$diff_cre
enrichment_path <- opt$enrichment_tsv
LOLA_db_path <- opt$LOLA_db
TF_list <- strsplit(opt$TF_list, " ")[[1]]

LOLA_db <- readRDS(LOLA_db_path)
enrichment <- read_delim(file = enrichment_path, delim = "\t", col_names = TRUE)
#diff_obj =  import(diff_cre_path)
df_diff <- read_delim(diff_cre_path,  delim = "\t", col_names =F)
colnames(df_diff) <- c("chrom", "start", "end")
diff_obj <- makeGRangesFromDataFrame(df_diff, keep.extra.columns = TRUE)

write_TFBS_bed <- function(TF, d_CRE, enrich = enrichment, GR_list = LOLA_db[[4]], out_dir = output_dir){
	db_set <- (enrich %>% filter(collection == TF))$dbSet
	gr <- GenomicRanges::reduce(unlist(GR_list[db_set])) 
	gr <- gr[nchar(as.character(seqnames(gr))) <= 5]
        TFBS_CREs <- subsetByOverlaps(d_CRE, gr)
	df <- data.frame(seqnames=seqnames(TFBS_CREs), starts=start(TFBS_CREs), ends=end(TFBS_CREs))
	write.table(df, file = file.path(out_dir, paste0(TF, ".bed")),  quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")
	#export(TFBS_CREs, con = file.path(out_dir, paste0(TF, ".bed")), format = "BED")
}

for(tf in TF_list){
	write_TFBS_bed(tf, diff_obj)
}



