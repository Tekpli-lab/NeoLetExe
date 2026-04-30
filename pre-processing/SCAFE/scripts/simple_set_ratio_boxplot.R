library(Seurat)
library(Matrix)
library(optparse)
library(tidyr)
library(tidyverse)
library(scales)
library(ggpubr)
library(ggrepel)

option_list = list(make_option(c("-o", "--output_dir"), type = "character", default = NULL, help = "Output file", metavar = "character"),
		   make_option(c("-t", "--tf_name"), type = "character", default = "TF", help = "Output file", metavar = "character"),
		   make_option(c("-p", "--pos_bed_files"), type = "character", help = "", metavar = "character"),
		   make_option(c("-n", "--neg_bed_files"), type = "character", default = character(), help = "", metavar = "character"),
		   make_option(c("-s", "--seurat_object"), type = "character", default = NULL, help = "Input seurat .RDS object", metavar = "character"),
		   make_option(c("-c", "--patient_split"), type = "character", default = "response", help = "", metavar = "character"),
		   make_option(c("-m", "--main_id"), type = "character", default = NULL, help = "Input seurat .RDS object", metavar = "character"))

opt = parse_args(OptionParser(option_list = option_list))

output_dir <- opt$output_dir
seurat_object_path <- opt$seurat_object
main_id <- opt$main_id
tf_name <- opt$tf_name
patient_split <- opt$patient_split
seurat_object <- readRDS(seurat_object_path)
pos_bed_files <- unlist(strsplit(opt$pos_bed_files, " "))

if(opt$neg_bed_files == ""){ neg_bed_files <- character()
} else { neg_bed_files <- unlist(strsplit(opt$neg_bed_files, " ")) }

seurat_object@assays$SCAFE@counts@x <- rep(1, length(seurat_object@assays$SCAFE@counts@x))

bed_files_to_CRE_set <- function(bed_files, union = FALSE) {
	if ( length(bed_files) == 0) {return( character() )}
	CRE_list <- lapply(bed_files, 
			   function(bed) {
				   df <- read.delim(bed, sep = "\t", header = FALSE)
				   CRE <- paste(df$V1, df$V2, df$V3, sep = "_")
				   return(CRE)
			   })
	if (union) {CRE_set <- Reduce(union, CRE_list)
	} else {CRE_set <- Reduce(intersect, CRE_list)}
	return(CRE_set)
}

CRE_pos <- bed_files_to_CRE_set(pos_bed_files)
CRE_neg <- bed_files_to_CRE_set(neg_bed_files, union = TRUE)


CRE_id_map <- tibble(CREID = rownames(seurat_object)) %>% mutate(CRE_ID = CREID %>% str_extract("^[^-]+-[^-]+-[^-]+") %>% str_replace_all("-", "_"))

CRE_pos_mapped <- CRE_id_map %>% slice(match(CRE_pos, CRE_ID)) %>% pull(CREID)
CRE_neg_mapped <- CRE_id_map %>% slice(match(CRE_neg, CRE_ID)) %>% pull(CREID)

calculate_CRE_set_activity_per_cell <- function(seurat_obj, pos_set, neg_set = character()){
        all_cre <- rownames(seurat_obj)
	exclude_cre <- intersect(pos_set, neg_set)
        used_cre <- all_cre[!(all_cre %in% exclude_cre)]
	pos_cre <- pos_set[pos_set %in% used_cre]

	seurat_sub <- subset(seurat_obj, features = used_cre)
	seurat_pos <- subset(seurat_obj, features = pos_cre)

	active_used <- Matrix::colSums(seurat_sub@assays$SCAFE@counts)
	active_pos <- Matrix::colSums(seurat_pos@assays$SCAFE@counts)
	df <- cbind(seurat_sub@meta.data, active_set = active_pos, active_CRE = active_used)
	
	return(df)
}


df_cell_activity <- calculate_CRE_set_activity_per_cell(seurat_object, CRE_pos_mapped, CRE_neg_mapped)

df_cell_activity$patient_number <- factor(df_cell_activity$patient_number)
df_cell_activity$time_point <- factor(df_cell_activity$time_point, ordered = TRUE, levels = c("A", "B", "C"))

data.table::fwrite(df_cell_activity, file = file.path(output_dir, paste0(patient_split, "__", main_id, "__", tf_name, "__cell_activity.csv")), row.names = TRUE, sep = ",")


df_meta_group <- tibble(df_cell_activity) %>% filter(main.ids.cor == main_id) %>% group_by(orig.ident) %>% summarize(active_set = sum(active_set), active_CRE = sum(active_CRE), time_point = first(time_point), patient_number = first(patient_number), nCells = n(), group = first(!!sym(patient_split)))

saveRDS(as.data.frame(df_meta_group), file = file.path(output_dir, paste0(patient_split, "__", main_id, "__", tf_name, ".rds")))

groups <- unique(df_meta_group %>% pull(group))

pdf(file.path(output_dir, paste0(patient_split, "__", main_id, "__", tf_name, ".pdf")))

comparisons = list(c("A", "B"), c("B", "C"), c("A", "C"))

for (g in groups) {
	print(ggplot(df_meta_group %>% filter(group==g), aes(x = time_point, y =  active_set / active_CRE))  + geom_boxplot() + geom_point(size = 2, alpha = 0.7) + geom_line(aes(group = patient_number, color = patient_number), alpha = 0.7, size = 0.4) + stat_compare_means(method = "t.test", paired = FALSE, comparisons = comparisons) + geom_label_repel(data = df_meta_group %>% filter(group==g) %>% group_by(patient_number) %>% summarize(x = last(time_point), y = last(active_set / active_CRE), col = last(patient_number)), aes(x = x, y = y, label = patient_number, color = col), size = 2, hjust = 1, nudge_x = 0.2) + labs(color = g, y = paste0(tf_name, " / active enhancers"))) 
}

dev.off()