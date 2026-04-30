library(Seurat)
library(dplyr)
library(optparse)


### COMMAND LINE OPTIONS ###
option_list = list(make_option(c("-o", "--outDir"), type = "character", default = NULL, help = "Output directory", metavar = "character"),
		   make_option(c("-f", "--nCount_threshold"), type = "character", default = NULL, help = "Min nCount", metavar = "character"),
		   make_option(c("-c", "--input_cell_type"), type = "character", default = NULL, help = "Cell type name", metavar = "character"))

opt = parse_args(OptionParser(option_list = option_list))
cell_type <- opt$input_cell_type
proximity <- "distal_cre"
outDir <- opt$outDir
nCount_threshold <- opt$nCount_threshold

####

seurat_object <- readRDS(file.path(outDir, "seurat_object_split", proximity, paste0(cell_type, "__", proximity, ".rds")))
metadata <- seurat_object@meta.data

improved_ids <- read.csv(file.path(outDir,"..", "data", "improved_ids.csv"), header = T, sep = ";")
metadata[,"cell_ids"] <- rownames(metadata)
metadata <- left_join(metadata, improved_ids, by = "final.ids")
rownames(metadata) <- metadata[,"cell_ids"]

new_metadata <- read.csv2(file.path(outDir,"..", "data", "metadata_11_10_24_sample_names_fixed.csv"), header = T, row.names = 1)
metadata_join <- left_join(metadata, new_metadata[,c(colnames(new_metadata)[!(colnames(new_metadata) %in% colnames(metadata))], "cell_ids")], by = "cell_ids")
rownames(metadata_join) <- metadata_join$cell_ids

write.csv(metadata_join, file = file.path(outDir, "cell_type_metadata", paste0(cell_type, "_metadata.csv")), row.names = TRUE)

