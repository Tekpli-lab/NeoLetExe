library(Seurat)
library(dplyr)
library(tidyverse)
library(optparse)


split_proximity <- function(seurat_object, cre_annotation){
	metadata <- seurat_object@meta.data

        cre_annotation <- cre_annotation %>% mutate(proximal_distal = case_when(gene_promoter == "annotated" ~ "proximal",class == "distal" & gene_promoter != "unnanotated" ~ "distal",TRUE ~ "none"))	
	

	proximal_cre_rows <- rownames(seurat_object)[cre_annotation[,"proximal_distal"] == "proximal"]
	distal_cre_rows <- rownames(seurat_object)[cre_annotation[,"proximal_distal"] == "distal"]
	
	proximal_counts <- seurat_object@assays$RNA@counts[proximal_cre_rows,]
	distal_counts <- seurat_object@assays$RNA@counts[distal_cre_rows,]

	proximal_seurat <- CreateSeuratObject(counts = proximal_counts, assay = "SCAFE", meta.data = metadata)
	distal_seurat <- CreateSeuratObject(counts = distal_counts, assay = "SCAFE", meta.data = metadata)

	proximal_seurat[["SCAFE"]]@meta.features <- cre_annotation[proximal_cre_rows,]
	distal_seurat[["SCAFE"]]@meta.features <- cre_annotation[distal_cre_rows,]

	return(c(proximal_seurat, distal_seurat))
}	




### COMMAND LINE OPTIONS ###
option_list = list(make_option(c("-o", "--outDir"), type = "character", default = NULL, help = "Output directory", metavar = "character"),
		   make_option(c("-c", "--input_cell_type"), type = "character", default = NULL, help = "Cell type name", metavar = "character"),
		   make_option(c("-a", "--cre_annotation"), type = "character", default = "", help = "cre_info.tsv from SCAFE annotation", metavar = "character"))

opt = parse_args(OptionParser(option_list = option_list))
cell_type <- opt$input_cell_type
outDir <- opt$outDir
cre_annotation <- opt$cre_annotation

####

seurat_object <- readRDS(file.path(outDir, "seurat_object_unprocessed", paste0(cell_type, ".rds")))

cre_annotation <- as.data.frame(read.csv(cre_annotation, sep = "\t"))
cre_annotation$CREID <- gsub("_", "-", cre_annotation$CREID)
rownames(cre_annotation) <- cre_annotation$CREID
cre_annotation <- cre_annotation[rownames(seurat_object), ]


out_list <- split_proximity(seurat_object, cre_annotation)
## Sharing the nCount / nFeature information across the two objects
out_list[[1]]@meta.data["nCount_distal"] <- out_list[[2]]@meta.data["nCount_SCAFE"]
out_list[[1]]@meta.data["nFeature_distal"] <- out_list[[2]]@meta.data["nFeature_SCAFE"]
out_list[[2]]@meta.data["nCount_distal"] <- out_list[[2]]@meta.data["nCount_SCAFE"]
out_list[[2]]@meta.data["nFeature_distal"] <- out_list[[2]]@meta.data["nFeature_SCAFE"]

out_list[[1]]@meta.data["nCount_proximal"] <- out_list[[1]]@meta.data["nCount_SCAFE"]
out_list[[1]]@meta.data["nFeature_proximal"] <- out_list[[1]]@meta.data["nFeature_SCAFE"]
out_list[[2]]@meta.data["nCount_proximal"] <- out_list[[1]]@meta.data["nCount_SCAFE"]
out_list[[2]]@meta.data["nFeature_proximal"] <- out_list[[1]]@meta.data["nFeature_SCAFE"]

saveRDS(out_list[[1]], file = file.path(outDir, "seurat_object_split", "proximal", paste0(cell_type, "__proximal.rds")))
saveRDS(out_list[[2]], file = file.path(outDir, "seurat_object_split", "distal_cre", paste0(cell_type, "__distal_cre.rds")))	

