library(optparse)
library(Seurat)

option_list = list(make_option(c("-m", "--sc_metadata"), type = "character", default = NULL, help = "Metadata file path", metavar = "character"),
		   make_option(c("-o", "--outDir"), type = "character", default = NULL, help = "Output directory", metavar = "character"),
		   make_option(c("-d", "--dataDir"), type = "character", default = NULL, help = "Directory with sample counts", metavar = "character"),
		   make_option(c("-s", "--sample_names_list"), type = "character", default = NULL, help = "String of sample names, each sample separated by space", metavar = "character"))

opt = parse_args(OptionParser(option_list = option_list))

sample_names <- strsplit(opt$sample_names_list, " ")[[1]]
outDir <- opt$outDir
dataDir <- opt$dataDir
metadata <- read.csv(opt$sc_metadata, header = T, row.names = 1)

sample_objects <- unlist(lapply(sample_names, function(x){
					mat <- Read10X(file.path(dataDir, "counts", x, "matrix"));
					seurat_object <- CreateSeuratObject(counts = mat);
					seurat_object[["orig.ident"]] <- x;
					return(seurat_object)
}))

merged <- merge(x = sample_objects[[1]], y = sample_objects[-1], add.cell.ids = sample_names)

### Make sure that the same cells are in both the single cell metadata and the scafe counts
filtered_cells <- colnames(merged)[(colnames(merged) %in% rownames(metadata))]
print("Cells removed from Scafe counts, because they werent in the single cell metadata. (Potentially because of QC on the single cell data, or because of different sample naming)")
print(ncol(merged) - length(filtered_cells))
merged <- subset(merged, cells = filtered_cells)

merged@meta.data <- metadata[colnames(merged),]
saveRDS(merged, file = file.path(outDir, "seurat_object_unprocessed", "all_cells.rds"))

