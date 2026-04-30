library(Seurat)
library(Matrix)
library(optparse)


option_list = list(make_option(c("-o", "--outDir"), type = "character", default = NULL, help = "Output directory", metavar = "character"),
		   make_option(c("-i", "--input"), type = "character", default = NULL, help = "Cell type name", metavar = "character"))

opt = parse_args(OptionParser(option_list = option_list))
input_path <- opt$input
outDir <- opt$outDir

seurat_object <- readRDS(input_path)

assay <- DefaultAssay(seurat_object)

write(rownames(seurat_object[[assay]]@counts), file = file.path(outDir, "rownames.txt"))
write(colnames(seurat_object[[assay]]@counts), file = file.path(outDir, "colnames.txt"))
writeMM(seurat_object[[assay]]@counts, file = file.path(outDir, "count_mat.mtx"))

