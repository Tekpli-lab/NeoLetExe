library(Seurat)
library(Matrix)
library(optparse)

### COMMAND LINE OPTIONS ###
option_list = list(
    make_option(c("--input_mat"), type = "character", default = NULL, help = "Path to input matrix (.mtx file)", metavar = "character"),
    make_option(c("--input_cell_names"), type = "character", default = NULL, help = "Path to cell names file", metavar = "character"),
    make_option(c("--input_dim_names"), type = "character", default = NULL, help = "Path to dimension names file", metavar = "character"),
    make_option(c("--input_seurat"), type = "character", default = NULL, help = "Path to input Seurat object (.rds file)", metavar = "character"),
    make_option(c("--output_seurat"), type = "character", default = NULL, help = "Path to output Seurat object (.rds file)", metavar = "character")
)

opt = parse_args(OptionParser(option_list = option_list))

# Load inputs
seurat_object <- readRDS(opt$input_seurat)
cell_names <- readLines(opt$input_cell_names)
dim_names <- readLines(opt$input_dim_names)
dim_reduc <- Matrix::readMM(opt$input_mat)

colnames(dim_reduc) <- cell_names
rownames(dim_reduc) <- dim_names

method <- "LDA_harmony"
seurat_object[[method]] <- CreateDimReducObject(embeddings = t(as.matrix(dim_reduc)), key = paste0(method, "_"), assay = DefaultAssay(seurat_object))
saveRDS(seurat_object, opt$output_seurat)
