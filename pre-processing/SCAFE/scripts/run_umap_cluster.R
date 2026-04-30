library(Seurat)
library(optparse)


option_list = list(make_option(c("-o", "--output"), type = "character", default = NULL, help = "Path to output pdf", metavar = "character"),
                   make_option(c("-i", "--input"), type = "character", default = NULL, help = "Input seurat object .rds file path", metavar = "character"))

opt = parse_args(OptionParser(option_list = option_list))
input_path <- opt$input
output_path <- opt$output

seurat_object <- readRDS(input_path)

cluster_umap <- function(seurat_object, method){
        nDims = ncol(seurat_object[[method]])

        seurat_object <- RunUMAP(seurat_object, dims = 1:nDims, reduction = method)

        seurat_object <- FindNeighbors(seurat_object, dims = 1:nDims, reduction = method)
        seurat_object <- FindClusters(seurat_object)

        return(seurat_object)
}

seurat_object <- cluster_umap(seurat_object = seurat_object, method = names(seurat_object@reductions)[1])

saveRDS(seurat_object, file = output_path)
