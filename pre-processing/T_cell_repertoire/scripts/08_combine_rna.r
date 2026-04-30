library(Seurat)
library(tidyverse)
library(scRepertoire)
library(argparse)

combine_rna <- function(vdj, seurat, output, meta_out, clone_call, align_tool) {
  load(seurat)
  print("Seurat object loaded")
  vdj <- readRDS(vdj)
  
  # Handle clone call
  call <- clone_call
  call <- gsub("CT", "", call)


  ## Fix discrepancy in barcodes between two objects
  vdj <- lapply(vdj, function(x) {
    if (nrow(x) > 0) {
      if (align_tool == "mixcr") {
        x$barcode <- paste0(x$barcode, "-1")
      }
    }
    return(x)
  })

  ## Combine two objects
  seurat <- combineExpression(
    vdj,
    Scaled,
    cloneCall = call,
    chain = "both",
    group.by = "sample",
    proportion = TRUE
  )
  print("Objects are combined")

  # save new object and meta data
  metadata <- seurat@meta.data

  saveRDS(metadata, meta_out)
  saveRDS(seurat, output)
}

parser <- ArgumentParser(description = "Find shared clones across patients")
parser$add_argument("--vdj", help = "scRep object")
parser$add_argument("--seurat", help = "Seurat object, either T or B cells filtered")
parser$add_argument("--output", help = "Output combined Seurat and scRep object")
parser$add_argument("--meta_out", help = "Output meta data of combined Seurat and scRep object")
parser$add_argument("--call", help = "Clone call")
parser$add_argument("--align_tool", help = "mixcr or cell_ranger")


args <- parser$parse_args()

combine_rna(
  vdj = args$vdj,
  seurat = args$seurat,
  output = args$output,
  meta_out = args$meta_out,
  clone_call = args$call,
  align_tool = args$align_tool
)