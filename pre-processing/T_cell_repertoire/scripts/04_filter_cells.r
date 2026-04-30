library(scRepertoire)
library(argparse)
library(dplyr)

filter_airr <- function(object, seurat_metadata,
                        align_tool, cell_types,
                        cell_states, cluster_list,
                        plot, output) {
  obj <- readRDS(object) #scRepertoire object
  obj_df <- do.call(rbind, obj)
  seurat <- readRDS(seurat_metadata)

  cluster_list <- unlist(strsplit(cluster_list, split = ","))

  seurat <- subset(seurat, get(cell_types) %in% cluster_list)
  barcode_seurat <- rownames(seurat)

  cell_barcodes <- barcode_seurat
  # Fix barcodes for cells aligned with MiXCR
  if(align_tool == "mixcr") {
    cell_barcodes <- gsub("-1", "", cell_barcodes)
  }

  for (i in seq_along(obj)) { #seq_along to iterate over indices
    obj[[i]] <- obj[[i]][obj[[i]]$barcode %in% cell_barcodes, ]
  }
  saveRDS(obj, output)

  # Plotting TCR-seq ratios in each cell state per sample
  vdj_barcodes <- unique(obj_df$barcode)

  # Fix barcodes for cells aligned with MiXCR
  if(align_tool == "mixcr") {
    vdj_barcodes <- paste0(vdj_barcodes, "-1")
  }
  
  seurat %>% filter(rownames(.) %in% vdj_barcodes)

  seurat$vdj <- ifelse(rownames(seurat) %in% vdj_barcodes, "yes", "no")

  vdj_ratio <- seurat %>%
    group_by(orig.ident, .data[[cell_states]]) %>%
    count(vdj) %>%
    reframe(vdj = vdj, ratio = n / sum(n)) %>%
    filter(vdj == "yes") %>%
    select(-vdj)

  ids <- unique(vdj_ratio[[cell_states]])

  pdf(plot, height = 5, width = 9)
  for (id in ids) {
    ratios <- vdj_ratio %>%
      filter(.data[[cell_states]] == id)
    print(ggplot(ratios, aes(x = orig.ident, y = ratio, fill = orig.ident)) +
            geom_col(color = "black") +
            ylim(0, 1) +
            theme(legend.position = "none",
                  axis.text.x = element_text(angle = 90, vjust = 0.5)) +
            ggtitle(paste0("Proportion of cells with TCR-seq in ", id)))
  }
  dev.off()
}

parser <- ArgumentParser(description = "Filtering scRepertoire list according to cell barcodes")
parser$add_argument("--object", "-i", help = "The product of combineTCR/BCR")
parser$add_argument("--seurat", "-s", help = "Input the annotated seurat object")
parser$add_argument("--align_tool", help = "The tool that was used for aligning cell_ranger or mixcr")
parser$add_argument("--col1", help = "Cell types column name in metadata of seurat clusters")
parser$add_argument("--col2", help = "Cell states column name in metadata of seurat clusters")
parser$add_argument("--clusters", "-c", nargs = "+", help = "List of clusters to keep")
parser$add_argument("--plot", help = "Plot for TCR-seq ratios in each cell state per sample")
parser$add_argument("--output", "-o", help = "Output file")

args <- parser$parse_args()


filter_airr(
  object = args$object,
  seurat_metadata = args$seurat,
  align_tool = args$align_tool,
  cell_types = args$col1,
  cell_states = args$col2,
  cluster_list = args$clusters,
  plot = args$plot,
  output = args$output
)