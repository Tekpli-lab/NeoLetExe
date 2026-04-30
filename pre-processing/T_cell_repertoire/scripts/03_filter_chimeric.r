library(scRepertoire)
library(argparse)
library(dplyr)

filter_chimeric <- function(tcr, bcr, tcr_out, bcr_out, plot) {
  tcr_obj <- readRDS(tcr)
  bcr_obj <- readRDS(bcr)
  
  # Extract all barcodes
  barcodes_t <- unique(unlist(lapply(tcr_obj, function(x) x$barcode)))
  barcodes_b <- unique(unlist(lapply(bcr_obj, function(x) x$barcode)))
  
  # Find chimeric cell barcodes
  chimeric <- intersect(barcodes_t, barcodes_b)
  print(paste0("Total number of chimeric cells: ", length(chimeric)))

  non_chimeric_t <- setdiff(barcodes_t, chimeric)
  non_chimeric_b <- setdiff(barcodes_b, chimeric)

  # Filter scRepertoire objects directly
  filtered_tcr_obj <- lapply(tcr_obj, function(x) x[x$barcode %in% non_chimeric_t, ])
  filtered_bcr_obj <- lapply(bcr_obj, function(x) x[x$barcode %in% non_chimeric_b, ])

  ### Check if there are any chimeric cells left
  filtered_tcr_df <- do.call(rbind, filtered_tcr_obj)
  filtered_bcr_df <- do.call(rbind, filtered_bcr_obj)

  chimeric_left <- intersect(filtered_tcr_df$barcode, filtered_bcr_df$barcode)
  # Print warning if chimeric cells are left
  if (length(chimeric_left) > 0) {
    warning(paste0("Warning: There are ", length(chimeric_left), " chimeric cells left in the filtered data."))
  }
  
  ### Save scRepertoire objects
  saveRDS(filtered_tcr_obj, tcr_out)
  saveRDS(filtered_bcr_obj, bcr_out)

  ### Make a plot that shows cell counts before and after filtering
  tcr_df <- do.call(rbind, tcr_obj)
  bcr_df <- do.call(rbind, bcr_obj)

  summary_tcr <- tcr_df %>%
    mutate(status = ifelse(barcode %in% chimeric, "Chimeric", "Non-Chimeric")) %>%
    group_by(sample, status) %>%
    summarize(count = n())

  summary_bcr <- bcr_df %>%
    mutate(status = ifelse(barcode %in% chimeric, "Chimeric", "Non-Chimeric")) %>%
    group_by(sample, status) %>%
    summarize(count = n())
  
  # Plot
  pdf(plot, height = 20, width = 35)
  print(ggplot(summary_tcr, aes(x = reorder(sample, count), y = count, fill = status)) +
    geom_bar(stat = "identity", position = "stack") +
    labs(y = "Number of Cells", fill = "Cell Status") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)))
  
  print(ggplot(summary_bcr, aes(x = reorder(sample, count), y = count, fill = status)) +
    geom_bar(stat = "identity", position = "stack") +
    labs(y = "Number of Cells", fill = "Cell Status") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)))
  dev.off()
}

parser <- ArgumentParser(description = "Filtering scRepertoire list according to cell barcodes")
parser$add_argument("--tcr", "-t", help = "The product of combineTCR")
parser$add_argument("--bcr", "-b", help = "The product of combineBCR")
parser$add_argument("--tcr_out", "-ot", help = "Output TCR file")
parser$add_argument("--bcr_out", "-ob", help = "Output BCR file")
parser$add_argument("--plot", help = "Output pdf file")

args <- parser$parse_args()

filter_chimeric(tcr = args$tcr,
                bcr = args$bcr,
                tcr_out = args$tcr_out,
                bcr_out = args$bcr_out,
                plot = args$plot)