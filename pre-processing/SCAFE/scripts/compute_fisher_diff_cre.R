library(Seurat)
library(DescTools)
library(Matrix)
library(tidyr)
library(optparse)

option_list = list(make_option(c("-o", "--outDir"), type = "character", default = NULL, help = "Output directory", metavar = "character"),
                           make_option(c("-g", "--contrast"), type = "character", default = NULL, help = "Name of the contrast", metavar = "character"),
			   make_option(c("-p", "--proximity"), type = "character", default = NULL, help = "", metavar = "character"),
			   make_option(c("-t", "--min_proportion"), type = "character", default = 0.0001, help = "min. threshold for CRE activity in atleast 1 of the groups", metavar = "character"),
			   make_option(c("-r", "--min_count_active_cells_threshold"), type = "character", default = 0.0001, help = "min. threshold for CRE activity in atleast 1 of the groups", metavar = "character"))

opt = parse_args(OptionParser(option_list = option_list))


proximity <- opt$proximity
outDir <- opt$outDir
contrast <- opt$contrast
min_proportion <- as.numeric(opt$min_proportion)
min_count <-  as.numeric(opt$min_count_active_cells_threshold)

seurat_object <- readRDS(file.path(outDir, "seurat_object_split", proximity, paste0("all_cells", "__", proximity, ".rds")))

contrast_vector <- read.csv(file = file.path(outDir, "contrasts", paste0("all_cells", "_model"), contrast, "cell_groups.csv"), row.names=1, header=T)
seurat_object@meta.data[["contrast"]] <- contrast_vector[rownames(seurat_object@meta.data), 1]

seurat_object@assays$SCAFE@counts@x <- rep(1, length(seurat_object@assays$SCAFE@counts@x))

foreground_seurat_object <- subset(seurat_object, subset = contrast == "foreground")
background_seurat_object <- subset(seurat_object, subset = contrast == "background")

cre_df <- data.frame(row.names = rownames(seurat_object))
cre_df[["foreground"]] <- 0
cre_df[["background"]] <- 0


cre_df[["background"]] <- Matrix::rowSums(background_seurat_object@assays$SCAFE@counts)
cre_df[["foreground"]] <- Matrix::rowSums(foreground_seurat_object@assays$SCAFE@counts)
cre_df[["combined"]] <- cre_df[["foreground"]] + cre_df[["background"]]

nCells_foreground <- ncol(foreground_seurat_object)
nCells_background <- ncol(background_seurat_object)

cre_df[["prop_foreground"]] <- cre_df[["foreground"]]/nCells_foreground
cre_df[["prop_background"]] <- cre_df[["background"]]/nCells_background
cre_df[["prop_combined"]] <- cre_df[["combined"]]/(nCells_foreground + nCells_background)


cre_df_test <- cre_df[((cre_df[["prop_foreground"]] > min_proportion & cre_df[["foreground"]]>min_count) | (cre_df[["prop_background"]] > min_proportion & cre_df[["background"]]>min_count)),]

cre_df_test_universe <- data.frame(cre_names = rownames(cre_df_test))
cre_df_test_universe <- separate(cre_df_test_universe, col = cre_names, into = c("C1", "C2", "C3", "C4"), sep = "-", extra = "merge")
write.table(cbind(cre_df_test_universe$C1, cre_df_test_universe$C2, cre_df_test_universe$C3), file = file.path(outDir, "contrasts", "all_cells_model", contrast,paste0(proximity, "_contig_tests"), "cre_unibind_universe.bed"), col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)


cre_df_test[["p_val_fisher"]] <- 1

for (cre in rownames(cre_df_test)){
	contig_table <- rbind(cre_df_test[cre, c("background", "foreground")], c(nCells_background -  cre_df_test[cre, "background"], nCells_foreground -  cre_df_test[cre, "foreground"]))	
	cre_df_test[cre, "p_val_fisher"] <- fisher.test(x = contig_table)$p.value
}

cre_df_test[["p_val_fisher"]][is.na(cre_df_test[["p_val_fisher"]])] <- 1


pdf(file.path(outDir, "contrasts", paste0("all_cells", "_model"), contrast,paste0(proximity, "_contig_tests"), "hist_p_val_tested_fisher.pdf"))
hist(as.numeric(cre_df_test[["p_val_fisher"]]))
dev.off()

cre_df_test[["fdr_fisher"]] <-  p.adjust(cre_df_test[["p_val_fisher"]], method="BH")


# Define paths for file outputs
base_output_path <- file.path(outDir, "contrasts", paste0("all_cells", "_model"), contrast, paste0(proximity, "_contig_tests"))

# Function for writing data frames to CSV and BED files
write_data_files <- function(data, prefix) {
  write.csv(data, row.names=TRUE, file=file.path(base_output_path, paste0(prefix, ".csv")))
  separated <- separate(data, col=cre_names, into=c("C1", "C2", "C3", "C4"), sep="-", extra="merge")
  write.table(cbind(separated$C1, separated$C2, separated$C3), file=file.path(base_output_path, paste0(prefix, ".bed")), col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE)
}

# Filter CREs for foreground and background enrichment
foreground_filter <- cre_df_test[["fdr_fisher"]] < 0.05 & (cre_df_test[["prop_foreground"]] - cre_df_test[["prop_background"]]) > 0 & cre_df_test[["prop_foreground"]] > min_proportion & cre_df_test[["foreground"]] > min_count
background_filter <- cre_df_test[["fdr_fisher"]] < 0.05 & (cre_df_test[["prop_foreground"]] - cre_df_test[["prop_background"]]) < 0 & cre_df_test[["prop_background"]] > min_proportion & cre_df_test[["background"]] > min_count

cre_fisher_foreground <- data.frame(cre_names=rownames(cre_df_test)[foreground_filter])
cre_fisher_background <- data.frame(cre_names=rownames(cre_df_test)[background_filter])

# Write foreground and background data files
write_data_files(cre_fisher_foreground, "diff_cre_fisher__foreground")
write_data_files(cre_fisher_background, "diff_cre_fisher__background")




# Function to plot histograms
plot_histogram <- function(column_name, title, output_file, cre_names_filter=NULL) {
  data <- as.numeric(log2(cre_df_test[[column_name]]))
  breaks <- hist(data, plot=FALSE)$breaks
  pdf(output_file)
  
  hist(data, col="blue", xaxt='n', xlab=paste("log2 of number of active cells", column_name), main=title)
  
  if (!is.null(cre_names_filter)) {
    enriched_data <- as.numeric(log2(cre_df_test[cre_names_filter, column_name]))
    hist(enriched_data, breaks=breaks, col=rgb(1,0,0, alpha=0.5), add=TRUE)
  }
  
  axis(side=1, at=breaks, labels=breaks)
  dev.off()
}


# Plot histograms for foreground, background, and combined active cells
histogram_files_path <- function(file_suffix) {
  file.path(base_output_path, paste0("hist_CRE_fisher_active_cells_", file_suffix, ".pdf"))
}

plot_histogram("foreground", "Active cells in the foreground, per CRE. Blue: Tested, Red: Enriched", histogram_files_path("foreground"), cre_fisher_foreground$cre_names)
plot_histogram("background", "Active cells in the background, per CRE. Blue: Tested, Red: Enriched in foreground", histogram_files_path("background"), cre_fisher_foreground$cre_names)
plot_histogram("combined", "Active cells per CRE. Blue: Tested, Red: Enriched in foreground", histogram_files_path("combined"), cre_fisher_foreground$cre_names)

# Plot foreground vs background histogram
pdf(histogram_files_path("fore_vs_back"))
breaks <- (-24:0)
max_y <- length(cre_fisher_foreground$cre_names) / 2

prop_background_data <- as.numeric(log2(cre_df_test[cre_fisher_foreground$cre_names, "prop_background"]))
prop_foreground_data <- as.numeric(log2(cre_df_test[cre_fisher_foreground$cre_names, "prop_foreground"]))

hist(prop_background_data, breaks=breaks, ylim=range(0:max_y), col="blue", xaxt='n', xlab="log2 of number of active cells", main="active cell prop of CREs enriched in foreground. Blue: Background cells, Red: Foreground cells")
hist(prop_foreground_data, breaks=breaks, ylim=range(0:max_y), col=rgb(1, 0, 0, alpha=0.5), add=TRUE)
axis(side=1, at=breaks, labels=breaks)
dev.off()


