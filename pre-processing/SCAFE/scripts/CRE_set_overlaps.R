library(dplyr)
library(optparse)
library(tidyverse)


### COMMAND LINE OPTIONS ###
option_list = list(make_option(c("--output"), type = "character", default = NULL, help = "Output file", metavar = "character"),
		   make_option(c("--CRE_sets"), type = "character", default = NULL, help = "CRE_set", metavar = "character")
)

opt = parse_args(OptionParser(option_list = option_list))
output <- opt$output
CRE_bed_files <- strsplit(opt$CRE_sets, " ")[[1]]

path_dir <- dirname(CRE_bed_files[1])
file_prefixes <- sapply(CRE_bed_files, function(x) unlist(strsplit(basename(x), "__"))[1])
CRE_SETS <- unique(file_prefixes)

## Initialize an empty matrix to store overlaps
#n_CRE_sets <- length(CRE_SETS)
#overlap_matrix <- matrix(nrow = n_CRE_sets, ncol = n_CRE_sets)
#rownames(overlap_matrix) <- colnames(overlap_matrix) <- CRE_SETS
#
#
## Function to count lines in a file (i.e., number of overlaps)
#count_overlaps <- function(file) {
#  if (file.exists(file)) {
#    return(length(readLines(file)))
#  }
#}
#
#
#
## Fill the matrix with overlap counts
#for (fore_CRE in CRE_SETS) {
#  for (back_CRE in CRE_SETS) {
#    overlap_file <- file.path(path_dir, paste0(fore_CRE, "__", back_CRE, ".bed"))
#    overlap_count <- count_overlaps(overlap_file)
#    overlap_matrix[fore_CRE, back_CRE] <- overlap_count
#  }
#}
#
## Output the overlap matrix to a CSV file
#write.csv(overlap_matrix, file = output)
#



# Initialize the output matrix with dimensions based on the number of CRE sets and possible counts
overlap_count_matrix <- matrix(0, nrow = length(CRE_SETS), ncol = length(CRE_SETS))
rownames(overlap_count_matrix) <- CRE_SETS
colnames(overlap_count_matrix) <- 1:length(CRE_SETS)  # Columns representing the count categories

# Fill the matrix with overlap counts
for (fore_CRE in CRE_SETS) {
  # Read the full foreground file
  full_foreground_file <- file.path(path_dir, paste0(fore_CRE, "__", fore_CRE, ".bed"))
  foreground_lines <- readLines(full_foreground_file)

  # Create a dataframe to hold overlap presence for each back_CRE
  overlap_df <- data.frame(matrix(0, nrow = length(foreground_lines), ncol = length(CRE_SETS) - 1))
  colnames(overlap_df) <- setdiff(CRE_SETS, fore_CRE)

  # For each back_CRE, calculate overlap presence
  for (back_CRE in setdiff(CRE_SETS, fore_CRE)) {
    overlap_file <- file.path(path_dir, paste0(fore_CRE, "__", back_CRE, ".bed"))
    back_lines <- readLines(overlap_file)

    # Check overlap presence
    overlap_df[[back_CRE]] <- as.integer(foreground_lines %in% back_lines)
  }

  # Calculate row sums (the number of overlap occurrences for each line)
  overlap_counts <- rowSums(overlap_df)

  # Use table to summarize how many times lines were seen in 0, 1, ..., length(CRE_SETS) - 1 files
  overlap_summary <- table(overlap_counts)

  # Fill the matrix with the summary, ensuring proper indexing
  for (count in names(overlap_summary)) {
    overlap_count_matrix[fore_CRE, as.integer(count) + 1] <- overlap_summary[count]
  }
}

# Output the overlap count matrix to a CSV file
write.csv(overlap_count_matrix, file = output)

print("Overlap count matrix calculated.")

