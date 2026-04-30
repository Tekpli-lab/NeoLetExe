library(scRepertoire)
library(tidyverse)
library(argparse)

load_tcr <- function(contigs, sample_names, out, remove_na, remove_multi,
                     filter_multi, timepoint) {
  # read contigs into a list
  data_list <- lapply(contigs, function(x){read_tsv(x, show_col_types = FALSE)})

  # Convert string arguments to logical
  remove_na <- tolower(remove_na) == "true"
  remove_multi <- tolower(remove_multi) == "true"
  filter_multi <- tolower(filter_multi) == "true"

  for (i in seq_along(data_list)) {
    if (dim(table(data_list[[i]]$topChains)) > 2) {
      print(table(data_list[[i]]$topChains))
    }
  }

  # prepare mixcr output for scRep
  filtered_data_list <- lapply(data_list, function(df) {
    df %>%
      filter(cellGroup != "contamination") %>%
      filter(topChains == "TRA" | topChains == "TRB") %>%
      select(tagValueCELL, topChains, readCount, bestVGene, bestDGene,
             bestJGene, bestCGene, nSeqCDR3, aaSeqCDR3) %>%
      rename(allVHitsWithScore = bestVGene, allDHitsWithScore = bestDGene,
             allJHitsWithScore = bestJGene, allCHitsWithScore = bestCGene)
  })

  filtered_data_list <- lapply(data_list, function(df) {
    # Step 1: Identify cells having chains other than "TRA" or "TRB"
    unwanted_cells <- df %>%
      filter(!(topChains %in% c("TRA", "TRB"))) %>%
      pull(tagValueCELL)

    # Step 2: Filter out unwanted cells and keep only "TRA" or "TRB" chains
    df %>%
      filter(!(tagValueCELL %in% unwanted_cells)) %>%
      filter(cellGroup != "contamination") %>%
      select(tagValueCELL, topChains, readCount, bestVGene, bestDGene,
             bestJGene, bestCGene, nSeqCDR3, aaSeqCDR3) %>%
      rename(allVHitsWithScore = bestVGene, allDHitsWithScore = bestDGene,
             allJHitsWithScore = bestJGene, allCHitsWithScore = bestCGene)
  })

  contig_list <- loadContigs(filtered_data_list, format = "MiXCR")

  # get sample names
  samples <- unlist(strsplit(sample_names, split = ","))

  for (i in seq_along(filtered_data_list)) {
    print("after filtering")
    print(table(filtered_data_list[[i]]$topChains))
  }

  # combine samples
  combined <- combineTCR(contig_list,
                         samples = samples,
                         removeNA = remove_na,
                         removeMulti = remove_multi,
                         filterMulti = filter_multi)

  combined <- lapply(combined, function(x) {
    x$cdr3s_pat <- paste(x$CTaa, x$sample, sep = "_")
    x
  })

  if (timepoint) {
    timepoints <- str_sub(samples, start = -1)
    print(timepoints)
    patients <- str_sub(samples, end = -2)
    print(patients)

    combined <- addVariable(combined,
      variable.name = "time_point",
      variables = timepoints
    )
    combined <- addVariable(combined,
      variable.name = "patient",
      variables = patients
    )
  }

  saveRDS(combined, out)
}

parser <- ArgumentParser(description = "Load contigs from TCR sequencing data")
parser$add_argument("--contigs", "-c", nargs = "+",
                    help = "input list TCR contigs")
parser$add_argument("--names", "-n", nargs = "+",
                    help = "input sample names list")
parser$add_argument("--out", "-o",
                    help = "output rds file to save combined clones list")
parser$add_argument("--na", "-a",
                    help = "T or F. Remove cell with an NA value in one chain")
parser$add_argument("--rem_multi", "-r",
                    help = "T or F. Remove cell with >2 receptor chains")
parser$add_argument("--fil_multi", "-f",
                    help = "T or F. Isolate the top 2 expressed chains")
parser$add_argument("--times", "-t",
                    help = "T or F. Add a timepoints variable")

args <- parser$parse_args()

load_tcr(
  contigs = args$contigs,
  sample_names = args$names,
  out = args$out,
  remove_na = args$na,
  remove_multi = args$rem_multi,
  filter_multi = args$fil_multi,
  timepoint = args$times
)
