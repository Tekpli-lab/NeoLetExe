library(scRepertoire)
library(argparse)
library(dplyr)

exclusion <- function (object, output, samples) {
    obj <- readRDS(object)
    
    print(paste0("Number of samples before exclusion ", length(obj)))

    samples <- unlist(strsplit(samples, split = ","))
    #print(paste("Samples to exclude", samples))

    all <- names(obj)
    print(all)
    keep <- all[!(all %in% samples)]

    #obj <- obj[!sapply(obj, function(df) any(df$sample %in% samples))]
    obj <- subsetClones(obj, name = "sample", variables = keep)

    print(paste0("Number of samples after exclusion ", length(obj)))
    saveRDS(obj, output)
}

parser <- ArgumentParser(description = "Find shared clones across patients")
parser$add_argument("--object", "-i", help = "scRep object that is assesed for QC")
parser$add_argument("--output", "-o", help = "Output scRep object")
parser$add_argument("--samples", nargs = "+", help = "List of samples to exclude")


args <- parser$parse_args()

exclusion(object = args$object,
          output = args$output,
          samples = args$samples)