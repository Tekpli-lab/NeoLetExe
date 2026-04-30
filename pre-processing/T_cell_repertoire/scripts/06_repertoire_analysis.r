library(scRepertoire)
library(ggplot2)
library(argparse)
library(tidyverse)
library(patchwork)
library(ggrepel)
library(pheatmap)

repertoire_analysis <- function(object, out, vdj, chain) {
  obj <- readRDS(object) #scRepertoire object
  obj_df <- do.call(rbind, obj)

  # Fix column names to match chain names
  if (vdj == "tcr") {
    obj_df <- obj_df %>%
      rename(TRA = TCR1, TRB = TCR2)
  } else {
    obj_df <- obj_df %>%
      rename(IGL = IGLC)
  }

  # Correct for D gene after splitting scRep genes column
  if (chain %in% c("TRB", "IGH")) {
    j_column <- "genes_3"
  } else {
    j_column <- "genes_2"
  }

  # Calculate gene usage correlation between samples
  chain_df <- data.frame("sample" = obj_df$sample,
                         "clone" = obj_df$CTstrict,
                         "genes" = obj_df[, chain])

  chain_df <- chain_df %>%
    separate_wider_delim(col = "genes", delim = ".", names_sep = "_") %>%
    distinct(.keep_all = TRUE)

  # How many clones use which V gene
  v_genes <- table(chain_df$sample, chain_df$genes_1)
  total_v <- rowSums(v_genes) # total use of each gene

  # Take percentage of gene usage
  percentage_v <- v_genes / total_v * 100
  correlation_v <- cor(t(percentage_v))

  # How many clones use which J gene
  j_genes <- table(chain_df$sample, chain_df[[j_column]])

  total_j <- rowSums(j_genes) # total use of each gene

  # Take percentage of gene usage
  percentage_j <- j_genes / total_j * 100
  correlation_j <- cor(t(percentage_j))

  #function to perform and plot pca
  plot_pca <- function(df_genes) {
    #Performing PCA
    pc <- prcomp(df_genes)

    #Getting data frame to plot from
    df <- as.data.frame(cbind(pc$x[, 1:2], rownames(df_genes)))
    df$PC1 <- as.numeric(df$PC1)
    df$PC2 <- as.numeric(df$PC2)

    #Plotting
    return(ggplot(df, aes(x = PC1, y = PC2, label = rownames(df))) +
             geom_point(aes(fill = df[, 3]), shape = 21, size = 5) +
             scale_fill_manual(values = hcl.colors(nrow(df), "inferno")) +
             theme_classic() +
             theme(legend.position = "none") +
             geom_label_repel())
  }

  ### Plotting repertoire analysis ###

  pdf(out, height = 13, width = 14)
  print(pheatmap(correlation_v, main = "V gene usage"))
  print(pheatmap(correlation_j, main = "J gene usage"))
  # Create plots for each chain
  # quantify the level of diversity across amino acid residues along cdr3
  print(positionalEntropy(obj, chain = chain, group.by = "patient",
                          method = "shannon") +
          ggtitle(chain))
  print(paste0("positional entropy done for ", chain))

  ### Gene usage ###

  # Concatenate gene names to chain
  v_gene <- paste0(chain, "V")
  j_gene <- paste0(chain, "J")

  genes <- c(v_gene, j_gene)

  # specific gene usage
  for (y in genes) {
    print(y)
    print(vizGenes(obj, x.axis = y, y.axis = NULL,
                   order = "variance", plot = "heatmap", scale = TRUE))
  }

  print(percentGenes(obj, chain = chain, gene = "V"))
  print(percentGenes(obj, chain = chain, gene = "J"))

  # PCA by specific gene usage
  df_v <- percentGenes(obj, chain = chain, gene = "V", exportTable = TRUE)
  df_j <- percentGenes(obj, chain = chain, gene = "J", exportTable = TRUE)

  print(plot_pca(df_v) + ggtitle(paste0("PCA by ", chain, "V gene usage")))
  print(paste0("PCA by ", chain, "V gene usage"))
  print(plot_pca(df_j) + ggtitle(paste0("PCA by ", chain, "J gene usage")))
  print(paste0("PCA by ", chain, "J gene usage"))


  print(percentKmer(obj, cloneCall = "aa", chain = chain))

  dev.off()
}

parser <- ArgumentParser(description = "Create repertoire analysis plots")
parser$add_argument("--object", "-i", help = "The product of combineTCR/BCR")
parser$add_argument("--output", "-o", help = "Output pdf file")
parser$add_argument("--vdj", "-v", help = "TCR or BCR wildcard")
parser$add_argument("--chain", "-c", help = "which chain to plot")
args <- parser$parse_args()

repertoire_analysis(object = args$object,
                    out = args$output,
                    vdj = args$vdj,
                    chain = args$chain)