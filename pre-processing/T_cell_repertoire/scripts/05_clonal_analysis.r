library(scRepertoire)
library(ggplot2)
library(argparse)
library(patchwork)
library(dplyr)
library(tidyr)
library(ggrepel)


clonal_analysis <- function(airr, output, clone_call, vdj) {
  obj <- readRDS(airr)

  pdf(output, height = 13, width = 17)
  #Basic QC
  #cell numbers across samples
  obj_unlist <- do.call(rbind, obj)
  cell_count <- obj_unlist %>%
    count(sample, name = "cell_number")
  print(ggplot(cell_count, aes(x = reorder(sample, -cell_number), y = cell_number)) +
          geom_bar(stat = "identity", fill = "#ffa6de", color = "maroon2") +
          theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
          xlab("samples"))
  #clone numbers across samples
  print(lapply(obj, dim))
  clones_df <- clonalAbundance(obj, cloneCall = clone_call, exportTable = TRUE)
  clone_count <- clones_df %>%
    count(values, name = "clone_number") %>%
    rename(sample = values)
  print(ggplot(clone_count, aes(x = reorder(sample, -clone_number), y = clone_number)) +
          geom_bar(stat = "identity", fill = "orange", color = "maroon2") +
          theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
          xlab("samples"))

  #cell numbers on clone numbers
  cell_on_clone <- merge(cell_count, clone_count, by = "sample")
  cell_on_clone$cells <- cell_on_clone$cell_number - cell_on_clone$clone_number
  df_plot <- cell_on_clone %>%
    select(!cell_number) %>%
    pivot_longer(cols = !sample,
                 names_to = "counts",
                 values_to = "number")
  print(ggplot(df_plot, aes(x = reorder(sample, -number), y = number, fill = counts)) +
          geom_bar(position = "stack", stat = "identity", color = "purple3") +
          scale_fill_manual(values = c("purple3", "mediumorchid1")) +
          theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
          xlab("samples"))

  #scRepertoire graphs
  print(clonalQuant(obj, cloneCall = clone_call, palette = "dynamic", scale = TRUE) +
          theme(legend.position = "none",
                axis.text.x = element_text(angle = 45, hjust = 1)))
  print(clonalQuant(obj, cloneCall = clone_call, group.by = "time_point", scale = TRUE))

  #relative space occupied by clones at specific proportions
  print(clonalHomeostasis(obj, cloneCall = clone_call) +
          theme(axis.text.x = element_text(angle = 45, hjust = 1)))
  print(clonalHomeostasis(obj, cloneCall = clone_call, group.by = "patient"))
  print(clonalHomeostasis(obj, cloneCall = clone_call, group.by = "time_point"))

  #rank the clones by total number and place them into bins
  print(clonalProportion(obj, cloneCall = clone_call) +
          theme(axis.text.x = element_text(angle = 45, hjust = 1)))
  print(clonalProportion(obj, cloneCall = clone_call, group.by = "time_point"))

  #Find the number of shared clones between samples
  print(clonalOverlap(obj, cloneCall = clone_call, method = "raw",
                      palette = "sunset") +
          theme(axis.text.x = element_text(angle = 45, hjust = 1)))

  ### Comparisons per patient ###

  #get the patient list
  patients <- list()
  for (i in seq_along(obj)) {
    patients <- append(patients, unique(obj[[i]]$patient))
  }
  patients <- unique(patients)

  #empty lists for plots to be plotted in the same page
  patch_quant <- list()
  patch_compare <- list()
  patch_homeostasis <- list()
  patch_diversity <- list()

  for (patient in patients) {
    sub_obj <- subsetClones(obj, name = "patient", variables = patient)

    # Percent of unique clones
    viz_quant <- clonalQuant(sub_obj, cloneCall = clone_call,
                             palette = "Red-Blue",
                             scale = TRUE, exportTable = FALSE) +
      theme(legend.position = "none",
            axis.text.x = element_text(angle = 45, hjust = 1)) + ylim(0, 100)
    # Debugging because clonalQuant returns not only ggplot objects
    if (!inherits(viz_quant, "ggplot")) {
      message("Non-ggplot object found for patient: ", patient)
      message("Object classes: ", paste(class(viz_quant), collapse = ", "))
      next
    }
    patch_quant[[length(patch_quant) + 1]] <- viz_quant

    # Change in proportion of top clones across timepoints
    viz_compare <- clonalCompare(sub_obj,
                                 cloneCall = clone_call,
                                 top.clones = 10,
                                 palette = "Zissou 1",
                                 graph = "alluvial") +
      theme(legend.position = "none")
    patch_compare[[length(patch_compare) + 1]] <- viz_compare

    # relative space occupied by clones at specific proportions
    viz_homeostasis <- clonalHomeostasis(sub_obj, cloneCall = clone_call) +
      theme(legend.position = "none")
    patch_homeostasis[[length(patch_homeostasis) + 1]] <- viz_homeostasis
    # Shannon index for clonal diversity
    viz_diversity <- clonalDiversity(sub_obj,
                                     cloneCall = clone_call,
                                     exportTable = FALSE,
                                     metrics = c("shannon", "ACE"),
                                     x.axis = "time_point") +
      theme(legend.position = "none") +
      ggtitle(patient)
    patch_diversity[[length(patch_diversity) + 1]] <- viz_diversity
  }

  print(wrap_plots(patch_quant))
  print(wrap_plots(patch_compare))
  print(wrap_plots(patch_homeostasis))
  print(wrap_plots(patch_diversity))

  ### Shannon index box plot ###
  diversity_result <- clonalDiversity(obj, cloneCall = clone_call,
                                      exportTable = TRUE)
  df <- data.frame(samples = diversity_result$Group,
                   shannon = diversity_result$shannon)
  df$patient_id <- substr(df$samples, 1, 5)
  df$time_point <- factor(substr(df$samples, 6, 6), levels = c("A", "B", "C"))

  df <- df %>%
    arrange(patient_id, time_point)

  box_diversity <- ggplot(df, aes(x = time_point,
                                  y = shannon, label = patient_id)) +
    geom_boxplot(color = c("#008033", "#00b000", "#43b800"),
                 fill = c("#00b749", "#24ff24", "#c3ffa1"), alpha = 0.2) +
    geom_line(aes(group = patient_id), size = 1,
              color = ifelse(is.na(sapply(1:nrow(df),
                                          function(i)ifelse(df$shannon[i] > df$shannon[i+1],
                                                            "blue", "red"))),
                             "orange",
                             sapply(1:nrow(df),
                                    function(i)ifelse(df$shannon[i] > df$shannon[i + 1], 
                                                      "blue", "red"))),
              alpha = 0.6) +
    geom_point(aes(colour = factor(time_point), fill = factor(time_point)),
               shape = 21, size = 5) +
    scale_fill_manual(values = c("white", "white", "white")) +
    scale_colour_manual(values = c('mediumpurple', '#e593f5' , '#ffbcdf')) +
    ylab("Shannon Index") +
    geom_text_repel(min.segment.length = 0.5, max.overlaps = Inf) +
    scale_x_discrete(labels = c("A" = "Baseline", "B" = "Crosssover", "C" = "End of treatment")) +
    theme(axis.line = element_line(colour = "black", linewidth = 1, linetype = "solid")) +
    theme(axis.text.x = element_text(size = 14),
          axis.text.y = element_text(size = 14)) +
    theme(axis.title = element_text(size = 18, face = "bold")) +
    theme(axis.ticks = element_line(linewidth = 1)) +
    theme(panel.background = element_rect(fill = 'white', color = 'black')) +
    labs(x = "Time point") +
    theme(legend.position = "none") +
    ggtitle("Shannon index change with treatment")
  print(box_diversity)

  dev.off()
}

parser <- ArgumentParser(description = "Create clonal analysis plots")
parser$add_argument("--object", "-i", help = "The product of combineTCR/BCR")
parser$add_argument("--output", "-o", help = "Output pdf file")
parser$add_argument("--clone", "-c", help = "How to call clones")
parser$add_argument("--vdj", "-v", help = "TCR or BCR wildcard")

args <- parser$parse_args()


clonal_analysis(airr = args$object,
                output = args$output,
                clone_call = args$clone,
                vdj = args$vdj)

