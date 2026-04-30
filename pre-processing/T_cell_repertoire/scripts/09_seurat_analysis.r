library(Seurat)
library(dplyr)
library(scRepertoire)
library(ggraph)
library(circlize)
library(scales)
library(APackOfTheClones)
library(argparse)

seurat_analysis <- function(object, plot, clone_call, cluster_column) {
  obj <- readRDS(object)

  # handle clone call
  clone_call <- gsub("CT", "", clone_call)

  pdf(plot)
  # Plot clonal size proportions on UMAP
  # define color palette 
  colorblind_vector <- hcl.colors(n=7, palette = "inferno", fixup = TRUE)

  DimPlot(obj, group.by = "cloneSize") +
      scale_color_manual(values=rev(colorblind_vector[c(1,3,5,7)]))
  
  # Overlays frequency (filters clones smaller than 3 cells)
  print(clonalOverlay(obj, 
            reduction = "umap", 
            cut.category = "clonalFrequency",
            cutpoint = 3, 
            bins = 10) + 
            guides(color = "none"))
  
  print(clonalOverlay(obj, 
            reduction = "umap", 
            cut.category = "clonalFrequency",
            cutpoint = 3, 
            bins = 10, 
            facet.by = "orig.ident") + 
            guides(color = "none"))

  # Clones that are shared between clusters
  print(clonalNetwork(obj, 
                reduction = "umap", 
                group.by = cluster_column,
                filter.clones = NULL,
                filter.identity = NULL,
                cloneCall = clone_call))
  
  print(clonalNetwork(obj, 
                reduction = "umap", 
                group.by = cluster_column,
                filter.clones = NULL,
                filter.identity = "CD8_EM",
                cloneCall = clone_call))
  
  # clone sizes in clusters
  print(clonalOccupy(obj, 
               x.axis = cluster_column, 
               proportion = FALSE, 
               label = FALSE) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)))
  print(clonalOccupy(obj, 
               x.axis = cluster_column, 
               proportion = TRUE, 
               label = FALSE) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)))
  
  # circle plot to see clone sharing between clusters
  circles <- getCirclize(obj, group.by = cluster_column)

  #Just assigning the normal colors to each cluster
  grid.cols <- hue_pal()(length(unique(obj$final.ids)))
  names(grid.cols) <- unique(obj$final.ids)

  #Graphing the chord diagram
  chordDiagram(circles, self.link = 1, grid.col = grid.cols)

  # produce the ball-packing plot
  print(vizAPOTC(obj))

  dev.off()
}

parser <- ArgumentParser(description = "Find shared clones across patients")
parser$add_argument("--object", help = "scRep object")
parser$add_argument("--plot", help = "Pdf output")
parser$add_argument("--call", help = "Clone call")
parser$add_argument("--clusters", help = "Name of the column that has cluster info")


args <- parser$parse_args()

seurat_analysis(
  obj = args$object,
  plot = args$plot,
  clone_call = args$call,
  cluster_column = args$clusters
)