###################################################
####### FIGURE 1 ##################################
###################################################

###################################################
###### Figure 1B ########## MAIN UMAP#############
###################################################

library(Seurat)
library(plyr)
library(dplyr)
library(future)
library (ggplot2)
library (ggrepel)
library (stringi)
library(SeuratData)
library(Palo)
library (RColorBrewer)

load("Scaled_NeoLetExe_final_git.RData")
#### Load Seurat object obtained from PP_Seurat_object
load("palopal.RData")
#### Load an object containing color informations

Idents(Scaled) <- "main.ids"  ### link betwwen barcodes and annotation is available in Supplementary Tables
u <- Scaled[["umap"]]@cell.embeddings
cl <- Idents(Scaled)

gg_color_hue <- function(n) {
  hues = seq(30, 500, length = n + 1)
  hues = seq(30, 500, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
pal <- gg_color_hue(length(unique(cl)))

getPalette = colorRampPalette(brewer.pal(12, "Paired"))
pal1 <- getPalette(length(unique(cl)))
palopal <- Palo(u,cl,pal1, rgb_weight=c(10,10,10))


palopal [1] <-"#A6CEE3" # Epith
palopal [2] <-"#F06C45" # CAF
palopal [3] <-"#B294C7" # Tcell
palopal [4] <- "#AB9AFF" # Bcells
palopal [5] <- "#FF9933"# Plasma
palopal [6] <- "#52AF43" # Endoth
palopal [7] <- "#FF99CC" # Myeloid
palopal [8] <- "#DAA70F" # PVFIBRO
palopal [9] <- "#666666"# Mast
palopal [10] <- "#006633" # Neutro
palopal [11] <- "#996600" # Adipocy
palopal [12] <- "#9A4D73" # DC

pdf ('Figure1B_main_id_UMAP.pdf')
DimPlot(Scaled, label = F, group.by='main.ids', cols=palopal, repel=F, label.box=F, label.size=4.5, raster.dpi = c(2048, 2048), raster=F, label.color= palopal) +NoLegend()
dev.off()


#########################################################
###### FIGURE 1C ########## HEATMAP OF MARKER genes######
#########################################################

library(Seurat)
library(plyr)
library(dplyr)
library(future)
library (ggplot2)
library (ggrepel)
library (stringi)
library(SeuratData)
library(Palo)
library (RColorBrewer)

load("Scaled_NeoLetExe_final_git.RData")
#### Load Seurat object obtained from PP_Seurat_object
load("palopal.RData")
#### Load an object containing color information

all_marker_genes <- c("EPCAM","FOXA1","ESR1", "KRT17", "KRT14",
"COL1A1","C1R","SERPINE1","CFD","COL1A2",
"CD3D", "IL7R", "CD8A","GZMK", "FOXP3","XCL1",
"CD79A","MS4A1" ,"CD83", "CD19", 
"IGHG1", "CD38","MZB1","IGHA2", 
"CA4", "CD36","VWF","ACKR1", "TIMP4",
"CD68", "CCL13", "CCL18", "CD163", "CD14", "IL1B" ,
"MYH11","RGS5","NOTCH3", "ACTA2", "MYL9",
"CPA3", "TPSAB1","TPSB2","GATA2",
"NAMPT","IL1RN","CAMP", "CSF3R",
"LEP","SLC7A10", "ADIPOQ","CIDEA",
"TBX1", "TFPI", "LYVE1", 
"LILRA4","PLAC8","IRF7")

Idents(Scaled) <- "main.ids" ### link betwwen barcodes and annotation is available in Supplementary Tables
cluster.averages <- AverageExpression(Scaled, return.seurat=TRUE)
vect <- c('Epithelial_cells', 'T_cells', 'Fibroblasts',  'Endothelial_cells', 'Myeloid_cells', 'PV_Fibroblasts', 'Plasmablasts', 'B_cells' ,'Mast_cells', 'Neutrophils', 'Adipocytes', 'Dendritic_cells')

u <- Scaled[["umap"]]@cell.embeddings
cl <- Idents(Scaled)

gg_color_hue <- function(n) {
  hues = seq(30, 500, length = n + 1)
  hues = seq(30, 500, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
pal <- gg_color_hue(length(unique(cl)))

getPalette = colorRampPalette(brewer.pal(12, "Paired"))
pal1 <- getPalette(length(unique(cl)))
palopal <- Palo(u,cl,pal1, rgb_weight=c(10,10,10))


palopal [1] <-"#A6CEE3" # Epith
palopal [2] <-"#F06C45" # CAF
palopal [3] <-"#B294C7" # Tcell
palopal [4] <- "#AB9AFF" # Bcells
palopal [5] <- "#FF9933"# Plasma
palopal [6] <- "#52AF43" # Endoth
palopal [7] <- "#FF99CC" # Myeloid
palopal [8] <- "#DAA70F" # PVFIBRO
palopal [9] <- "#666666"# Mast
palopal [10] <- "#006633" # Neutro
palopal [11] <- "#996600" # Adipocy
palopal [12] <- "#9A4D73" # DC

pdf ('Figure 1C marker_genes_heatmap.pdf')
DoHeatmap(cluster.averages,features = all_marker_genes, draw.lines = F, group.colors=palopal, size=4)+ scale_fill_gradientn(colors = c("dodgerblue4","white","brown3")) +  ggtitle('') +
  theme(axis.text.y = element_text(size = 7, face = "bold"),  
        legend.text = element_text(face = "bold"),  
        legend.title = element_text(face = "bold"))  
dev.off()


########################################################################
###### FIGURE 1D ########## BARPLOTS of celltypes per time point #######
########################################################################
##### from new metadata with new numbers all good

load("metadata_git.RData")
#### load the metadata available at EGA
load("palopal.RData")
#### Load an object containing color information

dta <- as.data.frame.matrix(prop.table(table (metadata$Time_point, metadata$main.ids),1)*100)
dta <- dta[, rev(order(dta[1,]))]
dta <- dta [c(3,1,2),]

pdf('Figure 1D_barplots_main.ids_proportions_time_points.pdf')
par(las = 2)
par(mar = c(5, 5, 5, 10))
barplot(t(dta), legend.text = TRUE ,  col =palopal[colnames(dta)], args.legend = list(x = "topright",inset = c(- 0.42, 0)))
dev.off()



###################################################
#######  FIGURE 1E  ###INFERCNV UMAP###############
###################################################


library(Seurat)
library(plyr)
library(dplyr)
library(future)
library (ggplot2)
library (ggrepel)
library (stringi)
library(SeuratData)
library(Palo)
library (RColorBrewer)

load("annot_umap_git.RData")
#### Need InfercNV UMAP and annotation object obtained from PP_InferCNV

color_map <- c(
  "1" = "darkgreen",
  "2" = "darkorange",
  "3" = "blue",
  "4" = "darkgoldenrod",
  "5" = "pink",
  "6" = "coral",
  "7" = "deeppink4",
  "8" = "brown",
  "9" = "cyan",
  "10" = "chartreuse",
  "11" = "purple",
  "12" = "black",
  "13" = "gray",
  "14" = "orange",
  "15" = "deep pink",
  "16" = "red",
  "17" = "chocolate",
  "18" = "yellow",
  "19" = "yellow3",
  "20" = "cadetblue",
  "21" = "darkblue",
  "22" = "khaki",
  "23" = "violet",
  "24" = "turquoise",
  "25" = "limegreen"
)

annot_umap2$patient_col <- color_map[as.character(annot_umap2$patient)]



pdf ('Figure 1E inferCNV_UMAP_centered.pdf')
#with text
ordered_patients <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22",  "23", "24", "25")
ordered_colors <- annot_umap2$patient_col[match(ordered_patients, annot_umap2$patient)]

plot(mat.umap2$layout, col = annot_umap2$patient_col, cex = .5, xlim = c(-12, 18), 
     xlab = expression(bold("UMAP 1")), ylab = expression(bold("UMAP 2")))
legend("topright", legend = ordered_patients, pch = 16, col = ordered_colors, cex = 0.9, bty = "n")

df <- data.frame(cbind(mat.umap2$layout, as.numeric(annot_umap2$patient)))
colnames(df) <- c('a', 'b', 'c')
result <- data.frame(df %>%group_by(c) %>%summarize(Average_Column1 = mean(a),Average_Column2 = mean(b)))
text(result[, 2], result[, 3], labels = result[, 1], pos = 3, col = ordered_colors, cex = 0.9)

#without text
plot(mat.umap2$layout, col = annot_umap2$patient_col, cex = .5, xlim = c(-12, 18), 
     xlab = expression(bold("UMAP 1")), ylab = expression(bold("UMAP 2")))
legend("topright", legend = ordered_patients, pch = 16, col = ordered_colors, cex = 0.9, bty = "n")

dev.off()


##################################################################################
#######  FIGURE 1F  ### Heatmap of chromosome alterations in patients ############
##################################################################################

library(Seurat)
library(plyr)
library(dplyr)
library(future)
library (ggplot2)
library (ggrepel)
library (stringi)
library(SeuratData)
library(Palo)
library (RColorBrewer)
library (pheatmap)

load("annot_umap_git.Rdata")
#### Need InfercNV UMAP and annotation object obtained from PP_InferCNV

annot_umap2$orig.ident <- gsub("_(.*)", "", rownames (annot_umap2))
dat<- data.frame(mat.umap2$data)

result <- data.frame (dat %>%
  group_by(c = annot_umap2$patient) %>%
  summarize(across(everything(), mean)))
rownames (result) <- result$c
mean_mat <- t (result[,-1])

gene <- data.frame(read.table("gene_location_hg38.txt", header=TRUE, sep="\t", stringsAsFactors=F))
###### need a bed file with gene hg38 locations

inters <-  intersect(rownames (mean_mat), gene$geneName)
gene_list <- gene[gene$geneName%in%inters,]
gene_list <- gene_list[,colnames (gene_list) %in% c('chrom','chromStart','chromEnd','geneName')]
gene_list <- gene_list[!duplicated (gene_list$geneName),]
mean_mat <- mean_mat[rownames (mean_mat)%in% inters,]
mean_mat <- mean_mat[match (gene_list$geneName, rownames (mean_mat)),]
table (rownames (mean_mat)==gene_list$geneName)
gene_list$chr_num <- gsub ("chr", "", gene_list$chrom)

chr<-gene_list$chr_num
cl <- (colnames (mean_mat))


anno_colors = list(chr = c("1" = "gray", "2" = "black","3" = "gray", "4" = "black","5" = "gray", "6" = "black","7" = "gray", "8" = "black","9" = "gray", "10" = "black","11" = "gray", "12" = "black","13" = "gray", "14" = "black",
"15" = "gray", "16" = "black","17" = "gray", "18" = "black","19" = "gray", "20" = "black","21" = "gray", "22" = "black","M" = "gray"),
cl = c(
  "1" = "darkgreen",
  "2" = "darkorange",
  "3" = "blue",
  "4" = "darkgoldenrod",
  "5" = "pink",
  "6" = "coral",
  "7" = "deeppink4",
  "8" = "brown",
  "9" = "cyan",
  "10" = "chartreuse",
  "11" = "purple",
  "12" = "black",
  "13" = "gray",
  "14" = "orange",
  "15" = "deep pink",
  "16" = "red",
  "17" = "chocolate",
  "18" = "yellow",
  "19" = "yellow3",
  "20" = "cadetblue",
  "21" = "darkblue",
  "22" = "khaki",
  "23" = "violet",
  "24" = "turquoise",
  "25" = "limegreen"))

rowcols<-data.frame (cbind(gene_list$chr_num))
annotation_row= data.frame (rowcols)
rownames(annotation_row) <- gene_list$geneName
colnames(annotation_row) <- c("chr")

colcols<-cbind( cl)
annotation_col= data.frame (colcols)
rownames(annotation_col) <- colnames (mean_mat)
colnames(annotation_col) <- c("cl")
library (pheatmap)
library (amap)
library (gwasrapidd)

dcols1 <- hcluster(t(mean_mat), method = 'correlation', link = 'ward')
breaksList = seq(-1, 1, by = 0.05)
col.pal <- colorRampPalette(c("dodgerblue4","white","brown3"))(length(breaksList))

pdf ('Figure 1F InferCNV_chromosomes_patiens.pdf')
pheatmap ((mean_mat-1)*15,  raster=F, color = col.pal,show_rownames = F, show_colnames = T, main = "",  
cluster_rows = F, cluster_cols = dcols1, annotation_col=annotation_col, annotation_row=annotation_row, annotation_colors = anno_colors, 
breaks = breaksList,treeheight_row=30 , treeheight_col=30, border_color = NA)
dev.off()


##################################################################################
#######  FIGURE 1G  ### Clone numbers across time boxplot #######################
##################################################################################
library(Seurat)
library(SeuratData)
library(Palo)
library(ggplot2)
library (RColorBrewer)
library (tabula)
library (stringi)
library (ggrepel)
library(ggsignif)

load("metadata_git.RData")
#### load the metadata available at EGA

meta <- metadata[metadata$clone_cat%in%c('res', 'sens', 'neut') & !is.na(metadata$clone_cat),]

a<- as.data.frame.matrix(table (meta$orig.ident, meta$clone_id))
a$nclone <- rowSums(a > 10)

a$category <- stri_sub(rownames (a),-1)
a$paired <- substr (rownames (a),2,3)
probe <- 'nclone'

pAB <- wilcox.test(a[a$category=="A", colnames (a)%in% probe], a[a$category=="B", colnames (a)%in% probe], paired = TRUE)[3]
pBC <- wilcox.test(a[a$category=="B", colnames (a)%in% probe], a[a$category=="C", colnames (a)%in% probe], paired = TRUE)[3]
pAC <- wilcox.test(a[a$category=="A", colnames (a)%in% probe], a[a$category=="C", colnames (a)%in% probe], paired = TRUE)[3]

pvalues <- c(pAB, pBC, pAC)

outline_colors <- c("#AA3830", '#9D2E99', '#350D7F') # Box outline colors
box_fill_colors <- c('#FFE2E0', '#FFD1FF', '#DAC3FF') # Box fill colors
point_fill_colors <- c('#FF8077', '#E66CE1', '#7E33FF') # Point fill colors

pdf ('Figure 1G clone numbers across time.pdf')

ggplot(a, aes(category,a[, colnames(a)%in%probe], label=paired)) + 
geom_boxplot(color=outline_colors,  fill=box_fill_colors,  alpha = 0.5, size = 1.5, stroke = 1.2)+ 
geom_point(aes(colour=factor(category), fill = factor(category)), shape=21, size = 5) + 
scale_fill_manual(values=point_fill_colors) + 
scale_colour_manual(values=outline_colors)+
ylab('Number of clones' )+ 
geom_text_repel(min.segment.length = 1, max.overlaps = Inf, size = 3)  + 
scale_x_discrete(labels=c("A" = "Baseline", "B" = "Crosssover", "C" = "End of treatment")) +
  theme(
    axis.line = element_line(colour = "black", size = 1, linetype = "solid"),
    axis.text.x = element_text(size = 14, color = "black"), 
    axis.text.y = element_text(size = 14, color = "black"),  
    axis.title = element_text(size = 18, face = "bold", color = "black"),  
    axis.ticks = element_line(size = 1, color = "black"),  
    panel.background = element_rect(fill = 'white', color = 'black'),
    legend.position = "none",
    plot.title = element_text(color = "black") 
  ) +
  labs(x = "") +
  
  geom_signif(comparisons = list(c("A", "B"), c("B", "C"), c("A", "C")),
              annotations = sprintf("p = %.3f", pvalues),
              y_position = 9,  
              map_signif_level = FALSE,
              test = "wilcox.test", 
              step_increase = 0.07)

dev.off()

