###################################################
#######  SUPPLEMENTARY FIGURE 1 ###################
###################################################

#########################################################################
###### Supp Figure 1A ########## VIOLIN PLOTS OF UMI per cell type ######
#########################################################################
library(Seurat)
library(plyr)
library(dplyr)
library(future)
library (ggplot2)
library (ggrepel)
library (stringi)
library(tidyverse)

#### load the metadata available at EGAload("metadata_git.RData")
pdf ('SuppFigure1A violin_plot_nCount_main_ids.pdf')
df <- metadata
ggplot(df, aes(x = main.ids, y =as.numeric(nCount_RNA), fill = main.ids)) +
  geom_violin() +
  labs(x = "Sample", y = "nCount_RNA", fill = 'main.ids') +
  scale_y_log10(limits = c(1, 120000), breaks = c(1000, 2500, 5000, 10000, 20000, 40000, 80000, 120000)) +
  scale_fill_manual(values = c("#996600", "#AB9AFF", "#9A4D73", "#52AF43","#A6CEE3", "#F06C45", "#666666", "#FF99CC" , "#006633", "#FF9933" ,  "#DAA70F" , "#B294C7")) + 
  coord_cartesian(ylim = c(1, 120000)) +
  stat_summary(fun=mean, geom="point", shape=23, size=2, color="black", fill="black") +
theme( axis.line = element_line(colour = "black", size = 1, linetype = "solid")) +
theme(axis.text.x = element_text(size=14), axis.text.y = element_text(size=14)) +
theme(axis.title=element_text(size=18,face="bold")) +
theme(axis.ticks = element_line(size = 1))+
theme(panel.background = element_rect(fill = 'white', color = 'black'))+
theme(legend.position = "none")+
theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

####################################################################
###### Supp Figure 1B ########## VIOLIN PLOTS OF genes per cell type
####################################################################

#### load the metadata available at EGAload("metadata_git.RData")
pdf ('SuppFigure1B violin_plot_nFeature_main_ids.pdf')
df <- metadata
ggplot(df, aes(x = main.ids, y =as.numeric(nFeature_RNA), fill = main.ids)) +
  geom_violin() +
  labs(x = "Sample", y = "nCount_RNA", fill = 'main.ids') +
  scale_fill_manual(values = c("#996600", "#AB9AFF", "#9A4D73", "#52AF43","#A6CEE3", "#F06C45", "#666666", "#FF99CC" , "#006633", "#FF9933" ,  "#DAA70F" , "#B294C7")) + 
  stat_summary(fun=mean, geom="point", shape=23, size=2, color="black", fill="black") +
theme( axis.line = element_line(colour = "black", size = 1, linetype = "solid")) +
theme(axis.text.x = element_text(size=14), axis.text.y = element_text(size=14)) +
theme(axis.title=element_text(size=18,face="bold")) +
theme(axis.ticks = element_line(size = 1))+
theme(panel.background = element_rect(fill = 'white', color = 'black'))+
theme(legend.position = "none")+
theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

####################################################################
###### Supp Figure 1C ########## VIOLIN PLOTS OF UMI per time point
####################################################################

#### load the metadata available at EGAload("metadata_git.RData")
pdf ('SuppFigure1C violin_plot_nCount_time_point.pdf')
df <- metadata
df$Time_point <- gsub ( 'Baseline', '1mnd' ,df$Time_point )
ggplot(df, aes(x = Time_point, y =as.numeric(nCount_RNA), fill = Time_point)) +
  geom_violin() +
  labs(x = "Sample", y = "nCount_RNA", fill = 'Time_point') +
  scale_y_log10(limits = c(1, 120000), breaks = c(1000, 2500, 5000, 10000, 20000, 40000, 80000, 120000)) +
  scale_fill_manual(values = c('#FFE2E0', '#FFD1FF' , '#DAC3FF')) + 
  coord_cartesian(ylim = c(1, 120000)) +
  scale_x_discrete(labels=c("1mnd" = "Baseline", "2mnd" = "Crosssover", "4mnd" = "End of treatment")) +
  stat_summary(fun=mean, geom="point", shape=23, size=4, color="black", fill="black") +
theme( axis.line = element_line(colour = "black", size = 1, linetype = "solid")) +
theme(axis.text.x = element_text(size=14), axis.text.y = element_text(size=14)) +
theme(axis.title=element_text(size=18,face="bold")) +
theme(axis.ticks = element_line(size = 1))+
theme(panel.background = element_rect(fill = 'white', color = 'black'))+
theme(legend.position = "none")
dev.off()

####################################################################
###### Supp Figure 1D ########## VIOLIN PLOTS OF genes per time point
####################################################################

#### load the metadata available at EGAload("metadata_git.RData")
pdf ('SuppFigure1D violin_plot_nFeature_time_point.pdf')
df <- metadata
df$Time_point <- gsub ( 'Baseline', '1mnd' ,df$Time_point )
ggplot(df, aes(x = Time_point, y =as.numeric(nFeature_RNA), fill = Time_point)) +
  geom_violin() +
  labs(x = "Sample", y = "nCount_RNA", fill = 'Time_point') +
  scale_fill_manual(values = c('#FFE2E0', '#FFD1FF' , '#DAC3FF')) + 
  scale_x_discrete(labels=c("1mnd" = "Baseline", "2mnd" = "Crosssover", "4mnd" = "End of treatment")) +
  stat_summary(fun=mean, geom="point", shape=23, size=4, color="black", fill="black") +
theme( axis.line = element_line(colour = "black", size = 1, linetype = "solid")) +
theme(axis.text.x = element_text(size=14), axis.text.y = element_text(size=14)) +
theme(axis.title=element_text(size=18,face="bold")) +
theme(axis.ticks = element_line(size = 1))+
theme(panel.background = element_rect(fill = 'white', color = 'black'))+
theme(legend.position = "none")
dev.off()

###################################################
#######  SUPPLEMENTARY FIGURE 2 ###################
###################################################

##############################################################
###### Supp Figure 2A ##########" FEATURE PLOT of MARKER genes
##############################################################

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

#### Load Seurat object obtained from PP_Seurat_object
load("Scaled_NeoLetExe_final_git.RData")

Idents(Scaled) <- "main.ids"
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

pdf('SuppFigure2A colored_feature_plots.pdf')
FeaturePlot(Scaled, features = c("EPCAM"), order = T,  cols =   c("grey89", "darkcyan"), min.cutoff = 0.5, repel=F, label.size=4.5, raster.dpi = c(2048, 2048), raster=T)
FeaturePlot(Scaled, features = c("CD3D"), order = T,  cols =   c("grey89", "#a085b3"), repel=F, label.size=4.5, raster.dpi = c(2048, 2048), raster=T)
FeaturePlot(Scaled, features = c("C1R"), order = T,  cols =   c("grey89", "#F06C45"), min.cutoff = 0.5, repel=F, label.size=4.5, raster.dpi = c(2048, 2048), raster=T)
FeaturePlot(Scaled, features = c("CD68"), order = F,  cols =   c("grey89", "#e589b7"), repel=F, label.size=4.5, raster.dpi = c(2048, 2048), raster=T)
FeaturePlot(Scaled, features = c("MYH11"), order = T,  cols =   c("grey89", "#DAA70F"), min.cutoff=0, repel=F, label.size=4.5, raster.dpi = c(2048, 2048), raster=T)
FeaturePlot(Scaled, features = c("VWF"), order = T,  cols =   c("grey89", "#52AF43"), min.cutoff=0, repel=F, label.size=4.5, raster.dpi = c(2048, 2048), raster=T)
FeaturePlot(Scaled, features = c("CPA3"), order = T,  cols =   c("grey89", "#5b5b5b"), min.cutoff=0, repel=F, label.size=4.5, raster.dpi = c(2048, 2048), raster=T)
FeaturePlot(Scaled, features = c("CD19"), order = T,  cols =   c("grey89", "#998ae5"), min.cutoff=0, repel=F, label.size=4.5, raster.dpi = c(2048, 2048), raster=T)
FeaturePlot(Scaled, features = c("MZB1"), order = F,  cols =   c("grey89", "#FF9933"), min.cutoff=0, repel=F, label.size=4.5, raster.dpi = c(2048, 2048), raster=T)
dev.off()

#################################################################
###### Supp Figure 2B ########## BARPLOTS of celltypes per samples
#################################################################

load("metadata_git.RData")
#### load the metadata available at EGA
load("palopal.RData")
#### Load an object containing color informations

dta <- as.data.frame.matrix(prop.table(table (metadata$orig.ident, metadata$main.ids),1)*100)
dta <- dta[, rev(order(dta[1,]))]


blank <- rep (0,12)
dta <- rbind(dta[1:3, ], blank, dta[4:6, ], blank, dta[7:9, ], blank,  dta[10, ], blank, dta[11:13, ], blank,dta[14:16, ], blank, dta[17:19, ], blank,
dta[20:22, ], blank, dta[23:25, ], blank, dta[26:28, ], blank, dta[29:31, ], blank,dta[32:34, ], blank, dta[35:37, ], blank,dta[38:40, ], blank, dta[41:43, ], blank,
dta[44:46, ], blank, dta[47:49, ], blank,dta[50:52, ], blank, dta[53:55, ], blank,dta[56:58, ], blank, dta[59:61, ], blank,dta[62:64, ], blank, dta[65:67, ], blank,
dta[68:70, ], blank, dta[71:73, ])

pdf('SuppFigure2B_barplots_main.ids_proportions.pdf', width = 24, height = 6)  
par(mar=c(4, 3, 1, 10))  
barplot(t(dta), legend.text = TRUE, col = palopal[colnames(dta)], 
        args.legend = list(x = "topright", inset = c(-0.42, 0)), 
        cey.axis=0.3, las=2)
dev.off()




#################################################################
###### Supp Figure 2C ########## BOXPLOTS of NORMAL EPITHELIAL
#################################################################
library(Seurat)
library(plyr)
library(dplyr)
library(future)
library (ggplot2)
library (ggrepel)
library (stringi)
library(tidyverse)
library(ggsignif)

load("metadata_git.RData")
#### load the metadata available at EGA

dta <- as.data.frame.matrix(prop.table(table (metadata$orig.ident, metadata$main.ids.cor),1)*100)
tabix <- as.data.frame.matrix(table (metadata$orig.ident, metadata$main.ids.cor))
pro_norm <- tabix[, 10]/( tabix[, 10] +tabix[, 6])*100

dta$category <- stri_sub(rownames (dta),-1)
dta$paired <- substr (rownames (dta),3,5)
dta$paired <- gsub ("A", "",dta$paired )
dta$paired <- gsub ("B", "",dta$paired )
dta$paired <- gsub ("C", "",dta$paired )
annot <- metadata[, colnames (metadata)%in% c("patient_number", "orig.ident", "resp", "met", "histo", "her2", "PGR", "KI67", "Arm", "resp_ext", "sc_resp")]
annot <- annot[!duplicated (metadata$orig.ident), ]
table (annot$orig.ident== rownames (dta))
dta <- cbind (dta , annot, pro_norm)

dta <- dta[!rownames(dta)=='P04A',]


probe <- 'pro_norm'
pAB <- wilcox.test(dta[dta$category=="A", colnames (dta)%in% probe], dta[dta$category=="B", colnames (dta)%in% probe], paired = TRUE)[3]
pBC <- wilcox.test(dta[dta$category=="B", colnames (dta)%in% probe], dta[dta$category=="C", colnames (dta)%in% probe], paired = TRUE)[3]
pAC <- wilcox.test(dta[dta$category=="A", colnames (dta)%in% probe], dta[dta$category=="C", colnames (dta)%in% probe], paired = TRUE)[3]
pvalues <- c(pAB, pBC, pAC)

outline_colors <- c("#AA3830", '#9D2E99', '#350D7F') 
box_fill_colors <- c('#FFE2E0', '#FFD1FF', '#DAC3FF') 
point_fill_colors <- c('#FF8077', '#E66CE1', '#7E33FF') 


pdf('SuppFigure2C norm_epith_boxplot_time_point.pdf')
a <- ggplot(dta, aes(category, dta[, colnames(dta) %in% probe], label = paired)) + 
 geom_boxplot(color = outline_colors, 
               fill = box_fill_colors, 
               alpha = 0.5, size = 1.5) + 
  geom_point(aes(fill = factor(category), colour = factor(category)), 
             shape = 21, size = 5, stroke = 1.2) +  
  scale_fill_manual(values = point_fill_colors) +
  scale_colour_manual(values = outline_colors) + 
  ylab('Percentage of normal epithelial cells') + 
  geom_text_repel(min.segment.length = 1, max.overlaps = Inf) + 
  scale_x_discrete(labels = c("A" = "A", "B" = "B", "C" = "C")) +
  theme(
    axis.line = element_line(colour = "black", size = 1, linetype = "solid"), 
    axis.text.x = element_text(size = 14, color = "black"), 
    axis.text.y = element_text(size = 14, color = "black"), 
    axis.title = element_text(size = 18, face = "bold", color = "black"), 
    axis.ticks = element_line(size = 1, color = "black"),
    panel.background = element_rect(fill = 'white', color = 'black'),
    plot.title = element_text(size = 18, face = "bold", color = "black"),
    legend.text = element_text(color = "black"),
    legend.title = element_text(color = "black"),
    legend.position = "none") +
  labs(x = "") +
  geom_signif(comparisons = list(c("A", "B"), c("B", "C"), c("A", "C")),
              annotations = sprintf("p = %.3f", pvalues), 
              y_position = 90,  
              map_signif_level = FALSE,
              test = "wilcox.test", 
              step_increase = 0.05)
print(a)
dev.off()


###################################################
#######  SUPPLEMENTARY FIGURE 3 ###################
###################################################

########################################################################
##### Supplementary Figure 3A ###### phylogenetic trees and fishplots ##
########################################################################
### need the input file from PP_InferCNV_per_sample
cd /Uphyloplot2/
module load Python/3.9.5-GCCcore-10.3.0
python uphyloplot2.py
sed '/^norm_/d' <  NL67.cell_groupings > NL67i.cell_groupings
cd ..
module load Python/3.9.5-GCCcore-10.3.0
python uphyloplot2.py


###################################################
##### Supplementary Figure 3B ###### volcano plot
###################################################

library(EnhancedVolcano)
test <- data.frame(read.table("markergenes_all_sensi_vs_all_resis_clone_cat2.txt", header=T, sep="\t", stringsAsFactors=F))
#correlation_matrix output of PP_differential_expression_clones

test$p_val_adj[test$p_val_adj == 0] <- 1.624931e-302
test$p_val_adj_transformed <- -log10(test$p_val_adj)
test <- test[test$cluster == 'res',]

logFC_cutoff <- 0.5
p_val_cutoff <- -log10(10e-32)

volcano <- EnhancedVolcano(
  test,
  lab = test$gene,
  x = "avg_log2FC",
  y = "p_val_adj",
  pointSize = 3.0,
  labSize = 3,
  parseLabels = FALSE, 
  pCutoff = p_val_cutoff, 
  FCcutoff = logFC_cutoff,
  max.overlaps = 0,  
  drawConnectors = TRUE, 
  widthConnectors = 0.75
) +
  geom_point(data = subset(test, abs(avg_log2FC) > logFC_cutoff & p_val_adj_transformed > p_val_cutoff),
             aes(x = avg_log2FC, y = p_val_adj_transformed, color = avg_log2FC > 0),
             size = 3) +
  scale_color_manual(values = c("TRUE" = "#A02629", "FALSE" = "#7BA58C")) + 
  geom_point(data = subset(test, gene %in% c('TFF1', 'SERPINA3', 'CCND1', 'IFI30', 'B2M', 'CD74', 'MAPK13', 'ERBB2', 'CLU')),
             aes(x = avg_log2FC, y = p_val_adj_transformed),
             color = c("purple", "purple", "purple", 'dark green', 'dark green', 'dark green', 'orange', 'orange', 'orange'), size = 3) +
  ggrepel::geom_text_repel(data = subset(test, gene %in% c('TFF1', 'SERPINA3', 'CCND1', 'IFI30', 'B2M', 'CD74', 'MAPK13', 'ERBB2', 'CLU')),
                           aes(x = avg_log2FC, y = p_val_adj_transformed, label = gene),
                           color = c("purple", "purple", "purple", 'dark green', 'dark green', 'dark green', 'orange', 'orange', 'orange'), size = 3) +
  labs(
    x = "log FC (resistant vs sensitive)",   
    y = "-log (adj-p-value)"                
  ) + 
  theme(
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),  
    plot.title = element_blank(),        
    plot.subtitle = element_blank(),     
    plot.caption = element_blank()       
  )
  
print(volcano)

pdf('SuppFigure3B volcano plots.pdf')
print(volcano)
dev.off()

#########################################################
##### Supplementary Figure 3C ###### Barplot sensitive ##
#########################################################

#correlation_matrix output of PP_differential_expression_clones followed by GSEA analysis
pdf('SuppFigure3C barplots GSEA sensitive.pdf')
sens <- data.frame(read.table("sensitive_hallmark_plot.txt", header=T, sep="\t", stringsAsFactors=F))
library(ggplot2)
library(scales)  

sens$Gene.set <- factor(sens$Gene.set, levels = sens$Gene.set[order(sens$X.log..p.value.)])

ggplot(sens, aes(x = Gene.set, y = X.log..p.value.)) +
  geom_bar(stat = "identity", fill = alpha("dark green", 0.3), color = "dark green", size = 1.5) + 
  coord_flip() + 
  labs(x = "Hallmark Gene Set", y = "-log (p-value)", title = "sensitive clones") +
  theme_minimal() + 
  theme(
    panel.background = element_rect(fill = "white"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.text = element_text(color = "black"),  
    axis.title = element_text(color = "black"), 
    axis.ticks.x = element_line(color = "black"), 
    plot.title = element_text(color = "dark green", size = 14, hjust = 0.5),  
    axis.text.x = element_text(size = 12), 
    axis.title.x = element_text(size = 14),  
    plot.subtitle = element_blank(),
    plot.caption = element_blank()
  )
dev.off()

#########################################################
##### Supplementary Figure 3D ###### Barplot resistant ##
#########################################################

#correlation_matrix output of PP_differential_expression_clones followed by GSEA analysis

pdf('SuppFigure3D barplots GSEA resistant.pdf')
res <-data.frame(read.table("resistant_hallmark_plot.txt", header=T, sep="\t", stringsAsFactors=F))

library(ggplot2)
library(scales)  

res$Gene.set <- factor(res$Gene.set, levels = res$Gene.set[order(res$X.log..p.value.)])

ggplot(res, aes(x = Gene.set, y = X.log..p.value.)) +
  geom_bar(stat = "identity", fill = alpha("dark red", 0.3), color = "dark red", size = 1.5) +  # More transparency
  coord_flip() + 
  labs(x = "Hallmark Gene Set", y = "-log (p-value)", title = "resistant clones") +
  theme_minimal() + 
  theme(
    panel.background = element_rect(fill = "white"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.text = element_text(color = "black"),  
    axis.title = element_text(color = "black"), 
    axis.ticks.x = element_line(color = "black"), 
    plot.title = element_text(color = "dark red", size = 14, hjust = 0.5),
    axis.text.x = element_text(size = 12),
    axis.title.x = element_text(size = 14),
    plot.subtitle = element_blank(),
    plot.caption = element_blank()
  )

dev.off()

###################################################
#######  SUPPLEMENTARY FIGURE 4 ###################
###################################################

#########################################################
##### Supplementary Figure 4 ###### BOXPLOTS CLONES and SCORES ##
#########################################################

library(Seurat)
library(ggplot2)
library(ggrepel)
library (ggsignif)

load("metadata_git.RData")
#### load the metadata available at EGA
metadata <- data.frame(read.table("maligmant_cels_metadata_with_modules.txt", header=T, sep="\t", stringsAsFactors=F))
#scores generated  in PP_gene_cluster_correlation

clones_ord <-c( unique (metadata$clone_id[metadata$clone_cat2=="res" &!is.na(metadata$clone_cat2)]), unique (metadata$clone_id[metadata$clone_cat2=="sens" &!is.na(metadata$clone_cat2)]))
clones_cat <-c( rep ('res' , length (unique (metadata$clone_id[metadata$clone_cat2=="res" &!is.na(metadata$clone_cat2)]))), rep ('sens', length (unique (metadata$clone_id[metadata$clone_cat2=="sens" &!is.na(metadata$clone_cat2)]))))


mapping_df <- data.frame( clone_id = clones_ord, clone_cat2 = clones_cat)
library(dplyr)
metadata <- metadata %>%left_join(mapping_df, by = "clone_id")

cluster_name <- 'cluster10'

  average <- aggregate(get(cluster_name) ~ clone_id, data = metadata, FUN = mean)
  colnames (average)[2] <- paste(cluster_name)
  status <- metadata[!duplicated(metadata$clone_id), c("clone_id", "clone_cat2.x")]
  
  average <- average[order(average$clone_id), ]
  status <- status[order(status$clone_id), ]
  print(all(average$clone_id == status$clone_id))
  average$clone_cat2 <- status$clone_cat2.x
  average <- average[!average$clone_cat2%in%c('excl', 'neut'),]
  
  pval <- t.test(get(cluster_name) ~ clone_cat2, data = average)[3]

outline_colors <- c("#A02629", '#7BA58C') 
box_fill_colors <- c('#E5A1A4', '#C3EAD6') 
point_fill_colors <- c('#C32E32', '#9ED3B9') 

  plot <- ggplot(average, aes(x = clone_cat2, y = average[, colnames(average)%in%cluster_name])) + 
  geom_boxplot(aes(fill = clone_cat2), color = outline_colors, fill = box_fill_colors, alpha = 0.2, size = 1.5) + 
  geom_point(aes(colour = clone_cat2, fill = clone_cat2), , shape = 21, size = 5, stroke = 1.2) +  # Fixed parentheses here
  scale_fill_manual(values= point_fill_colors) + 
  scale_colour_manual(values=outline_colors) +  
    ylab(paste (cluster_name, ' AddModuleScore')) + 
    scale_x_discrete(labels = c(res = 'resistant', sens = 'sensitive')) +
    theme(axis.line = element_line(colour = "black", size = 1, linetype = "solid")) +
    theme(axis.text.x = element_text(size = 14, color = "black"), axis.text.y = element_text(size = 14, color = "black")) +
    theme(axis.title = element_text(size = 18, face = "bold")) +
    theme(axis.ticks = element_line(size = 1)) +
    theme(panel.background = element_rect(fill = 'white', color = 'black')) +
	labs(x = '') +
    theme(legend.position = "none") +
	  geom_signif(comparisons = list(c('res', 'sens')),
               annotations = sprintf("p = %.3f", pval),
               map_signif_level = FALSE,
               test = "t.test", 
               step_increase = 0.07) 
			   
pdf('SuppFigure4 module scores in response.pdf')
    print (plot)
dev.off()

###################################################
####### SUPPLEMENTARY FIGURE 5 ####################
###################################################
#### PP_subclustering.r  ##########################
