
###################################################
#######  EXTENDED DATA FIGURE 1 ###################
###################################################
##########################################################################
###### Extended data Figure 1A ########## FEATURE PLOT of MARKER genes ###
##########################################################################

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

pdf('ExDaFigure1A colored_feature_plots.pdf')
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

################################################################################
###### Extended data Figure 1B ########## BARPLOTS of celltypes per samples ####
################################################################################

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

pdf('ExDaFigure1B_barplots_main.ids_proportions.pdf', width = 24, height = 6)  
par(mar=c(4, 3, 1, 10))  
barplot(t(dta), legend.text = TRUE, col = palopal[colnames(dta)], 
        args.legend = list(x = "topright", inset = c(-0.42, 0)), 
        cey.axis=0.3, las=2)
dev.off()


##########################################################################
###### Extended data Figure 1C ########## BOXPLOTS of NORMAL EPITHELIAL ###
##########################################################################
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


pdf('ExdaFigure1C norm_epith_boxplot_time_point.pdf')
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
#######  EXTENDED DATA FIGURE 2 ###################
###################################################

########################################################################
##### Extended data Figure 2A ###### phylogenetic trees and fishplots ##
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
##### Extended data Figure 2B ###### volcano plot
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

pdf('ExDaFigure2B volcano plots.pdf')
print(volcano)
dev.off()

#########################################################
##### Extended data Figure 2C ###### Barplot sensitive ##
#########################################################

#correlation_matrix output of PP_differential_expression_clones followed by GSEA analysis
pdf('ExDaFigure2C barplots GSEA sensitive.pdf')
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
##### Extended data Figure 2D ###### Barplot resistant ##
#########################################################

#correlation_matrix output of PP_differential_expression_clones followed by GSEA analysis

pdf('ExDaFigure2D barplots GSEA resistant.pdf')
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
####### EXTENDED DATA FIGURE 3 ####################
###################################################
#### PP_subclustering.r  ##########################


###################################################
####### EXTENDED DATA FIGURE 4 ####################
###################################################
#################################################################################################
###### EXTENDED DATA Figure 4 ########## Milo minor IDs SWARMPLOTS ##########
#################################################################################################
library(Seurat)
library(plyr)
library(dplyr)
library(future)
library (ggplot2)
library (ggrepel)
library (stringi)
library(tidyverse)
library(miloR)

###############################
### need outputs of PP_milo_R## 
###############################  
## final.id plot A vs C #######  
###############################

A= "#DA3A09"
B= "#9C2C98"
C= "#320E77"
da_results_B_4_rev <- as.data.frame(data.frame(read.table("finalID_AC.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)))

vector_r <- rev(data.frame(da_results_B_4_rev %>% group_by(final.ids) %>% summarize(avg_logFC = mean(logFC)) %>% arrange(avg_logFC)))
da_results_B_4_rev$final.ids <- factor(da_results_B_4_rev$final.ids, levels = rev(vector_r$final.ids))
da_results_B_4_rev <- da_results_B_4_rev[order(match(da_results_B_4_rev$final.ids, vector_r$final.ids)),]

p <- plotDAbeeswarm(da_results_B_4_rev, group.by = "final.ids", alpha=0.05) + 
scale_colour_gradient2(high = C, low = A, na.value = "grey85", space = "Lab") +  
  xlab("") + 
  scale_y_continuous(breaks = c(-5, -2.5, 0, 2.5, 5), limits = c(-6, 6.2)) +
  theme(axis.text.y = element_text(colour = (c("orange","#4f4fc5","#4f4fc5",
"#4f4fc5","#4f4fc5","orange",
"#4f4fc5","#8ba831","#4f4fc5",
"orange","#4f4fc5","#4f4fc5",
"orange","#4f4fc5","orange",
"#4f4fc5","#4f4fc5", "#8ba831",
"#8ba831","orange","#8ba831",
"#8ba831","#8ba831","#4f4fc5",
"#4f4fc5","#8ba831","#8ba831",
"#8ba831","orange","#8ba831",
"#8ba831","orange","#8ba831",
"#8ba831","#8ba831","orange",
"orange","#8ba831","#8ba831",
"orange","orange","orange",
"#8ba831","orange","#8ba831"))),
axis.text.x = element_text(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank())

da_results_B_4_sig_1 <- da_results_B_4_rev[da_results_B_4_rev$FDR <0.05 & da_results_B_4_rev$logFC <0,]
da_results_B_4_sig_2 <- da_results_B_4_rev[da_results_B_4_rev$FDR <0.05 & da_results_B_4_rev$logFC >0,]
count_sig_1 <- as.data.frame(table(da_results_B_4_sig_1$final.ids))
colnames(count_sig_1) <- c("final.ids", "count_sig_1")
count_sig_2 <- as.data.frame(table(da_results_B_4_sig_2$final.ids))
colnames(count_sig_2) <- c("final.ids", "count_sig_2")
count_data <- merge(count_sig_1, count_sig_2, by = "final.ids", all = TRUE)
count_data[is.na(count_data)] <- 0


p <- p + mapply(function(x, plus_label, minus_label) { 
    list(annotate("text", x = x, y = 6, label = paste0("(", plus_label, ")"), hjust = 1, size = 3),
      annotate("text", x = x, y = -5.5, label = paste0("(", minus_label, ")"), hjust = 1, size = 3))}, 
  seq_along(count_data$final.ids), 
  count_data$count_sig_2[rev(match(vector_r$final.ids, count_data$final.ids))], 
  count_data$count_sig_1[rev(match(vector_r$final.ids, count_data$final.ids))], 
  SIMPLIFY = FALSE)

pdf('ExDaFigure4a swarmplot final ids A C.pdf' , width = 7, height = 10)
print(p)
dev.off()


########################### 
## final.id plot A vs B ### 
###########################

da_results_B_4_rev <- as.data.frame(data.frame(read.table("finalID_AB.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)))
da_results_B_4_rev$final.ids <- factor(da_results_B_4_rev$final.ids, levels = rev(vector_r$final.ids))
da_results_B_4_rev <- da_results_B_4_rev[order(match(da_results_B_4_rev$final.ids, vector_r$final.ids)),]

p <- plotDAbeeswarm(da_results_B_4_rev, group.by = "final.ids", alpha=0.05) + 
scale_colour_gradient2(high = B, low = A, na.value = "grey85", space = "Lab") +  
  xlab("") + 
  scale_y_continuous(breaks = c(-5, -2.5, 0, 2.5, 5), limits = c(-6, 6.2)) +
  theme(axis.text.y = element_text(colour = (c("orange","#4f4fc5","#4f4fc5",
"#4f4fc5","#4f4fc5","orange",
"#4f4fc5","#8ba831","#4f4fc5",
"orange","#4f4fc5","#4f4fc5",
"orange","#4f4fc5","orange",
"#4f4fc5","#4f4fc5", "#8ba831",
"#8ba831","orange","#8ba831",
"#8ba831","#8ba831","#4f4fc5",
"#4f4fc5","#8ba831","#8ba831",
"#8ba831","orange","#8ba831",
"#8ba831","orange","#8ba831",
"#8ba831","#8ba831","orange",
"orange","#8ba831","#8ba831",
"orange","orange","orange",
"#8ba831","orange","#8ba831"))),
axis.text.x = element_text(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank())

da_results_B_4_sig_1 <- da_results_B_4_rev[da_results_B_4_rev$FDR <0.05 & da_results_B_4_rev$logFC <0,]
da_results_B_4_sig_2 <- da_results_B_4_rev[da_results_B_4_rev$FDR <0.05 & da_results_B_4_rev$logFC >0,]
count_sig_1 <- as.data.frame(table(da_results_B_4_sig_1$final.ids))
colnames(count_sig_1) <- c("final.ids", "count_sig_1")
count_sig_2 <- as.data.frame(table(da_results_B_4_sig_2$final.ids))
colnames(count_sig_2) <- c("final.ids", "count_sig_2")
count_data <- merge(count_sig_1, count_sig_2, by = "final.ids", all = TRUE)
count_data[is.na(count_data)] <- 0


p <- p + mapply(function(x, plus_label, minus_label) { 
    list(annotate("text", x = x, y = 6, label = paste0("(", plus_label, ")"), hjust = 1, size = 3),
      annotate("text", x = x, y = -5.5, label = paste0("(", minus_label, ")"), hjust = 1, size = 3))}, 
  seq_along(count_data$final.ids), 
  count_data$count_sig_2[rev(match(vector_r$final.ids, count_data$final.ids))], 
  count_data$count_sig_1[rev(match(vector_r$final.ids, count_data$final.ids))], 
  SIMPLIFY = FALSE)

pdf('ExDaFigure4b swarmplot final ids A B.pdf', width = 7, height = 10)
print(p)
dev.off()


########################### 
## final.id plot B vs C ### 
########################### 

da_results_B_4_rev <- as.data.frame(data.frame(read.table("finalID_BC.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)))
da_results_B_4_rev$final.ids <- factor(da_results_B_4_rev$final.ids, levels = rev(vector_r$final.ids))
da_results_B_4_rev <- da_results_B_4_rev[order(match(da_results_B_4_rev$final.ids, vector_r$final.ids)),]

p <- plotDAbeeswarm(da_results_B_4_rev, group.by = "final.ids", alpha=0.05) + 
scale_colour_gradient2(high = C, low = B, na.value = "grey85", space = "Lab") +  
  xlab("") + 
  scale_y_continuous(breaks = c(-5, -2.5, 0, 2.5, 5), limits = c(-6, 6.2)) +
  theme(axis.text.y = element_text(colour = (c("orange","#4f4fc5","#4f4fc5",
"#4f4fc5","#4f4fc5","orange",
"#4f4fc5","#8ba831","#4f4fc5",
"orange","#4f4fc5","#4f4fc5",
"orange","#4f4fc5","orange",
"#4f4fc5","#4f4fc5", "#8ba831",
"#8ba831","orange","#8ba831",
"#8ba831","#8ba831","#4f4fc5",
"#4f4fc5","#8ba831","#8ba831",
"#8ba831","orange","#8ba831",
"#8ba831","orange","#8ba831",
"#8ba831","#8ba831","orange",
"orange","#8ba831","#8ba831",
"orange","orange","orange",
"#8ba831","orange","#8ba831"))),
axis.text.x = element_text(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank())

da_results_B_4_sig_1 <- da_results_B_4_rev[da_results_B_4_rev$FDR <0.05 & da_results_B_4_rev$logFC <0,]
da_results_B_4_sig_2 <- da_results_B_4_rev[da_results_B_4_rev$FDR <0.05 & da_results_B_4_rev$logFC >0,]
count_sig_1 <- as.data.frame(table(da_results_B_4_sig_1$final.ids))
colnames(count_sig_1) <- c("final.ids", "count_sig_1")
count_sig_2 <- as.data.frame(table(da_results_B_4_sig_2$final.ids))
colnames(count_sig_2) <- c("final.ids", "count_sig_2")
count_data <- merge(count_sig_1, count_sig_2, by = "final.ids", all = TRUE)
count_data[is.na(count_data)] <- 0


p <- p + mapply(function(x, plus_label, minus_label) { 
    list(annotate("text", x = x, y = 6, label = paste0("(", plus_label, ")"), hjust = 1, size = 3),
      annotate("text", x = x, y = -5.5, label = paste0("(", minus_label, ")"), hjust = 1, size = 3))}, 
  seq_along(count_data$final.ids), 
  count_data$count_sig_2[rev(match(vector_r$final.ids, count_data$final.ids))], 
  count_data$count_sig_1[rev(match(vector_r$final.ids, count_data$final.ids))], 
  SIMPLIFY = FALSE)
  
pdf('ExDaFigure4c  swarmplot minor ids B C.pdf', width = 7, height = 10)
print(p)
dev.off()


################################################
###### EXTENDED DATA FIGURE Figure 5 ##########
################################################
###########################################################################
# Extended data Figure 5 A CORRELATION BETWEEN CLINICAL VARIABLES ########
###########################################################################

library(pheatmap)
library (tidyr)
library (scRepertoire)
library (stringi)
library (stringr)
library (ggplot2)
library (ggrepel)
library(ggsignif)
library  (dplyr)
library(patchwork)
library(RColorBrewer)
library(ggalluvial)
library(pheatmap)
library(amap)
library (gwasrapidd)
library (ggplot2)
library (corrplot)

df <- data.frame(read.table('raw_binary response.txt', header=T, sep="\t", stringsAsFactors=F))
#### load cinical data from Supplementary Table 1
mat <- as.matrix(df)

cor_sco <- cor (mat, method='pearson', use ='na.or.complete')
dcols <- hcluster(t(cor_sco), method = 'correlation', link = 'ward')
drows <- hcluster(cor_sco, method = 'correlation', link = 'ward')

breaksList = seq(-1, 1, by = 0.05)
col.pal <- colorRampPalette(c("dodgerblue4","white","brown3"))(length(breaksList))

pheatmap(cor_sco,  cutree_cols=4, color = col.pal, show_rownames = T, show_colnames = T, main = "",  cluster_rows = drows, cluster_cols = dcols, breaks = breaksList)

p_mat <- matrix(NA_real_,
                nrow = ncol(mat),
                ncol = ncol(mat),
                dimnames = list(colnames(mat), colnames(mat)))

for (i in seq_len(ncol(mat))) {
  for (j in seq_len(ncol(mat))) {
    if (j <= i) next  

    x <- mat[, i]
    y <- mat[, j]
    ok <- !is.na(x) & !is.na(y)
    if (sum(ok) < 3) next  

    ct <- suppressWarnings(
      cor.test(x[ok], y[ok],
               method = "pearson",
               exact  = TRUE) 
    )

    p_mat[i, j] <- p_mat[j, i] <- ct$p.value
  }
}
diag(p_mat) <- 0

cor_sco <- cor_sco[drows$order, dcols$order]
p_mat <- p_mat[drows$order, dcols$order]
corrplot(cor_sco, is.corr=F, col=col.pal, col.lim=c(-1, 1), tl.col = 'black')
corrplot(cor_sco, is.corr=F, col=col.pal, col.lim=c(-1, 1), tl.col = 'black', p.mat =p_mat , sig.level = 0.2, insig='blank')

pdf ('ExDaFigure5A corrplot between clinical variables.pdf')
pheatmap(cor_sco,  cutree_cols=4, color = col.pal, show_rownames = T, show_colnames = T, main = "",  cluster_rows = drows, cluster_cols = dcols, breaks = breaksList)
corrplot(cor_sco, is.corr=F, col=col.pal, col.lim=c(-1, 1), tl.col = 'black')
corrplot(cor_sco, is.corr=F, col=col.pal, col.lim=c(-1, 1), tl.col = 'black', p.mat =p_mat , sig.level = 0.2, insig='blank')
dev.off()

####################################################################################
# Extended data Figure 5 B ALL CLINICAL RESPONSE MEASURE WITH IMMUNE CLUSTERS ######
####################################################################################

library(Seurat)
library(plyr)
library(dplyr)
library(future)
library(ggplot2)
library(ggrepel)
library(stringi)
library(mixtools)
library(pheatmap)
library(amap)
library(gwasrapidd)
library(RColorBrewer)

load("metadata_git.RData")
##### needs the Seurat_metadata
df <- data.frame(read.table('raw_binary response.txt', header=T, sep="\t", stringsAsFactors=F))
#### load cinical data from Supplementary Table 1

mer <- merge(metadata, df, by.x='patient_number', by.y="Patient.no.", all.x=T)
metadata <- mer 

metadata1<- metadata[!metadata$main.ids=="Epithelial_cells" & !metadata$final.ids=="UNA",]
dta <- as.data.frame.matrix(prop.table(table (metadata1$orig.ident, metadata1$final.ids),1)*100)
dta1 <- dta

dta$category <- stri_sub(rownames (dta),-1)
dta$paired <- substr (rownames (dta),2,4)
dta$paired <- gsub ("A", "",dta$paired )
dta$paired <- gsub ("B", "",dta$paired )
dta$paired <- gsub ("C", "",dta$paired )

dta$category <- stri_sub(rownames (dta),-1)
annot <- metadata[, colnames (metadata) %in% c("orig.ident",
  'size_bin_20',        	
  'Prosigna_bin',                  
  'Ki67_bin_10',               	
  'Node.bin',                  	
  'scRNA_percentage_bin',          
  'scRNA_KI_bin',                 
  'perc_decrease_bin',    			
  'Metastasis_bin',          		
  'scRNA_decrease_bin85',        
  'grade_bin',                    
  'PR_bin',
  'residual_tumor_bin')]
annot <- annot[!duplicated (metadata$orig.ident), ]
table (annot$orig.ident == rownames (dta))
dta <- cbind (dta , annot)
a <- as.data.frame.matrix(prop.table (table (metadata$orig.ident, metadata$clone_cat),1)*100)
a$category <- stri_sub(rownames (a),-1)
a1 <- a[a$category == "A",]

dta1 <- dta1[dta$category == "A",]
dta <- dta[dta$category == "A",]

  x <- as.character(x)
annotation_row <- dta[,colnames (dta) %in% c(
  'residual_tumor_bin',
  'size_bin_20',        	
  'Prosigna_bin',                  
  'Ki67_bin_10',               	
  'Node.bin',                  	
  'scRNA_percentage_bin',          
  'scRNA_KI_bin',                 
  'perc_decrease_bin',    			
  'Metastasis_bin',          		
  'scRNA_decrease_bin85',        
  'grade_bin',                    
  'PR_bin' 					
)]


anno_colors <- list(
  residual_tumor_bin              = c("1" = "#326A43", "2" = "#A02629" ),
  size_bin_20                	= c("1" = "#E59EDD", "2" = "#78206E"),
  Prosigna_bin                   = c("1" = "#CC9900", "2" = "#996633"),
  Ki67_bin_10               	 = c("1" = "#F6C6AD", "2" = "#C04F15"),
  Node.bin                  	= c("1" = "#5CA46A", "2" = "#C55A5A"),
  scRNA_percentage_bin           = c("1" = "#326A43", "2" = "#A02629"),
  scRNA_KI_bin                  = c("1" = "#F6C6AD", "2" = "#C04F15"),
  perc_decrease_bin    			= c("1" = "#92D050", "2" = "#FF7C80"),
  Metastasis_bin          		= c("1" = "#5CA46A", "2" = "#C55A5A"),
  scRNA_decrease_bin85        	 = c("1" = "#92D050", "2" = "#FF7C80"),
  grade_bin                     = c("1" = "#5CA46A", "2" = "#C55A5A"),
  PR_bin 						= c("1" = "#5CA46A", "2" = "#C55A5A"),
  baseline_subtype 				= c("LumA" = "blue", "LumB" = "light blue", "Her2" = "deep pink"), 
  eot_subtype 				    = c("LumA" = "blue", "LumB" = "light blue", "Her2" = "deep pink")
)

# Adjusting color palette and breaks
max_val <- 20
breaks_vals <- c(seq(0, max_val, length.out = 101))
breaksList = seq(0, 35, by = 0.05)
colors <- colorRampPalette(c("white", "indianred1", "brown3", "brown4"))(length(breaksList))
col.pal <- colorRampPalette(c("dodgerblue4", "white", "white", "brown3", "brown3", "brown3"))(length(breaksList))

dcols <- hcluster(t(dta1), method = 'spearman', link = 'complete')
drows <- hcluster(dta1, method = 'spearman', link = 'complete')

# Save heatmap to a file with specific size to avoid cropping the legend
p<- pheatmap(
  t(dta1), 
  show_rownames = TRUE, show_colnames = TRUE,  
  cluster_rows = dcols, cluster_cols = drows, fontsize = 9,
  annotation_col = annotation_row,
  annotation_colors = anno_colors,  
  col = colors, 
   width = 10, height = 8, cutree_cols = 3 , border_color = NA
)

pdf ('ExDaFigure5B immune clusters and clinical variables.pdf')
print (p)
dev.off()



########################################
###### EXTENDED DATA Figure 6 ##########
########################################
######################## ##############################################
#### Extended data Figure 6 A Naive B cells BOXPLOTS ##################################################
######################################################################

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
met_clin <- data.frame(read.table("clinical_data.txt", header=T, sep="\t", stringsAsFactors=F))
##### load clinical data in Supp Table 1
mer <- merge(metadata, met_clin, by.x='patient_number', by.y="patient_number", all.x=T)
metadata <- mer

ids <- "final.ids"
subs <- "ME"
gory <- "perc_malig_surg10"
cate <- "A"
probe <- "B_cells_naive"

metadata1<- metadata[!metadata$main.ids=="Epithelial_cells" & !metadata$final.ids=="UNA",]
dta <- as.data.frame.matrix(prop.table(table (metadata1$orig.ident, metadata1[, colnames(metadata1)%in%ids]),1)*100)
dta1 <- dta

dta$category <- stri_sub(rownames (dta),-1)
dta$paired <- substr (rownames (dta),2,4)
dta$paired <- gsub ("A", "",dta$paired )
dta$paired <- gsub ("B", "",dta$paired )
dta$paired <- gsub ("C", "",dta$paired )

dta$category <- stri_sub(rownames (dta),-1)
annot <- metadata[, colnames (metadata)%in% c("orig.ident", "perc_malig_surg10")]

annot <- annot[!duplicated (metadata$orig.ident), ]
table (annot$orig.ident== rownames (dta))
dta <- cbind (dta , annot)

var1 <- rownames (table (dta[, colnames(dta)%in%gory]))[1]
var2 <- rownames (table (dta[, colnames(dta)%in%gory]))[2]

lis <- unique (metadata[, colnames (metadata)%in%ids])
lisME <- lis[!grepl("malignant", lis, ignore.case = TRUE) & !grepl("UNA", lis, ignore.case = TRUE) & !grepl("Normal_epithelial", lis, ignore.case = TRUE)]
selec <- get(paste ('lis', subs, sep=""))
metadata1<- metadata[metadata[,colnames (metadata)%in%ids]%in%selec,]

dta_bas<- dta[dta$category%in%cate & !is.na(dta[, colnames(dta)%in%gory]),]

pAB <- wilcox.test(dta_bas[dta_bas$category=="A" & dta_bas$perc_malig_surg10=="BR", colnames (dta)%in% probe], dta_bas[dta_bas$category=="A" & dta_bas$perc_malig_surg10=="GR", colnames (dta_bas)%in% probe], paired = F)[3]
pvalues <- c(pAB)

outline_colors <- c("#A02629", '#7BA58C') 
box_fill_colors <- c('#E5A1A4', '#C3EAD6') 
point_fill_colors <- c('#C32E32', '#9ED3B9') 

pdf ('ExDaFigure6A Naive B cells GR BR.pdf')

a <- ggplot(dta_bas, aes(x = !!sym(gory), y = !!sym(probe), label = paired)) + 
  geom_boxplot(aes(fill = factor(!!sym(gory))), color = outline_colors, fill = box_fill_colors, alpha = 0.2, size = 1.5) + 
  geom_point(aes(colour = factor(!!sym(gory)), fill = factor(!!sym(gory))), shape = 21, size = 5, stroke = 1.2) + 
  scale_fill_manual(values= point_fill_colors) + 
  scale_colour_manual(values=outline_colors) + 

  ylab('Percentage of naive B cells in the ME') + 
  geom_text_repel(aes(label = paired),  max.overlaps = Inf, box.padding = 0.5, point.padding = 0.3, segment.size = 0.2, segment.color = 'grey50') + 
  scale_x_discrete(labels = c(var1 = var1, var2 = var2)) +
  theme(axis.line = element_line(colour = "black", size = 1, linetype = "solid")) +
  theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14)) +
  theme(axis.title = element_text(size = 18, face = "bold")) +
  theme(axis.ticks = element_line(size = 1)) +
  theme(panel.background = element_rect(fill = 'white', color = 'black')) +
  labs(x = '') +
  theme(legend.position = "none")+  
  geom_signif(comparisons = list(c('BR', 'GR')),
			  y_position = 8.1,  
               annotations = sprintf("p = %.3f", pvalues),
               map_signif_level = FALSE,
               test = "wilcox.test", 
               step_increase = 0.07) 
			   
plot(a)
dev.off()

###################################################################################
#### Extended data Figure 6 B Boxplot immune cluster vs monocytes
###################################################################################
library(Seurat)
library(plyr)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(stringi)
library(mixtools)
library(ggsignif)
library (ggpattern)

load("metadata_git.RData")
#### load the metadata available at EGA
met_clin <- data.frame(read.table("clinical_data.txt", header=T, sep="\t", stringsAsFactors=F))
##### load clinical data in Supp Table 1
scores <- data.frame(read.table("all_maligmant_cells_metadata_with_modules.txt", header=T, sep="\t", stringsAsFactors=F))
##### module scores from  PP_gene_cluster_correlation


mer <- merge(metadata, met_clin, by.x='patient_number', by.y="patient_number", all.x=T)
metadata <- mer 

metadata1<- metadata[!metadata$main.ids=="Epithelial_cells" & !metadata$final.ids=="UNA",]
dta <- as.data.frame.matrix(prop.table(table (metadata1$orig.ident, metadata1$final.ids),1)*100)

dta$category <- stri_sub(rownames (dta),-1)
dta$paired <- substr (rownames (dta),2,4)
dta$paired <- gsub ("A", "",dta$paired )
dta$paired <- gsub ("B", "",dta$paired )
dta$paired <- gsub ("C", "",dta$paired )

dta$category <- stri_sub(rownames (dta),-1)
annot <- metadata[, colnames (metadata) %in% c("orig.ident", "perc_malig_surg10", "clone_cat", "BAS_subtype", "immune_cluster")]
annot <- annot[!duplicated (metadata$orig.ident), ]
table (annot$orig.ident == rownames (dta))
dta <- cbind (dta , annot)

a <- as.data.frame.matrix(prop.table (table (metadata$orig.ident, metadata$clone_cat),1)*100)
a$category <- stri_sub(rownames (a),-1)
a1 <- a[a$category == "A",]

dta <- dta[dta$category == "A",]
dta1 <- dta

dta1$immune_cluster <- factor(dta1$immune_cluster, levels = c('LumA_R', 'LumB_R', 'LumB_NR'))


pval_LumAGR_v_LumBGR <- wilcox.test(dta1$Mono_CD14[dta1$immune_cluster == 'LumA_R'],
                       dta1$Mono_CD14[dta1$immune_cluster == 'LumB_R'], 
                       paired = FALSE)$p.value
pval_LumBGR_v_LumBBR <- wilcox.test(dta1$Mono_CD14[dta1$immune_cluster == 'LumB_R'],
                       dta1$Mono_CD14[dta1$immune_cluster == 'LumB_NR'], 
                       paired = FALSE)$p.value
pval_LumAGR_v_LumBBR <- wilcox.test(dta1$Mono_CD14[dta1$immune_cluster == 'LumA_R'],
                       dta1$Mono_CD14[dta1$immune_cluster == 'LumB_NR'], 
                       paired = FALSE)$p.value
pvalues <- c(pval_LumAGR_v_LumBGR, pval_LumBGR_v_LumBBR, pval_LumAGR_v_LumBBR)

pdf ('ExDaFigure6B boxplot_monocytes_TME_clusters.pdf')

ggplot(dta1, aes(x = immune_cluster, y = Mono_CD14)) + 
  geom_boxplot_pattern(
    aes(fill = immune_cluster, pattern = immune_cluster), 
    pattern = c("none", "stripe", "stripe"),  
    pattern_density = c(0, 0.1, 0.1),          
    pattern_color = c(NA, "darkgreen", "red"),
    pattern_fill = c(NA, "darkgreen", "red"), 	
    alpha = 0.8, outlier.colour = NA, notch = FALSE, size = 1
  ) + 
  geom_point(aes(colour = immune_cluster), shape = 21, size = 5) + 
  scale_fill_manual(values = c('#0000FF', '#9ED2E2', '#9ED2E2')) +  
  scale_colour_manual(values = c('#3333CC', '#5EB3CE', '#5EB3CE')) + 
  ylab('Mono_CD14') + 
  geom_text_repel(aes(label = paired), 
                  max.overlaps = Inf, 
                  box.padding = 0.5, 
                  point.padding = 0.3, 
                  segment.size = 0.2, 
                  segment.color = 'grey50') +  
  theme(
    axis.line = element_line(colour = "black", size = 1, linetype = "solid"),
    axis.text.x = element_text(size = 14, color='black'), 
    axis.text.y = element_text(size = 14, color='black'),
    axis.title = element_text(size = 18, face = "bold"),
    axis.ticks = element_line(size = 1),
    panel.background = element_rect(fill = 'white', color = 'black'),
    legend.position = "none",
    strip.background = element_blank(),  
    strip.text = element_text(size = 15, face = "bold")  
  ) +
  labs(x = 'Cluster_response_TME') +  
  geom_signif(comparisons = list(c("LumA_R", "LumB_R"),  
                                  c("LumB_R", "LumB_NR"), 
                                  c("LumA_R", "LumB_NR")),
              annotations = sprintf("p = %.3f", pvalues), 
              map_signif_level = FALSE,
              test = "wilcox.test", 
              step_increase = 0.07)

dev.off()

############################################################################################
#### Extended data Figure 6 C Correlation between interferon and Myeloid cells #############
############################################################################################

library(Seurat)
library(plyr)
library(dplyr)
library(future)
library(ggplot2)
library(ggrepel)
library(stringi)
library(mixtools)
library(pheatmap)
library(amap)
library(gwasrapidd)
library(RColorBrewer)
library (corrplot)

load("metadata_git.RData")
##### needs the Seurat_metadata
met_clin <- data.frame(read.table("clinical_data.txt", header=T, sep="\t", stringsAsFactors=F))
##### needs clinical data in supp Table 1
scores <- data.frame(read.table("all_maligmant_cells_metadata_with_modules.txt", header=T, sep="\t", stringsAsFactors=F))
##### module scores from  PP_gene_cluster_correlation

mer <- merge(metadata, met_clin, by.x='patient_number', by.y="patient_number", all.x=T)
metadata <- mer 

metadata1<- metadata[!metadata$main.ids=="Epithelial_cells" & !metadata$final.ids=="UNA",]
dta <- as.data.frame.matrix(prop.table(table (metadata1$orig.ident, metadata1$final.ids),1)*100)
dta1 <- dta

dta$category <- stri_sub(rownames (dta),-1)
dta$paired <- substr (rownames (dta),2,4)
dta$paired <- gsub ("A", "",dta$paired )
dta$paired <- gsub ("B", "",dta$paired )
dta$paired <- gsub ("C", "",dta$paired )

dta$category <- stri_sub(rownames (dta),-1)
annot <- metadata[, colnames (metadata) %in% c("orig.ident")]
annot <- annot[!duplicated (metadata$orig.ident), ]
table (annot$orig.ident == rownames (dta))
dta <- cbind (dta , annot)

dtaA <- dta1[dta$category == "A",]
average_scores <- data.frame(scores %>%
  select(orig.ident, starts_with("cluster")) %>%  
  group_by(orig.ident) %>%                         
  summarize(across(starts_with("cluster"), mean, na.rm = TRUE)))  

score_bas <- average_scores[average_scores$orig.ident%in% rownames (dtaA),]
Inter_A <- score_bas$cluster10

A<- cor (Inter_A, dtaA, method='spearman')


pA <- sapply(dtaA, function(col) {    
  test <- cor.test(Inter_A, col, method = "spearman", exact = T)           
  test$p.value
})

dtaB <- dta1[dta$category == "B",]
average_scores <- data.frame(scores %>%
  select(orig.ident, starts_with("cluster")) %>%  
  group_by(orig.ident) %>%                         
  summarize(across(starts_with("cluster"), mean, na.rm = TRUE)))  

score_bas <- average_scores[average_scores$orig.ident%in% rownames (dtaB),]
Inter_B <- score_bas$cluster10

B<- cor (Inter_B, dtaB, method='spearman')

pB <- sapply(dtaB, function(col) {    
  test <- cor.test(Inter_B, col, method = "spearman", exact = T)           
  test$p.value
})

dtaC <- dta1[dta$category == "C",]
average_scores <- data.frame(scores %>%
  select(orig.ident, starts_with("cluster")) %>%  
  group_by(orig.ident) %>%                         
  summarize(across(starts_with("cluster"), mean, na.rm = TRUE)))  
score_bas <- average_scores[average_scores$orig.ident%in% rownames (dtaC),]
Inter_C <- score_bas$cluster10
C<- cor (Inter_C, dtaC, method='spearman')

pC <- sapply(dtaC, function(col) {    
  test <- cor.test(Inter_C, col, method = "spearman", exact = T)          
  test$p.value
})

cor_M <- data.frame (rbind (A, B, C))
rownames (cor_M) <- c('A', 'B', 'C')
cor_M <- t(cor_M)
cor_M <- cor_M[ rownames (cor_M) %in% c("DC", "DC_CCR7", "DC_CD1c", "DC_CLEC9A", "Granulocyte", "Mast_cells", "MF_angio", "MF_BNC2", "MF_IFN", "MF_LA_M2", "MF_LYVE1", "MF_REG_M1", "Mono_CD14", "Mono_CD16", "Neutrophils"),]

p_M <- data.frame (rbind (pA, pB, pC))
rownames (p_M) <- c('A', 'B', 'C')
p_M <- t(p_M)
p_M <- p_M[ rownames (p_M) %in% c("DC", "DC_CCR7", "DC_CD1c", "DC_CLEC9A", "Granulocyte", "Mast_cells", "MF_angio", "MF_BNC2", "MF_IFN", "MF_LA_M2", "MF_LYVE1", "MF_REG_M1", "Mono_CD14", "Mono_CD16", "Neutrophils"),]

breaksList = seq(-1, 1, by = 0.05)
col.pal <- colorRampPalette(c("dodgerblue4","white","brown3"))(length(breaksList))

dcols <- hcluster(t(cor_M), method = 'euclidean', link = 'ward')
drows <- hcluster(cor_M, method = 'euclidean', link = 'ward')
					
a<- pheatmap(cor_M, color = col.pal, show_rownames = T, show_colnames = T, main = "",  cluster_rows = drows, cluster_cols = dcols, 
breaks = breaksList, fontsize = 10, border_color = NA)
print (a) 

cor_M <- cor_M[drows$order, ]
p_M <- p_M[drows$order, ]
corrplot(cor_M, is.corr=F, col=col.pal, col.lim=c(-1, 1), tl.col = 'black')
corrplot(cor_M, is.corr=F, col=col.pal, col.lim=c(-1, 1), tl.col = 'black', p.mat =p_M , sig.level = 0.2, insig='blank')


pdf ('ExDaFigure6C corrplot interferon vs myeloid cells.pdf')
print (a)
corrplot(cor_M, is.corr=F, col=col.pal, col.lim=c(-1, 1), tl.col = 'black')
corrplot(cor_M, is.corr=F, col=col.pal, col.lim=c(-1, 1), tl.col = 'black', p.mat =p_M , sig.level = 0.2, insig='blank')
dev.off()




######################## ##############################################
#### Extended data Figure 6D CAF_ECM_BR_GR time ##############
######################################################################

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
met_clin <- data.frame(read.table("clinical_data.txt", header=T, sep="\t", stringsAsFactors=F))
##### load clinical data in Supp Table 1
mer <- merge(metadata, met_clin, by.x='patient_number', by.y="patient_number", all.x=T)
metadata <- mer 

ids <- "final.ids"
subs <- "ME"
gory <- "perc_malig_surg10"
cate <- "A"
probe <- "CAF_ECM"

metadata1<- metadata[!metadata$main.ids=="Epithelial_cells" & !metadata$final.ids=="UNA",]

lis <- unique (metadata[, colnames (metadata)%in%ids])
lisME <- lis[!grepl("malignant", lis, ignore.case = TRUE) & !grepl("UNA", lis, ignore.case = TRUE) & !grepl("Normal_epithelial", lis, ignore.case = TRUE)]
selec <- get(paste ('lis', subs, sep=""))
metadata1<- metadata[metadata[,colnames (metadata)%in%ids]%in%selec,]
dta <- as.data.frame.matrix(prop.table(table (metadata1$orig.ident, metadata1[, colnames(metadata1)%in%ids]),1)*100)


dta$category <- stri_sub(rownames (dta),-1)
dta$paired <- substr (rownames (dta),2,4)
dta$paired <- gsub ("A", "",dta$paired )
dta$paired <- gsub ("B", "",dta$paired )
dta$paired <- gsub ("C", "",dta$paired )

annot <- metadata[, colnames (metadata)%in% c("orig.ident", "perc_malig_surg10" )]
annot <- annot[!duplicated (metadata$orig.ident), ]
table (annot$orig.ident== rownames (dta))
dta <- cbind (dta , annot)

dta <- dta[!dta$orig.ident=='P04A',]

var1 <- rownames (table (dta[, colnames(dta)%in%gory]))[1]
var2 <- rownames (table (dta[, colnames(dta)%in%gory]))[2]
var3 <- rownames (table (dta[, colnames(dta)%in%gory]))[3]
var4 <- rownames (table (dta[, colnames(dta)%in%gory]))[4]

dta_var1 <- dta[!is.na(dta[, colnames(dta) %in% gory]) & dta[, colnames(dta) %in% gory] %in% var1, ]
dta_var2 <- dta[!is.na(dta[, colnames(dta) %in% gory]) & dta[, colnames(dta) %in% gory] %in% var2, ]

outline_colors <- c("#AA3830", '#9D2E99', '#350D7F') 
box_fill_colors <- c('#FFE2E0', '#FFD1FF', '#DAC3FF') 
point_fill_colors <- c('#FF8077', '#E66CE1', '#7E33FF') 


a<- ggplot(dta_var1, aes(category,dta_var1[, colnames(dta_var1)%in%probe], label=paired)) + 
geom_boxplot(color = outline_colors, fill = box_fill_colors, alpha = 0.2, size = 1.5)+ 
geom_point(aes(colour=factor(category), fill = factor(category)), shape = 21, size = 5, stroke = 1.2) +  
  scale_fill_manual(values= point_fill_colors) + 
  scale_colour_manual(values=outline_colors) + 
ylab(probe)+ 
geom_text_repel(min.segment.length = 1, max.overlaps = Inf)  + 
scale_x_discrete(labels=c("A" = "A", "B" = "B", "C" = "C")) +
theme( axis.line = element_line(colour = "black", size = 1, linetype = "solid")) +
theme(axis.text.x = element_text(size=14), axis.text.y = element_text(size=14)) +
theme(axis.title=element_text(size=18,face="bold")) +
theme(axis.ticks = element_line(size = 1))+
theme(panel.background = element_rect(fill = 'white', color = 'black'))+
labs( x = "Time point")+
theme(legend.position = "none")+
ggtitle (paste (gory, cate, sep="_"))
plot (a)

b<- ggplot(dta_var2, aes(category,dta_var2[, colnames(dta_var2)%in%probe], label=paired)) + 
geom_boxplot(color = outline_colors, fill = box_fill_colors, alpha = 0.2, size = 1.5)+ 
geom_point(aes(colour=factor(category), fill = factor(category)), shape = 21, size = 5, stroke = 1.2) +  # Fixed parentheses here
  scale_fill_manual(values= point_fill_colors) + 
  scale_colour_manual(values=outline_colors) + 
ylab(probe)+ 
geom_text_repel(min.segment.length = 1, max.overlaps = Inf)  + 
scale_x_discrete(labels=c("A" = "A", "B" = "B", "C" = "C")) +
theme( axis.line = element_line(colour = "black", size = 1, linetype = "solid")) +
theme(axis.text.x = element_text(size=14), axis.text.y = element_text(size=14)) +
theme(axis.title=element_text(size=18,face="bold")) +
theme(axis.ticks = element_line(size = 1))+
theme(panel.background = element_rect(fill = 'white', color = 'black'))+
labs( x = "Time point")+
theme(legend.position = "none")+
ggtitle (paste (gory, cate, sep="_"))
plot (b)


dta_var1$group <- var1
dta_var2$group <- var2
combined_data <- rbind(dta_var1, dta_var2)


outline_colors <- c("#A02629", '#7BA58C') 
box_fill_colors <- c('#E5A1A4', '#C3EAD6') 
point_fill_colors <- c('#C32E32', '#9ED3B9') 


combined_plot <- ggplot(combined_data, aes(x = category, y = !!sym(probe), fill = group, label = paired)) +
  geom_boxplot(aes(color = group), position = position_dodge(width = 0.75), alpha = 0.2, outlier.shape = NA) +
  geom_point(aes(colour = group), position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), shape = 21, size = 3) +
  scale_fill_manual(values= point_fill_colors) + 
  scale_colour_manual(values=outline_colors) + 
  ylab(paste("Percentage of", probe, "in the ME")) +
  scale_x_discrete(labels = c("A" = "A", "B" = "B", "C" = "C")) +
  geom_text_repel(aes(label = paired), position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), min.segment.length = 1) + 
  theme(axis.line = element_line(colour = "black", size = 1, linetype = "solid"),
        axis.text.x = element_text(size = 14),
        axis.title = element_text(size = 18, face = "bold"),
        axis.ticks = element_line(size = 1),
        panel.background = element_rect(fill = 'white', color = 'black'),
        legend.position = "none") +
  labs(x = "Time point") +
  scale_y_continuous(limits = c(0, 20))  

pdf ('ExDaFigure6D CAF_ECM BR GR time_pdf')
print(combined_plot)
dev.off()

########################################
###### EXTENDED DATA Figure 7 ##########
########################################  
########################################################################
############# ExDA FIG 7A - barplot of the number of clones per cell state ##
########################################################################
library(Seurat)
library(plyr)
library(dplyr)
library(future)
library (tidyr)
library (scRepertoire)
library (stringi)
library (stringr)
library (ggplot2)
library (ggrepel)

TCR<- readRDS("combined_metadata_repertoire_tcr_mix1.rds")
#### load the scRepertoire output from TCR snakemake pipeline
load("metadata_git.RData")
#### load the metadata available at EGA


TCR$TCR <- ifelse (is.na(TCR$clonalFrequency), 0,1)

count_data <- TCR %>%
  filter(TCR == 1) %>%
  group_by(detail.ids, CTnt) %>%
  summarise(clonalFrequency = n(), .groups = 'drop') %>%
  ungroup()

count_data <- count_data %>%
  mutate(frequency_category = case_when(
    clonalFrequency == 1 ~ "1",
    clonalFrequency >= 2 & clonalFrequency < 5 ~ "2-5",
    clonalFrequency >= 5 ~ ">5"
  ))



plot_data <- count_data %>% group_by(detail.ids, frequency_category) %>% summarise(clone_count = n(), .groups = 'drop')
plot_data <- plot_data %>% mutate(frequency_category = factor(frequency_category, levels = c(">5", "2-5", "1"))) 


plot_data <- plot_data %>%
  group_by(detail.ids) %>%
  mutate(total_clones = sum(clone_count)) %>%
  ungroup() %>%
  arrange(desc(total_clones)) %>%
  mutate(detail.ids = factor(detail.ids, levels = unique(detail.ids)))

pdf('ExDaFigure7A barplot of _clone_frequency.pdf')
ggplot(plot_data, aes(x = detail.ids, y = clone_count, fill = frequency_category)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("1" = "blue", "2-5" = "purple", ">5" = "red")) +
  labs(title = "Counts of TCR Clones by Cell State", 
       x = "Cell States (detail.ids)", 
       y = "Count of TCR Clones") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_blank(),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),  
        axis.line = element_line(colour = "black"),  
        axis.ticks = element_line(colour = "black"))		
dev.off()


########################################################################
############# ExDA FIG 7B - barplot of the number of cells in clonotypes ##
########################################################################

library(Seurat)
library(plyr)
library(dplyr)
library(future)
library (tidyr)
library (scRepertoire)
library (stringi)
library (stringr)
library (ggplot2)
library (ggrepel)

TCR<- readRDS("combined_metadata_repertoire_tcr_mix1.rds")
#### load the scRepertoire output from TCR snakemake pipeline
load("metadata_git.RData")
#### load the metadata available at EGA

TCR$TCR <- ifelse (is.na(TCR$clonalFrequency), 0,1)

cell_count_data <- TCR %>%
  filter(TCR == 1) %>%
  group_by(detail.ids, CTnt) %>%
  summarise(clonalFrequency = n(), .groups = 'drop') %>%
  ungroup()

cell_count_data <- cell_count_data %>%
  mutate(frequency_category = case_when(
    clonalFrequency == 1 ~ "1",
    clonalFrequency >= 2 & clonalFrequency < 5 ~ "2-5",
    clonalFrequency >= 5 ~ ">5"
  ))

cell_count_data <- cell_count_data %>% filter(detail.ids != "UNA")
cell_plot_data <- cell_count_data %>% 
  group_by(detail.ids, frequency_category) %>% 
  summarise(cell_count = sum(clonalFrequency), .groups = 'drop')

cell_plot_data <- cell_plot_data %>% 
  mutate(frequency_category = factor(frequency_category, levels = c(">5", "2-5", "1"))) 

cell_plot_data <- cell_plot_data %>%
  group_by(detail.ids) %>%
  mutate(total_cells = sum(cell_count)) %>%
  ungroup() %>%
  arrange(desc(total_cells)) %>%
  mutate(detail.ids = factor(detail.ids, levels = unique(detail.ids)))

pdf('ExDaFigure7B barplot of cell frequency.pdf')

ggplot(cell_plot_data, aes(x = detail.ids, y = cell_count, fill = frequency_category)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("1" = "blue", "2-5" = "purple", ">5" = "red")) +
  labs(title = "Counts of TCR Clones by Cell State", 
       x = "Cell States (detail.ids)", 
       y = "Count of TCR Clones") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_blank(),
        panel.grid.major = element_blank(),  #
        panel.grid.minor = element_blank(),  
        axis.line = element_line(colour = "black"),  
        axis.ticks = element_line(colour = "black")) 
dev.off()


#########################################################
## Extended data figure 7C Circos plot shared clone across states ##
#########################################################

library(Seurat)
library(dplyr)
library(scRepertoire)
library(ggraph)
library(circlize)
library(scales)
library(APackOfTheClones)
library(RColorBrewer)
library(SCpubr)
library(Nebulosa)

obj <- readRDS("combined_seurat_repertoire_tcr.rds")
#### load the scRepertoire output from TCR snakemake pipeline

obj_f <- obj[, !obj@meta.data$orig.ident %in% c('P24A', 'P06A', 'P11B', 'P01', 'P02C', 'P01B', 'P02B')]

palopal <- c(
  "#977899",  # NK_cyto
  "#E0DEDE",  # UNA
  "#5FB0E8",  # CD8_EM
  "#00C957",  # CD4_EM
  "#F57C7C",  # CD4_low
  "#AFDEA2",  # CD4_Naive
  "#3385BB",  # CD8_EMRA
  "#E2CAE3",  # NK_rest
  "#DE9E83",  # Treg
  "#7F9D55",  # CD4_Ex
  "#99CAF7",  # CD8_Naive
  "#FE8D19",  # GDT
  "#3F5361",  # CD8_Ex
  "#D40F0F",  # MAIT
  "#8A5A46",  # ILC3
  "#FBB268"   # CRTAM_GDT
)

circles <- getCirclize(obj_f, group.by = "final.ids")

unique_clusters <- unique(circles$from)
grid.cols <- palopal[match(unique_clusters, c("NK_cyto", "UNA", "CD8_EM", "CD4_EM", "CD4_low",
                                              "CD4_Naive", "CD8_EMRA", "NK_rest", "Treg", 
                                              "CD4_Ex", "CD8_Naive", "GDT", "CD8_Ex", 
                                              "MAIT", "ILC3", "CRTAM_GDT"))]

names(grid.cols) <- unique_clusters

circles <- getCirclize(obj_f, group.by = "final.ids")

pdf ('ExDaFigure7C circos_plot_clone state sharing.pdf')
chordDiagram(circles, self.link = 1, grid.col = grid.cols)
dev.off()

#########################################################################
## Extended data figure 7D Shannon index in between response and time ##
#########################################################################

library(ggplot2)
library(ggrepel)
library(ggpubr)
library(stringi)
library(rlang)
library(ggsignif)
library(ggpubr)

TCR<- readRDS("combined_metadata_repertoire_tcr_mix1.rds")
TCR <- TCR [!TCR$patient_number%in%c(4, 24, 6, 11, 2, 1),]

tab <- table(TCR$orig.ident, TCR$CTnt)
clone_sizes <- as.data.frame(tab, stringsAsFactors = FALSE)
colnames(clone_sizes) <- c("sample_id", "clone_id", "clone_size")
clone_sizes <- clone_sizes[clone_sizes$clone_size > 0, ]

sample_ids <- unique(clone_sizes$sample_id)

metrics <- data.frame(
  sample_id = sample_ids,
  shannon   = NA_real_,
  stringsAsFactors = FALSE
)

for (i in seq_along(sample_ids)) {
  sid <- sample_ids[i]
  cs  <- clone_sizes$clone_size[clone_sizes$sample_id == sid]
  p   <- cs / sum(cs)

  H <- -sum(p * log(p))        
  metrics$shannon[i] <- H
}

rownames(metrics) <- metrics$sample_id
metrics$paired   <- substr(rownames(metrics), 2, 4)
metrics$paired   <- gsub("A|B|C", "", metrics$paired)
metrics$category <- stri_sub(rownames(metrics), -1)
metrics$response <- ifelse(metrics$paired %in% c("14", "11", "16", "17", "13", "23", "10", "21"), "BR", metrics$paired)
metrics$response <- ifelse(metrics$response %in% c("01", "07", "15", "18", "19","05", "08", "24", "06", "09","20", "12", "02", "22", "25", "03"), "GR", metrics$response)

probe <- c('shannon')
metrics_f <- metrics

pdf ('ExDaFigure7D Sahnnon index TCR response time.pdf')

df_BR <- metrics_f [metrics_f$response=="BR",  ]
df_BR <- df_BR [order(df_BR$paired),  ]
df_GR <-  metrics_f [metrics_f$response=="GR" ,]
df_GR <- df_GR [order(df_GR$paired),  ]

outline_colors <- c("#AA3830", '#9D2E99', '#350D7F') 
box_fill_colors <- c('#FFE2E0', '#FFD1FF', '#DAC3FF') 
point_fill_colors <- c('#FF8077', '#E66CE1', '#7E33FF') 

BR_A <- df_BR[df_BR$category=='A', colnames(df_BR)%in%probe]
BR_B <- df_BR[df_BR$category=='B', colnames(df_BR)%in%probe]
BR_C <- df_BR[df_BR$category=='C', colnames(df_BR)%in%probe]

GR_A <- df_GR[df_GR$category=='A', colnames(df_GR)%in%probe]
GR_B <- df_GR[df_GR$category=='B', colnames(df_GR)%in%probe]
GR_C <- df_GR[df_GR$category=='C', colnames(df_GR)%in%probe]
  
test_A <- wilcox.test(BR_A, GR_A, paired = FALSE)$p.value
test_B <- wilcox.test(BR_B, GR_B, paired = FALSE)$p.value
test_C <- wilcox.test(BR_C, GR_C, paired = FALSE)$p.value

test_BR_AB <- wilcox.test(BR_A, BR_B, paired = TRUE)$p.value
test_BR_BC <- wilcox.test(BR_B, BR_C, paired = TRUE)$p.value
test_BR_AC <- wilcox.test(BR_A, BR_C, paired = TRUE)$p.value

test_GR_AB <- wilcox.test(GR_A, GR_B, paired = TRUE)$p.value
test_GR_BC <- wilcox.test(GR_B, GR_C, paired = TRUE)$p.value
test_GR_AC <- wilcox.test(GR_A, GR_C, paired = TRUE)$p.value

outline_colors <- c("#AA3830", '#9D2E99', '#350D7F') 
box_fill_colors <- c('#FFE2E0', '#FFD1FF', '#DAC3FF') 
point_fill_colors <- c('#FF8077', '#E66CE1', '#7E33FF') 

pvalues_a <- c(test_BR_AB, test_BR_BC, test_BR_AC)

a <- ggplot(df_BR, aes(category, df_BR[, colnames(df_BR) %in% probe], label = '')) +
  geom_boxplot(aes(color = factor(category), fill = factor(category)),
               alpha = 0.2, size = 1.5) +
  geom_point(aes(colour = factor(category), fill = factor(category)),
             shape = 21, size = 5, stroke = 1.2) +
  scale_fill_manual(values = point_fill_colors) +
  scale_colour_manual(values = outline_colors) +
  ylab(paste(probe, 'BR', sep="_")) +
  geom_text_repel(min.segment.length = 1, max.overlaps = Inf) +
  scale_x_discrete(labels = c("A" = "A", "B" = "B", "C" = "C")) +
  theme(axis.line = element_line(colour = "black", size = 1, linetype = "solid")) +
  theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14)) +
  theme(axis.title = element_text(size = 18, face = "bold")) +
  theme(axis.ticks = element_line(size = 1)) +
  theme(panel.background = element_rect(fill = 'white', color = 'black')) +
  labs(x = "Time point") +
  theme(legend.position = "none") +
  ggtitle("") +
  geom_signif(comparisons = list(c("A", "B"),
                                 c("B", "C"),
                                 c("A", "C")),
              annotations = sprintf("p = %.3f", pvalues_a),
              map_signif_level = FALSE,
              test = "wilcox.test",
              step_increase = 0.07)

plot(a)

            
pvalues_b <- c(test_GR_AB, test_GR_BC, test_GR_AC)
b<- ggplot(df_GR, aes(category, df_GR[, colnames(df_GR) %in% probe], label = '')) +
  geom_boxplot(aes(color = factor(category), fill = factor(category)),
               alpha = 0.2, size = 1.5) +
  geom_point(aes(colour = factor(category), fill = factor(category)),
             shape = 21, size = 5, stroke = 1.2) +
  scale_fill_manual(values = point_fill_colors) +
  scale_colour_manual(values = outline_colors) +
  ylab(paste(probe, 'GR', sep="_")) +
  geom_text_repel(min.segment.length = 1, max.overlaps = Inf) +
  scale_x_discrete(labels = c("A" = "A", "B" = "B", "C" = "C")) +
  theme(axis.line = element_line(colour = "black", size = 1, linetype = "solid")) +
  theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14)) +
  theme(axis.title = element_text(size = 18, face = "bold")) +
  theme(axis.ticks = element_line(size = 1)) +
  theme(panel.background = element_rect(fill = 'white', color = 'black')) +
  labs(x = "Time point") +
  theme(legend.position = "none") +
  ggtitle("") +
  geom_signif(comparisons = list(c("A", "B"),
                                 c("B", "C"),
                                 c("A", "C")),
              annotations = sprintf("p = %.3f", pvalues_b),
              map_signif_level = FALSE,
              test = "wilcox.test",
              step_increase = 0.07)

plot (b)

df_BR$group <- "BR"
df_GR$group <- "GR"
combined_data <- rbind(df_BR, df_GR)
combined_data$group_time <- paste0(combined_data$group, "_", combined_data$category)



outline_colors <- c("#A02629", '#7BA58C', '#4B4B4B') 
box_fill_colors <- c('#E5A1A4', '#C3EAD6', '#7E7E7E' ) 
point_fill_colors <- c('#C32E32', '#9ED3B9', '#B0B0B0') 

p_between <- ggplot(combined_data, aes(x = group, y = !!rlang::sym(probe), fill = group)) +
  geom_boxplot(aes(color = group), alpha = 0.2, outlier.shape = NA) +
  geom_jitter(width = 0.15, height = 0, aes(colour = group), shape = 21, size = 3) +
  scale_fill_manual(values = point_fill_colors) +
  scale_colour_manual(values = outline_colors) +
  facet_wrap(~ category) +
  ylab(paste( probe)) +
  theme_minimal() +
  theme(legend.position = "none")

p_between <- p_between +
  stat_compare_means(method = "wilcox.test", label = "p.format",
                     comparisons = list(c("BR", "GR")),
                     hide.ns = FALSE, size = 4, label.y = max(combined_data[[probe]], na.rm = TRUE) * 1.05)
print(p_between)




combined_plot <- ggplot(combined_data, aes(x = category, y = !!sym(probe), fill = group, label = paired)) +
  geom_boxplot(aes(color = group), position = position_dodge(width = 0.75), alpha = 0.2, outlier.shape = NA) +
  geom_point(aes(colour = group), position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), shape = 21, size = 3) +
  scale_fill_manual(values= point_fill_colors) + 
  scale_colour_manual(values=outline_colors) + 
  ylab(paste("Percentage of", probe, "in the ME")) +
  scale_x_discrete(labels = c("A" = "A", "B" = "B", "C" = "C")) +
  geom_text_repel(aes(label = paired), position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), min.segment.length = 1) + 
  theme(axis.line = element_line(colour = "black", size = 1, linetype = "solid"),
        axis.text.x = element_text(size = 14),
        axis.title = element_text(size = 18, face = "bold"),
        axis.ticks = element_line(size = 1),
        panel.background = element_rect(fill = 'white', color = 'black'),
        legend.position = "none") +
  labs(x = "Time point")

combined_plot <- combined_plot +
  stat_compare_means(method = "wilcox.test", label = "p.format",
                     comparisons = list(c("BR", "GR")),
                     hide.ns = FALSE, size = 4, label.y = max(combined_data[[probe]], na.rm = TRUE) * 1.05)
print(combined_plot)

dev.off()



########################################
###### EXTENDED DATA Figure 8 ##########
########################################
####################################################################################
############# EXTENDED DATA FIG 8A - BOXPLOT BCR CONTRACTION    ##################
####################################################################################
library (tidyr)
library (scRepertoire)
library (stringi)
library (stringr)
library (ggplot2)
library (ggrepel)
library(ggsignif)
library  (dplyr)
library(patchwork)
library(RColorBrewer)
library(ggalluvial)

load ('BCR_norm_cell_bind.RData')
#### load the summarized expansion_contraction data from PP_BCR_Expansion_all_cell_normalization.r
load("metadata_git.RData")
#### load the metadata available at EGA
merA <- data.frame(read.table("clinical_data.txt", header=T, sep='\t', row.names=))
merA <- merA [!merA$patient_number=="4",]

statusAB1_n_EC <- statusAB1_n_EC[order(rownames (statusAB1_n_EC)),]
statusAC1_n_EC <- statusAC1_n_EC[order(rownames (statusAC1_n_EC)),]

norm_EC_AB <- cbind (statusAB1_n_EC, merA$response, merA$patient_numb)
norm_EC_AC <- cbind (statusAC1_n_EC , merA$response, merA$patient_numb)
colnames (norm_EC_AB) <- c('lit_expansion', 'str_expansion','contraction',  'response', 'patient_number')
colnames (norm_EC_AC) <- c('lit_expansion', 'str_expansion', 'contraction',  'response', 'patient_number')
norm_EC_AB$category <- 'AB'
norm_EC_AC$category <- 'AC'


combined_data <- rbind(norm_EC_AB, norm_EC_AC)
combined_data <- combined_data %>% filter(!is.na(response) & !is.na(contraction))
pval_AB <- wilcox.test(norm_EC_AB$contraction[norm_EC_AB$response == 'BR'],
                       norm_EC_AB$contraction[norm_EC_AB$response == 'GR'], 
                       paired = FALSE)$p.value

pval_AC <- wilcox.test(norm_EC_AC$contraction[norm_EC_AC$response == 'BR'],
                       norm_EC_AC$contraction[norm_EC_AC$response == 'GR'], 
                       paired = FALSE)$p.value
					   
outline_colors <- c("#A02629", '#7BA58C') # Box outline colors
box_fill_colors <- c('#E5A1A4', '#C3EAD6') # Box fill colors
point_fill_colors <- c('#C32E32', '#9ED3B9') # Point fill colors

					   
pdf	('ExDaFigure8A boxplot_contraction_AB_AC_BCR.pdf')	   
ggplot(combined_data, aes(x = response, y = contraction)) + 
  geom_boxplot(aes(fill = response, colour = response), 
               alpha = 0.2, outlier.colour = NA, notch = FALSE, size = 1.5) + 
  geom_point(aes(fill = response, colour = response), shape = 21, size = 5, , stroke = 1.2) + 
  scale_fill_manual(values = setNames(c(box_fill_colors, point_fill_colors), 
                                      levels(combined_data$response))) +
  scale_colour_manual(values = setNames(c(outline_colors), 
                                         levels(combined_data$response))) +  
  ylab('Fraction of contracting clones') + 
  geom_text_repel(aes(label = patient_number), max.overlaps = Inf, 
                  box.padding = 0.5, point.padding = 0.3, 
                  segment.size = 0.2, segment.color = 'grey50') +  
  theme(
    axis.line = element_line(colour = "black", size = 1, linetype = "solid"),
    axis.text.x = element_text(size = 14, color = 'black'), 
    axis.text.y = element_text(size = 14, color = 'black'),
    axis.title = element_text(size = 18, face = "bold"),
    axis.ticks = element_line(size = 1),
    panel.background = element_rect(fill = 'white', color = 'black'),
    legend.position = "none",
    strip.background = element_blank(), 
    strip.text = element_text(size = 15, face = "bold") 
  ) +
  facet_wrap(~category) +  
  labs(x = 'Response') +  
  geom_signif(data = subset(combined_data, category == 'AB'), comparisons = list(c('BR', 'GR')),
               annotations = sprintf("p = %.3f", pval_AB),
               map_signif_level = FALSE,
               test = "wilcox.test", 
               step_increase = 0.1) +
  geom_signif(data = subset(combined_data, category == 'AC'), comparisons = list(c('BR', 'GR')),
               annotations = sprintf("p = %.3f", pval_AC),
               map_signif_level = FALSE,
               test = "wilcox.test", 
               step_increase = 0.1)  
 
dev.off()

#######################################################################################
############# EXTENDED DATA FIG 8B - BOXPLOT BCR LIT EXPANSION       ##################
#######################################################################################

combined_data <- combined_data %>% filter(!is.na(response) & !is.na(lit_expansion))
pval_AB <- wilcox.test(norm_EC_AB$lit_expansion[norm_EC_AB$response == 'BR'],
                       norm_EC_AB$lit_expansion[norm_EC_AB$response == 'GR'], 
                       paired = FALSE)$p.value
pval_AC <- wilcox.test(norm_EC_AC$lit_expansion[norm_EC_AC$response == 'BR'],
                       norm_EC_AC$lit_expansion[norm_EC_AC$response == 'GR'], 
                       paired = FALSE)$p.value

outline_colors <- c("#A02629", '#7BA58C') # Box outline colors
box_fill_colors <- c('#E5A1A4', '#C3EAD6') # Box fill colors
point_fill_colors <- c('#C32E32', '#9ED3B9') # Point fill colors


pdf	('ExDaFigure8B boxplot_lit_expansion_AB_AC_BCR.pdf')	 
ggplot(combined_data, aes(x = response, y = lit_expansion)) + 
  geom_boxplot(aes(fill = response, colour = response), 
               alpha = 0.2, outlier.colour = NA, notch = FALSE, size = 1.5) + 
  geom_point(aes(fill = response, colour = response), shape = 21, size = 5, , stroke = 1.2) + 
  scale_fill_manual(values = setNames(c(box_fill_colors, point_fill_colors), 
                                      levels(combined_data$response))) +
  scale_colour_manual(values = setNames(c(outline_colors), 
                                         levels(combined_data$response))) + 
  ylab('Fraction of expanding clones') + 
  geom_text_repel(aes(label = patient_number), max.overlaps = Inf, box.padding = 0.5, point.padding = 0.3, segment.size = 0.2, segment.color = 'grey50') +  
  theme(
    axis.line = element_line(colour = "black", size = 1, linetype = "solid"),
    axis.text.x = element_text(size = 14, color='black'), 
    axis.text.y = element_text(size = 14, color='black'),
    axis.title = element_text(size = 18, face = "bold"),
    axis.ticks = element_line(size = 1),
    panel.background = element_rect(fill = 'white', color = 'black'),
    legend.position = "none",
    strip.background = element_blank(),  
    strip.text = element_text(size = 15, face = "bold")  
  ) +
  facet_wrap(~category) +  
  labs(x = 'Response') +  
  geom_signif(data = subset(combined_data, category == 'AB'), comparisons = list(c('BR', 'GR')),
               annotations = sprintf("p = %.3f", pval_AB),
               map_signif_level = FALSE,
               test = "wilcox.test", 
               step_increase = 0.1) +
  geom_signif(data = subset(combined_data, category == 'AC'), comparisons = list(c('BR', 'GR')),
               annotations = sprintf("p = %.3f", pval_AC),
               map_signif_level = FALSE,
               test = "wilcox.test", 
               step_increase = 0.1)
dev.off()

#######################################################################################
############# EXTENDED DATA FIG 8C - BOXPLOT BCR LIT EXPANSION       ##################
#######################################################################################

combined_data <- combined_data %>% filter(!is.na(response) & !is.na(str_expansion))
pval_AB <- wilcox.test(norm_EC_AB$str_expansion[norm_EC_AB$response == 'BR'],
                       norm_EC_AB$str_expansion[norm_EC_AB$response == 'GR'], 
                       paired = FALSE)$p.value
pval_AC <- wilcox.test(norm_EC_AC$str_expansion[norm_EC_AC$response == 'BR'],
                       norm_EC_AC$str_expansion[norm_EC_AC$response == 'GR'], 
                       paired = FALSE)$p.value

outline_colors <- c("#A02629", '#7BA58C') # Box outline colors
box_fill_colors <- c('#E5A1A4', '#C3EAD6') # Box fill colors
point_fill_colors <- c('#C32E32', '#9ED3B9') # Point fill colors


pdf	('ExDaFigure8C boxplot_str_expansion_AB_AC_BCR.pdf')	 

ggplot(combined_data, aes(x = response, y = str_expansion)) + 
  geom_boxplot(aes(fill = response, colour = response), 
               alpha = 0.2, outlier.colour = NA, notch = FALSE, size = 1.5) + 
  geom_point(aes(fill = response, colour = response), shape = 21, size = 5, , stroke = 1.2) + 
  scale_fill_manual(values = setNames(c(box_fill_colors, point_fill_colors), 
                                      levels(combined_data$response))) +
  scale_colour_manual(values = setNames(c(outline_colors), 
                                         levels(combined_data$response))) + 
  ylab('Fraction of strongly expanding clones') + 
  geom_text_repel(aes(label = patient_number), max.overlaps = Inf, box.padding = 0.5, point.padding = 0.3, segment.size = 0.2, segment.color = 'grey50') +  
  theme(
    axis.line = element_line(colour = "black", size = 1, linetype = "solid"),
    axis.text.x = element_text(size = 14, color='black'), 
    axis.text.y = element_text(size = 14, color='black'),
    axis.title = element_text(size = 18, face = "bold"),
    axis.ticks = element_line(size = 1),
    panel.background = element_rect(fill = 'white', color = 'black'),
    legend.position = "none",
    strip.background = element_blank(),  
    strip.text = element_text(size = 15, face = "bold")  
  ) +
  facet_wrap(~category) +  
  geom_signif(data = subset(combined_data, category == 'AB'), comparisons = list(c('BR', 'GR')),
               annotations = sprintf("p = %.3f", pval_AB),
               map_signif_level = FALSE,
               test = "wilcox.test", 
               step_increase = 0.1) +
  geom_signif(data = subset(combined_data, category == 'AC'), comparisons = list(c('BR', 'GR')),
               annotations = sprintf("p = %.3f", pval_AC),
               map_signif_level = FALSE,
               test = "wilcox.test", 
			   step_increase = 0.1)
dev.off()




########################################################################
############# EXTENDED DATA 8D - BOXPLOT OF PATIENT behaviour thee time points       ##
########################################################################

library (tidyr)
library (scRepertoire)
library (stringi)
library (stringr)
library (ggplot2)
library (ggrepel)
library(ggsignif)
library  (dplyr)
library(patchwork)
library(RColorBrewer)
library(ggalluvial)
library(pheatmap)
library(amap)

load ('TCR_norm_cell_bind.RData')
#### load the summarized expansion_contraction data from PP_Expansion
load("metadata_git.RData")
#### load the metadata available at EGA
selected_cells <- statusABC1_n_EC[,grepl("ABC", colnames(statusABC1_n_EC)) & !colnames(statusABC1_n_EC) %in% c('ABC_late_contraction', 'ABC_lit_expansion')]
selected_cells <- selected_cells[complete.cases(selected_cells),]
 
colnames (selected_cells) <- gsub ("WOs_ABC_contraction_contraction", "contraction",   colnames (selected_cells))
colnames (selected_cells) <- gsub ("WOs_ABC_contraction_singleton", "contraction",   colnames (selected_cells))
colnames (selected_cells) <- gsub ("WOs_ABC_contraction_lit_expansion" , "contraction_expansion",   colnames (selected_cells))
colnames (selected_cells) <- gsub ("WOs_ABC_contraction_str_expansion" , "contraction_expansion",   colnames (selected_cells))
colnames (selected_cells) <- gsub ("WOs_ABC_lit_expansion_contraction" , "expansion_contraction",   colnames (selected_cells))
colnames (selected_cells) <- gsub ("WOs_ABC_str_expansion_contraction"  , "expansion_contraction",   colnames (selected_cells))
colnames (selected_cells) <- gsub ("WOs_ABC_lit_expansion_lit_expansion" , "expansion",   colnames (selected_cells))
colnames (selected_cells) <- gsub ("WOs_ABC_str_expansion_str_expansion" , "expansion",   colnames (selected_cells))
colnames (selected_cells) <- gsub ("WOs_ABC_lit_expansion_str_expansion" , "expansion",   colnames (selected_cells))
colnames (selected_cells) <- gsub ("WOs_ABC_str_expansion_lit_expansion" , "expansion",   colnames (selected_cells))
colnames (selected_cells) <- gsub ("WOs_ABC_new_late_expansion" , "new_late_expansion",   colnames (selected_cells))
colnames (selected_cells) <- gsub ("WOs_ABC_singleton_lit_expansion" , "singleton_expansion",   colnames (selected_cells))
colnames (selected_cells) <- gsub ("WOs_ABC_singleton_str_expansion", "singleton_expansion",   colnames (selected_cells))
 
agg <- res <- as.data.frame(lapply(split.default(selected_cells, names(selected_cells)),
function(group) rowSums(group, na.rm = TRUE)))

merA <- data.frame(read.table("clinical_data.txt", header=T, sep='\t', row.names=))
merA <- merA [!merA$patient_number=="4",]
merA <- merA[merA$patients %in% rownames(agg), ]
merA <- merA[match(rownames(agg), merA$patients), ]
 
all (rownames (agg) ==  merA$patients)
response <- merA$response
patient_number <- merA$patient_numb

norm_EC <- agg 
norm_EC$response <- response
norm_EC$patient_number <- patient_number
colnames (norm_EC) <- gsub ('WOs_ABC_' , '', colnames (norm_EC))


outline_colors <- c("#A02629", '#7BA58C') # Box outline colors
box_fill_colors <- c('#E5A1A4', '#C3EAD6') # Box fill colors
point_fill_colors <- c('#C32E32', '#9ED3B9') # Point fill colors
			  

pdf('ExDaFigure8D dynamic_TCR_contraction_expansion.pdf')
var1 <- 'BR'
var2 <- 'GR'
pval <- wilcox.test(norm_EC$contraction_expansion[norm_EC$response == var1],
                    norm_EC$contraction_expansion[norm_EC$response == var2], 
                    paired = FALSE)$p.value  
pvalues <- c(pval)
ggplot(norm_EC, aes(x = response, y = contraction_expansion)) + 
  geom_boxplot(aes(fill = response), color = outline_colors, fill = box_fill_colors, alpha = 0.2, size = 1.5) + 
  geom_point(aes(colour = response, fill = response), shape = 21, size = 5) +  
  scale_fill_manual(values= point_fill_colors) + 
  scale_colour_manual(values=outline_colors) + 
  ylab('Fraction of contractingAB - expansionBC clones') + 
  geom_text_repel(aes(label = patient_number), max.overlaps = Inf, box.padding = 0.5, point.padding = 0.3, segment.size = 0.2, segment.color = 'grey50') +  
  scale_x_discrete(labels = setNames(c(var1, var2), c(var1, var2))) +  
  theme(
    axis.line = element_line(colour = "black", size = 1, linetype = "solid"),
    axis.text.x = element_text(size = 14, color='black'), 
    axis.text.y = element_text(size = 14, color='black'),
    axis.title = element_text(size = 18, face = "bold"),
    axis.ticks = element_line(size = 1),
    panel.background = element_rect(fill = 'white', color = 'black'),
    legend.position = "none"
  ) +
  labs(x = 'Response') +  
  geom_signif(comparisons = list(c(var1, var2)),
               annotations = sprintf("p = %.3f", pvalues),
               map_signif_level = FALSE,
               test = "wilcox.test", 
               step_increase = 0.1)  		   
dev.off()			

pdf('ExDaFigure8E dynamic_TCR_expansion_contraction.pdf')
var1 <- 'BR'
var2 <- 'GR'
pval <- wilcox.test(norm_EC$expansion_contraction[norm_EC$response == var1],
                    norm_EC$expansion_contraction[norm_EC$response == var2], 
                    paired = FALSE)$p.value  
pvalues <- c(pval)
ggplot(norm_EC, aes(x = response, y = expansion_contraction)) + 
  geom_boxplot(aes(fill = response), color = outline_colors, fill = box_fill_colors, alpha = 0.2, size = 1.5) + 
  geom_point(aes(colour = response, fill = response), shape = 21, size = 5) +  
  scale_fill_manual(values= point_fill_colors) + 
  scale_colour_manual(values=outline_colors) + 
  ylab('Fraction of expansionAB - contractionBC clones') + 
  geom_text_repel(aes(label = patient_number), max.overlaps = Inf, box.padding = 0.5, point.padding = 0.3, segment.size = 0.2, segment.color = 'grey50') +  
  scale_x_discrete(labels = setNames(c(var1, var2), c(var1, var2))) +  
  theme(
    axis.line = element_line(colour = "black", size = 1, linetype = "solid"),
    axis.text.x = element_text(size = 14, color='black'), 
    axis.text.y = element_text(size = 14, color='black'),
    axis.title = element_text(size = 18, face = "bold"),
    axis.ticks = element_line(size = 1),
    panel.background = element_rect(fill = 'white', color = 'black'),
    legend.position = "none"
  ) +
  labs(x = 'Response') +  
  geom_signif(comparisons = list(c(var1, var2)),
               annotations = sprintf("p = %.3f", pvalues),
               map_signif_level = FALSE,
               test = "wilcox.test", 
               step_increase = 0.1)  		   
dev.off()	



########################################################################
############# EXTENDED DATA Figure 8 F - singletons across time      ##
########################################################################

library (tidyr)
library (scRepertoire)
library (stringi)
library (stringr)
library (ggplot2)
library (ggrepel)
library(ggsignif)
library  (dplyr)
library(patchwork)
library(RColorBrewer)
library(ggalluvial)
library(pheatmap)
library(amap)


load ('TCR_norm_cell_bind.RData')
#### load the summarized expansion_contraction data from PP_Expansion
load("metadata_git.RData")
#### load the metadata available at EGA
merA <- data.frame(read.table("clinical_data.txt", header=T, sep='\t', row.names=))
##### load clinical data in Supp Table 1

merA <- merA [!merA$patient_number=="4",]

colnames (merA) <- gsub ("ME", "response",colnames (merA))
all (rownames (statusABC1_n_all) == paste ("NL", merA$patient_numb, sep=""))
response <- merA$response
patient_number <- merA$patient_numb

norm_all <- statusABC1_n_all 
norm_all$response <- response
norm_all$patient_number <- patient_number
colnames (norm_all) <- gsub ('all_ABC_' , '', colnames (norm_all))


normBR <- norm_all [norm_all$response=="BR",]
normGR <- norm_all [norm_all$response=="GR",]


normGR_long <- as.data.frame (normGR %>%
  pivot_longer(cols = c(singletonA, singletonB, singletonC), 
               names_to = "singleton_type", 
               values_to = "value"))
pAB <- wilcox.test(normGR_long$value[normGR_long$singleton_type == "singletonA"],
                    normGR_long$value[normGR_long$singleton_type == "singletonB"], 
                    paired = F)$p.value
pBC <- wilcox.test(normGR_long$value[normGR_long$singleton_type == "singletonB"],
                    normGR_long$value[normGR_long$singleton_type == "singletonC"], 
                    paired = F)$p.value 
pAC <- wilcox.test(normGR_long$value[normGR_long$singleton_type == "singletonA"],
                    normGR_long$value[normGR_long$singleton_type == "singletonC"], 
                    paired = F)$p.value 
pvalues <- c(pAB, pBC, pAC)

outline_colors <- c("#AA3830", '#9D2E99', '#350D7F') # Box outline colors
box_fill_colors <- c('#FFE2E0', '#FFD1FF', '#DAC3FF') # Box fill colors
point_fill_colors <- c('#FF8077', '#E66CE1', '#7E33FF') # Point fill colors

pdf ('ExDaFigure8F singletons_GR.pdf')

ggplot(normGR_long, aes(singleton_type, value)) +
  geom_boxplot(color = outline_colors, 
               fill = box_fill_colors, 
               alpha = 0.5, size = 1.5) +
  geom_point(aes(colour = factor(singleton_type), fill = factor(singleton_type)), shape = 21, size = 5, stroke = 1.2) +  # Increase stroke to make outline visible
  scale_fill_manual(values = point_fill_colors) +
  scale_colour_manual(values = outline_colors) +  # Same colors for outline and fill
  ylab('Percentage of Singleton') + 
  geom_text_repel(aes(label = patient_number), min.segment.length = 1, max.overlaps = Inf) + 
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
  theme(legend.position = "none") +
        geom_signif(comparisons = list(c("singletonA", "singletonB"), c("singletonB", "singletonC"), c("singletonA", "singletonC")),
                annotations = sprintf("p = %.3f", pvalues),
                map_signif_level = FALSE,
                test = "t.test", 
                step_increase = 0.1)
dev.off()

########################################################################
############# EXTENDED DATA Figure 8 G - singletons across time      ##
########################################################################

library (tidyr)
library (scRepertoire)
library (stringi)
library (stringr)
library (ggplot2)
library (ggrepel)
library(ggsignif)
library  (dplyr)
library(patchwork)
library(RColorBrewer)
library(ggalluvial)
library(pheatmap)
library(amap)


load ('TCR_norm_cell_bind.RData')
#### load the summarized expansion_contraction data from PP_Expansion
load("metadata_git.RData")
#### load the metadata available at EGA
merA <- data.frame(read.table("clinical_data.txt", header=T, sep='\t', row.names=))
##### load clinical data in Supp Table 1

merA <- data.frame(read.table("clinical_data.txt", header=T, sep='\t', row.names=))
merA <- merA [!merA$patient_number=="4",]

colnames (merA) <- gsub ("ME", "response",colnames (merA))
all (rownames (statusABC1_n_all) == paste ("NL", merA$patient_numb, sep=""))
response <- merA$response
patient_number <- merA$patient_numb

norm_all <- statusABC1_n_all 
norm_all$response <- response
norm_all$patient_number <- patient_number
colnames (norm_all) <- gsub ('all_ABC_' , '', colnames (norm_all))


normBR <- norm_all [norm_all$response=="BR",]
normGR <- norm_all [norm_all$response=="GR",]

normBR_long <- as.data.frame (normBR %>%
  pivot_longer(cols = c(singletonA, singletonB, singletonC), 
               names_to = "singleton_type", 
               values_to = "value"))
pAB <- wilcox.test(normBR_long$value[normBR_long$singleton_type == "singletonA"],
                    normBR_long$value[normBR_long$singleton_type == "singletonB"], 
                    paired = F)$p.value
pBC <- wilcox.test(normBR_long$value[normBR_long$singleton_type == "singletonB"],
                    normBR_long$value[normBR_long$singleton_type == "singletonC"], 
                    paired = F)$p.value 
pAC <- wilcox.test(normBR_long$value[normBR_long$singleton_type == "singletonA"],
                    normBR_long$value[normBR_long$singleton_type == "singletonC"], 
                    paired = F)$p.value 
pvalues <- c(pAB, pBC, pAC)

pdf ('ExDaFigure8G singletons_BR.pdf')
ggplot(normBR_long, aes(singleton_type, value)) +
  geom_boxplot(color = outline_colors, 
               fill = box_fill_colors, 
               alpha = 0.5, size = 1.5) +
  geom_point(aes(colour = factor(singleton_type), fill = factor(singleton_type)), shape = 21, size = 5, stroke = 1.2) +  # Increase stroke to make outline visible
  scale_fill_manual(values = point_fill_colors) +
  scale_colour_manual(values = outline_colors) +  # Use the color mapping for lines
  ylab('Percentage of Singleton') + 
  geom_text_repel(aes(label = patient_number), min.segment.length = 1, max.overlaps = Inf) + 
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
  theme(legend.position = "none") +
  geom_signif(comparisons = list(c("singletonA", "singletonB"), c("singletonB", "singletonC"), c("singletonA", "singletonC")),
     annotations = sprintf("p = %.3f", pvalues),
     map_signif_level = FALSE,
     test = "t.test", 
     step_increase = 0.1)

dev.off() 


########################################
###### EXTENDED DATA Figure 9 #########
########################################

################################################################################################################
############# EXTENDED DATA FIGURE 9A- DINUCLEOTIDES    #######################################################
################################################################################################################
################################################################################################################
############# EXTENDED DATA FIGURE 9B- DINUCLEOTIDES     ######################################################
################################################################################################################

################################################################################################################
############# EXTENDED DATA FIGURE 9C- BARPLOT INTERSECT ATAC chromHMM     ####################################
################################################################################################################
library(ggplot2)
library(tidyr)

data <- data.frame(
  Category = c("Iglesia", "dCRE", "Regner", "chromHMM"),
  Iglesia = c(9.61, 37.46, 56.1, 65.5),  
  dCRE = c(45.02, 43.9, 58.1, 71.4), 
  Regner = c(64.2, 17, 73.1, 89.3)
)

data_long <- pivot_longer(data, cols = -Category, names_to = "Source", values_to = "Value")

row_colors <- c("blue", "red", "green", "orange")  

pdf("ExDatFigure9C Intersection with scATAC-seq.pdf", width=14, height=8, onefile=FALSE, useDingbats=FALSE)

ggplot(data_long, aes(x = Source, y = Value, fill = Category)) +  
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_manual(values = row_colors) +
  labs(title = "",  
       y = "% of intersection",  
       x = "Distal Regulatory Elements") +  
  theme_minimal(base_size = 15) + 
  theme(
    legend.position = "none",       
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),  
    axis.title.x = element_text(size = 20),  
    axis.title.y = element_text(size = 20),  
    axis.text.x = element_text(angle = 45, hjust = 1, size = 15), 
    axis.text.y = element_text(size = 15),  
    axis.ticks.length = unit(0.25, "cm"), 
    axis.ticks = element_line(size = 1),   
    axis.line = element_line(size = 1)   
  )

dev.off()

################################################################################################################
############# EXTENDED DATA FIGURE 9D-E - UPSET PLOTS - proximal distal     ################################
################################################################################################################
library (UpSetR)
pdf("ExDaFigure9D Intervene_upset_distal.pdf", width=14, height=8, onefile=FALSE, useDingbats=FALSE)
expressionInput <- c('T_cells_0'=5952,
'Normal_epithelial_0'=169,
'Normal_epithelial_0&T_cells_0'=4,
'Myeloid_0'=13465,
'Myeloid_0&T_cells_0'=927,'Myeloid_0&Normal_epithelial_0'=119,'Myeloid_0&Normal_epithelial_0&T_cells_0'=0,'Malignant_cells_0'=10953,'Malignant_cells_0&T_cells_0'=164,'Malignant_cells_0&Normal_epithelial_0'=1582,'Malignant_cells_0&Normal_epithelial_0&T_cells_0'=7,'Malignant_cells_0&Myeloid_0'=382,'Malignant_cells_0&Myeloid_0&T_cells_0'=10,'Malignant_cells_0&Myeloid_0&Normal_epithelial_0'=23,'Malignant_cells_0&Myeloid_0&Normal_epithelial_0&T_cells_0'=1,'Fibroblasts_0'=8407,'Fibroblasts_0&T_cells_0'=121,'Fibroblasts_0&Normal_epithelial_0'=37,'Fibroblasts_0&Normal_epithelial_0&T_cells_0'=0,'Fibroblasts_0&Myeloid_0'=578,'Fibroblasts_0&Myeloid_0&T_cells_0'=11,'Fibroblasts_0&Myeloid_0&Normal_epithelial_0'=1,'Fibroblasts_0&Myeloid_0&Normal_epithelial_0&T_cells_0'=0,'Fibroblasts_0&Malignant_cells_0'=487,'Fibroblasts_0&Malignant_cells_0&T_cells_0'=2,'Fibroblasts_0&Malignant_cells_0&Normal_epithelial_0'=31,'Fibroblasts_0&Malignant_cells_0&Normal_epithelial_0&T_cells_0'=0,'Fibroblasts_0&Malignant_cells_0&Myeloid_0'=22,'Fibroblasts_0&Malignant_cells_0&Myeloid_0&T_cells_0'=0,'Fibroblasts_0&Malignant_cells_0&Myeloid_0&Normal_epithelial_0'=0,'Fibroblasts_0&Malignant_cells_0&Myeloid_0&Normal_epithelial_0&T_cells_0'=0,'Endothelial_cells_0'=3811,'Endothelial_cells_0&T_cells_0'=185,'Endothelial_cells_0&Normal_epithelial_0'=14,'Endothelial_cells_0&Normal_epithelial_0&T_cells_0'=0,'Endothelial_cells_0&Myeloid_0'=564,'Endothelial_cells_0&Myeloid_0&T_cells_0'=30,'Endothelial_cells_0&Myeloid_0&Normal_epithelial_0'=3,'Endothelial_cells_0&Myeloid_0&Normal_epithelial_0&T_cells_0'=0,'Endothelial_cells_0&Malignant_cells_0'=289,'Endothelial_cells_0&Malignant_cells_0&T_cells_0'=4,'Endothelial_cells_0&Malignant_cells_0&Normal_epithelial_0'=16,'Endothelial_cells_0&Malignant_cells_0&Normal_epithelial_0&T_cells_0'=0,'Endothelial_cells_0&Malignant_cells_0&Myeloid_0'=14,'Endothelial_cells_0&Malignant_cells_0&Myeloid_0&T_cells_0'=0,'Endothelial_cells_0&Malignant_cells_0&Myeloid_0&Normal_epithelial_0'=1,'Endothelial_cells_0&Malignant_cells_0&Myeloid_0&Normal_epithelial_0&T_cells_0'=0,'Endothelial_cells_0&Fibroblasts_0'=786,'Endothelial_cells_0&Fibroblasts_0&T_cells_0'=11,'Endothelial_cells_0&Fibroblasts_0&Normal_epithelial_0'=2,'Endothelial_cells_0&Fibroblasts_0&Normal_epithelial_0&T_cells_0'=0,'Endothelial_cells_0&Fibroblasts_0&Myeloid_0'=58,'Endothelial_cells_0&Fibroblasts_0&Myeloid_0&T_cells_0'=3,'Endothelial_cells_0&Fibroblasts_0&Myeloid_0&Normal_epithelial_0'=0,'Endothelial_cells_0&Fibroblasts_0&Myeloid_0&Normal_epithelial_0&T_cells_0'=0,'Endothelial_cells_0&Fibroblasts_0&Malignant_cells_0'=20,'Endothelial_cells_0&Fibroblasts_0&Malignant_cells_0&T_cells_0'=0,'Endothelial_cells_0&Fibroblasts_0&Malignant_cells_0&Normal_epithelial_0'=3,'Endothelial_cells_0&Fibroblasts_0&Malignant_cells_0&Normal_epithelial_0&T_cells_0'=0,'Endothelial_cells_0&Fibroblasts_0&Malignant_cells_0&Myeloid_0'=0,'Endothelial_cells_0&Fibroblasts_0&Malignant_cells_0&Myeloid_0&T_cells_0'=0,'Endothelial_cells_0&Fibroblasts_0&Malignant_cells_0&Myeloid_0&Normal_epithelial_0'=0,'Endothelial_cells_0&Fibroblasts_0&Malignant_cells_0&Myeloid_0&Normal_epithelial_0&T_cells_0'=0,'BP_cells_0'=2238,'BP_cells_0&T_cells_0'=809,'BP_cells_0&Normal_epithelial_0'=3,'BP_cells_0&Normal_epithelial_0&T_cells_0'=0,'BP_cells_0&Myeloid_0'=932,'BP_cells_0&Myeloid_0&T_cells_0'=274,'BP_cells_0&Myeloid_0&Normal_epithelial_0'=1,'BP_cells_0&Myeloid_0&Normal_epithelial_0&T_cells_0'=1,'BP_cells_0&Malignant_cells_0'=770,'BP_cells_0&Malignant_cells_0&T_cells_0'=38,'BP_cells_0&Malignant_cells_0&Normal_epithelial_0'=15,'BP_cells_0&Malignant_cells_0&Normal_epithelial_0&T_cells_0'=0,'BP_cells_0&Malignant_cells_0&Myeloid_0'=47,'BP_cells_0&Malignant_cells_0&Myeloid_0&T_cells_0'=4,'BP_cells_0&Malignant_cells_0&Myeloid_0&Normal_epithelial_0'=0,'BP_cells_0&Malignant_cells_0&Myeloid_0&Normal_epithelial_0&T_cells_0'=0,'BP_cells_0&Fibroblasts_0'=79,'BP_cells_0&Fibroblasts_0&T_cells_0'=27,'BP_cells_0&Fibroblasts_0&Normal_epithelial_0'=0,'BP_cells_0&Fibroblasts_0&Normal_epithelial_0&T_cells_0'=0,'BP_cells_0&Fibroblasts_0&Myeloid_0'=31,'BP_cells_0&Fibroblasts_0&Myeloid_0&T_cells_0'=9,'BP_cells_0&Fibroblasts_0&Myeloid_0&Normal_epithelial_0'=0,'BP_cells_0&Fibroblasts_0&Myeloid_0&Normal_epithelial_0&T_cells_0'=0,'BP_cells_0&Fibroblasts_0&Malignant_cells_0'=16,'BP_cells_0&Fibroblasts_0&Malignant_cells_0&T_cells_0'=3,'BP_cells_0&Fibroblasts_0&Malignant_cells_0&Normal_epithelial_0'=0,'BP_cells_0&Fibroblasts_0&Malignant_cells_0&Normal_epithelial_0&T_cells_0'=0,'BP_cells_0&Fibroblasts_0&Malignant_cells_0&Myeloid_0'=1,'BP_cells_0&Fibroblasts_0&Malignant_cells_0&Myeloid_0&T_cells_0'=0,'BP_cells_0&Fibroblasts_0&Malignant_cells_0&Myeloid_0&Normal_epithelial_0'=0,'BP_cells_0&Fibroblasts_0&Malignant_cells_0&Myeloid_0&Normal_epithelial_0&T_cells_0'=0,'BP_cells_0&Endothelial_cells_0'=61,'BP_cells_0&Endothelial_cells_0&T_cells_0'=15,'BP_cells_0&Endothelial_cells_0&Normal_epithelial_0'=1,'BP_cells_0&Endothelial_cells_0&Normal_epithelial_0&T_cells_0'=0,'BP_cells_0&Endothelial_cells_0&Myeloid_0'=27,'BP_cells_0&Endothelial_cells_0&Myeloid_0&T_cells_0'=5,'BP_cells_0&Endothelial_cells_0&Myeloid_0&Normal_epithelial_0'=0,'BP_cells_0&Endothelial_cells_0&Myeloid_0&Normal_epithelial_0&T_cells_0'=0,'BP_cells_0&Endothelial_cells_0&Malignant_cells_0'=6,'BP_cells_0&Endothelial_cells_0&Malignant_cells_0&T_cells_0'=0,'BP_cells_0&Endothelial_cells_0&Malignant_cells_0&Normal_epithelial_0'=0,'BP_cells_0&Endothelial_cells_0&Malignant_cells_0&Normal_epithelial_0&T_cells_0'=0,'BP_cells_0&Endothelial_cells_0&Malignant_cells_0&Myeloid_0'=1,'BP_cells_0&Endothelial_cells_0&Malignant_cells_0&Myeloid_0&T_cells_0'=0,'BP_cells_0&Endothelial_cells_0&Malignant_cells_0&Myeloid_0&Normal_epithelial_0'=0,'BP_cells_0&Endothelial_cells_0&Malignant_cells_0&Myeloid_0&Normal_epithelial_0&T_cells_0'=0,'BP_cells_0&Endothelial_cells_0&Fibroblasts_0'=16,'BP_cells_0&Endothelial_cells_0&Fibroblasts_0&T_cells_0'=3,'BP_cells_0&Endothelial_cells_0&Fibroblasts_0&Normal_epithelial_0'=0,'BP_cells_0&Endothelial_cells_0&Fibroblasts_0&Normal_epithelial_0&T_cells_0'=0,'BP_cells_0&Endothelial_cells_0&Fibroblasts_0&Myeloid_0'=2,'BP_cells_0&Endothelial_cells_0&Fibroblasts_0&Myeloid_0&T_cells_0'=1,'BP_cells_0&Endothelial_cells_0&Fibroblasts_0&Myeloid_0&Normal_epithelial_0'=0,'BP_cells_0&Endothelial_cells_0&Fibroblasts_0&Myeloid_0&Normal_epithelial_0&T_cells_0'=0,'BP_cells_0&Endothelial_cells_0&Fibroblasts_0&Malignant_cells_0'=0,'BP_cells_0&Endothelial_cells_0&Fibroblasts_0&Malignant_cells_0&T_cells_0'=1,'BP_cells_0&Endothelial_cells_0&Fibroblasts_0&Malignant_cells_0&Normal_epithelial_0'=0,'BP_cells_0&Endothelial_cells_0&Fibroblasts_0&Malignant_cells_0&Normal_epithelial_0&T_cells_0'=0,'BP_cells_0&Endothelial_cells_0&Fibroblasts_0&Malignant_cells_0&Myeloid_0'=0,'BP_cells_0&Endothelial_cells_0&Fibroblasts_0&Malignant_cells_0&Myeloid_0&T_cells_0'=0,'BP_cells_0&Endothelial_cells_0&Fibroblasts_0&Malignant_cells_0&Myeloid_0&Normal_epithelial_0'=0,'BP_cells_0&Endothelial_cells_0&Fibroblasts_0&Malignant_cells_0&Myeloid_0&Normal_epithelial_0&T_cells_0'=0,'Adipocytes_0'=165,'Adipocytes_0&T_cells_0'=9,'Adipocytes_0&Normal_epithelial_0'=5,'Adipocytes_0&Normal_epithelial_0&T_cells_0'=0,'Adipocytes_0&Myeloid_0'=29,'Adipocytes_0&Myeloid_0&T_cells_0'=0,'Adipocytes_0&Myeloid_0&Normal_epithelial_0'=0,'Adipocytes_0&Myeloid_0&Normal_epithelial_0&T_cells_0'=0,'Adipocytes_0&Malignant_cells_0'=31,'Adipocytes_0&Malignant_cells_0&T_cells_0'=0,'Adipocytes_0&Malignant_cells_0&Normal_epithelial_0'=5,'Adipocytes_0&Malignant_cells_0&Normal_epithelial_0&T_cells_0'=0,'Adipocytes_0&Malignant_cells_0&Myeloid_0'=2,'Adipocytes_0&Malignant_cells_0&Myeloid_0&T_cells_0'=0,'Adipocytes_0&Malignant_cells_0&Myeloid_0&Normal_epithelial_0'=0,'Adipocytes_0&Malignant_cells_0&Myeloid_0&Normal_epithelial_0&T_cells_0'=0,'Adipocytes_0&Fibroblasts_0'=129,'Adipocytes_0&Fibroblasts_0&T_cells_0'=2,'Adipocytes_0&Fibroblasts_0&Normal_epithelial_0'=1,'Adipocytes_0&Fibroblasts_0&Normal_epithelial_0&T_cells_0'=0,'Adipocytes_0&Fibroblasts_0&Myeloid_0'=8,'Adipocytes_0&Fibroblasts_0&Myeloid_0&T_cells_0'=0,'Adipocytes_0&Fibroblasts_0&Myeloid_0&Normal_epithelial_0'=0,'Adipocytes_0&Fibroblasts_0&Myeloid_0&Normal_epithelial_0&T_cells_0'=0,'Adipocytes_0&Fibroblasts_0&Malignant_cells_0'=8,'Adipocytes_0&Fibroblasts_0&Malignant_cells_0&T_cells_0'=0,'Adipocytes_0&Fibroblasts_0&Malignant_cells_0&Normal_epithelial_0'=2,'Adipocytes_0&Fibroblasts_0&Malignant_cells_0&Normal_epithelial_0&T_cells_0'=0,'Adipocytes_0&Fibroblasts_0&Malignant_cells_0&Myeloid_0'=1,'Adipocytes_0&Fibroblasts_0&Malignant_cells_0&Myeloid_0&T_cells_0'=0,'Adipocytes_0&Fibroblasts_0&Malignant_cells_0&Myeloid_0&Normal_epithelial_0'=0,'Adipocytes_0&Fibroblasts_0&Malignant_cells_0&Myeloid_0&Normal_epithelial_0&T_cells_0'=0,'Adipocytes_0&Endothelial_cells_0'=183,'Adipocytes_0&Endothelial_cells_0&T_cells_0'=7,'Adipocytes_0&Endothelial_cells_0&Normal_epithelial_0'=0,'Adipocytes_0&Endothelial_cells_0&Normal_epithelial_0&T_cells_0'=0,'Adipocytes_0&Endothelial_cells_0&Myeloid_0'=18,'Adipocytes_0&Endothelial_cells_0&Myeloid_0&T_cells_0'=0,'Adipocytes_0&Endothelial_cells_0&Myeloid_0&Normal_epithelial_0'=0,'Adipocytes_0&Endothelial_cells_0&Myeloid_0&Normal_epithelial_0&T_cells_0'=0,'Adipocytes_0&Endothelial_cells_0&Malignant_cells_0'=12,'Adipocytes_0&Endothelial_cells_0&Malignant_cells_0&T_cells_0'=0,'Adipocytes_0&Endothelial_cells_0&Malignant_cells_0&Normal_epithelial_0'=2,'Adipocytes_0&Endothelial_cells_0&Malignant_cells_0&Normal_epithelial_0&T_cells_0'=0,'Adipocytes_0&Endothelial_cells_0&Malignant_cells_0&Myeloid_0'=1,'Adipocytes_0&Endothelial_cells_0&Malignant_cells_0&Myeloid_0&T_cells_0'=0,'Adipocytes_0&Endothelial_cells_0&Malignant_cells_0&Myeloid_0&Normal_epithelial_0'=0,'Adipocytes_0&Endothelial_cells_0&Malignant_cells_0&Myeloid_0&Normal_epithelial_0&T_cells_0'=0,'Adipocytes_0&Endothelial_cells_0&Fibroblasts_0'=53,'Adipocytes_0&Endothelial_cells_0&Fibroblasts_0&T_cells_0'=0,'Adipocytes_0&Endothelial_cells_0&Fibroblasts_0&Normal_epithelial_0'=1,'Adipocytes_0&Endothelial_cells_0&Fibroblasts_0&Normal_epithelial_0&T_cells_0'=0,'Adipocytes_0&Endothelial_cells_0&Fibroblasts_0&Myeloid_0'=6,'Adipocytes_0&Endothelial_cells_0&Fibroblasts_0&Myeloid_0&T_cells_0'=0,'Adipocytes_0&Endothelial_cells_0&Fibroblasts_0&Myeloid_0&Normal_epithelial_0'=0,'Adipocytes_0&Endothelial_cells_0&Fibroblasts_0&Myeloid_0&Normal_epithelial_0&T_cells_0'=0,'Adipocytes_0&Endothelial_cells_0&Fibroblasts_0&Malignant_cells_0'=5,'Adipocytes_0&Endothelial_cells_0&Fibroblasts_0&Malignant_cells_0&T_cells_0'=0,'Adipocytes_0&Endothelial_cells_0&Fibroblasts_0&Malignant_cells_0&Normal_epithelial_0'=0,'Adipocytes_0&Endothelial_cells_0&Fibroblasts_0&Malignant_cells_0&Normal_epithelial_0&T_cells_0'=0,'Adipocytes_0&Endothelial_cells_0&Fibroblasts_0&Malignant_cells_0&Myeloid_0'=1,'Adipocytes_0&Endothelial_cells_0&Fibroblasts_0&Malignant_cells_0&Myeloid_0&T_cells_0'=0,'Adipocytes_0&Endothelial_cells_0&Fibroblasts_0&Malignant_cells_0&Myeloid_0&Normal_epithelial_0'=0,'Adipocytes_0&Endothelial_cells_0&Fibroblasts_0&Malignant_cells_0&Myeloid_0&Normal_epithelial_0&T_cells_0'=0,'Adipocytes_0&BP_cells_0'=4,'Adipocytes_0&BP_cells_0&T_cells_0'=4,'Adipocytes_0&BP_cells_0&Normal_epithelial_0'=0,'Adipocytes_0&BP_cells_0&Normal_epithelial_0&T_cells_0'=0,'Adipocytes_0&BP_cells_0&Myeloid_0'=3,'Adipocytes_0&BP_cells_0&Myeloid_0&T_cells_0'=1,'Adipocytes_0&BP_cells_0&Myeloid_0&Normal_epithelial_0'=0,'Adipocytes_0&BP_cells_0&Myeloid_0&Normal_epithelial_0&T_cells_0'=0,'Adipocytes_0&BP_cells_0&Malignant_cells_0'=1,'Adipocytes_0&BP_cells_0&Malignant_cells_0&T_cells_0'=1,'Adipocytes_0&BP_cells_0&Malignant_cells_0&Normal_epithelial_0'=0,'Adipocytes_0&BP_cells_0&Malignant_cells_0&Normal_epithelial_0&T_cells_0'=0,'Adipocytes_0&BP_cells_0&Malignant_cells_0&Myeloid_0'=1,'Adipocytes_0&BP_cells_0&Malignant_cells_0&Myeloid_0&T_cells_0'=0,'Adipocytes_0&BP_cells_0&Malignant_cells_0&Myeloid_0&Normal_epithelial_0'=0,'Adipocytes_0&BP_cells_0&Malignant_cells_0&Myeloid_0&Normal_epithelial_0&T_cells_0'=0,'Adipocytes_0&BP_cells_0&Fibroblasts_0'=0,'Adipocytes_0&BP_cells_0&Fibroblasts_0&T_cells_0'=1,'Adipocytes_0&BP_cells_0&Fibroblasts_0&Normal_epithelial_0'=0,'Adipocytes_0&BP_cells_0&Fibroblasts_0&Normal_epithelial_0&T_cells_0'=0,'Adipocytes_0&BP_cells_0&Fibroblasts_0&Myeloid_0'=3,'Adipocytes_0&BP_cells_0&Fibroblasts_0&Myeloid_0&T_cells_0'=0,'Adipocytes_0&BP_cells_0&Fibroblasts_0&Myeloid_0&Normal_epithelial_0'=0,'Adipocytes_0&BP_cells_0&Fibroblasts_0&Myeloid_0&Normal_epithelial_0&T_cells_0'=0,'Adipocytes_0&BP_cells_0&Fibroblasts_0&Malignant_cells_0'=0,'Adipocytes_0&BP_cells_0&Fibroblasts_0&Malignant_cells_0&T_cells_0'=0,'Adipocytes_0&BP_cells_0&Fibroblasts_0&Malignant_cells_0&Normal_epithelial_0'=0,'Adipocytes_0&BP_cells_0&Fibroblasts_0&Malignant_cells_0&Normal_epithelial_0&T_cells_0'=0,'Adipocytes_0&BP_cells_0&Fibroblasts_0&Malignant_cells_0&Myeloid_0'=0,'Adipocytes_0&BP_cells_0&Fibroblasts_0&Malignant_cells_0&Myeloid_0&T_cells_0'=0,'Adipocytes_0&BP_cells_0&Fibroblasts_0&Malignant_cells_0&Myeloid_0&Normal_epithelial_0'=0,'Adipocytes_0&BP_cells_0&Fibroblasts_0&Malignant_cells_0&Myeloid_0&Normal_epithelial_0&T_cells_0'=0,'Adipocytes_0&BP_cells_0&Endothelial_cells_0'=4,'Adipocytes_0&BP_cells_0&Endothelial_cells_0&T_cells_0'=1,'Adipocytes_0&BP_cells_0&Endothelial_cells_0&Normal_epithelial_0'=0,'Adipocytes_0&BP_cells_0&Endothelial_cells_0&Normal_epithelial_0&T_cells_0'=0,'Adipocytes_0&BP_cells_0&Endothelial_cells_0&Myeloid_0'=0,'Adipocytes_0&BP_cells_0&Endothelial_cells_0&Myeloid_0&T_cells_0'=2,'Adipocytes_0&BP_cells_0&Endothelial_cells_0&Myeloid_0&Normal_epithelial_0'=0,'Adipocytes_0&BP_cells_0&Endothelial_cells_0&Myeloid_0&Normal_epithelial_0&T_cells_0'=0,'Adipocytes_0&BP_cells_0&Endothelial_cells_0&Malignant_cells_0'=0,'Adipocytes_0&BP_cells_0&Endothelial_cells_0&Malignant_cells_0&T_cells_0'=0,'Adipocytes_0&BP_cells_0&Endothelial_cells_0&Malignant_cells_0&Normal_epithelial_0'=0,'Adipocytes_0&BP_cells_0&Endothelial_cells_0&Malignant_cells_0&Normal_epithelial_0&T_cells_0'=0,'Adipocytes_0&BP_cells_0&Endothelial_cells_0&Malignant_cells_0&Myeloid_0'=0,'Adipocytes_0&BP_cells_0&Endothelial_cells_0&Malignant_cells_0&Myeloid_0&T_cells_0'=0,'Adipocytes_0&BP_cells_0&Endothelial_cells_0&Malignant_cells_0&Myeloid_0&Normal_epithelial_0'=0,'Adipocytes_0&BP_cells_0&Endothelial_cells_0&Malignant_cells_0&Myeloid_0&Normal_epithelial_0&T_cells_0'=0,'Adipocytes_0&BP_cells_0&Endothelial_cells_0&Fibroblasts_0'=6,'Adipocytes_0&BP_cells_0&Endothelial_cells_0&Fibroblasts_0&T_cells_0'=0,'Adipocytes_0&BP_cells_0&Endothelial_cells_0&Fibroblasts_0&Normal_epithelial_0'=0,'Adipocytes_0&BP_cells_0&Endothelial_cells_0&Fibroblasts_0&Normal_epithelial_0&T_cells_0'=0,'Adipocytes_0&BP_cells_0&Endothelial_cells_0&Fibroblasts_0&Myeloid_0'=1,'Adipocytes_0&BP_cells_0&Endothelial_cells_0&Fibroblasts_0&Myeloid_0&T_cells_0'=0,'Adipocytes_0&BP_cells_0&Endothelial_cells_0&Fibroblasts_0&Myeloid_0&Normal_epithelial_0'=0,'Adipocytes_0&BP_cells_0&Endothelial_cells_0&Fibroblasts_0&Myeloid_0&Normal_epithelial_0&T_cells_0'=0,'Adipocytes_0&BP_cells_0&Endothelial_cells_0&Fibroblasts_0&Malignant_cells_0'=0,'Adipocytes_0&BP_cells_0&Endothelial_cells_0&Fibroblasts_0&Malignant_cells_0&T_cells_0'=0,'Adipocytes_0&BP_cells_0&Endothelial_cells_0&Fibroblasts_0&Malignant_cells_0&Normal_epithelial_0'=0,'Adipocytes_0&BP_cells_0&Endothelial_cells_0&Fibroblasts_0&Malignant_cells_0&Normal_epithelial_0&T_cells_0'=0,'Adipocytes_0&BP_cells_0&Endothelial_cells_0&Fibroblasts_0&Malignant_cells_0&Myeloid_0'=0,'Adipocytes_0&BP_cells_0&Endothelial_cells_0&Fibroblasts_0&Malignant_cells_0&Myeloid_0&T_cells_0'=0,'Adipocytes_0&BP_cells_0&Endothelial_cells_0&Fibroblasts_0&Malignant_cells_0&Myeloid_0&Normal_epithelial_0'=0,'Adipocytes_0&BP_cells_0&Endothelial_cells_0&Fibroblasts_0&Malignant_cells_0&Myeloid_0&Normal_epithelial_0&T_cells_0'=0)
upset(fromExpression(expressionInput), nsets=8, nintersects=30, show.numbers="yes", main.bar.color="#ea5d4e", sets.bar.color="#317eab", empty.intersections=NULL, order.by = "freq", number.angles = 0, mainbar.y.label ="No. of Intersections", sets.x.label ="Set size")
invisible(dev.off())

library (UpSetR)
pdf("ExDaFigure9E Intervene_upset_proximal.pdf", width=14, height=8, onefile=FALSE, useDingbats=FALSE)
expressionInput <- c('T_cells_0'=1348,'Normal_epithelial_0'=135,'Normal_epithelial_0&T_cells_0'=18,'Myeloid_0'=1567,'Myeloid_0&T_cells_0'=337,'Myeloid_0&Normal_epithelial_0'=77,'Myeloid_0&Normal_epithelial_0&T_cells_0'=4,'Malignant_cells_0'=4338,'Malignant_cells_0&T_cells_0'=347,'Malignant_cells_0&Normal_epithelial_0'=776,'Malignant_cells_0&Normal_epithelial_0&T_cells_0'=27,'Malignant_cells_0&Myeloid_0'=410,'Malignant_cells_0&Myeloid_0&T_cells_0'=116,'Malignant_cells_0&Myeloid_0&Normal_epithelial_0'=48,'Malignant_cells_0&Myeloid_0&Normal_epithelial_0&T_cells_0'=3,'Fibroblasts_0'=3745,'Fibroblasts_0&T_cells_0'=423,'Fibroblasts_0&Normal_epithelial_0'=124,'Fibroblasts_0&Normal_epithelial_0&T_cells_0'=5,'Fibroblasts_0&Myeloid_0'=581,'Fibroblasts_0&Myeloid_0&T_cells_0'=175,'Fibroblasts_0&Myeloid_0&Normal_epithelial_0'=24,'Fibroblasts_0&Myeloid_0&Normal_epithelial_0&T_cells_0'=1,'Fibroblasts_0&Malignant_cells_0'=2476,'Fibroblasts_0&Malignant_cells_0&T_cells_0'=445,'Fibroblasts_0&Malignant_cells_0&Normal_epithelial_0'=120,'Fibroblasts_0&Malignant_cells_0&Normal_epithelial_0&T_cells_0'=14,'Fibroblasts_0&Malignant_cells_0&Myeloid_0'=487,'Fibroblasts_0&Malignant_cells_0&Myeloid_0&T_cells_0'=136,'Fibroblasts_0&Malignant_cells_0&Myeloid_0&Normal_epithelial_0'=11,'Fibroblasts_0&Malignant_cells_0&Myeloid_0&Normal_epithelial_0&T_cells_0'=1,'Endothelial_cells_0'=695,'Endothelial_cells_0&T_cells_0'=80,'Endothelial_cells_0&Normal_epithelial_0'=18,'Endothelial_cells_0&Normal_epithelial_0&T_cells_0'=1,'Endothelial_cells_0&Myeloid_0'=147,'Endothelial_cells_0&Myeloid_0&T_cells_0'=28,'Endothelial_cells_0&Myeloid_0&Normal_epithelial_0'=9,'Endothelial_cells_0&Myeloid_0&Normal_epithelial_0&T_cells_0'=0,'Endothelial_cells_0&Malignant_cells_0'=229,'Endothelial_cells_0&Malignant_cells_0&T_cells_0'=16,'Endothelial_cells_0&Malignant_cells_0&Normal_epithelial_0'=29,'Endothelial_cells_0&Malignant_cells_0&Normal_epithelial_0&T_cells_0'=0,'Endothelial_cells_0&Malignant_cells_0&Myeloid_0'=30,'Endothelial_cells_0&Malignant_cells_0&Myeloid_0&T_cells_0'=5,'Endothelial_cells_0&Malignant_cells_0&Myeloid_0&Normal_epithelial_0'=5,'Endothelial_cells_0&Malignant_cells_0&Myeloid_0&Normal_epithelial_0&T_cells_0'=0,'Endothelial_cells_0&Fibroblasts_0'=625,'Endothelial_cells_0&Fibroblasts_0&T_cells_0'=91,'Endothelial_cells_0&Fibroblasts_0&Normal_epithelial_0'=16,'Endothelial_cells_0&Fibroblasts_0&Normal_epithelial_0&T_cells_0'=1,'Endothelial_cells_0&Fibroblasts_0&Myeloid_0'=165,'Endothelial_cells_0&Fibroblasts_0&Myeloid_0&T_cells_0'=36,'Endothelial_cells_0&Fibroblasts_0&Myeloid_0&Normal_epithelial_0'=2,'Endothelial_cells_0&Fibroblasts_0&Myeloid_0&Normal_epithelial_0&T_cells_0'=0,'Endothelial_cells_0&Fibroblasts_0&Malignant_cells_0'=266,'Endothelial_cells_0&Fibroblasts_0&Malignant_cells_0&T_cells_0'=26,'Endothelial_cells_0&Fibroblasts_0&Malignant_cells_0&Normal_epithelial_0'=12,'Endothelial_cells_0&Fibroblasts_0&Malignant_cells_0&Normal_epithelial_0&T_cells_0'=0,'Endothelial_cells_0&Fibroblasts_0&Malignant_cells_0&Myeloid_0'=39,'Endothelial_cells_0&Fibroblasts_0&Malignant_cells_0&Myeloid_0&T_cells_0'=7,'Endothelial_cells_0&Fibroblasts_0&Malignant_cells_0&Myeloid_0&Normal_epithelial_0'=1,'Endothelial_cells_0&Fibroblasts_0&Malignant_cells_0&Myeloid_0&Normal_epithelial_0&T_cells_0'=0,'BP_cells_0'=507,'BP_cells_0&T_cells_0'=627,'BP_cells_0&Normal_epithelial_0'=6,'BP_cells_0&Normal_epithelial_0&T_cells_0'=4,'BP_cells_0&Myeloid_0'=368,'BP_cells_0&Myeloid_0&T_cells_0'=628,'BP_cells_0&Myeloid_0&Normal_epithelial_0'=5,'BP_cells_0&Myeloid_0&Normal_epithelial_0&T_cells_0'=11,'BP_cells_0&Malignant_cells_0'=475,'BP_cells_0&Malignant_cells_0&T_cells_0'=466,'BP_cells_0&Malignant_cells_0&Normal_epithelial_0'=36,'BP_cells_0&Malignant_cells_0&Normal_epithelial_0&T_cells_0'=10,'BP_cells_0&Malignant_cells_0&Myeloid_0'=132,'BP_cells_0&Malignant_cells_0&Myeloid_0&T_cells_0'=224,'BP_cells_0&Malignant_cells_0&Myeloid_0&Normal_epithelial_0'=14,'BP_cells_0&Malignant_cells_0&Myeloid_0&Normal_epithelial_0&T_cells_0'=17,'BP_cells_0&Fibroblasts_0'=207,'BP_cells_0&Fibroblasts_0&T_cells_0'=616,'BP_cells_0&Fibroblasts_0&Normal_epithelial_0'=3,'BP_cells_0&Fibroblasts_0&Normal_epithelial_0&T_cells_0'=1,'BP_cells_0&Fibroblasts_0&Myeloid_0'=140,'BP_cells_0&Fibroblasts_0&Myeloid_0&T_cells_0'=347,'BP_cells_0&Fibroblasts_0&Myeloid_0&Normal_epithelial_0'=0,'BP_cells_0&Fibroblasts_0&Myeloid_0&Normal_epithelial_0&T_cells_0'=2,'BP_cells_0&Fibroblasts_0&Malignant_cells_0'=450,'BP_cells_0&Fibroblasts_0&Malignant_cells_0&T_cells_0'=669,'BP_cells_0&Fibroblasts_0&Malignant_cells_0&Normal_epithelial_0'=11,'BP_cells_0&Fibroblasts_0&Malignant_cells_0&Normal_epithelial_0&T_cells_0'=4,'BP_cells_0&Fibroblasts_0&Malignant_cells_0&Myeloid_0'=164,'BP_cells_0&Fibroblasts_0&Malignant_cells_0&Myeloid_0&T_cells_0'=244,'BP_cells_0&Fibroblasts_0&Malignant_cells_0&Myeloid_0&Normal_epithelial_0'=7,'BP_cells_0&Fibroblasts_0&Malignant_cells_0&Myeloid_0&Normal_epithelial_0&T_cells_0'=10,'BP_cells_0&Endothelial_cells_0'=29,'BP_cells_0&Endothelial_cells_0&T_cells_0'=49,'BP_cells_0&Endothelial_cells_0&Normal_epithelial_0'=1,'BP_cells_0&Endothelial_cells_0&Normal_epithelial_0&T_cells_0'=0,'BP_cells_0&Endothelial_cells_0&Myeloid_0'=29,'BP_cells_0&Endothelial_cells_0&Myeloid_0&T_cells_0'=73,'BP_cells_0&Endothelial_cells_0&Myeloid_0&Normal_epithelial_0'=2,'BP_cells_0&Endothelial_cells_0&Myeloid_0&Normal_epithelial_0&T_cells_0'=0,'BP_cells_0&Endothelial_cells_0&Malignant_cells_0'=22,'BP_cells_0&Endothelial_cells_0&Malignant_cells_0&T_cells_0'=19,'BP_cells_0&Endothelial_cells_0&Malignant_cells_0&Normal_epithelial_0'=1,'BP_cells_0&Endothelial_cells_0&Malignant_cells_0&Normal_epithelial_0&T_cells_0'=0,'BP_cells_0&Endothelial_cells_0&Malignant_cells_0&Myeloid_0'=8,'BP_cells_0&Endothelial_cells_0&Malignant_cells_0&Myeloid_0&T_cells_0'=7,'BP_cells_0&Endothelial_cells_0&Malignant_cells_0&Myeloid_0&Normal_epithelial_0'=0,'BP_cells_0&Endothelial_cells_0&Malignant_cells_0&Myeloid_0&Normal_epithelial_0&T_cells_0'=2,'BP_cells_0&Endothelial_cells_0&Fibroblasts_0'=53,'BP_cells_0&Endothelial_cells_0&Fibroblasts_0&T_cells_0'=58,'BP_cells_0&Endothelial_cells_0&Fibroblasts_0&Normal_epithelial_0'=0,'BP_cells_0&Endothelial_cells_0&Fibroblasts_0&Normal_epithelial_0&T_cells_0'=0,'BP_cells_0&Endothelial_cells_0&Fibroblasts_0&Myeloid_0'=46,'BP_cells_0&Endothelial_cells_0&Fibroblasts_0&Myeloid_0&T_cells_0'=36,'BP_cells_0&Endothelial_cells_0&Fibroblasts_0&Myeloid_0&Normal_epithelial_0'=0,'BP_cells_0&Endothelial_cells_0&Fibroblasts_0&Myeloid_0&Normal_epithelial_0&T_cells_0'=0,'BP_cells_0&Endothelial_cells_0&Fibroblasts_0&Malignant_cells_0'=23,'BP_cells_0&Endothelial_cells_0&Fibroblasts_0&Malignant_cells_0&T_cells_0'=25,'BP_cells_0&Endothelial_cells_0&Fibroblasts_0&Malignant_cells_0&Normal_epithelial_0'=0,'BP_cells_0&Endothelial_cells_0&Fibroblasts_0&Malignant_cells_0&Normal_epithelial_0&T_cells_0'=0,'BP_cells_0&Endothelial_cells_0&Fibroblasts_0&Malignant_cells_0&Myeloid_0'=11,'BP_cells_0&Endothelial_cells_0&Fibroblasts_0&Malignant_cells_0&Myeloid_0&T_cells_0'=4,'BP_cells_0&Endothelial_cells_0&Fibroblasts_0&Malignant_cells_0&Myeloid_0&Normal_epithelial_0'=0,'BP_cells_0&Endothelial_cells_0&Fibroblasts_0&Malignant_cells_0&Myeloid_0&Normal_epithelial_0&T_cells_0'=0,'Adipocytes_0'=82,'Adipocytes_0&T_cells_0'=11,'Adipocytes_0&Normal_epithelial_0'=11,'Adipocytes_0&Normal_epithelial_0&T_cells_0'=0,'Adipocytes_0&Myeloid_0'=37,'Adipocytes_0&Myeloid_0&T_cells_0'=8,'Adipocytes_0&Myeloid_0&Normal_epithelial_0'=3,'Adipocytes_0&Myeloid_0&Normal_epithelial_0&T_cells_0'=2,'Adipocytes_0&Malignant_cells_0'=95,'Adipocytes_0&Malignant_cells_0&T_cells_0'=8,'Adipocytes_0&Malignant_cells_0&Normal_epithelial_0'=32,'Adipocytes_0&Malignant_cells_0&Normal_epithelial_0&T_cells_0'=0,'Adipocytes_0&Malignant_cells_0&Myeloid_0'=16,'Adipocytes_0&Malignant_cells_0&Myeloid_0&T_cells_0'=1,'Adipocytes_0&Malignant_cells_0&Myeloid_0&Normal_epithelial_0'=3,'Adipocytes_0&Malignant_cells_0&Myeloid_0&Normal_epithelial_0&T_cells_0'=0,'Adipocytes_0&Fibroblasts_0'=393,'Adipocytes_0&Fibroblasts_0&T_cells_0'=19,'Adipocytes_0&Fibroblasts_0&Normal_epithelial_0'=25,'Adipocytes_0&Fibroblasts_0&Normal_epithelial_0&T_cells_0'=0,'Adipocytes_0&Fibroblasts_0&Myeloid_0'=68,'Adipocytes_0&Fibroblasts_0&Myeloid_0&T_cells_0'=11,'Adipocytes_0&Fibroblasts_0&Myeloid_0&Normal_epithelial_0'=5,'Adipocytes_0&Fibroblasts_0&Myeloid_0&Normal_epithelial_0&T_cells_0'=1,'Adipocytes_0&Fibroblasts_0&Malignant_cells_0'=108,'Adipocytes_0&Fibroblasts_0&Malignant_cells_0&T_cells_0'=9,'Adipocytes_0&Fibroblasts_0&Malignant_cells_0&Normal_epithelial_0'=19,'Adipocytes_0&Fibroblasts_0&Malignant_cells_0&Normal_epithelial_0&T_cells_0'=0,'Adipocytes_0&Fibroblasts_0&Malignant_cells_0&Myeloid_0'=32,'Adipocytes_0&Fibroblasts_0&Malignant_cells_0&Myeloid_0&T_cells_0'=3,'Adipocytes_0&Fibroblasts_0&Malignant_cells_0&Myeloid_0&Normal_epithelial_0'=2,'Adipocytes_0&Fibroblasts_0&Malignant_cells_0&Myeloid_0&Normal_epithelial_0&T_cells_0'=0,'Adipocytes_0&Endothelial_cells_0'=247,'Adipocytes_0&Endothelial_cells_0&T_cells_0'=16,'Adipocytes_0&Endothelial_cells_0&Normal_epithelial_0'=9,'Adipocytes_0&Endothelial_cells_0&Normal_epithelial_0&T_cells_0'=0,'Adipocytes_0&Endothelial_cells_0&Myeloid_0'=59,'Adipocytes_0&Endothelial_cells_0&Myeloid_0&T_cells_0'=10,'Adipocytes_0&Endothelial_cells_0&Myeloid_0&Normal_epithelial_0'=0,'Adipocytes_0&Endothelial_cells_0&Myeloid_0&Normal_epithelial_0&T_cells_0'=0,'Adipocytes_0&Endothelial_cells_0&Malignant_cells_0'=62,'Adipocytes_0&Endothelial_cells_0&Malignant_cells_0&T_cells_0'=2,'Adipocytes_0&Endothelial_cells_0&Malignant_cells_0&Normal_epithelial_0'=9,'Adipocytes_0&Endothelial_cells_0&Malignant_cells_0&Normal_epithelial_0&T_cells_0'=1,'Adipocytes_0&Endothelial_cells_0&Malignant_cells_0&Myeloid_0'=8,'Adipocytes_0&Endothelial_cells_0&Malignant_cells_0&Myeloid_0&T_cells_0'=0,'Adipocytes_0&Endothelial_cells_0&Malignant_cells_0&Myeloid_0&Normal_epithelial_0'=0,'Adipocytes_0&Endothelial_cells_0&Malignant_cells_0&Myeloid_0&Normal_epithelial_0&T_cells_0'=0,'Adipocytes_0&Endothelial_cells_0&Fibroblasts_0'=402,'Adipocytes_0&Endothelial_cells_0&Fibroblasts_0&T_cells_0'=14,'Adipocytes_0&Endothelial_cells_0&Fibroblasts_0&Normal_epithelial_0'=26,'Adipocytes_0&Endothelial_cells_0&Fibroblasts_0&Normal_epithelial_0&T_cells_0'=0,'Adipocytes_0&Endothelial_cells_0&Fibroblasts_0&Myeloid_0'=102,'Adipocytes_0&Endothelial_cells_0&Fibroblasts_0&Myeloid_0&T_cells_0'=10,'Adipocytes_0&Endothelial_cells_0&Fibroblasts_0&Myeloid_0&Normal_epithelial_0'=7,'Adipocytes_0&Endothelial_cells_0&Fibroblasts_0&Myeloid_0&Normal_epithelial_0&T_cells_0'=0,'Adipocytes_0&Endothelial_cells_0&Fibroblasts_0&Malignant_cells_0'=75,'Adipocytes_0&Endothelial_cells_0&Fibroblasts_0&Malignant_cells_0&T_cells_0'=5,'Adipocytes_0&Endothelial_cells_0&Fibroblasts_0&Malignant_cells_0&Normal_epithelial_0'=10,'Adipocytes_0&Endothelial_cells_0&Fibroblasts_0&Malignant_cells_0&Normal_epithelial_0&T_cells_0'=0,'Adipocytes_0&Endothelial_cells_0&Fibroblasts_0&Malignant_cells_0&Myeloid_0'=16,'Adipocytes_0&Endothelial_cells_0&Fibroblasts_0&Malignant_cells_0&Myeloid_0&T_cells_0'=1,'Adipocytes_0&Endothelial_cells_0&Fibroblasts_0&Malignant_cells_0&Myeloid_0&Normal_epithelial_0'=2,'Adipocytes_0&Endothelial_cells_0&Fibroblasts_0&Malignant_cells_0&Myeloid_0&Normal_epithelial_0&T_cells_0'=0,'Adipocytes_0&BP_cells_0'=8,'Adipocytes_0&BP_cells_0&T_cells_0'=10,'Adipocytes_0&BP_cells_0&Normal_epithelial_0'=1,'Adipocytes_0&BP_cells_0&Normal_epithelial_0&T_cells_0'=0,'Adipocytes_0&BP_cells_0&Myeloid_0'=9,'Adipocytes_0&BP_cells_0&Myeloid_0&T_cells_0'=6,'Adipocytes_0&BP_cells_0&Myeloid_0&Normal_epithelial_0'=1,'Adipocytes_0&BP_cells_0&Myeloid_0&Normal_epithelial_0&T_cells_0'=1,'Adipocytes_0&BP_cells_0&Malignant_cells_0'=9,'Adipocytes_0&BP_cells_0&Malignant_cells_0&T_cells_0'=4,'Adipocytes_0&BP_cells_0&Malignant_cells_0&Normal_epithelial_0'=2,'Adipocytes_0&BP_cells_0&Malignant_cells_0&Normal_epithelial_0&T_cells_0'=0,'Adipocytes_0&BP_cells_0&Malignant_cells_0&Myeloid_0'=6,'Adipocytes_0&BP_cells_0&Malignant_cells_0&Myeloid_0&T_cells_0'=4,'Adipocytes_0&BP_cells_0&Malignant_cells_0&Myeloid_0&Normal_epithelial_0'=0,'Adipocytes_0&BP_cells_0&Malignant_cells_0&Myeloid_0&Normal_epithelial_0&T_cells_0'=0,'Adipocytes_0&BP_cells_0&Fibroblasts_0'=21,'Adipocytes_0&BP_cells_0&Fibroblasts_0&T_cells_0'=5,'Adipocytes_0&BP_cells_0&Fibroblasts_0&Normal_epithelial_0'=0,'Adipocytes_0&BP_cells_0&Fibroblasts_0&Normal_epithelial_0&T_cells_0'=0,'Adipocytes_0&BP_cells_0&Fibroblasts_0&Myeloid_0'=10,'Adipocytes_0&BP_cells_0&Fibroblasts_0&Myeloid_0&T_cells_0'=6,'Adipocytes_0&BP_cells_0&Fibroblasts_0&Myeloid_0&Normal_epithelial_0'=0,'Adipocytes_0&BP_cells_0&Fibroblasts_0&Myeloid_0&Normal_epithelial_0&T_cells_0'=0,'Adipocytes_0&BP_cells_0&Fibroblasts_0&Malignant_cells_0'=6,'Adipocytes_0&BP_cells_0&Fibroblasts_0&Malignant_cells_0&T_cells_0'=8,'Adipocytes_0&BP_cells_0&Fibroblasts_0&Malignant_cells_0&Normal_epithelial_0'=1,'Adipocytes_0&BP_cells_0&Fibroblasts_0&Malignant_cells_0&Normal_epithelial_0&T_cells_0'=0,'Adipocytes_0&BP_cells_0&Fibroblasts_0&Malignant_cells_0&Myeloid_0'=4,'Adipocytes_0&BP_cells_0&Fibroblasts_0&Malignant_cells_0&Myeloid_0&T_cells_0'=2,'Adipocytes_0&BP_cells_0&Fibroblasts_0&Malignant_cells_0&Myeloid_0&Normal_epithelial_0'=1,'Adipocytes_0&BP_cells_0&Fibroblasts_0&Malignant_cells_0&Myeloid_0&Normal_epithelial_0&T_cells_0'=0,'Adipocytes_0&BP_cells_0&Endothelial_cells_0'=8,'Adipocytes_0&BP_cells_0&Endothelial_cells_0&T_cells_0'=9,'Adipocytes_0&BP_cells_0&Endothelial_cells_0&Normal_epithelial_0'=0,'Adipocytes_0&BP_cells_0&Endothelial_cells_0&Normal_epithelial_0&T_cells_0'=1,'Adipocytes_0&BP_cells_0&Endothelial_cells_0&Myeloid_0'=10,'Adipocytes_0&BP_cells_0&Endothelial_cells_0&Myeloid_0&T_cells_0'=10,'Adipocytes_0&BP_cells_0&Endothelial_cells_0&Myeloid_0&Normal_epithelial_0'=0,'Adipocytes_0&BP_cells_0&Endothelial_cells_0&Myeloid_0&Normal_epithelial_0&T_cells_0'=0,'Adipocytes_0&BP_cells_0&Endothelial_cells_0&Malignant_cells_0'=4,'Adipocytes_0&BP_cells_0&Endothelial_cells_0&Malignant_cells_0&T_cells_0'=1,'Adipocytes_0&BP_cells_0&Endothelial_cells_0&Malignant_cells_0&Normal_epithelial_0'=0,'Adipocytes_0&BP_cells_0&Endothelial_cells_0&Malignant_cells_0&Normal_epithelial_0&T_cells_0'=0,'Adipocytes_0&BP_cells_0&Endothelial_cells_0&Malignant_cells_0&Myeloid_0'=2,'Adipocytes_0&BP_cells_0&Endothelial_cells_0&Malignant_cells_0&Myeloid_0&T_cells_0'=2,'Adipocytes_0&BP_cells_0&Endothelial_cells_0&Malignant_cells_0&Myeloid_0&Normal_epithelial_0'=0,'Adipocytes_0&BP_cells_0&Endothelial_cells_0&Malignant_cells_0&Myeloid_0&Normal_epithelial_0&T_cells_0'=0,'Adipocytes_0&BP_cells_0&Endothelial_cells_0&Fibroblasts_0'=28,'Adipocytes_0&BP_cells_0&Endothelial_cells_0&Fibroblasts_0&T_cells_0'=12,'Adipocytes_0&BP_cells_0&Endothelial_cells_0&Fibroblasts_0&Normal_epithelial_0'=1,'Adipocytes_0&BP_cells_0&Endothelial_cells_0&Fibroblasts_0&Normal_epithelial_0&T_cells_0'=0,'Adipocytes_0&BP_cells_0&Endothelial_cells_0&Fibroblasts_0&Myeloid_0'=5,'Adipocytes_0&BP_cells_0&Endothelial_cells_0&Fibroblasts_0&Myeloid_0&T_cells_0'=8,'Adipocytes_0&BP_cells_0&Endothelial_cells_0&Fibroblasts_0&Myeloid_0&Normal_epithelial_0'=0,'Adipocytes_0&BP_cells_0&Endothelial_cells_0&Fibroblasts_0&Myeloid_0&Normal_epithelial_0&T_cells_0'=1,'Adipocytes_0&BP_cells_0&Endothelial_cells_0&Fibroblasts_0&Malignant_cells_0'=7,'Adipocytes_0&BP_cells_0&Endothelial_cells_0&Fibroblasts_0&Malignant_cells_0&T_cells_0'=8,'Adipocytes_0&BP_cells_0&Endothelial_cells_0&Fibroblasts_0&Malignant_cells_0&Normal_epithelial_0'=0,'Adipocytes_0&BP_cells_0&Endothelial_cells_0&Fibroblasts_0&Malignant_cells_0&Normal_epithelial_0&T_cells_0'=0,'Adipocytes_0&BP_cells_0&Endothelial_cells_0&Fibroblasts_0&Malignant_cells_0&Myeloid_0'=1,'Adipocytes_0&BP_cells_0&Endothelial_cells_0&Fibroblasts_0&Malignant_cells_0&Myeloid_0&T_cells_0'=1,'Adipocytes_0&BP_cells_0&Endothelial_cells_0&Fibroblasts_0&Malignant_cells_0&Myeloid_0&Normal_epithelial_0'=0,'Adipocytes_0&BP_cells_0&Endothelial_cells_0&Fibroblasts_0&Malignant_cells_0&Myeloid_0&Normal_epithelial_0&T_cells_0'=0)
upset(fromExpression(expressionInput), nsets=8, nintersects=30, show.numbers="yes", main.bar.color="#ea5d4e", sets.bar.color="#317eab", empty.intersections=NULL, order.by = "freq", number.angles = 0, mainbar.y.label ="No. of Intersections", sets.x.label ="Set size")
invisible(dev.off())

########################################
###### EXTENDED DATA Figure 10 #########
########################################
###############################################################################################################################
############# EXTENDED DATA Figure 10 A-D - SWARMPLOT - Baseline - Crossover  - Sensitive - Resistant   ########################
###############################################################################################################################

library(EnhancedVolcano)
library(ggrepel)
library(ggplot2)
library(miloR)
library(plyr)
library(dplyr)
library(future)
library(stringi)
library(tidyverse)

### available from PP_snakemake_scafe
## files for resistant clones ExDat fig 10 A 
#d2 <- read.table("clones_res_A_vs_B__distal_cre.tsv", sep = '\t', header = TRUE, fill = TRUE)
#d1 <- read.table("clones_res_B_vs_A__distal_cre.tsv", sep = '\t', header = TRUE, fill = TRUE)

## files for sensitive clones ExDat fig 10 C 
d2 <- read.table("clones_sens_A_vs_B__distal_cre.tsv", sep = '\t', header = TRUE, fill = TRUE)
d1 <- read.table("clones_sens_B_vs_A__distal_cre.tsv", sep = '\t', header = TRUE, fill = TRUE)

d2$oddsRatio <- (-d2$oddsRatio)

  d <- rbind(d1[,colnames(d1) %in% c('collection', 'pValueLog', 'oddsRatio')], 
             d2[,colnames(d2) %in% c('collection', 'pValueLog', 'oddsRatio')])
  
  max <- max(d$oddsRatio[!d$oddsRatio %in% c('Inf')])
  min <- min(d$oddsRatio[!d$oddsRatio %in% c('-Inf')])
  
  d$pValue <- 10^(-d$pValueLog)

filtered_data <- subset(d, (pValue < 0.05 & oddsRatio > 2) | pValue < 0.05 & oddsRatio < -2 )
collection_count <- table(filtered_data$collection)


print(head (colnames(d)))
collection_count_plus <- as.data.frame(table (d$collection[d$oddsRatio > 2 & d$pValue < 0.05])) ; colnames (collection_count_plus)[2] <- 'plus'
collection_count_minus <- as.data.frame(table (d$collection[!d$oddsRatio > -2 & d$pValue < 0.05])) ;colnames (collection_count_minus)[2] <- 'minus'

selected_collections <- names(collection_count[collection_count > 5])

selected_collections_f <- collection_count_plus %>% filter(plus >= 10) %>% arrange(desc(plus)) %>% head(5)                     
selected_collections_plus <- selected_collections_f$Var1
selected_collections_f <- collection_count_minus %>% filter(minus >= 10) %>% arrange(desc(minus)) %>% head(5)                     
selected_collections_minus <- selected_collections_f$Var1
selected_collections <- c(selected_collections_plus, selected_collections_minus)
print(selected_collections)

d1 <-  d[d$collection%in%selected_collections,]

colnames(d1)[4] <- 'SpatialFDR'
colnames(d1)[3] <- 'logFC'
d1$Nhood <- rownames(d1)

vector_r <- rev(data.frame(d1 %>% group_by(collection) %>% summarize(avg_logFC = mean(logFC)) %>% arrange(avg_logFC)))
d1$collection <- factor(d1$collection, levels = rev(vector_r$collection))


collection_count_tot <- as.data.frame(table (d$collection)/2)
  d <- d[!d$oddsRatio %in% c("Inf", "-Inf"), ]

print(colnames(collection_count_tot))

mer<- merge(collection_count_tot, collection_count_plus, by.x='Var1', by.y='Var1', all.x=T)
mer1<- merge(mer, collection_count_minus, by.x='Var1', by.y='Var1', all.x=T)
mer1$plus <- ifelse (is.na(mer1$plus), 0, mer1$plus)
mer1$minus <- ifelse (is.na(mer1$minus), 0, mer1$minus)


mer1_filtered <- mer1[mer1$Var1 %in% vector_r$collection, ]
mer1_filtered$Var1 <- factor(mer1_filtered$Var1, levels = vector_r$collection)
mer1_filtered <- mer1_filtered[order(mer1_filtered$Var1), ] 

  plot_title <- ''
  
 p <- plotDAbeeswarm(d1, group.by = "collection", alpha = 0.05) +
    ggtitle(plot_title) +
    scale_colour_gradient2(
      high = "#9D2E99",  
      low = "#AA3830",    
      na.value = "gray88", 
      space = "Lab"
    ) + 
    xlab("") +
    theme(
	  plot.title = element_text(size = 10, face = "bold"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.y = element_text(size = 8)
    )
  
  p$data$logFC_color <- ifelse(p$data$logFC_color >= 2, 2, p$data$logFC_color)
  p$data$logFC_color <- ifelse(p$data$logFC_color < -2, -2, p$data$logFC_color)
  p$data$logFC_color <- ifelse((p$data$logFC_color == -2 | p$data$logFC_color == 2), p$data$logFC_color, NA)
  p$data <- p$data[order(!is.na(p$data$logFC_color), p$data$logFC_color), ]
  
  p <- p + mapply(function(x, plus_label, minus_label) { 
    list(
      annotate("text", x = x, y = max(d1$logFC) + 5, label = plus_label, hjust = 1, size = 3),
      annotate("text", x = x, y = min(d1$logFC) - 5, label = minus_label, hjust = 1, size = 3)
    ) 
  }, seq_along(mer1_filtered$Var1), paste ("(", rev(mer1_filtered$plus), ")"), paste ("(",rev(mer1_filtered$minus), ")"))



p1 <- plotDAbeeswarm(d1, group.by = "collection", alpha = 0.05) +
  ggtitle('Plot with max logFC') +
  scale_colour_gradient2(
      high = "#9D2E99",  
      low = "#AA3830",    
      na.value = "gray88", 
    space = "Lab"
  ) + 
  ylim(0, max + 3) + 
  xlab("") +
  theme(
    plot.title = element_text(size = 10, face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 20),
    axis.text.x = element_text(size = 20)
  )

  p1$data$logFC_color <- ifelse(p1$data$logFC_color >= 2, 2, p1$data$logFC_color)
  p1$data$logFC_color <- ifelse(p1$data$logFC_color < -2, -2, p1$data$logFC_color)
  p1$data$logFC_color <- ifelse((p1$data$logFC_color == -2 | p1$data$logFC_color == 2), p1$data$logFC_color, NA)
  p1$data <- p1$data[order(!is.na(p1$data$logFC_color), p1$data$logFC_color), ]

p1 <- p1 + mapply(function(x, plus_label, minus_label) { 
    list(
      annotate("text", x = x, y = max + 1, label = plus_label, hjust = 1, size = 3),
      annotate("text", x = x, y = min - 1, label = minus_label, hjust = 1, size = 3)
    ) 
  }, seq_along(mer1_filtered$Var1), paste ("(", rev(mer1_filtered$plus), ")"), paste ("(",rev(mer1_filtered$minus), ")"))

  
p2 <- plotDAbeeswarm(d1, group.by = "collection", alpha = 0.05) +
  ggtitle('Plot with min logFC') +
  scale_colour_gradient2(
      high = "#9D2E99",  
      low = "#AA3830",    
      na.value = "gray88", 
      space = "Lab"
  ) + 
  ylim(0, min - 6) + 
  xlab("") +
  theme(
    plot.title = element_text(size = 10, face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 20),
    axis.text.x = element_text(size = 20)
  )

  p2$data$logFC_color <- ifelse(p2$data$logFC_color >= 2, 2, p2$data$logFC_color)
  p2$data$logFC_color <- ifelse(p2$data$logFC_color < -2, -2, p2$data$logFC_color)
  p2$data$logFC_color <- ifelse((p2$data$logFC_color == -2 | p2$data$logFC_color == 2), p2$data$logFC_color, NA)
  p2$data <- p2$data[order(!is.na(p2$data$logFC_color), p2$data$logFC_color), ]

p2 <- p2 + mapply(function(x, plus_label, minus_label) { 
    list(
      annotate("text", x = x, y = max(d1$logFC) + 1, label = plus_label, hjust = 1, size = 3),
      annotate("text", x = x, y = min(d1$logFC) - 1, label = minus_label, hjust = 1, size = 3)
    ) 
  }, seq_along(mer1_filtered$Var1), paste ("(", rev(mer1_filtered$plus), ")"), paste ("(",rev(mer1_filtered$minus), ")"))

#pdf ("ExDaFigure10A swarmplot resistant_clones A_v_B.pdf")
pdf ("ExDaFigure10C swarmplot sensitive_clones A_v_B.pdf")
p
p1
p2  
dev.off()


###############################################################################################################################
############# EXTENDED DATA Figure 10 B-D - SWARMPLOT - Baseline - EOT  - Sensitive - Resistant   ########################
###############################################################################################################################

library(EnhancedVolcano)
library(ggrepel)
library(ggplot2)
library(miloR)
library(plyr)
library(dplyr)
library(future)
library(stringi)
library(tidyverse)

### available from PP_snakemake_scafe
## files for resistant clones ExDa fig 10 B 
#d2 <- read.table("clones_res_A_vs_C__distal_cre.tsv", sep = '\t', header = TRUE, fill = TRUE)
#d1 <- read.table("clones_res_C_vs_A__distal_cre.tsv", sep = '\t', header = TRUE, fill = TRUE)
## files for resistant clones ExDa fig 10 D 
d2 <- read.table("clones_sens_A_vs_C__distal_cre.tsv", sep = '\t', header = TRUE, fill = TRUE)
d1 <- read.table("clones_sens_C_vs_A__distal_cre.tsv", sep = '\t', header = TRUE, fill = TRUE)

d2$oddsRatio <- (-d2$oddsRatio)

  d <- rbind(d1[,colnames(d1) %in% c('collection', 'pValueLog', 'oddsRatio')], 
             d2[,colnames(d2) %in% c('collection', 'pValueLog', 'oddsRatio')])
  
  max <- max(d$oddsRatio[!d$oddsRatio %in% c('Inf')])
  min <- min(d$oddsRatio[!d$oddsRatio %in% c('-Inf')])
  
  d$pValue <- 10^(-d$pValueLog)

filtered_data <- subset(d, (pValue < 0.05 & oddsRatio > 2) | pValue < 0.05 & oddsRatio < -2 )
collection_count <- table(filtered_data$collection)


print(head (colnames(d)))
collection_count_plus <- as.data.frame(table (d$collection[d$oddsRatio > 2 & d$pValue < 0.05])) ; colnames (collection_count_plus)[2] <- 'plus'
collection_count_minus <- as.data.frame(table (d$collection[!d$oddsRatio > -2 & d$pValue < 0.05])) ;colnames (collection_count_minus)[2] <- 'minus'

selected_collections <- names(collection_count[collection_count > 5])

selected_collections_f <- collection_count_plus %>% filter(plus >= 10) %>% arrange(desc(plus)) %>% head(5)                     
selected_collections_plus <- selected_collections_f$Var1
selected_collections_f <- collection_count_minus %>% filter(minus >= 10) %>% arrange(desc(minus)) %>% head(5)                     
selected_collections_minus <- selected_collections_f$Var1
selected_collections <- c(selected_collections_plus, selected_collections_minus)
print(selected_collections)

d1 <-  d[d$collection%in%selected_collections,]

colnames(d1)[4] <- 'SpatialFDR'
colnames(d1)[3] <- 'logFC'
d1$Nhood <- rownames(d1)

vector_r <- rev(data.frame(d1 %>% group_by(collection) %>% summarize(avg_logFC = mean(logFC)) %>% arrange(avg_logFC)))
d1$collection <- factor(d1$collection, levels = rev(vector_r$collection))


collection_count_tot <- as.data.frame(table (d$collection)/2)
  d <- d[!d$oddsRatio %in% c("Inf", "-Inf"), ]

print(colnames(collection_count_tot))

mer<- merge(collection_count_tot, collection_count_plus, by.x='Var1', by.y='Var1', all.x=T)
mer1<- merge(mer, collection_count_minus, by.x='Var1', by.y='Var1', all.x=T)
mer1$plus <- ifelse (is.na(mer1$plus), 0, mer1$plus)
mer1$minus <- ifelse (is.na(mer1$minus), 0, mer1$minus)


mer1_filtered <- mer1[mer1$Var1 %in% vector_r$collection, ]
mer1_filtered$Var1 <- factor(mer1_filtered$Var1, levels = vector_r$collection)
mer1_filtered <- mer1_filtered[order(mer1_filtered$Var1), ] 

  plot_title <- ''
  
 p <- plotDAbeeswarm(d1, group.by = "collection", alpha = 0.05) +
    ggtitle(plot_title) +
    scale_colour_gradient2(
      high = "#350D7F",  
      low = "#AA3830",    
      na.value = "gray88", 
      space = "Lab"
    ) + 
    xlab("") +
    theme(
	  plot.title = element_text(size = 10, face = "bold"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.y = element_text(size = 8)
    )
  
  p$data$logFC_color <- ifelse(p$data$logFC_color >= 2, 2, p$data$logFC_color)
  p$data$logFC_color <- ifelse(p$data$logFC_color < -2, -2, p$data$logFC_color)
  p$data$logFC_color <- ifelse((p$data$logFC_color == -2 | p$data$logFC_color == 2), p$data$logFC_color, NA)
  p$data <- p$data[order(!is.na(p$data$logFC_color), p$data$logFC_color), ]
  
  p <- p + mapply(function(x, plus_label, minus_label) { 
    list(
      annotate("text", x = x, y = max(d1$logFC) + 5, label = plus_label, hjust = 1, size = 3),
      annotate("text", x = x, y = min(d1$logFC) - 5, label = minus_label, hjust = 1, size = 3)
    ) 
  }, seq_along(mer1_filtered$Var1), paste ("(", rev(mer1_filtered$plus), ")"), paste ("(",rev(mer1_filtered$minus), ")"))



p1 <- plotDAbeeswarm(d1, group.by = "collection", alpha = 0.05) +
  ggtitle('Plot with max logFC') +
  scale_colour_gradient2(
      high = "#350D7F",  
      low = "#AA3830",    
      na.value = "gray88", 
      space = "Lab"
  ) + 
  ylim(0, max + 3) + 
  xlab("") +
  theme(
    plot.title = element_text(size = 10, face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 20),
    axis.text.x = element_text(size = 20)
  )

  p1$data$logFC_color <- ifelse(p1$data$logFC_color >= 2, 2, p1$data$logFC_color)
  p1$data$logFC_color <- ifelse(p1$data$logFC_color < -2, -2, p1$data$logFC_color)
  p1$data$logFC_color <- ifelse((p1$data$logFC_color == -2 | p1$data$logFC_color == 2), p1$data$logFC_color, NA)
  p1$data <- p1$data[order(!is.na(p1$data$logFC_color), p1$data$logFC_color), ]

p1 <- p1 + mapply(function(x, plus_label, minus_label) { 
    list(
      annotate("text", x = x, y = max + 1, label = plus_label, hjust = 1, size = 3),
      annotate("text", x = x, y = min - 1, label = minus_label, hjust = 1, size = 3)
    ) 
  }, seq_along(mer1_filtered$Var1), paste ("(", rev(mer1_filtered$plus), ")"), paste ("(",rev(mer1_filtered$minus), ")"))

  
p2 <- plotDAbeeswarm(d1, group.by = "collection", alpha = 0.05) +
  ggtitle('Plot with min logFC') +
  scale_colour_gradient2(
      high = "#350D7F",  
      low = "#AA3830",    
      na.value = "gray88", 
      space = "Lab"
  ) + 
  ylim(0, min - 6) + 
  xlab("") +
  theme(
    plot.title = element_text(size = 10, face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 20),
    axis.text.x = element_text(size = 20)
  )

  p2$data$logFC_color <- ifelse(p2$data$logFC_color >= 2, 2, p2$data$logFC_color)
  p2$data$logFC_color <- ifelse(p2$data$logFC_color < -2, -2, p2$data$logFC_color)
  p2$data$logFC_color <- ifelse((p2$data$logFC_color == -2 | p2$data$logFC_color == 2), p2$data$logFC_color, NA)
  p2$data <- p2$data[order(!is.na(p2$data$logFC_color), p2$data$logFC_color), ]

p2 <- p2 + mapply(function(x, plus_label, minus_label) { 
    list(
      annotate("text", x = x, y = max(d1$logFC) + 1, label = plus_label, hjust = 1, size = 3),
      annotate("text", x = x, y = min(d1$logFC) - 1, label = minus_label, hjust = 1, size = 3)
    ) 
  }, seq_along(mer1_filtered$Var1), paste ("(", rev(mer1_filtered$plus), ")"), paste ("(",rev(mer1_filtered$minus), ")"))

#pdf ("ExDaFigure10B swarmplot resistant_clones A_v_C.pdf")
pdf ("ExDaFigure10D swarmplot sensitive_clones A_v_C.pdf.pdf")
p
p1
p2  
dev.off()
