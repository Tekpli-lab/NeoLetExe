###################################################
####### FIGURE 3 ##################################
###################################################

################################################################
###### Figure 3A ########## Milo minor IDs SWARMPLOTS ##########
################################################################

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
## minor.id plot A vs C #######  
###############################   

A= "#DA3A09"
B= "#9C2C98"
C= "#320E77"

pdf('Figure 3Aa swarmplot minor ids A C.pdf', width = 7, height = 10)
da_results_B_4_rev <- as.data.frame(data.frame(read.table("minorID_AC.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)))

vector_r <- rev(data.frame(da_results_B_4_rev %>% group_by(minor.ids) %>% summarize(avg_logFC = mean(logFC)) %>% arrange(avg_logFC)))
da_results_B_4_rev$minor.ids <- factor(da_results_B_4_rev$minor.ids, levels = rev(vector_r$minor.ids))
da_results_B_4_rev <- da_results_B_4_rev[order(match(da_results_B_4_rev$minor.ids, vector_r$minor.ids)),]

p <- plotDAbeeswarm(da_results_B_4_rev, group.by = "minor.ids") + 
scale_colour_gradient2(high = C, low = A, na.value = "grey85", space = "Lab") +  
  xlab("") + 
  scale_y_continuous(breaks = c(-5, -2.5, 0, 2.5, 5), limits = c(-6, 6.2)) +
  theme(axis.text.y = element_text(colour = rev(c('orange', '#8ba831', 'orange', 'orange', '#8ba831', '#8ba831', '#8ba831', '#8ba831', 'orange', '#4f4fc5', '#4f4fc5', '#8ba831', '#4f4fc5', '#4f4fc5', 'orange', '#4f4fc5', '#4f4fc5', '#4f4fc5', 'orange'))),
axis.text.x = element_text(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank())

da_results_B_4_sig_1 <- da_results_B_4_rev[da_results_B_4_rev$FDR <0.05 & da_results_B_4_rev$logFC <0,]
da_results_B_4_sig_2 <- da_results_B_4_rev[da_results_B_4_rev$FDR <0.05 & da_results_B_4_rev$logFC >0,]
count_sig_1 <- as.data.frame(table(da_results_B_4_sig_1$minor.ids))
colnames(count_sig_1) <- c("minor.ids", "count_sig_1")
count_sig_2 <- as.data.frame(table(da_results_B_4_sig_2$minor.ids))
colnames(count_sig_2) <- c("minor.ids", "count_sig_2")
count_data <- merge(count_sig_1, count_sig_2, by = "minor.ids", all = TRUE)
count_data[is.na(count_data)] <- 0

  p$data <- p$data[order(!is.na(p$data$logFC_color), p$data$logFC_color), ]

p <- p + mapply(function(x, plus_label, minus_label) { 
    list(annotate("text", x = x, y = 6, label = paste0("(", plus_label, ")"), hjust = 1, size = 3),
      annotate("text", x = x, y = -5.5, label = paste0("(", minus_label, ")"), hjust = 1, size = 3))}, 
  seq_along(count_data$minor.ids), 
  count_data$count_sig_2[rev(match(vector_r$minor.ids, count_data$minor.ids))], 
  count_data$count_sig_1[rev(match(vector_r$minor.ids, count_data$minor.ids))], 
  SIMPLIFY = FALSE)
print(p)
dev.off()

########################### 
## minor.id plot A vs B ### 
########################### 
pdf('Figure 3Ab swarmplot minor ids A B.pdf', width = 7, height = 10)
da_results_B_4_rev <- as.data.frame(data.frame(read.table("minorID_AB.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)))
da_results_B_4_rev$minor.ids <- factor(da_results_B_4_rev$minor.ids, levels = rev(vector_r$minor.ids))
da_results_B_4_rev <- da_results_B_4_rev[order(match(da_results_B_4_rev$minor.ids, vector_r$minor.ids)),]

p <- plotDAbeeswarm(da_results_B_4_rev, group.by = "minor.ids") + 
scale_colour_gradient2(high = B, low = A, na.value = "grey85", space = "Lab") +  
  xlab("") + 
  scale_y_continuous(breaks = c(-5, -2.5, 0, 2.5, 5), limits = c(-6, 6.2)) +
  theme(axis.text.y = element_text(colour = rev(c('orange', '#8ba831', 'orange', 'orange', '#8ba831', '#8ba831', '#8ba831', '#8ba831', 'orange', '#4f4fc5', '#4f4fc5', '#8ba831', '#4f4fc5', '#4f4fc5', 'orange', '#4f4fc5', '#4f4fc5', '#4f4fc5', 'orange'))),
axis.text.x = element_text(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank())

da_results_B_4_sig_1 <- da_results_B_4_rev[da_results_B_4_rev$FDR <0.05 & da_results_B_4_rev$logFC <0,]
da_results_B_4_sig_2 <- da_results_B_4_rev[da_results_B_4_rev$FDR <0.05 & da_results_B_4_rev$logFC >0,]
count_sig_1 <- as.data.frame(table(da_results_B_4_sig_1$minor.ids))
colnames(count_sig_1) <- c("minor.ids", "count_sig_1")
count_sig_2 <- as.data.frame(table(da_results_B_4_sig_2$minor.ids))
colnames(count_sig_2) <- c("minor.ids", "count_sig_2")
count_data <- merge(count_sig_1, count_sig_2, by = "minor.ids", all = TRUE)
count_data[is.na(count_data)] <- 0

  p$data <- p$data[order(!is.na(p$data$logFC_color), p$data$logFC_color), ]

p <- p + mapply(function(x, plus_label, minus_label) { 
    list(annotate("text", x = x, y = 6, label = paste0("(", plus_label, ")"), hjust = 1, size = 3),
      annotate("text", x = x, y = -5.5, label = paste0("(", minus_label, ")"), hjust = 1, size = 3))}, 
  seq_along(count_data$minor.ids), 
  count_data$count_sig_2[rev(match(vector_r$minor.ids, count_data$minor.ids))], 
  count_data$count_sig_1[rev(match(vector_r$minor.ids, count_data$minor.ids))], 
  SIMPLIFY = FALSE)
print(p)

dev.off()

########################### 
## minor.id plot B vs C ### 
###########################
pdf('Figure 3Ac swarmplot minor ids B C.pdf', width = 7, height = 10)
da_results_B_4_rev <- as.data.frame(data.frame(read.table("minorID_BC.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)))
da_results_B_4_rev$minor.ids <- factor(da_results_B_4_rev$minor.ids, levels = rev(vector_r$minor.ids))
da_results_B_4_rev <- da_results_B_4_rev[order(match(da_results_B_4_rev$minor.ids, vector_r$minor.ids)),]

p <- plotDAbeeswarm(da_results_B_4_rev, group.by = "minor.ids") + 
scale_colour_gradient2(high = C, low = B, na.value = "grey85", space = "Lab") +  
  xlab("") + 
  scale_y_continuous(breaks = c(-5, -2.5, 0, 2.5, 5), limits = c(-6, 6.2)) +
  theme(axis.text.y = element_text(colour = rev(c('orange', '#8ba831', 'orange', 'orange', '#8ba831', '#8ba831', '#8ba831', '#8ba831', 'orange', '#4f4fc5', '#4f4fc5', '#8ba831', '#4f4fc5', '#4f4fc5', 'orange', '#4f4fc5', '#4f4fc5', '#4f4fc5', 'orange'))),
axis.text.x = element_text(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank())

da_results_B_4_sig_1 <- da_results_B_4_rev[da_results_B_4_rev$FDR <0.05 & da_results_B_4_rev$logFC <0,]
da_results_B_4_sig_2 <- da_results_B_4_rev[da_results_B_4_rev$FDR <0.05 & da_results_B_4_rev$logFC >0,]
count_sig_1 <- as.data.frame(table(da_results_B_4_sig_1$minor.ids))
colnames(count_sig_1) <- c("minor.ids", "count_sig_1")
count_sig_2 <- as.data.frame(table(da_results_B_4_sig_2$minor.ids))
colnames(count_sig_2) <- c("minor.ids", "count_sig_2")
count_data <- merge(count_sig_1, count_sig_2, by = "minor.ids", all = TRUE)
count_data[is.na(count_data)] <- 0

  p$data <- p$data[order(!is.na(p$data$logFC_color), p$data$logFC_color), ]

p <- p + mapply(function(x, plus_label, minus_label) { 
    list(annotate("text", x = x, y = 6, label = paste0("(", plus_label, ")"), hjust = 1, size = 3),
      annotate("text", x = x, y = -5.5, label = paste0("(", minus_label, ")"), hjust = 1, size = 3))}, 
  seq_along(count_data$minor.ids), 
  count_data$count_sig_2[rev(match(vector_r$minor.ids, count_data$minor.ids))], 
  count_data$count_sig_1[rev(match(vector_r$minor.ids, count_data$minor.ids))], 
  SIMPLIFY = FALSE)
print(p)
dev.off()


######################################################################
###### Figure 3B ########## BOXPLOT CD4 MINOR IDS in time   ##########
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

metadata <- metadata[!metadata$main.ids=="Epithelial_cells",]
metadata <- metadata[!metadata$minor.ids=="UNA",]
dta <- as.data.frame.matrix(prop.table(table (metadata$orig.ident, metadata$minor.ids),1)*100)
tabix <- as.data.frame.matrix(table (metadata$orig.ident, metadata$minor.ids))

dta$category <- stri_sub(rownames (dta),-1)
dta$paired <- substr (rownames (dta),2,4)
dta$paired <- gsub ("A", "",dta$paired )
dta$paired <- gsub ("B", "",dta$paired )
dta$paired <- gsub ("C", "",dta$paired )
annot <- metadata[, colnames (metadata)%in% c("patient_number", "orig.ident")]
annot <- annot[!duplicated (metadata$orig.ident), ]
table (annot$orig.ident== rownames (dta))
dta <- cbind (dta , annot)

dta <- dta[!rownames(dta)=='P04A',]
probe <- 'CD4'

pAB <- wilcox.test(dta[dta$category=="A", colnames (dta)%in% probe], dta[dta$category=="B", colnames (dta)%in% probe], paired = TRUE)[3]
pBC <- wilcox.test(dta[dta$category=="B", colnames (dta)%in% probe], dta[dta$category=="C", colnames (dta)%in% probe], paired = TRUE)[3]
pAC <- wilcox.test(dta[dta$category=="A", colnames (dta)%in% probe], dta[dta$category=="C", colnames (dta)%in% probe], paired = TRUE)[3]
pvalues <- c(pAB, pBC, pAC)

dta$category <- factor(dta$category, levels = c("A", "B", "C"))

outline_colors <- c("#AA3830", '#9D2E99', '#350D7F') 
box_fill_colors <- c('#FFE2E0', '#FFD1FF', '#DAC3FF')
point_fill_colors <- c('#FF8077', '#E66CE1', '#7E33FF') 

pdf ('Figure 3B CD4_boxplot.pdf')
a <- ggplot(dta, aes(x = category, y = dta[, colnames(dta) %in% probe])) + 
  geom_boxplot(color = outline_colors, 
               fill = box_fill_colors, 
               alpha = 0.5, size = 1.5) + 
  geom_point(aes(fill = factor(category), colour = factor(category)), 
             shape = 21, size = 5, stroke = 1.2) +  
  scale_fill_manual(values = point_fill_colors) +
  scale_colour_manual(values = outline_colors) +  
  ylab('Percentage of CD4 cells in the TME') + 
  geom_text_repel(aes(label = paired), min.segment.length = 1, max.overlaps = Inf) + 
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
    legend.position = "none"
  ) +
  labs(x = "") +

  geom_signif(comparisons = list(c("A", "B"), c("B", "C"), c("A", "C")),
              annotations = sprintf("p = %.3f", pvalues),
              y_position = 37,  
              map_signif_level = FALSE,
              test = "wilcox.test", 
              step_increase = 0.07)
print(a)
dev.off()



######################################################################
###### Figure 3C ########## BOXPLOT Malignant_cells  in time   #######
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

dta <- as.data.frame.matrix(prop.table(table (metadata$orig.ident, metadata$main.ids.cor),1)*100)
tabix <- as.data.frame.matrix(table (metadata$orig.ident, metadata$minor.ids))

dta$category <- stri_sub(rownames (dta),-1)
dta$paired <- substr (rownames (dta),2,4)
dta$paired <- gsub ("A", "",dta$paired )
dta$paired <- gsub ("B", "",dta$paired )
dta$paired <- gsub ("C", "",dta$paired )
annot <- metadata[, colnames (metadata)%in% c("patient_number", "orig.ident")]
annot <- annot[!duplicated (metadata$orig.ident), ]
table (annot$orig.ident== rownames (dta))
dta <- cbind (dta , annot)

dta <- dta[!rownames(dta)=='P04A',]
probe <- 'Malignant_cells'

pAB <- wilcox.test(dta[dta$category=="A", colnames (dta)%in% probe], dta[dta$category=="B", colnames (dta)%in% probe], paired = TRUE)[3]
pBC <- wilcox.test(dta[dta$category=="B", colnames (dta)%in% probe], dta[dta$category=="C", colnames (dta)%in% probe], paired = TRUE)[3]
pAC <- wilcox.test(dta[dta$category=="A", colnames (dta)%in% probe], dta[dta$category=="C", colnames (dta)%in% probe], paired = TRUE)[3]
pvalues <- c(pAB, pBC, pAC)

dta$category <- factor(dta$category, levels = c("A", "B", "C"))


outline_colors <- c("#AA3830", '#9D2E99', '#350D7F') 
box_fill_colors <- c('#FFE2E0', '#FFD1FF', '#DAC3FF') 
point_fill_colors <- c('#FF8077', '#E66CE1', '#7E33FF') 


pdf (' Figure 3C Malignant_ cells boxplot.pdf')

a <- ggplot(dta, aes(x = category, y = dta[, colnames(dta) %in% probe])) + 
  geom_boxplot(color = outline_colors, 
               fill = box_fill_colors, 
               alpha = 0.5, size = 1.5) + 
  geom_point(aes(fill = factor(category), colour = factor(category)), 
             shape = 21, size = 5, stroke = 1.2) +  
  scale_fill_manual(values = point_fill_colors) +
  scale_colour_manual(values = outline_colors) +  
  ylab('Percentage of malignant cells') + 
  geom_text_repel(aes(label = paired), min.segment.length = 1, max.overlaps = Inf) + 
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
    legend.position = "none"
  ) +
  labs(x = "") +
theme(legend.position = "none")+
geom_signif(comparisons = list(c("A", "B"), c("B", "C"), c("A", "C")),
              annotations = sprintf("p = %.3f", pvalues),
              y_position = 85,  
              map_signif_level = FALSE,
              test = "wilcox.test", 
              step_increase = 0.07)
print(a)
dev.off()



#################################################################################################
###### Figure 3D ########## CLUSTERING OF BINARY CLINICAL VARIABLES #############################
#################################################################################################

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

df <- data.frame(read.table('raw_binary response.txt', header=T, sep="\t", stringsAsFactors=F))
#### load cinical data from Supplementary Table 1

rownames (df) <- df[,1]
df <- df[,-1]
df <- df [-4, ]
mat <- as.matrix(df)
bin_mat <- mat

# ── NA-aware distance function ────────────────────────────────────────────────
na_aware_dist <- function(x, method = "euclidean") {
  n <- nrow(x)
  d <- matrix(0, n, n)
  for (i in seq_len(n - 1)) {
    for (j in (i + 1):n) {
      valid <- !is.na(x[i, ]) & !is.na(x[j, ])
      n_valid <- sum(valid)
      if (n_valid == 0) {
        d[i, j] <- d[j, i] <- NA
      } else {
        d[i, j] <- d[j, i] <- dist(rbind(x[i, valid], x[j, valid]),
                                     method = method) / n_valid
      }
    }
  }
  rownames(d) <- colnames(d) <- rownames(x)
  as.dist(d)
}


dist_method <- "manhattan"
row_dist <- na_aware_dist(bin_mat, method = dist_method)
col_dist <- na_aware_dist(t(bin_mat), method = dist_method)

row_clust <- hclust(row_dist, method = "ward")
col_clust <- hclust(col_dist, method = "ward")


colours <- c("steelblue", "firebrick")


annotation_row <- data.frame(mat[, rev(col_clust$order)])
annotation_row$baseline_subtype <- df1$Baseline_subtype
annotation_row$eot_subtype <- df1$EOT_subtype

anno_colors <- list(
  residual_tumor_bin              = c("1" = "#326A43", "2" = "#A02629", "3" = "grey80"),
  size_bin_20                	= c("1" = "#E59EDD", "2" = "#78206E", "3" = "grey80"),
  Prosigna_bin                   = c("1" = "#CC9900", "2" = "#996633", "3" = "grey80"),
  Ki67_bin_10               	 = c("1" = "#F6C6AD", "2" = "#C04F15", "3" = "grey80"),
  Node.bin                  	= c("1" = "#5CA46A", "2" = "#C55A5A", "3" = "grey80"),
  scRNA_percentage_bin           = c("1" = "#326A43", "2" = "#A02629", "3" = "grey80"),
  scRNA_KI_bin                  = c("1" = "#F6C6AD", "2" = "#C04F15", "3" = "grey80"),
  perc_decrease_bin    			= c("1" = "#92D050", "2" = "#FF7C80", "3" = "grey80"),
  Metastasis_bin          		= c("1" = "#5CA46A", "2" = "#C55A5A", "3" = "grey80"),
  scRNA_decrease_bin85        	 = c("1" = "#92D050", "2" = "#FF7C80", "3" = "grey80"),
  grade_bin                     = c("1" = "#5CA46A", "2" = "#C55A5A", "3" = "grey80"),
  PR_bin 						= c("1" = "#5CA46A", "2" = "#C55A5A", "3" = "grey80"),
  baseline_subtype 				= c("LumA" = "blue", "LumB" = "light blue", "Her2" = "deep pink"), 
  eot_subtype 				    = c("LumA" = "blue", "LumB" = "light blue", "Her2" = "deep pink")
)

plot <- pheatmap(
  t(mat),
  cluster_rows  = col_clust,
  cluster_cols  = row_clust,
  col           = colours, 
  annotation_col    = annotation_row,
  annotation_colors = anno_colors,
  display_numbers = FALSE,
  border_color  = "white",
  cutree_cols   = 2
)
print (plot)

pdf ('Figure 3D Clustering of binary clinical and molecula data.pdf')
print (plot)
dev.off()



############################################################
###### Figure 3E ########## Milo BR_GR SWARMPLOTS ##########
############################################################

library(Seurat)
library(plyr)
library(dplyr)
library(future)
library (ggplot2)
library (ggrepel)
library (stringi)
library(tidyverse)
library(miloR)

##################################
### need outputs of PP_milo_R##### 
##################################  
################################## 
## final.id plot timeA R vs NR ### 
################################## 

pdf('Figure 3Ea swarmplot final ids BR GR A.pdf', width = 7, height = 10)
da_results_B_4_rev <- as.data.frame(data.frame(read.table("final_BR_GR_A.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)))
vector_r <- rev(data.frame(da_results_B_4_rev %>% group_by(final.ids) %>% summarize(avg_logFC = mean(logFC)) %>% arrange(avg_logFC)))
da_results_B_4_rev$final.ids <- factor(da_results_B_4_rev$final.ids, levels = rev(vector_r$final.ids))
da_results_B_4_rev <- da_results_B_4_rev[order(match(da_results_B_4_rev$final.ids, vector_r$final.ids)),]


p <- plotDAbeeswarm(da_results_B_4_rev, group.by = "final.ids", alpha=0.05) +
 scale_colour_gradient2(high = "#326A43", low = "#7F0002", na.value = "grey85", space = "Lab") + 

  xlab("") + 
  scale_y_continuous(breaks = c(-5, -2.5, 0, 2.5, 5 ), limits = c(-6.5, 6.5)) +
  theme(axis.text.y = element_text(colour = rev(c('#8ba831','#8ba831','#8ba831',  'orange',  'orange', '#8ba831','#8ba831','#8ba831',
  'orange', 'orange', '#8ba831','#8ba831',  '#4f4fc5', '#8ba831',  'orange', '#8ba831', 'orange',  'orange',  'orange', 
 '#8ba831','#8ba831', 'orange', '#8ba831','#8ba831','#8ba831','#8ba831','#8ba831',  'orange',  '#4f4fc5',  'orange',  '#4f4fc5', 
 '#8ba831', '#4f4fc5',  'orange',  '#4f4fc5',  '#4f4fc5',  '#4f4fc5',  '#4f4fc5',  '#4f4fc5',  'orange',  '#4f4fc5',  '#4f4fc5', 
 'orange',  '#4f4fc5',  '#4f4fc5'))),
axis.text.x = element_text(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank())

da_results_B_4_sig_1 <- da_results_B_4_rev[da_results_B_4_rev$FDR <0.05 & da_results_B_4_rev$logFC <0,]
da_results_B_4_sig_2 <- da_results_B_4_rev[da_results_B_4_rev$FDR <0.05 & da_results_B_4_rev$logFC >0,]
count_sig_1 <- as.data.frame(table(da_results_B_4_sig_1$final.ids))
colnames(count_sig_1) <- c("final.ids", "count_sig_1")
count_sig_2 <- as.data.frame(table(da_results_B_4_sig_2$final.ids))
colnames(count_sig_2) <- c("final.ids", "count_sig_2")
count_data <- merge(count_sig_1, count_sig_2, by = "final.ids", all = TRUE)
count_data[is.na(count_data)] <- 0

  p$data <- p$data[order(!is.na(p$data$logFC_color), p$data$logFC_color), ]

p <- p + mapply(function(x, plus_label, minus_label) { 
    list(annotate("text", x = x, y = 6.25, label = paste0("(", plus_label, ")"), hjust = 1, size = 3),
      annotate("text", x = x, y = -6.25, label = paste0("(", minus_label, ")"), hjust = 1, size = 3))}, 
  seq_along(count_data$final.ids), 
  count_data$count_sig_2[rev(match(vector_r$final.ids, count_data$final.ids))], 
  count_data$count_sig_1[rev(match(vector_r$final.ids, count_data$final.ids))], 
  SIMPLIFY = FALSE)
print(p)

dev.off()


################################## 
## final.id plot timeB R vs NR ### 
################################## 

pdf('Figure 3Eb swarmplot final ids BR GR B.pdf', width = 7, height = 10)
da_results_B_4_rev <- as.data.frame(data.frame(read.table("final_BR_GR_B.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)))

da_results_B_4_rev$final.ids <- factor(da_results_B_4_rev$final.ids, levels = rev(vector_r$final.ids))
da_results_B_4_rev <- da_results_B_4_rev[order(match(da_results_B_4_rev$final.ids, vector_r$final.ids)),]


p <- plotDAbeeswarm(da_results_B_4_rev, group.by = "final.ids", alpha=0.05) +
 scale_colour_gradient2(high = "#326A43", low = "#7F0002", na.value = "grey85", space = "Lab") + 
  xlab("") + 
  scale_y_continuous(breaks = c(-5, -2.5, 0, 2.5, 5 ), limits = c(-6.5, 6.5)) +
  theme(axis.text.y = element_text(colour = rev(c('#8ba831','#8ba831','#8ba831',  'orange',  'orange', '#8ba831','#8ba831','#8ba831',
  'orange', 'orange', '#8ba831','#8ba831',  '#4f4fc5', '#8ba831',  'orange', '#8ba831', 'orange',  'orange',  'orange', 
 '#8ba831','#8ba831', 'orange', '#8ba831','#8ba831','#8ba831','#8ba831','#8ba831',  'orange',  '#4f4fc5',  'orange',  '#4f4fc5', 
 '#8ba831', '#4f4fc5',  'orange',  '#4f4fc5',  '#4f4fc5',  '#4f4fc5',  '#4f4fc5',  '#4f4fc5',  'orange',  '#4f4fc5',  '#4f4fc5', 
 'orange',  '#4f4fc5',  '#4f4fc5'))),
axis.text.x = element_text(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank())

da_results_B_4_sig_1 <- da_results_B_4_rev[da_results_B_4_rev$FDR <0.05 & da_results_B_4_rev$logFC <0,]
da_results_B_4_sig_2 <- da_results_B_4_rev[da_results_B_4_rev$FDR <0.05 & da_results_B_4_rev$logFC >0,]
count_sig_1 <- as.data.frame(table(da_results_B_4_sig_1$final.ids))
colnames(count_sig_1) <- c("final.ids", "count_sig_1")
count_sig_2 <- as.data.frame(table(da_results_B_4_sig_2$final.ids))
colnames(count_sig_2) <- c("final.ids", "count_sig_2")
count_data <- merge(count_sig_1, count_sig_2, by = "final.ids", all = TRUE)
count_data[is.na(count_data)] <- 0

  p$data <- p$data[order(!is.na(p$data$logFC_color), p$data$logFC_color), ]

  
p <- p + mapply(function(x, plus_label, minus_label) { 
    list(annotate("text", x = x, y = 6.25, label = paste0("(", plus_label, ")"), hjust = 1, size = 3),
      annotate("text", x = x, y = -6.25, label = paste0("(", minus_label, ")"), hjust = 1, size = 3))}, 
  seq_along(count_data$final.ids), 
  count_data$count_sig_2[rev(match(vector_r$final.ids, count_data$final.ids))], 
  count_data$count_sig_1[rev(match(vector_r$final.ids, count_data$final.ids))], 
  SIMPLIFY = FALSE)
print(p)

dev.off()


################################## 
## final.id plot timeC R vs NR ### 
################################## 

pdf('Figure 3Ec swarmplot final ids BR GR C.pdf', width = 7, height = 10)

da_results_B_4_rev <- as.data.frame(data.frame(read.table("final_BR_GR_A.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)))
vector_r <- rev(data.frame(da_results_B_4_rev %>% group_by(final.ids) %>% summarize(avg_logFC = mean(logFC)) %>% arrange(avg_logFC)))

da_results_B_4_rev <- as.data.frame(data.frame(read.table("final_BR_GR_C.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)))
vector_r <- vector_r[-1,]
da_results_B_4_rev$final.ids <- factor(da_results_B_4_rev$final.ids, levels = rev(vector_r$final.ids))
da_results_B_4_rev <- da_results_B_4_rev[order(match(da_results_B_4_rev$final.ids, vector_r$final.ids)),]


p <- plotDAbeeswarm(da_results_B_4_rev, group.by = "final.ids", alpha=0.05) +
 scale_colour_gradient2(high = "#326A43", low = "#7F0002", na.value = "grey85", space = "Lab") + 
  xlab("") + 
  scale_y_continuous(breaks = c(-5, -2.5, 0, 2.5, 5 ), limits = c(-6.5, 6.5)) +
  theme(axis.text.y = element_text(colour = rev(c('#8ba831','#8ba831','#8ba831',  'orange',  'orange', '#8ba831','#8ba831','#8ba831',
  'orange', 'orange', '#8ba831','#8ba831',  '#4f4fc5', '#8ba831',  'orange', '#8ba831', 'orange',  'orange',  'orange', 
 '#8ba831','#8ba831', 'orange', '#8ba831','#8ba831','#8ba831','#8ba831','#8ba831',  'orange',  '#4f4fc5',  'orange',  '#4f4fc5', 
 '#8ba831', '#4f4fc5',  'orange',  '#4f4fc5',  '#4f4fc5',  '#4f4fc5',  '#4f4fc5',  '#4f4fc5',  'orange',  '#4f4fc5',  '#4f4fc5', 
 'orange',  '#4f4fc5',  '#4f4fc5'))),
axis.text.x = element_text(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank())

da_results_B_4_sig_1 <- da_results_B_4_rev[da_results_B_4_rev$FDR <0.05 & da_results_B_4_rev$logFC <0,]
da_results_B_4_sig_2 <- da_results_B_4_rev[da_results_B_4_rev$FDR <0.05 & da_results_B_4_rev$logFC >0,]
count_sig_1 <- as.data.frame(table(da_results_B_4_sig_1$final.ids))
colnames(count_sig_1) <- c("final.ids", "count_sig_1")
count_sig_2 <- as.data.frame(table(da_results_B_4_sig_2$final.ids))
colnames(count_sig_2) <- c("final.ids", "count_sig_2")
count_data <- merge(count_sig_1, count_sig_2, by = "final.ids", all = TRUE)
count_data[is.na(count_data)] <- 0

  p$data <- p$data[order(!is.na(p$data$logFC_color), p$data$logFC_color), ]

p <- p + mapply(function(x, plus_label, minus_label) { 
    list(annotate("text", x = x, y = 6.25, label = paste0("(", plus_label, ")"), hjust = 1, size = 3),
      annotate("text", x = x, y = -6.25, label = paste0("(", minus_label, ")"), hjust = 1, size = 3))}, 
  seq_along(count_data$final.ids), 
  count_data$count_sig_2[rev(match(vector_r$final.ids, count_data$final.ids))], 
  count_data$count_sig_1[rev(match(vector_r$final.ids, count_data$final.ids))], 
  SIMPLIFY = FALSE)
print(p)

dev.off()

######################## ##############################################
#### Figure 3F CD4 NAIVE BOXPLOTS #####################################
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
probe <- "CD4_Naive"

metadata1<- metadata[!metadata$main.ids=="Epithelial_cells" & !metadata$final.ids=="UNA",]
dta <- as.data.frame.matrix(prop.table(table (metadata1$orig.ident, metadata1[, colnames(metadata1)%in%ids]),1)*100)

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

pdf ('Figure 3F CD4 Naive GR BR.pdf')

a <- ggplot(dta_bas, aes(x = !!sym(gory), y = !!sym(probe), label = paired)) + 
  geom_boxplot(aes(fill = factor(!!sym(gory))), color = outline_colors, fill = box_fill_colors, alpha = 0.2, size = 1.5) + 
  geom_point(aes(colour = factor(!!sym(gory)), fill = factor(!!sym(gory))), shape = 21, size = 5, stroke = 1.2) + 
  scale_fill_manual(values= point_fill_colors) + 
  scale_colour_manual(values=outline_colors) + 
  ylab('Percentage of naive CD4 cells in the ME') + 

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


######################## ##############################################
#### Figure 3G CD14 MONOCYTES BOXPLOTS ################################
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
probe <- "Mono_CD14"

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

pdf ('Figure 3G CD14 Monocytes GR BR.pdf')

a <- ggplot(dta_bas, aes(x = !!sym(gory), y = !!sym(probe), label = paired)) + 
  geom_boxplot(aes(fill = factor(!!sym(gory))), color = outline_colors, fill = box_fill_colors, alpha = 0.2, size = 1.5) + 
  geom_point(aes(colour = factor(!!sym(gory)), fill = factor(!!sym(gory))), shape = 21, size = 5, stroke = 1.2) + 
  scale_fill_manual(values= point_fill_colors) + 
  scale_colour_manual(values=outline_colors) + 
  ylab('Percentage of CD14 monocytes in the ME') + 
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


###########################################################
#### Figure 3H HEATMAP BASELINE MICROENVIRONMENT ##########
###########################################################

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
#### load the metadata available at EGA
met_clin <- data.frame(read.table("clinical_data.txt", header=T, sep="\t", stringsAsFactors=F))
##### load clinical data in Supp Table 1
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
annot <- metadata[, colnames (metadata) %in% c("orig.ident", "perc_malig_surg10", "clone_cat", "BAS_subtype")]
annot <- annot[!duplicated (metadata$orig.ident), ]
table (annot$orig.ident == rownames (dta))
dta <- cbind (dta , annot)

a <- as.data.frame.matrix(prop.table (table (metadata$orig.ident, metadata$clone_cat),1)*100)
a$category <- stri_sub(rownames (a),-1)
a1 <- a[a$category == "A",]

a1$cat_res<- ifelse (a1$res <10, '[0-10]', a1$res)
a1$cat_res<- ifelse (a1$cat_res >10, '>10', a1$cat_res )
a1$cat_res<- ifelse (a1$excl==100, NA, a1$cat_res)

a1$cat_sens<- ifelse (a1$sens >50, '[50-100]', '[0-50]')
a1$cat_sens<- ifelse (a1$excl==100, NA, a1$cat_sens)

dta1 <- dta1[dta$category == "A",]
dta <- dta[dta$category == "A",]

annotation_row <- data.frame(cbind(dta$perc_malig_surg10, dta$BAS_subtype, a1$cat_res, a1$cat_sens))
rownames(annotation_row) <- rownames(dta)

average_scores <- data.frame(scores %>%
  select(orig.ident, starts_with("cluster")) %>%  
  group_by(orig.ident) %>%                         
  summarize(across(starts_with("cluster"), mean, na.rm = TRUE)))  

score_bas <- average_scores[average_scores$orig.ident%in% rownames (annotation_row),]
(score_bas$orig.ident==rownames (annotation_row))
annotation_row$interferon <- score_bas$cluster10

color_palette <- colorRampPalette(brewer.pal(9, "Purples"))(100)

colnames(annotation_row) <- c('Response', 'PAM50', 'res', 'sens', 'interferon')
annotation_row <- annotation_row[,c(1,2,5,4,3)]

max_val <- 20
breaks_vals <- c(seq(0, max_val, length.out = 101))
breaksList = seq(0, 35, by = 0.05)
colors <- colorRampPalette(c("white", "indianred1", "brown3", "brown4"))(length(breaksList))
col.pal <- colorRampPalette(c("dodgerblue4", "white", "white", "brown3", "brown3", "brown3"))(length(breaksList))

anno_colors = list(
  Response = c("GR" = "green", "BR" = "red"),
  PAM50 = c("LUM-A" = "blue", "LUM-B" = "light blue", "HER-2-enriched" = "deep pink"), 
  res = c("[0-10]" = "lightcoral", '>10' = "red"),  
  sens = c('[0-50]' = "lightgreen", '[50-100]' = "darkgreen"),
  interferon = color_palette
)

dcols <- hcluster(t(dta1), method = 'spearman', link = 'complete')
drows <- hcluster(dta1, method = 'spearman', link = 'complete')

pdf ('Figure 3H TME Baseline clustering.pdf')
pheatmap(
  t(dta1), 
  show_rownames = TRUE, show_colnames = TRUE,  
  cluster_rows = dcols, cluster_cols = drows, fontsize = 9,
  annotation_col = annotation_row, annotation_colors = anno_colors,
  col = colors, 
   width = 10, height = 8,  , border_color = NA
)
dev.off()

#########################################################
#### Figure 3I BOXPLOT INTERFERON BY IMMUNE CLUSTERS ####
#########################################################
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

annotation_row <- data.frame(cbind(dta$perc_malig_surg10, dta$immune_cluster))
rownames(annotation_row) <- rownames(dta)

average_scores <- data.frame(scores %>%
  select(orig.ident, starts_with("cluster")) %>%  
  group_by(orig.ident) %>%                         
  summarize(across(starts_with("cluster"), mean, na.rm = TRUE)))  

score_bas <- average_scores[average_scores$orig.ident%in% rownames (annotation_row),]
(score_bas$orig.ident==rownames (annotation_row))
annotation_row$interferon <- score_bas$cluster10

colnames(annotation_row) <- c('Response', 'cutree1', 'interferon')
annotation_row$cutree1 <- factor(annotation_row$cutree1, levels = c('LumA_R', 'LumB_R', 'LumB_NR'))

annotation_row$patient_number <- gsub ("A", "", rownames(annotation_row))
annotation_row$patient_number <- gsub ("P", "", annotation_row$patient_number)

pval_LumAGR_v_LumBGR <- wilcox.test(annotation_row$interferon[annotation_row$cutree1 == 'LumA_R'],
                       annotation_row$interferon[annotation_row$cutree1 == 'LumB_R'], 
                       paired = FALSE)$p.value
pval_LumBGR_v_LumBBR <- wilcox.test(annotation_row$interferon[annotation_row$cutree1 == 'LumB_R'],
                       annotation_row$interferon[annotation_row$cutree1 == 'LumB_NR'], 
                       paired = FALSE)$p.value
pval_LumAGR_v_LumBBR <- wilcox.test(annotation_row$interferon[annotation_row$cutree1 == 'LumA_R'],
                       annotation_row$interferon[annotation_row$cutree1 == 'LumB_NR'], 
                       paired = FALSE)$p.value
pvalues <- c(pval_LumAGR_v_LumBGR, pval_LumBGR_v_LumBBR, pval_LumAGR_v_LumBBR)

pdf ('Figure 3I boxplot_interferon_score_TME_clusters.pdf')

ggplot(annotation_row, aes(x = cutree1, y = interferon)) + 
geom_boxplot_pattern(
    aes(fill = cutree1, pattern = cutree1), 
    pattern = c("none", "stripe", "stripe"),  
    pattern_density = c(0, 0.1, 0.1),          
    pattern_color = c(NA, "darkgreen", "red"),
    pattern_fill = c(NA, "darkgreen", "red"), 	
    alpha = 0.8, outlier.colour = NA, notch = FALSE, size = 1
  ) + 
  geom_point(aes(colour = cutree1), shape = 21, size = 5) + 
  scale_fill_manual(values = c('#0000FF', '#9ED2E2', '#9ED2E2')) +  
  scale_colour_manual(values = c('#3333CC', '#5EB3CE', '#5EB3CE')) + 
  ylab('Interferon score') + 
  geom_text_repel(aes(label = patient_number), 
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

			  
#########################################################
#### Figure 3J Boxplot immune cluster vs M1 #############
#########################################################
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

pval_LumAGR_v_LumBGR <- wilcox.test(dta1$MF_REG_M1[dta1$immune_cluster == 'LumA_R'],
                       dta1$MF_REG_M1[dta1$immune_cluster == 'LumB_R'], 
                       paired = FALSE)$p.value
pval_LumBGR_v_LumBBR <- wilcox.test(dta1$MF_REG_M1[dta1$immune_cluster == 'LumB_R'],
                       dta1$MF_REG_M1[dta1$immune_cluster == 'LumB_NR'], 
                       paired = FALSE)$p.value
pval_LumAGR_v_LumBBR <- wilcox.test(dta1$MF_REG_M1[dta1$immune_cluster == 'LumA_R'],
                       dta1$MF_REG_M1[dta1$immune_cluster == 'LumB_NR'], 
                       paired = FALSE)$p.value
pvalues <- c(pval_LumAGR_v_LumBGR, pval_LumBGR_v_LumBBR, pval_LumAGR_v_LumBBR)

pdf ('Figure 3J boxplot_MF_proinf_TME_clusters.pdf')

ggplot(dta1, aes(x = immune_cluster, y = MF_REG_M1)) + 
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
  ylab('MF_REG_M1') + 
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