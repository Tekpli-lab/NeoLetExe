###################################################
####### FIGURE 5 ##################################
###################################################

########################################################################
############# FIG 5A - HEATMAP BEHAVIOURS VERSUS PATIENT       ##
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
met_clin <- data.frame(read.table("clinical_data.txt", header=T, sep='\t', row.names=))
##### load clinical data in Supp Table 1

annotation_row <- data.frame(cbind(met_clin$BAS_subtype, met_clin$response))
rownames (annotation_row) <- met_clin$patients
colnames(annotation_row) <- c('PAM50','Response')
anno_colors = list(
  PAM50 = c("LUM-A" = "blue", "LUM-B" = "light blue", "HER-2-enriched" = "deep pink"),
  Response = c("GR" = "#7BA58C", "BR" = "#A02629"))

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

  annotation_row <- annotation_row[rownames(annotation_row) %in% rownames(agg), ]
  annotation_row <- annotation_row[match(rownames(agg), rownames(annotation_row)), ]
 
  dcols <- hcluster(t(agg), method = 'correlation', link = 'complete')
  drows <- hcluster(agg, method = 'correlation', link = 'complete')
  breaksList = seq(-1, 1, by = 0.05)
  col.pal <- colorRampPalette(c("white","brown3","brown3"))(length(breaksList))

pdf	('Figure 5A Heatmap patient behavior ABC.pdf')	 
res <-  pheatmap(
    t(agg), 
  show_rownames = TRUE, show_colnames = TRUE,  
  cluster_rows = dcols, cluster_cols = drows, fontsize = 9,
  annotation_col = annotation_row, annotation_colors = anno_colors, color = col.pal, , border_color = NA)
  print (res)
dev.off()



########################################################################
############# FIG 5 BCD - BOXPLOT OF PATIENT behaviour thee time points       ##
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

outline_colors <- c("#A02629", '#7BA58C') 
box_fill_colors <- c('#E5A1A4', '#C3EAD6') 
point_fill_colors <- c('#C32E32', '#9ED3B9') 

pdf('Figure 5B dynamic_TCR_contraction_contraction.pdf')	
var1 <- 'BR'
var2 <- 'GR'
pval <- wilcox.test(norm_EC$contraction[norm_EC$response == var1],
                    norm_EC$contraction[norm_EC$response == var2], 
                    paired = FALSE)$p.value  
pvalues <- c(pval)


ggplot(norm_EC, aes(x = response, y = contraction)) + 
  geom_boxplot(aes(fill = response), color = outline_colors, fill = box_fill_colors, alpha = 0.2, size = 1.5) + 
  scale_fill_manual(values= point_fill_colors) + 
  scale_colour_manual(values=outline_colors) + 
  ylab('Fraction of contractionAB - contractionBC clones') + 
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

pdf('Figure 5C dynamic_TCR_expansion_expansion.pdf')	
var1 <- 'BR'
var2 <- 'GR'
pval <- wilcox.test(norm_EC$expansion[norm_EC$response == var1],
                    norm_EC$expansion[norm_EC$response == var2], 
                    paired = FALSE)$p.value  
pvalues <- c(pval)


ggplot(norm_EC, aes(x = response, y = expansion)) + 
  geom_boxplot(aes(fill = response), color = outline_colors, fill = box_fill_colors, alpha = 0.2, size = 1.5) + 
  scale_fill_manual(values= point_fill_colors) + 
  scale_colour_manual(values=outline_colors) + 
  ylab('Fraction of expansionAB - expansionBC clones') + 
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


pdf('Figure 5D dynamic_TCR_new_late_expansion.pdf')			   
var1 <- 'BR'
var2 <- 'GR'
pval <- wilcox.test(norm_EC$new_late_expansion[norm_EC$response == var1],
                    norm_EC$new_late_expansion[norm_EC$response == var2], 
                    paired = FALSE)$p.value  
pvalues <- c(pval)


ggplot(norm_EC, aes(x = response, y = new_late_expansion)) + 
  geom_boxplot(aes(fill = response), color = outline_colors, fill = box_fill_colors, alpha = 0.2, size = 1.5) + 
  geom_point(aes(colour = response, fill = response), shape = 21, size = 5) +  # Fixed parentheses here
  scale_fill_manual(values= point_fill_colors) + 
  scale_colour_manual(values=outline_colors) + 
  ylab('Fraction of new_late_expansion clones') + 
  geom_text_repel(aes(label = patient_number), max.overlaps = Inf, box.padding = 0.5, point.padding = 0.3, segment.size = 0.2, segment.color = 'grey50') +  
  scale_x_discrete(labels = setNames(c(var1, var2), c(var1, var2))) +  # Make sure var1 and var2 relate to response
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


		
###############################################################################################
############# FIG 5E - HEATMAP of TCR behavior correlation with cell states abundance #########
###############################################################################################
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

all_data_frames <- list()

for (ids in c("final.ids")) {
  for (subs in c("ME")) {

    lis <- unique(metadata[, colnames(metadata) %in% ids])
    lisME <- lis[!grepl("malignant", lis, ignore.case = TRUE) & !grepl("UNA", lis, ignore.case = TRUE) & !grepl("Normal_epithelial", lis, ignore.case = TRUE)]
    
    selec <- get(paste0('lis', subs))
    metadata1 <- metadata[metadata[, colnames(metadata) %in% ids] %in% selec,]
    dta <- as.data.frame.matrix(prop.table(table(metadata1$orig.ident, metadata1[, colnames(metadata) %in% ids]), 1) * 100)
    
    colnames(dta) <- paste(colnames(dta), ids, subs, sep="_")
    all_data_frames[[paste(ids, subs, sep = "_")]] <- dta
  }
}


abund <- Reduce(function(x, y) cbind(x, y), all_data_frames)

abund_A <- abund [gsub(".*(.)$", "\\1", rownames (abund))=="A",]
abund_B <- abund [gsub(".*(.)$", "\\1", rownames (abund))=="B",]
abund_C <- abund [gsub(".*(.)$", "\\1", rownames (abund))=="C",]

delta_abund_A_B <- abund_A[!rownames (abund_A)=="P04A",] - abund_B
delta_abund_A_C <- abund_A[!rownames (abund_A)=="P04A",] - abund_C
delta_abund_B_C <- abund_B - abund_C

colnames(abund_A) <- paste0('abund_A_', colnames(abund_A))
colnames(abund_B) <- paste0('abund_B_', colnames(abund_B))
colnames(abund_C) <- paste0('abund_C_', colnames(abund_C))
colnames(delta_abund_A_B) <- paste0('delta_abund_A_B_', colnames(delta_abund_A_B))
colnames(delta_abund_A_C) <- paste0('delta_abund_A_C_', colnames(delta_abund_A_C))
colnames(delta_abund_B_C) <- paste0('delta_abund_B_C_', colnames(delta_abund_B_C))

tot_abund <- cbind (abund_A[!rownames (abund_A)=="P04A",], abund_B, abund_C, delta_abund_A_B, delta_abund_A_C, delta_abund_B_C)
columns_to_remove <- names((colSums(tot_abund == 0) / nrow(tot_abund))[(colSums(tot_abund == 0) / nrow(tot_abund)) > 0.9])
tot_abund <- tot_abund[, !(names(tot_abund) %in% columns_to_remove)]


calculate_correlation <- function(df1, df2, method = c("spearman", "pearson")) {
  
  cor_matrix <- matrix(NA, nrow = ncol(df1), ncol = ncol(df2))
  p_matrix <- matrix(NA, nrow = ncol(df1), ncol = ncol(df2))
  
  rownames(cor_matrix) <- colnames(df1)
  colnames(cor_matrix) <- colnames(df2)
  rownames(p_matrix) <- colnames(df1)
  colnames(p_matrix) <- colnames(df2)


  for (i in 1:ncol(df1)) {
    for (j in 1:ncol(df2)) {
      test <- cor.test(df1[[i]], df2[[j]], method = method, use = "pairwise.complete.obs")
      cor_matrix[i, j] <- test$estimate  
      p_matrix[i, j] <- test$p.value     
    }
  }
  
  return(list(correlation_df = cor_matrix, p_values_df = p_matrix))
}

df1 = tot_abund[, grep("^abund_(A|B|C).*_final.ids_ME$", names(tot_abund))]
df2 = TCR_norm_cell_bind[, grep("(all_AC_|all_BC_|all_AB_)", names(TCR_norm_cell_bind), value = TRUE)]
df2 = df2[, !grepl("_amb_", names(df2))]
df2 = df2[, !grepl("(_singleton)", names(df2))]

results <- calculate_correlation(df1, df2, method = "spearman")
correlation_df <- results$correlation_df
p_values_df <- results$p_values_df

threshold <- 0.05
selected_rows <- rownames(p_values_df)[apply(p_values_df, 1, function(x) { any( x <= threshold) })]
cor_all_f <- correlation_df[rownames (correlation_df)%in%c(selected_rows) , ]

dcols <- hcluster(t(cor_all_f), method = 'correlation', link = 'ward')
drows <- hcluster(cor_all_f, method = 'correlation', link = 'ward')
breaksList = seq(-1, 1, by = 0.05)
col.pal <- colorRampPalette(c("dodgerblue4","white","brown3"))(length(breaksList))

annotation_row= data.frame ( sub("abund_([^_]*)_.*", "\\1", rownames(cor_all_f)))
colnames(annotation_row) <- c("time_point")
rownames(annotation_row) <- rownames(cor_all_f)

anno_colors = list(time_point = c("A" = "#DA3A09", "B" = "#9C2C98","C" = "#320E77"))

pdf ('Figure 5E correlation abundance versus TCR.pdf')	
map <- pheatmap (cor_all_f,  color = col.pal,show_rownames = T, show_colnames = T, main = "",  cluster_rows = drows, cluster_cols = dcols, breaks = breaksList, border_color = NA, annotation_row=annotation_row, annotation_colors = anno_colors)
dev.off()


###############################################################################################
############# FIG 5F - Expanded clones in cell states baseline immune clusters     ############
###############################################################################################

library(Seurat)
library(plyr)
library(dplyr)
library(future)
library (ggplot2)
library (ggrepel)
library (stringi)
library(tidyverse)
library(ggsignif)
library (ggpattern)

obj <- readRDS("combined_metadata_repertoire_tcr_mix1.rds")
#### load the Seurat_object_update with scRepertoire from @@@@@@@@@@@@@@@@

obj_f <- obj[!obj$orig.ident%in% c('P24A', 'P06A'),]
obj_f <- obj_f[obj_f$time_point=="A",]
obj_f <- obj_f[obj_f$clonalFrequency>=2 & !is.na (obj_f$clonalFrequency),]


tabix <- as.data.frame.matrix(prop.table(table (obj_f$orig.ident, obj_f$final.ids),1)*100)
rownames (tabix) <- substr(rownames(tabix), 1, nchar(rownames(tabix)) - 1)
met_clin <- data.frame(read.table("clinical_data.txt", header=T, sep='\t', row.names=))
met_clin <- met_clin[met_clin$patients%in%rownames(tabix),]
all (met_clin$patients%in%rownames(tabix))
tabix$response <- met_clin$immune_cluster
tabix$patient_number <- met_clin$patient_number

pdf ('Figure 5F boxplot expanded clones cell states.pdf')	

##############
#### CD8EM ###
##############

pval_LumAGR_v_LumBGR <- wilcox.test(tabix$ CD8_EM[tabix$response == 'LumA_R'],
                       tabix$ CD8_EM[tabix$response == 'LumB_R'], 
                       paired = FALSE)$p.value
pval_LumBGR_v_LumBBR <- wilcox.test(tabix$ CD8_EM[tabix$response == 'LumB_R'],
                       tabix$ CD8_EM[tabix$response == 'LumB_NR'], 
                       paired = FALSE)$p.value
pval_LumAGR_v_LumBBR <- wilcox.test(tabix$ CD8_EM[tabix$response == 'LumA_R'],
                       tabix$ CD8_EM[tabix$response == 'LumB_NR'], 
                       paired = FALSE)$p.value
pvalues <- c(pval_LumAGR_v_LumBGR, pval_LumBGR_v_LumBBR, pval_LumAGR_v_LumBBR)


tabix$response <- factor(tabix$response, levels = c('LumA_R', 'LumB_R', 'LumB_NR'))

ggplot(tabix, aes(x = response, y = CD8_EM)) + 
  geom_boxplot_pattern(
    aes(fill = response, pattern = response), 
    pattern = c("none", "stripe", "stripe"),  
    pattern_density = c(0, 0.1, 0.1),          
    pattern_color = c(NA, "darkgreen", "red"),
    pattern_fill = c(NA, "darkgreen", "red"), 	
    alpha = 0.8, outlier.colour = NA, notch = FALSE, size = 1
  ) + 
  geom_point(aes(colour = response), shape = 21, size = 5) + 
  scale_fill_manual(values = c('#0000FF', '#9ED2E2', '#9ED2E2')) +  
  scale_colour_manual(values = c('#3333CC', '#5EB3CE', '#5EB3CE')) + 
  ylab('% of expanded CD8EM') + 
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
              step_increase = 0.07)+
			  coord_cartesian(ylim = c(0, 80))

##############
#### CD4EM ###
##############

pval_LumAGR_v_LumBGR <- wilcox.test(tabix$ CD4_EM[tabix$response == 'LumA_R'],
                       tabix$ CD4_EM[tabix$response == 'LumB_R'], 
                       paired = FALSE)$p.value
pval_LumBGR_v_LumBBR <- wilcox.test(tabix$ CD4_EM[tabix$response == 'LumB_R'],
                       tabix$ CD4_EM[tabix$response == 'LumB_NR'], 
                       paired = FALSE)$p.value
pval_LumAGR_v_LumBBR <- wilcox.test(tabix$ CD4_EM[tabix$response == 'LumA_R'],
                       tabix$ CD4_EM[tabix$response == 'LumB_NR'], 
                       paired = FALSE)$p.value
pvalues <- c(pval_LumAGR_v_LumBGR, pval_LumBGR_v_LumBBR, pval_LumAGR_v_LumBBR)


tabix$response <- factor(tabix$response, levels = c('LumA_R', 'LumB_R', 'LumB_NR'))

ggplot(tabix, aes(x = response, y = CD4_EM)) + 
  geom_boxplot_pattern(
    aes(fill = response, pattern = response), 
    pattern = c("none", "stripe", "stripe"),  
    pattern_density = c(0, 0.1, 0.1),          
    pattern_color = c(NA, "darkgreen", "red"),
    pattern_fill = c(NA, "darkgreen", "red"), 	
    alpha = 0.8, outlier.colour = NA, notch = FALSE, size = 1
  ) + 
  geom_point(aes(colour = response), shape = 21, size = 5) + 
  scale_fill_manual(values = c('#0000FF', '#9ED2E2', '#9ED2E2')) +  
  scale_colour_manual(values = c('#3333CC', '#5EB3CE', '#5EB3CE')) + 
  ylab('% of expanded CD4EM') + 
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
              step_increase = 0.07)+
			  coord_cartesian(ylim = c(0, 80))
			  
			  

#################
#### CD8-EMRA ###
#################

pval_LumAGR_v_LumBGR <- wilcox.test(tabix$ CD8_EMRA[tabix$response == 'LumA_R'],
                       tabix$ CD8_EMRA[tabix$response == 'LumB_R'], 
                       paired = FALSE)$p.value
pval_LumBGR_v_LumBBR <- wilcox.test(tabix$ CD8_EMRA[tabix$response == 'LumB_R'],
                       tabix$ CD8_EMRA[tabix$response == 'LumB_NR'], 
                       paired = FALSE)$p.value
pval_LumAGR_v_LumBBR <- wilcox.test(tabix$ CD8_EMRA[tabix$response == 'LumA_R'],
                       tabix$ CD8_EMRA[tabix$response == 'LumB_NR'], 
                       paired = FALSE)$p.value
pvalues <- c(pval_LumAGR_v_LumBGR, pval_LumBGR_v_LumBBR, pval_LumAGR_v_LumBBR)


tabix$response <- factor(tabix$response, levels = c('LumA_R', 'LumB_R', 'LumB_NR'))

ggplot(tabix, aes(x = response, y = CD8_EMRA)) + 
  geom_boxplot_pattern(
    aes(fill = response, pattern = response), 
    pattern = c("none", "stripe", "stripe"),  
    pattern_density = c(0, 0.1, 0.1),          
    pattern_color = c(NA, "darkgreen", "red"),
    pattern_fill = c(NA, "darkgreen", "red"), 	
    alpha = 0.8, outlier.colour = NA, notch = FALSE, size = 1
  ) + 
  geom_point(aes(colour = response), shape = 21, size = 5) + 
  scale_fill_manual(values = c('#0000FF', '#9ED2E2', '#9ED2E2')) +  
  scale_colour_manual(values = c('#3333CC', '#5EB3CE', '#5EB3CE')) + 
  ylab('% of expanded CD8_EMRA') + 
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
              step_increase = 0.07)+
			  coord_cartesian(ylim = c(0, 80))

dev.off()