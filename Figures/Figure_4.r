###################################################
####### FIGURE 4 ##################################
###################################################

##################################################
############# FIG 4A - Reclustering of T cells  ##
##################################################

#output of PP-subclustering

#########################################################
## Figure 4B - bubble plots of clone size in UMAP form ##
#########################################################
library(Seurat)
library(dplyr)
library(scRepertoire)
library(ggraph)
library(circlize)
library(scales)
library(APackOfTheClones)
library(APackOfTheClones)
library(RColorBrewer)
library(SCpubr)
library(Nebulosa)

obj <- readRDS("combined_seurat_repertoire_tcr.rds")
#### load the Seurat_object_update with scRepertoire from @@@@@@@@@@@@@@@@

obj_f <- obj[,!obj@meta.data$orig.ident%in% c('P24A', 'P06A', 'P11B', 'P01C', 'P02C', 'P01B', 'P02B')]

palopal <- vector()
palopal [1] <-"#977899"  # NK_cyto
palopal [2] <-"#E0DEDE" # UNA
palopal [3] <-"#5FB0E8" # CD8_EM
palopal [4] <-"#00C957" # CD4_EM
palopal [5] <-"#F57C7C" # CD4_low
palopal [6] <-"#AFDEA2" # CD4_Naive
palopal [7] <-"#3385BB" # CD8_EMRA
palopal [8] <-"#E2CAE3" # NK_rest
palopal [9] <-"#DE9E83" # Treg
palopal [10] <-"#7F9D55" # CD4_Ex
palopal [11] <-"#99CAF7" # CD8_Naive
palopal [12] <- "#FE8D19" # GDT
palopal [13] <-"#3F5361" # CD8_Ex
palopal [14] <-"#D40F0F" # MAIT
palopal [15] <-"#8A5A46" # ILC3
palopal [16] <-"#FBB268" # CRTAM_GDT

desired_labels <- c("NK_cyto", "UNA", "CD8_EM", "CD4_EM", "CD4_low",
                    "CD4_Naive", "CD8_EMRA", "NK_rest", "Treg", "CD4_Ex",
                    "CD8_Naive", "GDT", "CD8_Ex", "MAIT", "ILC3", "CRTAM_GDT")
cluster_indices <- match(desired_labels, levels(Idents(obj_f)))

obj_apotc <- RunAPOTC(obj_f, verbose = FALSE)

obj_apotc <- AdjustAPOTC(
  obj_apotc,
  recolor_cluster = cluster_indices,
  new_color = as.vector(palopal),
  verbose = FALSE
)

pdf ('Figure 4B ballpacking_clone_size_states.pdf ')
viz <- APOTCPlot(obj_apotc)
viz
dev.off()



##################################################################################
############# FIG 4C - density plots only three time points homogene density   ##    
##################################################################################
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


load ('TCR_norm_cell_bind.RData')
#### load the summarized expansion_contraction data from PP_Expansion
TCR<- readRDS("combined_metadata_repertoire_tcr_mix1.rds")
#### load the scRepertoire output from  @@@@@@@@@@@@@@@@
obj <- readRDS("combined_seurat_repertoire_tcr.rds")
#### load the Seurat_object_update with scRepertoire from @@@@@@@@@@@@@@@@

TCR$umap_1 <- obj@reductions$umap@cell.embeddings[,1]
TCR$umap_2 <- obj@reductions$umap@cell.embeddings[,2]

concatenated_cols <- apply(statusAC1, 2, function(x) {
  paste(na.omit(x), collapse = ", ")
})
concAC <- data.frame(Concatenated = concatenated_cols)
colnames (concAC) <- 'AC'

concatenated_cols <- apply(statusAB1, 2, function(x) {
  paste(na.omit(x), collapse = ", ")
})
concAB <- data.frame(Concatenated = concatenated_cols)
colnames (concAB) <- 'AB'

concatenated_cols <- apply(statusBC1, 2, function(x) {
  paste(na.omit(x), collapse = ", ")
})
concBC <- data.frame(Concatenated = concatenated_cols)
colnames (concBC) <- 'BC'

mer <- merge (TCR, concAB, by.x='CTnt', by.y=0, all.x=T)
mer1 <- merge (mer, concBC, by.x='CTnt', by.y=0, all.x=T)
mer2 <- merge (mer1, concAC, by.x='CTnt', by.y=0, all.x=T)

mer2_A <- mer2[ mer2$time_point=="A" &!is.na(mer$CTnt),]
dim (mer2_A)
mer2_B <- mer2[ mer2$time_point=="B"&!is.na(mer$CTnt),]
dim (mer2_B)
mer2_C <- mer2[ mer2$time_point=="C"&!is.na(mer$CTnt),]
dim (mer2_C)

mer2_A_density <- ggplot() + 
  stat_density_2d(data = mer2_A, aes(x = umap_1, y = umap_2, fill = stat(level)), 
                  geom = "polygon", bins = 10, contour = TRUE)

mer2_B_density <- ggplot() + 
  stat_density_2d(data = mer2_B, aes(x = umap_1, y = umap_2, fill = stat(level)), 
                  geom = "polygon", bins = 10, contour = TRUE)

mer2_C_density <- ggplot() + 
  stat_density_2d(data = mer2_C, aes(x = umap_1, y = umap_2, fill = stat(level)), 
                  geom = "polygon", bins = 10, contour = TRUE)

density_range <- range(
  c(ggplot_build(mer2_A_density)$data[[1]]$level,
    ggplot_build(mer2_B_density)$data[[1]]$level,
    ggplot_build(mer2_C_density)$data[[1]]$level)
)

min_level <- density_range[1]
max_level <- density_range[2]

common_fill_scale <- scale_fill_gradientn(colors = c("white", rev(brewer.pal(6, "RdYlBu"))), 
                                           na.value = "white", limits = c(min_level, max_level), 
                                           name = "Density Level")
p1 <- ggplot() + 
  stat_density_2d(data = mer2_A[sample(nrow(mer2_A), 10000), ], aes(x = umap_1, y = umap_2, fill = stat(level)), 
                  geom = "polygon", bins = 10, contour = TRUE) +  
  common_fill_scale + 
  theme_void() + 
  ggtitle("Density Plot for A")

p2 <- ggplot() + 
  stat_density_2d(data = mer2_B[sample(nrow(mer2_B), 10000), ], aes(x = umap_1, y = umap_2, fill = stat(level)), 
                  geom = "polygon", bins = 10, contour = TRUE) +  
  common_fill_scale + 
  theme_void() + 
  ggtitle("Density Plot for B")

p3 <- ggplot() + 
  stat_density_2d(data = mer2_C[sample(nrow(mer2_C), 10000), ], aes(x = umap_1, y = umap_2, fill = stat(level)), 
                  geom = "polygon", bins = 10, contour = TRUE) +  
  common_fill_scale + 
  theme_void() + 
  ggtitle("Density Plot for C")


pdf('Figure 4C UMAP density TCR_three_time_points.pdf')
p1
p2
p3
dev.off()



########################################################################
############# FIG 4D - circle plot of average size and counts        ##
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

load ('TCR_norm_cell_bind.RData')
#### load the summarized expansion_contraction data from PP_Expansion
counts_df1 <- counts_A

counts_df1$paired <- substr (rownames (counts_df1),2,4)
counts_df1$paired <- gsub ("A", "",counts_df1$paired )
counts_df1$paired <- gsub ("B", "",counts_df1$paired )
counts_df1$paired <- gsub ("C", "",counts_df1$paired )

counts_df1$response <- ifelse (counts_df1$paired %in% c('14', '11','16', '17', '13', '23', '10', '21'), "BR", counts_df1$paired)
counts_df1$response <- ifelse (counts_df1$response%in% c('01', '07', '15', '18', '19',  '05', '08', '24', '06', '09', '20', '12', '02', '22','25', '03'), "GR", counts_df1$response)
counts_df1$response <- ifelse (counts_df1$response%in% c('04'), NA, counts_df1$response)

counts_df1A <- counts_df1

var1 <- 'BR'
var2 <- 'GR'
pval <- wilcox.test(counts_df1A$av_size_A[counts_df1A$response == var1],
                    counts_df1A$av_size_A[counts_df1A$response == var2], 
                    paired = FALSE)$p.value  
pvalues <- c(pval)


pdf('Figure 4D counts_average_size_BR_GR_baseline.pdf')
summary_stats <- counts_df1A %>%
  group_by(response) %>%
  summarise(median = median(av_size_A, na.rm = TRUE),
            se = sd(av_size_A, na.rm = TRUE) / sqrt(n()),  
            .groups = 'drop') 


ggplot(counts_df1A, aes(x = response, y = av_size_A, size = counts_A)) +
  geom_point(aes(fill = response, colour = response), 
             shape = 21, 
             alpha = 0.5, 
             position = position_jitter(width = 0.4), 
             stroke = 1.5) + 
  geom_errorbar(data = summary_stats %>% filter(response == "BR"),
                aes(x = response, y = median, ymin = median - se, ymax = median + se),
                width = 0.2, color = "darkred", size = 1) +  
  geom_point(data = summary_stats %>% filter(response == "BR"), 
             aes(x = response, y = median), 
             shape = 15, size = 4, color = "darkred") +  
  geom_errorbar(data = summary_stats %>% filter(response == "GR"),
                aes(x = response, y = median, ymin = median - se, ymax = median + se),
                width = 0.2, color = "darkgreen", size = 1) +  
  geom_point(data = summary_stats %>% filter(response == "GR"), 
             aes(x = response, y = median), 
             shape = 15, size = 4, color = "darkgreen") +  
  geom_text_repel(aes(label = paired), max.overlaps = Inf, box.padding = 0.5, point.padding = 0.3 , size = 3) +
  scale_fill_manual(values = c("BR" = "lightcoral", "GR" = "lightgreen")) +
  scale_colour_manual(values = c("BR" = "darkred", "GR" = "darkgreen")) +  
  scale_size(range = c(0, 20), name = "Counts") +  
  labs(x = "Response", y = "Average Size of Clones", title = "Average Size vs Response for Clones") +
  theme_minimal() +
  theme(legend.position = "right", 
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),  
        axis.ticks = element_line(color = "black"),  
        axis.line = element_line(color = "black")) +  
  scale_x_discrete(limits = c("BR", "GR")) +  
  scale_y_continuous(limits = c(0, 0.16)) +  
  geom_signif(comparisons = list(c(var1, var2)),
               annotations = sprintf("p = %.3f", pvalues),
               map_signif_level = FALSE,
               test = "wilcox.test", 
               step_increase = 0.1)  

dev.off()  



########################################################################
############# FIG 4E - ALLUVIAL PLOTS AB RESPONDERS ####################
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

load ('TCR_norm_cell_bind.RData')
#### load the summarized expansion_contraction data from PP_Expansion
load("metadata_git.RData")
#### load the metadata available at EGA
merA <- data.frame(read.table("clinical_data.txt", header=T, sep="\t", stringsAsFactors=F))
##### load clinical data in Supp Table 1

merA <- merA [!merA$patient_number=="4",]

nA_annot <- ifelse(nA == 1, 'singletonA', ifelse(nA >= 5, 'str_expansion', ifelse((nA >= 2) & (nA < 5), 'expansion',  'zero')))
rownames (nA_annot)
paste ("P", merA$patient_number, "A", sep="")
all(rownames (nA_annot)==paste ("P", merA$patient_number, "A", sep=""))

statusAB1[is.na(statusAB1)] <- 'zero'
statusAC1[is.na(statusAC1)] <- 'zero'

nA_GR <- nA_annot[merA$response=="GR",]
AB_GR <- statusAB1[merA$response=="GR",]
AC_GR <- statusAC1[merA$response=="GR",]

nA_BR <- nA_annot[merA$response=="BR",]
AB_BR <- statusAB1[merA$response=="BR",]
AC_BR <- statusAC1[merA$response=="BR",]

t_AB_GR <- table(nA_GR, AB_GR)
t_AB_BR <- table(nA_BR, AB_BR)

t_AC_GR <- table(nA_GR, AC_GR)
t_AC_BR <- table(nA_BR, AC_BR)
t_AC_all <- table(nA_annot, statusAC1)


pdf ('Figure 4E alluvial_plots_AB_GR.pdf')

two_way_table <- t_AB_GR
two_way_table <- two_way_table[, colnames (two_way_table)%in%c("contraction",   "lit_expansion",  "str_expansion")]
filtered_table <- as.data.frame(two_way_table)

colnames(filtered_table) <- c("nA_annot_Status", "statusAB1_Status", "Count")
colnames (filtered_table) <- c('StatusA', 'BehaviorAB', 'Counts')
filtered_table$StatusA <- gsub ('zero', 'absent', filtered_table$StatusA)
filtered_table$StatusA <- gsub ('str_expansion', 'strongly expanded', filtered_table$StatusA)
filtered_table$StatusA <- gsub ('singletonA', 'singleton', filtered_table$StatusA)
filtered_table$StatusA <- gsub ('expansion', 'expanded', filtered_table$StatusA)
filtered_table$BehaviorAB <- gsub ('str_expansion', 'strong expansion', filtered_table$BehaviorAB)
filtered_table$BehaviorAB <- gsub ('lit_expansion', 'expansion', filtered_table$BehaviorAB)

filtered_table$StatusA <- factor(filtered_table$StatusA, levels = c("strongly expanded", "expanded", "singleton", 'absent'))
filtered_table$BehaviorAB <- factor(filtered_table$BehaviorAB, levels = c("strong expansion", "expansion", "contraction"))


ggplot(filtered_table, aes(axis1 = StatusA, axis2 = BehaviorAB, y = Counts)) +
  geom_alluvium(aes(fill = BehaviorAB)) +
  geom_stratum(width = 0.2) + 
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3, 
            vjust = -0.5, hjust = 0.5) +  
  scale_x_discrete(limits = c("Status A", "Behavior AB"), expand = c(0.05, 0.05)) +
  scale_fill_manual(values = c("strong expansion" = "#0B6623",  "expansion" = "#00A86B",   "contraction" = "lightcoral")) +  
  theme_minimal() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 14, face = "bold", angle = 0, hjust = 0.5, margin = margin(t = 0, r = 5000)), 
        axis.text.y = element_text(size = 10, face = "bold"),  
        axis.title.x = element_text(size = 16, face = "bold", margin = margin(t = 10)),  
        axis.title.y = element_text(size = 16, face = "bold"),  
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5)) + 
  labs(title = "Clone dynamics - Good Responders", x = "Clone Statuses", y = "Count")

dev.off()

########################################################################
############# FIG 4F - ALLUVIAL PLOTS AB non responders ################
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

load ('TCR_norm_cell_bind.RData')
#### load the summarized expansion_contraction data from PP_Expansion
load("metadata_git.RData")
#### load the metadata available at EGA
merA <- data.frame(read.table("clinical_data.txt", header=T, sep="\t", stringsAsFactors=F))
##### load clinical data in Supp Table 1

merA <- merA [!merA$patient_number=="4",]

nA_annot <- ifelse(nA == 1, 'singletonA', ifelse(nA >= 5, 'str_expansion', ifelse((nA >= 2) & (nA < 5), 'expansion',  'zero')))
rownames (nA_annot)
paste ("P", merA$patient_number, "A", sep="")
all(rownames (nA_annot)==paste ("P", merA$patient_number, "A", sep=""))

statusAB1[is.na(statusAB1)] <- 'zero'
statusAC1[is.na(statusAC1)] <- 'zero'

nA_GR <- nA_annot[merA$response=="GR",]
AB_GR <- statusAB1[merA$response=="GR",]
AC_GR <- statusAC1[merA$response=="GR",]

nA_BR <- nA_annot[merA$response=="BR",]
AB_BR <- statusAB1[merA$response=="BR",]
AC_BR <- statusAC1[merA$response=="BR",]

t_AB_GR <- table(nA_GR, AB_GR)
t_AB_BR <- table(nA_BR, AB_BR)

t_AC_GR <- table(nA_GR, AC_GR)
t_AC_BR <- table(nA_BR, AC_BR)
t_AC_all <- table(nA_annot, statusAC1)

pdf ('Figure 4F alluvial_plots_AB_BR.pdf')
two_way_table <- t_AB_BR
two_way_table <- two_way_table[, colnames (two_way_table)%in%c("contraction",   "lit_expansion",  "str_expansion")]

filtered_table <- as.data.frame(two_way_table)

colnames(filtered_table) <- c("nA_annot_Status", "statusAB1_Status", "Count")
colnames (filtered_table) <- c('StatusA', 'BehaviorAB', 'Counts')
filtered_table$StatusA <- gsub ('zero', 'absent', filtered_table$StatusA)
filtered_table$StatusA <- gsub ('str_expansion', 'strongly expanded', filtered_table$StatusA)
filtered_table$StatusA <- gsub ('singletonA', 'singleton', filtered_table$StatusA)
filtered_table$StatusA <- gsub ('expansion', 'expanded', filtered_table$StatusA)
filtered_table$BehaviorAB <- gsub ('str_expansion', 'strong expansion', filtered_table$BehaviorAB)
filtered_table$BehaviorAB <- gsub ('lit_expansion', 'expansion', filtered_table$BehaviorAB)

filtered_table$StatusA <- factor(filtered_table$StatusA, levels = c("strongly expanded", "expanded", "singleton", 'absent'))
filtered_table$BehaviorAB <- factor(filtered_table$BehaviorAB, levels = c("strong expansion", "expansion", "contraction"))

library(ggplot2)
library(ggalluvial)  
ggplot(filtered_table, aes(axis1 = StatusA, axis2 = BehaviorAB, y = Counts)) +
  geom_alluvium(aes(fill = BehaviorAB)) +
  geom_stratum(width = 0.2) +  
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3, 
            vjust = -0.5, hjust = 0.5) + 
  scale_x_discrete(limits = c("Status A", "Behavior AB"), expand = c(0.05, 0.05)) +
  scale_fill_manual(values = c("strong expansion" = "#0B6623",  "expansion" = "#00A86B",   "contraction" = "lightcoral")) +  
  theme_minimal() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 14, face = "bold", angle = 0, hjust = 0.5, margin = margin(t = 0, r = 5000)), 
        axis.text.y = element_text(size = 10, face = "bold"),  
        axis.title.x = element_text(size = 16, face = "bold", margin = margin(t = 10)),  
        axis.title.y = element_text(size = 16, face = "bold"),  
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5)) +  
  labs(title = "Clone dynamics - Bad Responders", x = "Clone Statuses", y = "Count")

dev.off()




########################################################################
############# FIG 4 GHI - BOXPLOT CONTRACTION AB AC       ##
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

load ('TCR_norm_cell_bind.RData')
#### load the summarized expansion_contraction data from PP_Expansion
load("metadata_git.RData")
#### load the metadata available at EGA
merA <- data.frame(read.table("clinical_data.txt", header=T, sep='\t', row.names=))
##### load clinical data in Supp Table 1

merA <- merA [!merA$patient_number=="4",]

norm_EC_AB <- cbind (statusAB1_n_EC, merA$response, merA$patient_numb)
norm_EC_AC <- cbind (statusAC1_n_EC , merA$response, merA$patient_numb)
colnames (norm_EC_AB) <- c('contraction','lit_expansion', 'str_expansion',  'response', 'patient_number')
colnames (norm_EC_AC) <- c('contraction','lit_expansion', 'str_expansion',  'response', 'patient_number')
norm_EC_AB$category <- 'AB'
norm_EC_AC$category <- 'AC'


###### Figure 4G
combined_data <- rbind(norm_EC_AB, norm_EC_AC)
combined_data <- combined_data %>% filter(!is.na(response) & !is.na(contraction))
pval_AB <- wilcox.test(norm_EC_AB$contraction[norm_EC_AB$response == 'BR'],
                       norm_EC_AB$contraction[norm_EC_AB$response == 'GR'], 
                       paired = FALSE)$p.value

pval_AC <- wilcox.test(norm_EC_AC$contraction[norm_EC_AC$response == 'BR'],
                       norm_EC_AC$contraction[norm_EC_AC$response == 'GR'], 
                       paired = FALSE)$p.value
					   
outline_colors <- c("#A02629", '#7BA58C') 
box_fill_colors <- c('#E5A1A4', '#C3EAD6') 
point_fill_colors <- c('#C32E32', '#9ED3B9') 

					   
pdf	('Figure 4G boxplot_contraction_AB_AC.pdf')	   
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

###### Figure 4H
combined_data <- combined_data %>% filter(!is.na(response) & !is.na(lit_expansion))
pval_AB <- wilcox.test(norm_EC_AB$lit_expansion[norm_EC_AB$response == 'BR'],
                       norm_EC_AB$lit_expansion[norm_EC_AB$response == 'GR'], 
                       paired = FALSE)$p.value
pval_AC <- wilcox.test(norm_EC_AC$lit_expansion[norm_EC_AC$response == 'BR'],
                       norm_EC_AC$lit_expansion[norm_EC_AC$response == 'GR'], 
                       paired = FALSE)$p.value

outline_colors <- c("#A02629", '#7BA58C') 
box_fill_colors <- c('#E5A1A4', '#C3EAD6') 
point_fill_colors <- c('#C32E32', '#9ED3B9') 


pdf	('Figure 4H boxplot_lit_expansion_AB_AC.pdf')	 
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

###### Figure 4I
combined_data <- combined_data %>% filter(!is.na(response) & !is.na(str_expansion))
pval_AB <- wilcox.test(norm_EC_AB$str_expansion[norm_EC_AB$response == 'BR'],
                       norm_EC_AB$str_expansion[norm_EC_AB$response == 'GR'], 
                       paired = FALSE)$p.value
pval_AC <- wilcox.test(norm_EC_AC$str_expansion[norm_EC_AC$response == 'BR'],
                       norm_EC_AC$str_expansion[norm_EC_AC$response == 'GR'], 
                       paired = FALSE)$p.value

outline_colors <- c("#A02629", '#7BA58C') 
box_fill_colors <- c('#E5A1A4', '#C3EAD6') 
point_fill_colors <- c('#C32E32', '#9ED3B9') 


pdf	('Figure 4I boxplot_str_expansion_AB_AC.pdf')	 

ggplot(combined_data, aes(x = response, y = str_expansion)) + 
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
    strip.background = element_blank(),  # Remove the gray background
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

