
###################################################
####### SUPPLEMENTARY FIGURE 6 ####################
###################################################

#################################################################################################
###### Supplementary Figure 6 ########## Milo minor IDs SWARMPLOTS ##########
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

pdf('SuppFigure6a swarmplot final ids A C.pdf' , width = 7, height = 10)
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

pdf('SuppFigure6b swarmplot final ids A B.pdf', width = 7, height = 10)
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
  
pdf('SuppFigure6c  swarmplot minor ids B C.pdf', width = 7, height = 10)
print(p)
dev.off()

########################################
###### Supplementary Figure 7 ##########
########################################
###########################################################################
# Supplementary Figure 7 A Percentage Malignant cells versus clincal ######
###########################################################################

library(ggplot2)
library(ggrepel)
library (ggsignif)

met <- read.table("clinical_data.txt", header=T, sep='\t')
##### load clinical data in Supp Table 1


BR_GR <- t.test(met$perc_malig_cells[met$size_surgery2_cat == "BR"], 
                met$perc_malig_cells[met$size_surgery2_cat == "GR"])

BR_MET <- t.test(met$perc_malig_cells[met$size_surgery2_cat == "BR"], 
                 met$perc_malig_cells[met$size_surgery2_cat == "MET"])

GR_MET <- t.test(met$perc_malig_cells[met$size_surgery2_cat == "GR"], 
                 met$perc_malig_cells[met$size_surgery2_cat == "MET"])

pvalues <- c(BR_GR$p.value, GR_MET$p.value, BR_MET$p.value)

outline_colors <- c("#A02629", '#7BA58C', '#4B4B4B') 
box_fill_colors <- c('#E5A1A4', '#C3EAD6', '#7E7E7E' ) 
point_fill_colors <- c('#C32E32', '#9ED3B9', '#B0B0B0') 

pdf ('SuppFigure7A  clinico_pathological vs malignant cells.pdf')
ggplot(met, aes(x = size_surgery2_cat, y = perc_malig_cells)) +
  geom_boxplot(color = outline_colors, fill = box_fill_colors, alpha = 0.2, size = 1.5) +
  geom_point(aes(colour = factor(size_surgery2_cat), fill = factor(size_surgery2_cat)), 
            , shape = 21, size = 5, stroke = 1.2) +
  scale_fill_manual(values = point_fill_colors) +
  scale_colour_manual(values = outline_colors) +
  labs(title = "",
       x = "Clinicopathological Response",
       y = "Percentage Malignant Cells per Biopsy at EOT") +
  geom_text_repel(aes(label = patient_number), max.overlaps = Inf, box.padding = 0.5, point.padding = 0.3, segment.size = 0.2, segment.color = 'grey50') +
  theme(axis.line = element_line(colour = "black", size = 1, linetype = "solid"),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title = element_text(size = 18, face = "bold"),
        axis.ticks = element_line(size = 1),
        panel.background = element_rect(fill = 'white', color = 'black'),
        legend.position = "none") +
  annotate("text", x = Inf, y = Inf, 
           label = paste("p-values: BR vs GR =", sprintf("%.3f", BR_GR$p.value), 
                         ", GR vs MET =", sprintf("%.3f", GR_MET$p.value),
                         ", BR vs MET =", sprintf("%.3f", BR_MET$p.value)), 
           vjust = 2, hjust = 2, size = 5, color = "black", parse = FALSE) +
  geom_signif(comparisons = list(c("BR", "GR"), c("GR", "MET"), c("BR", "MET")),
              annotations = sprintf("p = %.3f", pvalues),
              y_position = 90,  
              map_signif_level = FALSE,
              test = "t.test", 
              step_increase = 0.07)
dev.off()


###########################################################################
# Supplementary Figure 7 B Percentage Malignant cells versus ROR     ######
###########################################################################
library(ggplot2)
library(ggrepel)
library (ggsignif)

met <- read.table("clinical_data.txt", header=T, sep='\t')
##### load clinical data in Supp Table 1

met <- met [!is.na(met$PROSIGNA_resp),]
ttest_result <- t.test(met$perc_malig_cells[met$PROSIGNA_resp == "BR"], met$perc_malig_cells[met$size_surgery2_cat == "GR"])[3]

outline_colors <- c("#A02629", '#7BA58C') 
box_fill_colors <- c('#E5A1A4', '#C3EAD6') 
point_fill_colors <- c('#C32E32', '#9ED3B9') 

pdf ('SuppFigure7B  ROR vs malignant cells.pdf')

ggplot(met, aes(x = PROSIGNA_resp, y = perc_malig_cells)) +
  geom_boxplot(color = outline_colors, fill = box_fill_colors, alpha = 0.2, size = 1.5) +
  geom_point(aes(colour=factor(PROSIGNA_resp), fill = factor(PROSIGNA_resp)), shape = 21, size = 5, stroke = 1.2) +
  scale_fill_manual(values= point_fill_colors) + 
  scale_colour_manual(values=outline_colors) + 
  labs(title = "",
       x = "Prosigna Response",
       y = "Percentage Malignant Cells per Biopsy at EOT") +
  geom_text_repel(aes(label = patient_number), max.overlaps = Inf, box.padding = 0.5, point.padding = 0.3, segment.size = 0.2, segment.color = 'grey50') +
  theme(axis.line = element_line(colour = "black", size = 1, linetype = "solid"),
        axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14),
        axis.title = element_text(size=18, face="bold"),
        axis.ticks = element_line(size = 1),
        panel.background = element_rect(fill = 'white', color = 'black'),
        legend.position = "none") +
  geom_signif(comparisons = list(c('BR', 'GR')),
			  y_position = 73,  
               annotations = sprintf("p = %.3f", ttest_result),
               map_signif_level = FALSE,
               test = "t.test", 
               step_increase = 0.07) 
dev.off()



######################## ##############################################
#### Supplementary Figure 7 C Naive B cells BOXPLOTS ##################################################
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

pdf ('SuppFigure7C Naive B cells GR BR.pdf')

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
#### Supplementary Figure 7 D Boxplot immune cluster vs monocytes #############
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

pdf ('SuppFigure7D boxplot_monocytes_TME_clusters.pdf')

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

######################## ##############################################
#### Supplementary Figure 7E CAF_ECM_BR_GR time ##############
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

setwd ("C:\\Xavier\\single_cells\\analysis\\Files3.0\\Results\\NeoLetEXe_final2\\supplementary_figures\\")
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

pdf ('SuppFigure7E CAF_ECM BR GR time.pdf')
print(combined_plot)
dev.off()

########################################
###### Supplementary Figure 8 ##########
########################################
###############################################
############# SUPP FIG 8A - Ratio TCR T cell ##
###############################################

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
#### load the scRepertoire output from  @@@@@@@@@@@@@@@@
load("metadata_git.RData")
#### load the metadata available at EGA
cell_n <- as.data.frame (table (metadata$orig.ident))
rownames (cell_n) <- cell_n$Var1
colnames (cell_n) <- c('orig.ident', 'cell_n')

TCRa <- as.data.frame (table (TCR$orig.ident[!is.na (TCR$CTnt)]))
rownames (TCRa) <- TCRa$Var1
colnames (TCRa) <- c('orig.ident', 'TCRa')

Tcell_n <- as.data.frame (table (TCR$orig.ident[!TCR$detail.ids%in%c('NK_cyto', 'CD4_low', 'GDT', 'UNA', 'NK_rest', 'Treg_min', 'CRTAM_GDT' )]))
rownames (Tcell_n) <- Tcell_n$Var1
colnames (Tcell_n) <- c('orig.ident', 'Tcell_n')

ME_n <- as.data.frame (table (metadata$orig.ident[!metadata$main.ids.cor%in%c('Malignant_cells','Normal_epithelial')]))
rownames (ME_n) <- ME_n$Var1
colnames (ME_n) <- c('orig.ident', 'ME_n')
mer <- merge (cell_n, TCRa, by.x='orig.ident', by.y='orig.ident', all.x=T)
mer <- merge (mer, Tcell_n, by.x='orig.ident', by.y='orig.ident', all.x=T)
mer <- merge (mer, ME_n, by.x='orig.ident', by.y='orig.ident', all.x=T)
mer $ratio <- mer$TCRa/mer$Tcell_n
mer $ratio_cell <- mer$TCRa/mer$cell_n
rownames (mer) <-  mer$orig.ident
mer$ratio[is.na(mer$ratio)] <- 0
mer <- mer [order(mer$ratio),]

mer$label_color <- ifelse(mer$ratio < 0.38, "red", "blue")
mer$color <- ifelse(mer$ratio < 0.38, "red", "blue") 

pdf ('SuppFigure8A ratio_TCR_MixCR.pdf')
ggplot(mer, aes(x = factor(rownames(mer), levels = rownames(mer)), y = ratio, fill = color)) +
  geom_bar(stat = "identity") +
  labs(title = "", x = "Samples", y = "Ratio of TCR / Tcell") +
  scale_fill_identity() +  
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6, 
                                   color = mer$label_color),  
        legend.position = "none",  
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank())  
  dev.off()
  
  
  
########################################################################
############# SUPP FIG 8B - barplot of the number of clones per cell state ##
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
#### load the scRepertoire output from  @@@@@@@@@@@@@@@@
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

pdf('SuppFigure8B barplot of _clone_frequency.pdf')
ggplot(plot_data, aes(x = detail.ids, y = clone_count, fill = frequency_category)) +
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


########################################################################
############# SUPP FIG 8C - barplot of the number of cells in clonotypes ##
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
#### load the scRepertoire output from  @@@@@@@@@@@@@@@@
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

pdf('SuppFigure8C barplot of cell frequency.pdf')

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
## Supplementary figure 8 Circos plot shared clone across states ##
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
#### load the Seurat_object_update with scRepertoire from @@@@@@@@@@@@@@@@

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

pdf ('SuppFigure8D circos_plot_clone state sharing.pdf')
chordDiagram(circles, self.link = 1, grid.col = grid.cols)
dev.off()


########################################
###### Supplementary Figure 9 ##########
########################################
########################################################################
############# SUPPLEMENTARY FIG 9B - ALLUVIAL PLOTS AC - responders   ##
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

pdf ('SuppFigure9B alluvial_plots_AC_GR.pdf')

two_way_table <- t_AC_GR
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
  scale_x_discrete(limits = c("Status A", "Behavior AC"), expand = c(0.05, 0.05)) +
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
############# SUPPLEMENTARY FIG 9C - ALLUVIAL PLOTS AC - non reponders #
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

pdf ('SuppFigure9C alluvial_plots_AC_BR.pdf')
two_way_table <- t_AC_BR
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
  scale_x_discrete(limits = c("Status A", "Behavior AC"), expand = c(0.05, 0.05)) +
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


###################################################################################
############# SUPPLEMANTARY FIG 9D - BOXPLOT CONTRACTION TCR normalization      ##
###################################################################################

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


load ('TCR_norm_TCR_bind.RData')
#### load the summarized expansion_contraction data from PP_Expansion_TCR
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

					   
pdf	('SuppFigure9D boxplot_contraction_AB_AC_TCRnorm.pdf')	   
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

###################################################################################
############# SUPPLEMANTARY FIG 9E - BOXPLOT EXPANSION TCR normalization      ####
###################################################################################

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


load ('TCR_norm_TCR_bind.RData')
#### load the summarized expansion_contraction data from PP_Expansion_TCR
load("metadata_git.RData")
#### load the metadata available at EGA
merA <- data.frame(read.table("clinical_data.txt", header=T, sep='\t', row.names=))
##### load clinical data in Supp Table 1

merA <- data.frame(read.table("clinical_data.txt", header=T, sep='\t', row.names=))
merA <- merA [!merA$patient_number=="4",]

norm_EC_AB <- cbind (statusAB1_n_EC, merA$response, merA$patient_numb)
norm_EC_AC <- cbind (statusAC1_n_EC , merA$response, merA$patient_numb)
colnames (norm_EC_AB) <- c('contraction','lit_expansion', 'str_expansion',  'response', 'patient_number')
colnames (norm_EC_AC) <- c('contraction','lit_expansion', 'str_expansion',  'response', 'patient_number')
norm_EC_AB$category <- 'AB'
norm_EC_AC$category <- 'AC'


combined_data <- rbind(norm_EC_AB, norm_EC_AC)
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


pdf	('SuppFigure9E boxplot_lit_expansion_AB_AC_TCRnorm.pdf')	 
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
    strip.background = element_blank(),  # Remove the gray background
    strip.text = element_text(size = 15, face = "bold")  # Optional: customize the text size and style
  ) +
  facet_wrap(~category) +  # Create separate panels for each category
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


###################################################################################
############# SUPPLEMANTARY FIG 9F - BOXPLOT STR TCR normalization      ####
###################################################################################

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

load ('TCR_norm_TCR_bind.RData')
#### load the summarized expansion_contraction data from PP_Expansion_TCR
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


combined_data <- rbind(norm_EC_AB, norm_EC_AC)
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


pdf	('SuppFigure9F boxplot_str_expansion_AB_AC_TCRnorm.pdf')	 

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


########################################
###### Supplementary Figure 10 ##########
########################################

########################################################################
############# FIG 10A - BOXPLOT OF PATIENT behaviour thee time points       ##
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

setwd ("C:\\Xavier\\single_cells\\analysis\\Files3.0\\Results\\NeoLetEXe_final2\\supplementary_figures\\")
load ('TCR_norm_cell_bind.RData')
#### load the summarized expansion_contraction data from PP_Expansion
load("metadata_git.RData")
#### load the metadata available at EGA
merA <- data.frame(read.table("clinical_data.txt", header=T, sep='\t', row.names=))
##### load clinical data in Supp Table 1


merA <- merA [!merA$patient_number=="4",]

all (rownames (statusABC1_n_EC) ==  merA$patients)
response <- merA$response
patient_number <- merA$patient_numb

norm_EC <- statusABC1_n_EC 
norm_EC$response <- response
norm_EC$patient_number <- patient_number
colnames (norm_EC) <- gsub ('WOs_ABC_' , '', colnames (norm_EC))

outline_colors <- c("#A02629", '#7BA58C') # Box outline colors
box_fill_colors <- c('#E5A1A4', '#C3EAD6') # Box fill colors
point_fill_colors <- c('#C32E32', '#9ED3B9') # Point fill colors


pdf('SuppFigure10A dynamic_TCR_contraction_contraction.pdf')
var1 <- 'BR'
var2 <- 'GR'
pval <- wilcox.test(norm_EC$contraction_contraction[norm_EC$response == var1],
                    norm_EC$contraction_contraction[norm_EC$response == var2], 
                    paired = FALSE)$p.value  
pvalues <- c(pval)
ggplot(norm_EC, aes(x = response, y = contraction_contraction)) + 
  geom_boxplot(aes(fill = response), color = outline_colors, fill = box_fill_colors, alpha = 0.2, size = 1.5) + 
  geom_point(aes(colour = response, fill = response), shape = 21, size = 5) +  
  scale_fill_manual(values= point_fill_colors) + 
  scale_colour_manual(values=outline_colors) + 
  ylab('Fraction of contractingAB - contractingBC clones') + 
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
############# SUPPLEMENTARY FIG 10B - BOXPLOT OF PATIENT behaviour thee time points       ##
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

setwd ("C:\\Xavier\\single_cells\\analysis\\Files3.0\\Results\\NeoLetEXe_final2\\supplementary_figures\\")
load ('TCR_norm_cell_bind.RData')
#### load the summarized expansion_contraction data from PP_Expansion
load("metadata_git.RData")
#### load the metadata available at EGA
merA <- data.frame(read.table("clinical_data.txt", header=T, sep='\t', row.names=))
##### load clinical data in Supp Table 1

merA <- merA [!merA$patient_number=="4",]

all (rownames (statusABC1_n_EC) ==  merA$patients)
response <- merA$response
patient_number <- merA$patient_numb

norm_EC <- statusABC1_n_EC 
norm_EC$response <- response
norm_EC$patient_number <- patient_number
colnames (norm_EC) <- gsub ('WOs_ABC_' , '', colnames (norm_EC))

outline_colors <- c("#A02629", '#7BA58C') # Box outline colors
box_fill_colors <- c('#E5A1A4', '#C3EAD6') # Box fill colors
point_fill_colors <- c('#C32E32', '#9ED3B9') # Point fill colors



pdf('SuppFigure10B dynamic_TCR_str_expansion_str_expansion.pdf')				   
var1 <- 'BR'
var2 <- 'GR'
pval <- wilcox.test(norm_EC$str_expansion_str_expansion[norm_EC$response == var1],
                    norm_EC$str_expansion_str_expansion[norm_EC$response == var2], 
                    paired = FALSE)$p.value  
pvalues <- c(pval)


ggplot(norm_EC, aes(x = response, y = str_expansion_str_expansion)) + 
  geom_boxplot(aes(fill = response), color = outline_colors, fill = box_fill_colors, alpha = 0.2, size = 1.5) + 
  geom_point(aes(colour = response, fill = response), shape = 21, size = 5) +  # Fixed parentheses here
  scale_fill_manual(values= point_fill_colors) + 
  scale_colour_manual(values=outline_colors) + 
  ylab('Fraction of strong_expnasionAB - strong_expnasionBC clones') + 
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



########################################################################
############# Supplementary Figure 10 C - singletons across time      ##
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


setwd ("C:\\Xavier\\single_cells\\analysis\\Files3.0\\Results\\NeoLetEXe_final2\\supplementary_figures\\")
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

pdf ('SuppFigure10C singletons_GR.pdf')

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
############# Supplementary Figure 10 D - singletons across time      ##
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


setwd ("C:\\Xavier\\single_cells\\analysis\\Files3.0\\Results\\NeoLetEXe_final2\\supplementary_figures\\")
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

pdf ('SuppFigure10D singletons_BR.pdf')
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
