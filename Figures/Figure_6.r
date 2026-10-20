###################################################
####### FIGURE 6 ##################################
###################################################

################################################################################
############# FIG 6A - UMAP with distla CRE     ################################
################################################################################

### PP_snakemake_scafe

###############################################################################################
############# FIG 6B - SWARMPLOT PROXIMAL DISTAL MALIGNANT     ################################
###############################################################################################
library(EnhancedVolcano)
library(ggrepel)
library(ggplot2)
library(miloR)
library(plyr)
library(dplyr)
library(future)
library(stringi)
library(tidyverse)

d1 <- read.table("SCAFE/enichment_proximal_ids/improved_main_ids__Malignant_cells__proximal.tsv", sep = '\t', header = TRUE, fill = TRUE)
d2 <- read.table("SCAFE/enrichment_distals_ids/improved_main_ids__Malignant_cells.tsv", sep = '\t', header = TRUE, fill = TRUE)
### available from PP_snakemake_scafe

d2$oddsRatio <- (-d2$oddsRatio)
min_significant_ratio= 0.5
min_datasets = 10
 
  d <- rbind(d1[,colnames(d1) %in% c('collection', 'pValueLog', 'oddsRatio')], 
             d2[,colnames(d2) %in% c('collection', 'pValueLog', 'oddsRatio')])
  
    
  collection_count_tot <- as.data.frame(table(d$collection) / 2)
  collection_count_tot <- collection_count_tot[collection_count_tot$Freq > min_datasets, ]
  
  d <- d[!d$oddsRatio %in% c("Inf", "-Inf"), ]
  d$pValue <- 10^(-d$pValueLog)
  
  collection_count_plus <- as.data.frame(table(d$collection[d$oddsRatio > 2 & d$pValue < 0.05]))
  colnames(collection_count_plus)[2] <- 'plus'
  
  collection_count_minus <- as.data.frame(table(d$collection[d$oddsRatio < -2 & d$pValue < 0.05]))
  colnames(collection_count_minus)[2] <- 'minus'

  mer <- merge(collection_count_tot, collection_count_plus, by.x = 'Var1', by.y = 'Var1', all.x = TRUE)
  mer1 <- merge(mer, collection_count_minus, by.x = 'Var1', by.y = 'Var1', all.x = TRUE)
  
  mer1$plus <- ifelse(is.na(mer1$plus), 0, mer1$plus)
  mer1$minus <- ifelse(is.na(mer1$minus), 0, mer1$minus)
  
  mer1$ratio_plus <- mer1$plus / mer1$Freq
  mer1$ratio_minus <- mer1$minus / mer1$Freq
 
selected <- mer1 %>%
  arrange(desc(ratio_minus)) %>%
  slice(1:5) %>%
  bind_rows(
    mer1 %>%
      arrange(desc(ratio_plus)) %>%
      slice(1:5))
selected_collections <- selected$Var1

  d1 <- d[d$collection %in% selected_collections, ]
  
  colnames(d1)[4] <- 'SpatialFDR'
  colnames(d1)[3] <- 'logFC'
  d1$Nhood <- rownames(d1)
  
  vector_r <- rev(d1 %>% group_by(collection) %>% summarize(avg_logFC = mean(logFC)) %>% arrange(avg_logFC))
  d1$collection <- factor(d1$collection, levels = rev(vector_r$collection))
  
  n <- length(mer1$Var1[mer1$ratio_minus > min_significant_ratio])
  total_levels <- length(unique(d1$collection))
  line_position = total_levels - n + 0.5
  
  mer1_filtered <- mer1[mer1$Var1 %in% vector_r$collection, ]
  mer1_filtered$Var1 <- factor(mer1_filtered$Var1, levels = vector_r$collection)
  mer1_filtered <- mer1_filtered[order(mer1_filtered$Var1), ]

  p <- plotDAbeeswarm(d1, group.by = "collection", alpha = 0.05) +
    ggtitle('plot') +
    scale_colour_gradient2(
      high = "dark blue",  
      low = "blue",    
      na.value = "#cbcbcb", 
      space = "Lab"
    ) + 
    xlab("") +
    theme(
	  plot.title = element_text(size = 10, face = "bold"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.y = element_text(size = 20),
	  axis.text.x = element_text(size = 20)
    )
  
  p$data$logFC_color <- ifelse(p$data$logFC_color >= 2, 2, p$data$logFC_color)
  p$data$logFC_color <- ifelse(p$data$logFC_color < -2, -2, p$data$logFC_color)
  p$data$logFC_color <- ifelse((p$data$logFC_color == -2 | p$data$logFC_color == 2), p$data$logFC_color, NA)
  p$data <- p$data[order(!is.na(p$data$logFC_color), p$data$logFC_color), ]
  
  
  p <- p + mapply(function(x, plus_label, minus_label) { 
    list(
      annotate("text", x = x, y = max(d1$logFC) + 6, label = plus_label, hjust = 1, size = 4),
      annotate("text", x = x, y = min(d1$logFC) - 6, label = minus_label, hjust = 1, size = 4)
    ) 
  }, seq_along(mer1_filtered$Var1), rev(paste("(", round(mer1_filtered$ratio_plus,3)*100, ")", sep="")), rev(paste("(",round(mer1_filtered$ratio_minus,3)*100, ")", sep="")))

p1 <- plotDAbeeswarm(d1, group.by = "collection", alpha = 0.05) +
  ggtitle('Plot with max logFC') +
  scale_colour_gradient2(
    high = "dark blue",  
    low = "blue",    
    na.value = "#cbcbcb", 
    space = "Lab"
  ) + 
  ylim(0, max(d1$logFC) + 3) + 
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

p2 <- plotDAbeeswarm(d1, group.by = "collection", alpha = 0.05) +
  ggtitle('Plot with min logFC') +
  scale_colour_gradient2(
    high = "dark blue",  
    low = "blue",    
    na.value = "#cbcbcb", 
    space = "Lab"
  ) + 
  ylim(0, min(d1$logFC) - 6) +  
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

pdf('Figure 6B Malignant_proximal_vs_distal_double.pdf')
print(p)
print(p1)
print(p2)
dev.off()

###############################################################################################
############# FIG 6C - SWARMPLOT PROXIMAL DISTAL Tcells     ################################
###############################################################################################

library(EnhancedVolcano)
library(ggrepel)
library(ggplot2)
library(miloR)
library(plyr)
library(dplyr)
library(future)
library(stringi)
library(tidyverse)

d1 <- read.table("SCAFE/enichment_proximal_ids/improved_main_ids__T_cells__proximal.tsv", sep = '\t', header = TRUE, fill = TRUE)
d2 <- read.table("SCAFE/enrichment_distals_ids/improved_main_ids__T_cells.tsv", sep = '\t', header = TRUE, fill = TRUE)
### available from PP_snakemake_scafe

d2$oddsRatio <- (-d2$oddsRatio)

min_significant_ratio= 0.5
min_datasets = 10
 
  d <- rbind(d1[,colnames(d1) %in% c('collection', 'pValueLog', 'oddsRatio')], 
             d2[,colnames(d2) %in% c('collection', 'pValueLog', 'oddsRatio')])
  
    
  collection_count_tot <- as.data.frame(table(d$collection) / 2)
  collection_count_tot <- collection_count_tot[collection_count_tot$Freq > min_datasets, ]
  
  d <- d[!d$oddsRatio %in% c("Inf", "-Inf"), ]
  d$pValue <- 10^(-d$pValueLog)
  
  collection_count_plus <- as.data.frame(table(d$collection[d$oddsRatio > 2 & d$pValue < 0.05]))
  colnames(collection_count_plus)[2] <- 'plus'
  
  collection_count_minus <- as.data.frame(table(d$collection[d$oddsRatio < -2 & d$pValue < 0.05]))
  colnames(collection_count_minus)[2] <- 'minus'

  mer <- merge(collection_count_tot, collection_count_plus, by.x = 'Var1', by.y = 'Var1', all.x = TRUE)
  mer1 <- merge(mer, collection_count_minus, by.x = 'Var1', by.y = 'Var1', all.x = TRUE)
  
  mer1$plus <- ifelse(is.na(mer1$plus), 0, mer1$plus)
  mer1$minus <- ifelse(is.na(mer1$minus), 0, mer1$minus)
  
  mer1$ratio_plus <- mer1$plus / mer1$Freq
  mer1$ratio_minus <- mer1$minus / mer1$Freq 

top_minus <- mer1 %>% arrange(desc(ratio_minus)) %>% slice_head(n = 5)
remaining_data <- mer1 %>% anti_join(top_minus, by = "Var1")
top_plus <- remaining_data %>% arrange(desc(ratio_plus)) %>% slice_head(n = 5)
selected <- bind_rows(top_minus, top_plus)

selected_collections <- selected$Var1

  d1 <- d[d$collection %in% selected_collections, ]
  
  colnames(d1)[4] <- 'SpatialFDR'
  colnames(d1)[3] <- 'logFC'
  d1$Nhood <- rownames(d1)
  
  vector_r <- rev(d1 %>% group_by(collection) %>% summarize(avg_logFC = mean(logFC)) %>% arrange(avg_logFC))
  d1$collection <- factor(d1$collection, levels = rev(vector_r$collection))
  
  n <- length(mer1$Var1[mer1$ratio_minus > min_significant_ratio])
  total_levels <- length(unique(d1$collection))
  line_position = total_levels - n + 0.5
  
  mer1_filtered <- mer1[mer1$Var1 %in% vector_r$collection, ]
  mer1_filtered$Var1 <- factor(mer1_filtered$Var1, levels = vector_r$collection)
  mer1_filtered <- mer1_filtered[order(mer1_filtered$Var1), ]

  p <- plotDAbeeswarm(d1, group.by = "collection", alpha = 0.05) +
    ggtitle('plot') +
    scale_colour_gradient2(
      high = "purple",  
      low = "orchid3", 
      na.value = "#cbcbcb", 
      space = "Lab"
    ) + 
    xlab("") +
    theme(
	  plot.title = element_text(size = 10, face = "bold"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.y = element_text(size = 20),
	  axis.text.x = element_text(size = 20)
    )
  
  p$data$logFC_color <- ifelse(p$data$logFC_color >= 2, 2, p$data$logFC_color)
  p$data$logFC_color <- ifelse(p$data$logFC_color < -2, -2, p$data$logFC_color)
  p$data$logFC_color <- ifelse((p$data$logFC_color == -2 | p$data$logFC_color == 2), p$data$logFC_color, NA)
  p$data <- p$data[order(!is.na(p$data$logFC_color), p$data$logFC_color), ]
  
  
  p <- p + mapply(function(x, plus_label, minus_label) { 
    list(
      annotate("text", x = x, y = max(d1$logFC) + 1, label = plus_label, hjust = 1, size = 4),
      annotate("text", x = x, y = min(d1$logFC) - 1, label = minus_label, hjust = 1, size = 4)
    ) 
  }, seq_along(mer1_filtered$Var1), rev(paste("(", round(mer1_filtered$ratio_plus,3)*100, ")", sep="")), rev(paste("(",round(mer1_filtered$ratio_minus,3)*100, ")", sep="")))


p1 <- plotDAbeeswarm(d1, group.by = "collection", alpha = 0.05) +
  ggtitle('Plot with max logFC') +
  scale_colour_gradient2(
      high = "purple",  
      low = "orchid3",  
    na.value = "#cbcbcb", 
    space = "Lab"
  ) + 
  ylim(0, max(d1$logFC) + 3) + 
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

p2 <- plotDAbeeswarm(d1, group.by = "collection", alpha = 0.05) +
  ggtitle('Plot with min logFC') +
  scale_colour_gradient2(
      high = "purple",  
      low = "orchid3",  
    na.value = "#cbcbcb", 
    space = "Lab"
  ) + 
  ylim(0, min(d1$logFC) - 6) +  
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

pdf('Figure 6C Tcell_proximal_vs_distal_double.pdf')
print (p)
print(p1)
print(p2)
dev.off()

 
 
###############################################################################################
############# FIG 6D - SWARMPLOT PROXIMAL DISTAL Myeloid     ################################
###############################################################################################

library(EnhancedVolcano)
library(ggrepel)
library(ggplot2)
library(miloR)
library(plyr)
library(dplyr)
library(future)
library(stringi)
library(tidyverse)

d1 <- read.table("SCAFE/enichment_proximal_ids/improved_main_ids__Myeloid__proximal.tsv", sep = '\t', header = TRUE, fill = TRUE)
d2 <- read.table("SCAFE/enrichment_distals_ids/improved_main_ids__Myeloid.tsv", sep = '\t', header = TRUE, fill = TRUE)
### available from PP_snakemake_scafe

d2$oddsRatio <- (-d2$oddsRatio)

min_significant_ratio= 0.5
min_datasets = 10
 
  d <- rbind(d1[,colnames(d1) %in% c('collection', 'pValueLog', 'oddsRatio')], 
             d2[,colnames(d2) %in% c('collection', 'pValueLog', 'oddsRatio')])
  
  collection_count_tot <- as.data.frame(table(d$collection) / 2)
  collection_count_tot <- collection_count_tot[collection_count_tot$Freq > min_datasets, ]
  
  d <- d[!d$oddsRatio %in% c("Inf", "-Inf"), ]
  d$pValue <- 10^(-d$pValueLog)
  
  collection_count_plus <- as.data.frame(table(d$collection[d$oddsRatio > 2 & d$pValue < 0.05]))
  colnames(collection_count_plus)[2] <- 'plus'
  
  collection_count_minus <- as.data.frame(table(d$collection[d$oddsRatio < -2 & d$pValue < 0.05]))
  colnames(collection_count_minus)[2] <- 'minus'

  mer <- merge(collection_count_tot, collection_count_plus, by.x = 'Var1', by.y = 'Var1', all.x = TRUE)
  mer1 <- merge(mer, collection_count_minus, by.x = 'Var1', by.y = 'Var1', all.x = TRUE)
  
  mer1$plus <- ifelse(is.na(mer1$plus), 0, mer1$plus)
  mer1$minus <- ifelse(is.na(mer1$minus), 0, mer1$minus)
  
  mer1$ratio_plus <- mer1$plus / mer1$Freq
  mer1$ratio_minus <- mer1$minus / mer1$Freq
    

top_minus <- mer1 %>% arrange(desc(ratio_minus)) %>% slice_head(n = 5)
remaining_data <- mer1 %>% anti_join(top_minus, by = "Var1")
top_plus <- remaining_data %>% arrange(desc(ratio_plus)) %>% slice_head(n = 5)
selected <- bind_rows(top_minus, top_plus)

selected_collections <- selected$Var1

  d1 <- d[d$collection %in% selected_collections, ]
  
  colnames(d1)[4] <- 'SpatialFDR'
  colnames(d1)[3] <- 'logFC'
  d1$Nhood <- rownames(d1)
  
  vector_r <- rev(d1 %>% group_by(collection) %>% summarize(avg_logFC = mean(logFC)) %>% arrange(avg_logFC))
  d1$collection <- factor(d1$collection, levels = rev(vector_r$collection))
  
  n <- length(mer1$Var1[mer1$ratio_minus > min_significant_ratio])
  total_levels <- length(unique(d1$collection))
  line_position = total_levels - n + 0.5
  
  mer1_filtered <- mer1[mer1$Var1 %in% vector_r$collection, ]
  mer1_filtered$Var1 <- factor(mer1_filtered$Var1, levels = vector_r$collection)
  mer1_filtered <- mer1_filtered[order(mer1_filtered$Var1), ]

  p <- plotDAbeeswarm(d1, group.by = "collection", alpha = 0.05) +
    ggtitle('plot') +
    scale_colour_gradient2(
      high = "deeppink3",  
      low = "hotpink",   
      na.value = "#cbcbcb", 
      space = "Lab"
    ) + 
    xlab("") +
    theme(
	  plot.title = element_text(size = 10, face = "bold"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.y = element_text(size = 20),
	  axis.text.x = element_text(size = 20)
    )
  
  p$data$logFC_color <- ifelse(p$data$logFC_color >= 2, 2, p$data$logFC_color)
  p$data$logFC_color <- ifelse(p$data$logFC_color < -2, -2, p$data$logFC_color)
  p$data$logFC_color <- ifelse((p$data$logFC_color == -2 | p$data$logFC_color == 2), p$data$logFC_color, NA)
  p$data <- p$data[order(!is.na(p$data$logFC_color), p$data$logFC_color), ]
  
  
  p <- p + mapply(function(x, plus_label, minus_label) { 
    list(
      annotate("text", x = x, y = max(d1$logFC) + 1, label = plus_label, hjust = 1, size = 4),
      annotate("text", x = x, y = min(d1$logFC) - 1, label = minus_label, hjust = 1, size = 4)
    ) 
  }, seq_along(mer1_filtered$Var1), rev(paste("(", round(mer1_filtered$ratio_plus,3)*100, ")", sep="")), rev(paste("(",round(mer1_filtered$ratio_minus,3)*100, ")", sep="")))

p1 <- plotDAbeeswarm(d1, group.by = "collection", alpha = 0.05) +
  ggtitle('Plot with max logFC') +
  scale_colour_gradient2(
      high = "deeppink3",  
      low = "hotpink",    
    na.value = "#cbcbcb", 
    space = "Lab"
  ) + 
  ylim(0, max(d1$logFC) + 3) + 
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

p2 <- plotDAbeeswarm(d1, group.by = "collection", alpha = 0.05) +
  ggtitle('Plot with min logFC') +
  scale_colour_gradient2(
      high = "deeppink3",  
      low = "hotpink",   ,  
    na.value = "#cbcbcb", 
    space = "Lab"
  ) + 
  ylim(0, min(d1$logFC) - 6) + 
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


pdf('Figure 6D Myeloid_proximal_vs_distal_double.pdf')
print (p)
print(p1)
print(p2)

dev.off()




###############################################################################################
############# FIG 6E - SWARMPLOT PROXIMAL DISTAL FIBROBLAST     ################################
###############################################################################################

library(EnhancedVolcano)
library(ggrepel)
library(ggplot2)
library(miloR)
library(plyr)
library(dplyr)
library(future)
library(stringi)
library(tidyverse)

d1 <- read.table("C:/Xavier/single_cells/analysis/Files3.0/Results/NeoLetEXe_final2/SCAFE/enichment_proximal_ids/improved_main_ids__Fibroblasts__proximal.tsv", sep = '\t', header = TRUE, fill = TRUE)
d2 <- read.table("C:/Xavier/single_cells/analysis/Files3.0/Results/NeoLetEXe_final2/SCAFE/enrichment_distals_ids/improved_main_ids__Fibroblasts.tsv", sep = '\t', header = TRUE, fill = TRUE)
### available from PP_snakemake_scafe

d2$oddsRatio <- (-d2$oddsRatio)

min_significant_ratio= 0.5
min_datasets = 10
 
  d <- rbind(d1[,colnames(d1) %in% c('collection', 'pValueLog', 'oddsRatio')], 
             d2[,colnames(d2) %in% c('collection', 'pValueLog', 'oddsRatio')])
  
    
  collection_count_tot <- as.data.frame(table(d$collection) / 2)
  collection_count_tot <- collection_count_tot[collection_count_tot$Freq > min_datasets, ]
  
  d <- d[!d$oddsRatio %in% c("Inf", "-Inf"), ]
  d$pValue <- 10^(-d$pValueLog)
  
  collection_count_plus <- as.data.frame(table(d$collection[d$oddsRatio > 2 & d$pValue < 0.05]))
  colnames(collection_count_plus)[2] <- 'plus'
  
  collection_count_minus <- as.data.frame(table(d$collection[d$oddsRatio < -2 & d$pValue < 0.05]))
  colnames(collection_count_minus)[2] <- 'minus'

  mer <- merge(collection_count_tot, collection_count_plus, by.x = 'Var1', by.y = 'Var1', all.x = TRUE)
  mer1 <- merge(mer, collection_count_minus, by.x = 'Var1', by.y = 'Var1', all.x = TRUE)
  
  mer1$plus <- ifelse(is.na(mer1$plus), 0, mer1$plus)
  mer1$minus <- ifelse(is.na(mer1$minus), 0, mer1$minus)
  
  mer1$ratio_plus <- mer1$plus / mer1$Freq
  mer1$ratio_minus <- mer1$minus / mer1$Freq
   
top_minus <- mer1 %>% arrange(desc(ratio_minus)) %>% slice_head(n = 5)
remaining_data <- mer1 %>% anti_join(top_minus, by = "Var1")
top_plus <- remaining_data %>% arrange(desc(ratio_plus)) %>% slice_head(n = 10)
selected <- bind_rows( top_plus)

selected_collections <- selected$Var1

  d1 <- d[d$collection %in% selected_collections, ]
  
  colnames(d1)[4] <- 'SpatialFDR'
  colnames(d1)[3] <- 'logFC'
  d1$Nhood <- rownames(d1)
  
  vector_r <- rev(d1 %>% group_by(collection) %>% summarize(avg_logFC = mean(logFC)) %>% arrange(avg_logFC))
  d1$collection <- factor(d1$collection, levels = rev(vector_r$collection))
  
  n <- length(mer1$Var1[mer1$ratio_minus > min_significant_ratio])
  total_levels <- length(unique(d1$collection))
  line_position = total_levels - n + 0.5
  
  mer1_filtered <- mer1[mer1$Var1 %in% vector_r$collection, ]
  mer1_filtered$Var1 <- factor(mer1_filtered$Var1, levels = vector_r$collection)
  mer1_filtered <- mer1_filtered[order(mer1_filtered$Var1), ]

  p <- plotDAbeeswarm(d1, group.by = "collection", alpha = 0.05) +
    ggtitle('plot') +
    scale_colour_gradient2(
      high = "darkorange3",  
      low = "orange", ,   
      na.value = "#cbcbcb", 
      space = "Lab"
    ) + 
    xlab("") +
    theme(
	  plot.title = element_text(size = 10, face = "bold"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.y = element_text(size = 20),
	  axis.text.x = element_text(size = 20)
    )
  
  p$data$logFC_color <- ifelse(p$data$logFC_color >= 2, 2, p$data$logFC_color)
  p$data$logFC_color <- ifelse(p$data$logFC_color < -1.4, -1.4, p$data$logFC_color)
  p$data$logFC_color <- ifelse((p$data$logFC_color == -1.4 | p$data$logFC_color == 2), p$data$logFC_color, NA)
  
  p <- p + mapply(function(x, plus_label, minus_label) { 
    list(
      annotate("text", x = x, y = max(d1$logFC) + 1, label = plus_label, hjust = 1, size = 4),
      annotate("text", x = x, y = min(d1$logFC) - 1, label = minus_label, hjust = 1, size = 4)
    ) 
  }, seq_along(mer1_filtered$Var1), rev(paste("(", round(mer1_filtered$ratio_plus,3)*100, ")", sep="")), rev(paste("(",round(mer1_filtered$ratio_minus,3)*100, ")", sep="")))


p1 <- plotDAbeeswarm(d1, group.by = "collection", alpha = 0.05) +
  ggtitle('Plot with max logFC') +
  scale_colour_gradient2(
      high = "darkorange3",  
      low = "orange",  
    na.value = "#cbcbcb", 
    space = "Lab"
  ) + 
  ylim(0, max(d1$logFC) + 3) + 
  xlab("") +
  theme(
    plot.title = element_text(size = 10, face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 20),
    axis.text.x = element_text(size = 20)
  )

  p$data$logFC_color <- ifelse(p$data$logFC_color >= 2, 2, p$data$logFC_color)
  p$data$logFC_color <- ifelse(p$data$logFC_color < -1.4, -1.4, p$data$logFC_color)
  p$data$logFC_color <- ifelse((p$data$logFC_color == -1.4 | p$data$logFC_color == 2), p$data$logFC_color, NA)
  
p2 <- plotDAbeeswarm(d1, group.by = "collection", alpha = 0.05) +
  ggtitle('Plot with min logFC') +
  scale_colour_gradient2(
      high = "darkorange3",  
      low = "orange", 
    na.value = "#cbcbcb", 
    space = "Lab"
  ) + 
  ylim(0, min(d1$logFC) - 6) +  
  theme(
    plot.title = element_text(size = 10, face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 20),
    axis.text.x = element_text(size = 20)
  )

  p$data$logFC_color <- ifelse(p$data$logFC_color >= 2, 2, p$data$logFC_color)
  p$data$logFC_color <- ifelse(p$data$logFC_color < -1.4, -1.4, p$data$logFC_color)
  p$data$logFC_color <- ifelse((p$data$logFC_color == -1.4 | p$data$logFC_color == 2), p$data$logFC_color, NA)
  
  
pdf('Figure 6E Fibroblast_proximal_vs_distal_double.pdf')
print (p)
print(p1)
print(p2)
dev.off()



###############################################################################################
############# FIG 6 GHJ - SWARMPLOT - A vs B     ################################
###############################################################################################

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

d2 <- read.table("SCAFE/malignant_timepoint_analysis/Malignant_baseline_vs_2mnd__distal_cre__unibind.tsv", sep = '\t', header = TRUE, fill = TRUE)
d1 <- read.table("SCAFE/malignant_timepoint_analysis/Malignant_2mnd_vs_baseline__distal_cre__unibind.tsv", sep = '\t', header = TRUE, fill = TRUE)

#d2 <- read.table("SCAFE/malignant_timepoint_analysis/response_BR_A_vs_B__distal_cre__unibind.tsv", sep = '\t', header = TRUE, fill = TRUE)
#d1 <- read.table("SCAFE/malignant_timepoint_analysis/response_BR_B_vs_A__distal_cre__unibind.tsv", sep = '\t', header = TRUE, fill = TRUE)

#d2 <- read.table("SCAFE/malignant_timepoint_analysis/response_GR_A_vs_B__distal_cre__unibind.tsv", sep = '\t', header = TRUE, fill = TRUE)
#d1 <- read.table("SCAFE/malignant_timepoint_analysis/response_GR_B_vs_A__distal_cre__unibind.tsv", sep = '\t', header = TRUE, fill = TRUE)

d2$oddsRatio <- (-d2$oddsRatio)

  d <- rbind(d1[,colnames(d1) %in% c('collection', 'pValueLog', 'oddsRatio')], 
             d2[,colnames(d2) %in% c('collection', 'pValueLog', 'oddsRatio')])
  
  max <- max(d$oddsRatio[!d$oddsRatio %in% c('Inf')])
  min <- min(d$oddsRatio[!d$oddsRatio %in% c('-Inf')])
  
  d$pValue <- 10^(-d$pValueLog)

filtered_data <- subset(d, (pValue < 0.05 & oddsRatio > 1.5) | pValue < 0.05 & oddsRatio < -1.5 )
collection_count <- table(filtered_data$collection)


print(head (colnames(d)))
collection_count_plus <- as.data.frame(table (d$collection[d$oddsRatio > 1.5 & d$pValue < 0.05])) ; colnames (collection_count_plus)[2] <- 'plus'
collection_count_minus <- as.data.frame(table (d$collection[!d$oddsRatio > -1.5 & d$pValue < 0.05])) ;colnames (collection_count_minus)[2] <- 'minus'

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
  
  p$data$logFC_color <- ifelse(p$data$logFC_color >= 1.5, 2, p$data$logFC_color)
  p$data$logFC_color <- ifelse(p$data$logFC_color < -1.5, -2, p$data$logFC_color)
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

  p1$data$logFC_color <- ifelse(p1$data$logFC_color >= 1.5, 2, p1$data$logFC_color)
  p1$data$logFC_color <- ifelse(p1$data$logFC_color < -1.5, -2, p1$data$logFC_color)
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

  p2$data$logFC_color <- ifelse(p2$data$logFC_color >= 1.5, 2, p2$data$logFC_color)
  p2$data$logFC_color <- ifelse(p2$data$logFC_color < -1.5, -2, p2$data$logFC_color)
  p2$data$logFC_color <- ifelse((p2$data$logFC_color == -2 | p2$data$logFC_color == 2), p2$data$logFC_color, NA)
  p2$data <- p2$data[order(!is.na(p2$data$logFC_color), p2$data$logFC_color), ]

p2 <- p2 + mapply(function(x, plus_label, minus_label) { 
    list(
      annotate("text", x = x, y = max(d1$logFC) + 1, label = plus_label, hjust = 1, size = 3),
      annotate("text", x = x, y = min(d1$logFC) - 1, label = minus_label, hjust = 1, size = 3)
    ) 
  }, seq_along(mer1_filtered$Var1), paste ("(", rev(mer1_filtered$plus), ")"), paste ("(",rev(mer1_filtered$minus), ")"))

pdf ("Figure 6F sarmplot all A_v_B.pdf")
#pdf ("Figure 6H sarmplot all A_v_B _BR.pdf")
#pdf ("Figure 6J sarmplot all A_v_B _GR.pdf")
p
p1
p2  
dev.off()




###############################################################################################
############# FIG 6 GIK - SWARMPLOT - A vs C    SWARMPLOT - A vs C    ################################
###############################################################################################

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
d2 <- read.table("SCAFE/malignant_timepoint_analysis/Malignant_baseline_vs_4mnd__distal_cre__unibind.tsv", sep = '\t', header = TRUE, fill = TRUE)
d1 <- read.table("SCAFE/malignant_timepoint_analysis/Malignant_4mnd_vs_baseline__distal_cre__unibind.tsv", sep = '\t', header = TRUE, fill = TRUE)

#d2 <- read.table("SCAFE/malignant_timepoint_analysis/response_BR_A_vs_C__distal_cre__unibind.tsv", sep = '\t', header = TRUE, fill = TRUE)
#d1 <- read.table("SCAFE/malignant_timepoint_analysis/response_BR_C_vs_A__distal_cre__unibind.tsv", sep = '\t', header = TRUE, fill = TRUE)

#d2 <- read.table("SCAFE/malignant_timepoint_analysis/response_GR_A_vs_C__distal_cre__unibind.tsv", sep = '\t', header = TRUE, fill = TRUE)
#d1 <- read.table("SCAFE/malignant_timepoint_analysis/response_GR_C_vs_A__distal_cre__unibind.tsv", sep = '\t', header = TRUE, fill = TRUE)

d2$oddsRatio <- (-d2$oddsRatio)

  d <- rbind(d1[,colnames(d1) %in% c('collection', 'pValueLog', 'oddsRatio')], 
             d2[,colnames(d2) %in% c('collection', 'pValueLog', 'oddsRatio')])
  
  max <- max(d$oddsRatio[!d$oddsRatio %in% c('Inf')])
  min <- min(d$oddsRatio[!d$oddsRatio %in% c('-Inf')])
  
  d$pValue <- 10^(-d$pValueLog)

filtered_data <- subset(d, (pValue < 0.05 & oddsRatio > 1.5) | pValue < 0.05 & oddsRatio < -1.5 )
collection_count <- table(filtered_data$collection)

print(head (colnames(d)))
collection_count_plus <- as.data.frame(table (d$collection[d$oddsRatio > 1.5 & d$pValue < 0.05])) ; colnames (collection_count_plus)[2] <- 'plus'
collection_count_minus <- as.data.frame(table (d$collection[!d$oddsRatio > -1.5 & d$pValue < 0.05])) ;colnames (collection_count_minus)[2] <- 'minus'

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
  
  p$data$logFC_color <- ifelse(p$data$logFC_color >= 1.5, 2, p$data$logFC_color)
  p$data$logFC_color <- ifelse(p$data$logFC_color < -1.5, -2, p$data$logFC_color)
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

  p1$data$logFC_color <- ifelse(p1$data$logFC_color >= 1.5, 2, p1$data$logFC_color)
  p1$data$logFC_color <- ifelse(p1$data$logFC_color < -1.5, -2, p1$data$logFC_color)
  p1$data$logFC_color <- ifelse((p1$data$logFC_color == -2 | p1$data$logFC_color == 2), p1$data$logFC_color, NA)
  p1$data <- p1$data[order(!is.na(p1$data$logFC_color), p1$data$logFC_color), ]

p1 <- p1 + mapply(function(x, plus_label, minus_label) { 
    list(
      annotate("text", x = x, y = max(d1$logFC) + 1, label = plus_label, hjust = 1, size = 3),
      annotate("text", x = x, y = min(d1$logFC) - 1, label = minus_label, hjust = 1, size = 3)
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

  p2$data$logFC_color <- ifelse(p2$data$logFC_color >= 1.5, 2, p2$data$logFC_color)
  p2$data$logFC_color <- ifelse(p2$data$logFC_color < -1.5, -2, p2$data$logFC_color)
  p2$data$logFC_color <- ifelse((p2$data$logFC_color == -2 | p2$data$logFC_color == 2), p2$data$logFC_color, NA)
  p2$data <- p2$data[order(!is.na(p2$data$logFC_color), p2$data$logFC_color), ]

p2 <- p2 + mapply(function(x, plus_label, minus_label) { 
    list(
      annotate("text", x = x, y = max(d1$logFC) + 1, label = plus_label, hjust = 1, size = 3),
      annotate("text", x = x, y = min(d1$logFC) - 1, label = minus_label, hjust = 1, size = 3)
    ) 
  }, seq_along(mer1_filtered$Var1), paste ("(", rev(mer1_filtered$plus), ")"), paste ("(",rev(mer1_filtered$minus), ")"))

pdf ("Figure 6G sarmplot all A_v_C.pdf")
#pdf ("Figure 6I sarmplot all A_v_C_BR.pdf")
#pdf ("Figure 6K sarmplot all A_v_C_GR.pdf")
p
p1
p2  
dev.off()


###############################################################################################
############# FIG 6L - BOXPLOTS ESR1 score across time in GR     ################################
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

dat <- readRDS("SCAFE/TF_time_boxplots/data/Malignant_cells__ESR1_AR_FOXA1.rds")
### available from PP_snakemake_scafe

dat$ESR1_nCells <- dat$TF / dat$nCells
dat_GR <- dat[dat$response=="GR",]

pAB <- wilcox.test(dat_GR[dat_GR$time_point=="A", colnames (dat_GR)=="ESR1_nCells"], dat_GR[dat_GR$time_point=="B", colnames (dat_GR)=="ESR1_nCells"], paired = F)[3]
pBC <- wilcox.test(dat_GR[dat_GR$time_point=="B", colnames (dat_GR)=="ESR1_nCells"], dat_GR[dat_GR$time_point=="C", colnames (dat_GR)=="ESR1_nCells"], paired = F)[3]
pAC <- wilcox.test(dat_GR[dat_GR$time_point=="A", colnames (dat_GR)=="ESR1_nCells"], dat_GR[dat_GR$time_point=="C", colnames (dat_GR)=="ESR1_nCells"], paired = F)[3]
pvalues <- c(pAB, pBC, pAC)

dat_GR_complete <- dat_GR %>% complete(patient_number, time_point)


outline_colors <- c("#AA3830", '#9D2E99', '#350D7F') 
box_fill_colors <- c('#FFE2E0', '#FFD1FF', '#DAC3FF') 
point_fill_colors <- c('#FF8077', '#E66CE1', '#7E33FF') 


pdf('Figure 6L Boxplot ESR1_FOXA1_AR score GR .pdf')

a <- ggplot(dat_GR, aes(time_point, ESR1_nCells, label = patient_number)) + 
 geom_boxplot(color = outline_colors, 
               fill = box_fill_colors, 
               alpha = 0.5, size = 1.5) + 
  geom_point(aes(fill = factor(time_point), colour = factor(time_point)), 
             shape = 21, size = 5, stroke = 1.2) +  
  scale_fill_manual(values = point_fill_colors) +
  scale_colour_manual(values = outline_colors) + 
  ylab('Percentage of malignant cells')+ 
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
              y_position = 35
			 ,  
              map_signif_level = FALSE,
              test = "wilcox.test", 
              step_increase = 0.07)
print(a)
dev.off()

###############################################################################################
############# FIG 6M - BOXPLOTS ESR1 score across time in GR     ################################
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

dat <- readRDS("SCAFE/TF_time_boxplots/data/Malignant_cells__ESR1_AR_FOXA1.rds")
### available from PP_snakemake_scafe

dat$ESR1_nCells <- dat$TF / dat$nCells
dat_BR <- dat[dat$response=="BR",]

pAB <- wilcox.test(dat_BR[dat_BR$time_point=="A", colnames (dat_BR)=="ESR1_nCells"], dat_BR[dat_BR$time_point=="B", colnames (dat_BR)=="ESR1_nCells"], paired = F)[3]
pBC <- wilcox.test(dat_BR[dat_BR$time_point=="B", colnames (dat_BR)=="ESR1_nCells"], dat_BR[dat_BR$time_point=="C", colnames (dat_BR)=="ESR1_nCells"], paired = F)[3]
pAC <- wilcox.test(dat_BR[dat_BR$time_point=="A", colnames (dat_BR)=="ESR1_nCells"], dat_BR[dat_BR$time_point=="C", colnames (dat_BR)=="ESR1_nCells"], paired = F)[3]
pvalues <- c(pAB, pBC, pAC)

dat_BR_complete <- dat_BR %>% complete(patient_number, time_point)


outline_colors <- c("#AA3830", '#9D2E99', '#350D7F') 
box_fill_colors <- c('#FFE2E0', '#FFD1FF', '#DAC3FF') 
point_fill_colors <- c('#FF8077', '#E66CE1', '#7E33FF')



pdf('Figure 6M Boxplot ESR1_FOXA1_AR score BR .pdf')

a <- ggplot(dat_BR, aes(time_point, ESR1_nCells, label = patient_number)) + 
 geom_boxplot(color = outline_colors, 
               fill = box_fill_colors, 
               alpha = 0.5, size = 1.5) + 
  geom_point(aes(fill = factor(time_point), colour = factor(time_point)), 
             shape = 21, size = 5, stroke = 1.2) +  
  scale_fill_manual(values = point_fill_colors) +
  scale_colour_manual(values = outline_colors) + 
  ylab('Percentage of malignant cells')+ 
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
              y_position = 25,  
              map_signif_level = FALSE,
              test = "wilcox.test", 
              step_increase = 0.07)
print(a)
 
dev.off()
