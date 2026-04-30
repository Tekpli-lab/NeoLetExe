########################################
###### Supplementary Figure 13 #########
########################################

################################################################################################################
############# SUPPLEMENTARY FIGURE 13A- DINUCLEOTIDES    #######################################################
################################################################################################################
################################################################################################################
############# SUPPLEMENTARY FIGURE 13B- DINUCLEOTIDES     ######################################################
################################################################################################################

################################################################################################################
############# SUPPLEMENTARY FIGURE 13C- BARPLOT INTERSECT ATAC chromHMM     ####################################
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

pdf("SuppFigure13C Intersection with scATAC-seq.pdf", width=14, height=8, onefile=FALSE, useDingbats=FALSE)

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
############# SUPPLEMENTARY FIGURE 13D-E - UPSET PLOTS - proximal distal     ################################
################################################################################################################
library (UpSetR)
pdf("SuppFigure11D Intervene_upset_distal.pdf", width=14, height=8, onefile=FALSE, useDingbats=FALSE)
expressionInput <- c('T_cells_0'=5952,
'Normal_epithelial_0'=169,
'Normal_epithelial_0&T_cells_0'=4,
'Myeloid_0'=13465,
'Myeloid_0&T_cells_0'=927,'Myeloid_0&Normal_epithelial_0'=119,'Myeloid_0&Normal_epithelial_0&T_cells_0'=0,'Malignant_cells_0'=10953,'Malignant_cells_0&T_cells_0'=164,'Malignant_cells_0&Normal_epithelial_0'=1582,'Malignant_cells_0&Normal_epithelial_0&T_cells_0'=7,'Malignant_cells_0&Myeloid_0'=382,'Malignant_cells_0&Myeloid_0&T_cells_0'=10,'Malignant_cells_0&Myeloid_0&Normal_epithelial_0'=23,'Malignant_cells_0&Myeloid_0&Normal_epithelial_0&T_cells_0'=1,'Fibroblasts_0'=8407,'Fibroblasts_0&T_cells_0'=121,'Fibroblasts_0&Normal_epithelial_0'=37,'Fibroblasts_0&Normal_epithelial_0&T_cells_0'=0,'Fibroblasts_0&Myeloid_0'=578,'Fibroblasts_0&Myeloid_0&T_cells_0'=11,'Fibroblasts_0&Myeloid_0&Normal_epithelial_0'=1,'Fibroblasts_0&Myeloid_0&Normal_epithelial_0&T_cells_0'=0,'Fibroblasts_0&Malignant_cells_0'=487,'Fibroblasts_0&Malignant_cells_0&T_cells_0'=2,'Fibroblasts_0&Malignant_cells_0&Normal_epithelial_0'=31,'Fibroblasts_0&Malignant_cells_0&Normal_epithelial_0&T_cells_0'=0,'Fibroblasts_0&Malignant_cells_0&Myeloid_0'=22,'Fibroblasts_0&Malignant_cells_0&Myeloid_0&T_cells_0'=0,'Fibroblasts_0&Malignant_cells_0&Myeloid_0&Normal_epithelial_0'=0,'Fibroblasts_0&Malignant_cells_0&Myeloid_0&Normal_epithelial_0&T_cells_0'=0,'Endothelial_cells_0'=3811,'Endothelial_cells_0&T_cells_0'=185,'Endothelial_cells_0&Normal_epithelial_0'=14,'Endothelial_cells_0&Normal_epithelial_0&T_cells_0'=0,'Endothelial_cells_0&Myeloid_0'=564,'Endothelial_cells_0&Myeloid_0&T_cells_0'=30,'Endothelial_cells_0&Myeloid_0&Normal_epithelial_0'=3,'Endothelial_cells_0&Myeloid_0&Normal_epithelial_0&T_cells_0'=0,'Endothelial_cells_0&Malignant_cells_0'=289,'Endothelial_cells_0&Malignant_cells_0&T_cells_0'=4,'Endothelial_cells_0&Malignant_cells_0&Normal_epithelial_0'=16,'Endothelial_cells_0&Malignant_cells_0&Normal_epithelial_0&T_cells_0'=0,'Endothelial_cells_0&Malignant_cells_0&Myeloid_0'=14,'Endothelial_cells_0&Malignant_cells_0&Myeloid_0&T_cells_0'=0,'Endothelial_cells_0&Malignant_cells_0&Myeloid_0&Normal_epithelial_0'=1,'Endothelial_cells_0&Malignant_cells_0&Myeloid_0&Normal_epithelial_0&T_cells_0'=0,'Endothelial_cells_0&Fibroblasts_0'=786,'Endothelial_cells_0&Fibroblasts_0&T_cells_0'=11,'Endothelial_cells_0&Fibroblasts_0&Normal_epithelial_0'=2,'Endothelial_cells_0&Fibroblasts_0&Normal_epithelial_0&T_cells_0'=0,'Endothelial_cells_0&Fibroblasts_0&Myeloid_0'=58,'Endothelial_cells_0&Fibroblasts_0&Myeloid_0&T_cells_0'=3,'Endothelial_cells_0&Fibroblasts_0&Myeloid_0&Normal_epithelial_0'=0,'Endothelial_cells_0&Fibroblasts_0&Myeloid_0&Normal_epithelial_0&T_cells_0'=0,'Endothelial_cells_0&Fibroblasts_0&Malignant_cells_0'=20,'Endothelial_cells_0&Fibroblasts_0&Malignant_cells_0&T_cells_0'=0,'Endothelial_cells_0&Fibroblasts_0&Malignant_cells_0&Normal_epithelial_0'=3,'Endothelial_cells_0&Fibroblasts_0&Malignant_cells_0&Normal_epithelial_0&T_cells_0'=0,'Endothelial_cells_0&Fibroblasts_0&Malignant_cells_0&Myeloid_0'=0,'Endothelial_cells_0&Fibroblasts_0&Malignant_cells_0&Myeloid_0&T_cells_0'=0,'Endothelial_cells_0&Fibroblasts_0&Malignant_cells_0&Myeloid_0&Normal_epithelial_0'=0,'Endothelial_cells_0&Fibroblasts_0&Malignant_cells_0&Myeloid_0&Normal_epithelial_0&T_cells_0'=0,'BP_cells_0'=2238,'BP_cells_0&T_cells_0'=809,'BP_cells_0&Normal_epithelial_0'=3,'BP_cells_0&Normal_epithelial_0&T_cells_0'=0,'BP_cells_0&Myeloid_0'=932,'BP_cells_0&Myeloid_0&T_cells_0'=274,'BP_cells_0&Myeloid_0&Normal_epithelial_0'=1,'BP_cells_0&Myeloid_0&Normal_epithelial_0&T_cells_0'=1,'BP_cells_0&Malignant_cells_0'=770,'BP_cells_0&Malignant_cells_0&T_cells_0'=38,'BP_cells_0&Malignant_cells_0&Normal_epithelial_0'=15,'BP_cells_0&Malignant_cells_0&Normal_epithelial_0&T_cells_0'=0,'BP_cells_0&Malignant_cells_0&Myeloid_0'=47,'BP_cells_0&Malignant_cells_0&Myeloid_0&T_cells_0'=4,'BP_cells_0&Malignant_cells_0&Myeloid_0&Normal_epithelial_0'=0,'BP_cells_0&Malignant_cells_0&Myeloid_0&Normal_epithelial_0&T_cells_0'=0,'BP_cells_0&Fibroblasts_0'=79,'BP_cells_0&Fibroblasts_0&T_cells_0'=27,'BP_cells_0&Fibroblasts_0&Normal_epithelial_0'=0,'BP_cells_0&Fibroblasts_0&Normal_epithelial_0&T_cells_0'=0,'BP_cells_0&Fibroblasts_0&Myeloid_0'=31,'BP_cells_0&Fibroblasts_0&Myeloid_0&T_cells_0'=9,'BP_cells_0&Fibroblasts_0&Myeloid_0&Normal_epithelial_0'=0,'BP_cells_0&Fibroblasts_0&Myeloid_0&Normal_epithelial_0&T_cells_0'=0,'BP_cells_0&Fibroblasts_0&Malignant_cells_0'=16,'BP_cells_0&Fibroblasts_0&Malignant_cells_0&T_cells_0'=3,'BP_cells_0&Fibroblasts_0&Malignant_cells_0&Normal_epithelial_0'=0,'BP_cells_0&Fibroblasts_0&Malignant_cells_0&Normal_epithelial_0&T_cells_0'=0,'BP_cells_0&Fibroblasts_0&Malignant_cells_0&Myeloid_0'=1,'BP_cells_0&Fibroblasts_0&Malignant_cells_0&Myeloid_0&T_cells_0'=0,'BP_cells_0&Fibroblasts_0&Malignant_cells_0&Myeloid_0&Normal_epithelial_0'=0,'BP_cells_0&Fibroblasts_0&Malignant_cells_0&Myeloid_0&Normal_epithelial_0&T_cells_0'=0,'BP_cells_0&Endothelial_cells_0'=61,'BP_cells_0&Endothelial_cells_0&T_cells_0'=15,'BP_cells_0&Endothelial_cells_0&Normal_epithelial_0'=1,'BP_cells_0&Endothelial_cells_0&Normal_epithelial_0&T_cells_0'=0,'BP_cells_0&Endothelial_cells_0&Myeloid_0'=27,'BP_cells_0&Endothelial_cells_0&Myeloid_0&T_cells_0'=5,'BP_cells_0&Endothelial_cells_0&Myeloid_0&Normal_epithelial_0'=0,'BP_cells_0&Endothelial_cells_0&Myeloid_0&Normal_epithelial_0&T_cells_0'=0,'BP_cells_0&Endothelial_cells_0&Malignant_cells_0'=6,'BP_cells_0&Endothelial_cells_0&Malignant_cells_0&T_cells_0'=0,'BP_cells_0&Endothelial_cells_0&Malignant_cells_0&Normal_epithelial_0'=0,'BP_cells_0&Endothelial_cells_0&Malignant_cells_0&Normal_epithelial_0&T_cells_0'=0,'BP_cells_0&Endothelial_cells_0&Malignant_cells_0&Myeloid_0'=1,'BP_cells_0&Endothelial_cells_0&Malignant_cells_0&Myeloid_0&T_cells_0'=0,'BP_cells_0&Endothelial_cells_0&Malignant_cells_0&Myeloid_0&Normal_epithelial_0'=0,'BP_cells_0&Endothelial_cells_0&Malignant_cells_0&Myeloid_0&Normal_epithelial_0&T_cells_0'=0,'BP_cells_0&Endothelial_cells_0&Fibroblasts_0'=16,'BP_cells_0&Endothelial_cells_0&Fibroblasts_0&T_cells_0'=3,'BP_cells_0&Endothelial_cells_0&Fibroblasts_0&Normal_epithelial_0'=0,'BP_cells_0&Endothelial_cells_0&Fibroblasts_0&Normal_epithelial_0&T_cells_0'=0,'BP_cells_0&Endothelial_cells_0&Fibroblasts_0&Myeloid_0'=2,'BP_cells_0&Endothelial_cells_0&Fibroblasts_0&Myeloid_0&T_cells_0'=1,'BP_cells_0&Endothelial_cells_0&Fibroblasts_0&Myeloid_0&Normal_epithelial_0'=0,'BP_cells_0&Endothelial_cells_0&Fibroblasts_0&Myeloid_0&Normal_epithelial_0&T_cells_0'=0,'BP_cells_0&Endothelial_cells_0&Fibroblasts_0&Malignant_cells_0'=0,'BP_cells_0&Endothelial_cells_0&Fibroblasts_0&Malignant_cells_0&T_cells_0'=1,'BP_cells_0&Endothelial_cells_0&Fibroblasts_0&Malignant_cells_0&Normal_epithelial_0'=0,'BP_cells_0&Endothelial_cells_0&Fibroblasts_0&Malignant_cells_0&Normal_epithelial_0&T_cells_0'=0,'BP_cells_0&Endothelial_cells_0&Fibroblasts_0&Malignant_cells_0&Myeloid_0'=0,'BP_cells_0&Endothelial_cells_0&Fibroblasts_0&Malignant_cells_0&Myeloid_0&T_cells_0'=0,'BP_cells_0&Endothelial_cells_0&Fibroblasts_0&Malignant_cells_0&Myeloid_0&Normal_epithelial_0'=0,'BP_cells_0&Endothelial_cells_0&Fibroblasts_0&Malignant_cells_0&Myeloid_0&Normal_epithelial_0&T_cells_0'=0,'Adipocytes_0'=165,'Adipocytes_0&T_cells_0'=9,'Adipocytes_0&Normal_epithelial_0'=5,'Adipocytes_0&Normal_epithelial_0&T_cells_0'=0,'Adipocytes_0&Myeloid_0'=29,'Adipocytes_0&Myeloid_0&T_cells_0'=0,'Adipocytes_0&Myeloid_0&Normal_epithelial_0'=0,'Adipocytes_0&Myeloid_0&Normal_epithelial_0&T_cells_0'=0,'Adipocytes_0&Malignant_cells_0'=31,'Adipocytes_0&Malignant_cells_0&T_cells_0'=0,'Adipocytes_0&Malignant_cells_0&Normal_epithelial_0'=5,'Adipocytes_0&Malignant_cells_0&Normal_epithelial_0&T_cells_0'=0,'Adipocytes_0&Malignant_cells_0&Myeloid_0'=2,'Adipocytes_0&Malignant_cells_0&Myeloid_0&T_cells_0'=0,'Adipocytes_0&Malignant_cells_0&Myeloid_0&Normal_epithelial_0'=0,'Adipocytes_0&Malignant_cells_0&Myeloid_0&Normal_epithelial_0&T_cells_0'=0,'Adipocytes_0&Fibroblasts_0'=129,'Adipocytes_0&Fibroblasts_0&T_cells_0'=2,'Adipocytes_0&Fibroblasts_0&Normal_epithelial_0'=1,'Adipocytes_0&Fibroblasts_0&Normal_epithelial_0&T_cells_0'=0,'Adipocytes_0&Fibroblasts_0&Myeloid_0'=8,'Adipocytes_0&Fibroblasts_0&Myeloid_0&T_cells_0'=0,'Adipocytes_0&Fibroblasts_0&Myeloid_0&Normal_epithelial_0'=0,'Adipocytes_0&Fibroblasts_0&Myeloid_0&Normal_epithelial_0&T_cells_0'=0,'Adipocytes_0&Fibroblasts_0&Malignant_cells_0'=8,'Adipocytes_0&Fibroblasts_0&Malignant_cells_0&T_cells_0'=0,'Adipocytes_0&Fibroblasts_0&Malignant_cells_0&Normal_epithelial_0'=2,'Adipocytes_0&Fibroblasts_0&Malignant_cells_0&Normal_epithelial_0&T_cells_0'=0,'Adipocytes_0&Fibroblasts_0&Malignant_cells_0&Myeloid_0'=1,'Adipocytes_0&Fibroblasts_0&Malignant_cells_0&Myeloid_0&T_cells_0'=0,'Adipocytes_0&Fibroblasts_0&Malignant_cells_0&Myeloid_0&Normal_epithelial_0'=0,'Adipocytes_0&Fibroblasts_0&Malignant_cells_0&Myeloid_0&Normal_epithelial_0&T_cells_0'=0,'Adipocytes_0&Endothelial_cells_0'=183,'Adipocytes_0&Endothelial_cells_0&T_cells_0'=7,'Adipocytes_0&Endothelial_cells_0&Normal_epithelial_0'=0,'Adipocytes_0&Endothelial_cells_0&Normal_epithelial_0&T_cells_0'=0,'Adipocytes_0&Endothelial_cells_0&Myeloid_0'=18,'Adipocytes_0&Endothelial_cells_0&Myeloid_0&T_cells_0'=0,'Adipocytes_0&Endothelial_cells_0&Myeloid_0&Normal_epithelial_0'=0,'Adipocytes_0&Endothelial_cells_0&Myeloid_0&Normal_epithelial_0&T_cells_0'=0,'Adipocytes_0&Endothelial_cells_0&Malignant_cells_0'=12,'Adipocytes_0&Endothelial_cells_0&Malignant_cells_0&T_cells_0'=0,'Adipocytes_0&Endothelial_cells_0&Malignant_cells_0&Normal_epithelial_0'=2,'Adipocytes_0&Endothelial_cells_0&Malignant_cells_0&Normal_epithelial_0&T_cells_0'=0,'Adipocytes_0&Endothelial_cells_0&Malignant_cells_0&Myeloid_0'=1,'Adipocytes_0&Endothelial_cells_0&Malignant_cells_0&Myeloid_0&T_cells_0'=0,'Adipocytes_0&Endothelial_cells_0&Malignant_cells_0&Myeloid_0&Normal_epithelial_0'=0,'Adipocytes_0&Endothelial_cells_0&Malignant_cells_0&Myeloid_0&Normal_epithelial_0&T_cells_0'=0,'Adipocytes_0&Endothelial_cells_0&Fibroblasts_0'=53,'Adipocytes_0&Endothelial_cells_0&Fibroblasts_0&T_cells_0'=0,'Adipocytes_0&Endothelial_cells_0&Fibroblasts_0&Normal_epithelial_0'=1,'Adipocytes_0&Endothelial_cells_0&Fibroblasts_0&Normal_epithelial_0&T_cells_0'=0,'Adipocytes_0&Endothelial_cells_0&Fibroblasts_0&Myeloid_0'=6,'Adipocytes_0&Endothelial_cells_0&Fibroblasts_0&Myeloid_0&T_cells_0'=0,'Adipocytes_0&Endothelial_cells_0&Fibroblasts_0&Myeloid_0&Normal_epithelial_0'=0,'Adipocytes_0&Endothelial_cells_0&Fibroblasts_0&Myeloid_0&Normal_epithelial_0&T_cells_0'=0,'Adipocytes_0&Endothelial_cells_0&Fibroblasts_0&Malignant_cells_0'=5,'Adipocytes_0&Endothelial_cells_0&Fibroblasts_0&Malignant_cells_0&T_cells_0'=0,'Adipocytes_0&Endothelial_cells_0&Fibroblasts_0&Malignant_cells_0&Normal_epithelial_0'=0,'Adipocytes_0&Endothelial_cells_0&Fibroblasts_0&Malignant_cells_0&Normal_epithelial_0&T_cells_0'=0,'Adipocytes_0&Endothelial_cells_0&Fibroblasts_0&Malignant_cells_0&Myeloid_0'=1,'Adipocytes_0&Endothelial_cells_0&Fibroblasts_0&Malignant_cells_0&Myeloid_0&T_cells_0'=0,'Adipocytes_0&Endothelial_cells_0&Fibroblasts_0&Malignant_cells_0&Myeloid_0&Normal_epithelial_0'=0,'Adipocytes_0&Endothelial_cells_0&Fibroblasts_0&Malignant_cells_0&Myeloid_0&Normal_epithelial_0&T_cells_0'=0,'Adipocytes_0&BP_cells_0'=4,'Adipocytes_0&BP_cells_0&T_cells_0'=4,'Adipocytes_0&BP_cells_0&Normal_epithelial_0'=0,'Adipocytes_0&BP_cells_0&Normal_epithelial_0&T_cells_0'=0,'Adipocytes_0&BP_cells_0&Myeloid_0'=3,'Adipocytes_0&BP_cells_0&Myeloid_0&T_cells_0'=1,'Adipocytes_0&BP_cells_0&Myeloid_0&Normal_epithelial_0'=0,'Adipocytes_0&BP_cells_0&Myeloid_0&Normal_epithelial_0&T_cells_0'=0,'Adipocytes_0&BP_cells_0&Malignant_cells_0'=1,'Adipocytes_0&BP_cells_0&Malignant_cells_0&T_cells_0'=1,'Adipocytes_0&BP_cells_0&Malignant_cells_0&Normal_epithelial_0'=0,'Adipocytes_0&BP_cells_0&Malignant_cells_0&Normal_epithelial_0&T_cells_0'=0,'Adipocytes_0&BP_cells_0&Malignant_cells_0&Myeloid_0'=1,'Adipocytes_0&BP_cells_0&Malignant_cells_0&Myeloid_0&T_cells_0'=0,'Adipocytes_0&BP_cells_0&Malignant_cells_0&Myeloid_0&Normal_epithelial_0'=0,'Adipocytes_0&BP_cells_0&Malignant_cells_0&Myeloid_0&Normal_epithelial_0&T_cells_0'=0,'Adipocytes_0&BP_cells_0&Fibroblasts_0'=0,'Adipocytes_0&BP_cells_0&Fibroblasts_0&T_cells_0'=1,'Adipocytes_0&BP_cells_0&Fibroblasts_0&Normal_epithelial_0'=0,'Adipocytes_0&BP_cells_0&Fibroblasts_0&Normal_epithelial_0&T_cells_0'=0,'Adipocytes_0&BP_cells_0&Fibroblasts_0&Myeloid_0'=3,'Adipocytes_0&BP_cells_0&Fibroblasts_0&Myeloid_0&T_cells_0'=0,'Adipocytes_0&BP_cells_0&Fibroblasts_0&Myeloid_0&Normal_epithelial_0'=0,'Adipocytes_0&BP_cells_0&Fibroblasts_0&Myeloid_0&Normal_epithelial_0&T_cells_0'=0,'Adipocytes_0&BP_cells_0&Fibroblasts_0&Malignant_cells_0'=0,'Adipocytes_0&BP_cells_0&Fibroblasts_0&Malignant_cells_0&T_cells_0'=0,'Adipocytes_0&BP_cells_0&Fibroblasts_0&Malignant_cells_0&Normal_epithelial_0'=0,'Adipocytes_0&BP_cells_0&Fibroblasts_0&Malignant_cells_0&Normal_epithelial_0&T_cells_0'=0,'Adipocytes_0&BP_cells_0&Fibroblasts_0&Malignant_cells_0&Myeloid_0'=0,'Adipocytes_0&BP_cells_0&Fibroblasts_0&Malignant_cells_0&Myeloid_0&T_cells_0'=0,'Adipocytes_0&BP_cells_0&Fibroblasts_0&Malignant_cells_0&Myeloid_0&Normal_epithelial_0'=0,'Adipocytes_0&BP_cells_0&Fibroblasts_0&Malignant_cells_0&Myeloid_0&Normal_epithelial_0&T_cells_0'=0,'Adipocytes_0&BP_cells_0&Endothelial_cells_0'=4,'Adipocytes_0&BP_cells_0&Endothelial_cells_0&T_cells_0'=1,'Adipocytes_0&BP_cells_0&Endothelial_cells_0&Normal_epithelial_0'=0,'Adipocytes_0&BP_cells_0&Endothelial_cells_0&Normal_epithelial_0&T_cells_0'=0,'Adipocytes_0&BP_cells_0&Endothelial_cells_0&Myeloid_0'=0,'Adipocytes_0&BP_cells_0&Endothelial_cells_0&Myeloid_0&T_cells_0'=2,'Adipocytes_0&BP_cells_0&Endothelial_cells_0&Myeloid_0&Normal_epithelial_0'=0,'Adipocytes_0&BP_cells_0&Endothelial_cells_0&Myeloid_0&Normal_epithelial_0&T_cells_0'=0,'Adipocytes_0&BP_cells_0&Endothelial_cells_0&Malignant_cells_0'=0,'Adipocytes_0&BP_cells_0&Endothelial_cells_0&Malignant_cells_0&T_cells_0'=0,'Adipocytes_0&BP_cells_0&Endothelial_cells_0&Malignant_cells_0&Normal_epithelial_0'=0,'Adipocytes_0&BP_cells_0&Endothelial_cells_0&Malignant_cells_0&Normal_epithelial_0&T_cells_0'=0,'Adipocytes_0&BP_cells_0&Endothelial_cells_0&Malignant_cells_0&Myeloid_0'=0,'Adipocytes_0&BP_cells_0&Endothelial_cells_0&Malignant_cells_0&Myeloid_0&T_cells_0'=0,'Adipocytes_0&BP_cells_0&Endothelial_cells_0&Malignant_cells_0&Myeloid_0&Normal_epithelial_0'=0,'Adipocytes_0&BP_cells_0&Endothelial_cells_0&Malignant_cells_0&Myeloid_0&Normal_epithelial_0&T_cells_0'=0,'Adipocytes_0&BP_cells_0&Endothelial_cells_0&Fibroblasts_0'=6,'Adipocytes_0&BP_cells_0&Endothelial_cells_0&Fibroblasts_0&T_cells_0'=0,'Adipocytes_0&BP_cells_0&Endothelial_cells_0&Fibroblasts_0&Normal_epithelial_0'=0,'Adipocytes_0&BP_cells_0&Endothelial_cells_0&Fibroblasts_0&Normal_epithelial_0&T_cells_0'=0,'Adipocytes_0&BP_cells_0&Endothelial_cells_0&Fibroblasts_0&Myeloid_0'=1,'Adipocytes_0&BP_cells_0&Endothelial_cells_0&Fibroblasts_0&Myeloid_0&T_cells_0'=0,'Adipocytes_0&BP_cells_0&Endothelial_cells_0&Fibroblasts_0&Myeloid_0&Normal_epithelial_0'=0,'Adipocytes_0&BP_cells_0&Endothelial_cells_0&Fibroblasts_0&Myeloid_0&Normal_epithelial_0&T_cells_0'=0,'Adipocytes_0&BP_cells_0&Endothelial_cells_0&Fibroblasts_0&Malignant_cells_0'=0,'Adipocytes_0&BP_cells_0&Endothelial_cells_0&Fibroblasts_0&Malignant_cells_0&T_cells_0'=0,'Adipocytes_0&BP_cells_0&Endothelial_cells_0&Fibroblasts_0&Malignant_cells_0&Normal_epithelial_0'=0,'Adipocytes_0&BP_cells_0&Endothelial_cells_0&Fibroblasts_0&Malignant_cells_0&Normal_epithelial_0&T_cells_0'=0,'Adipocytes_0&BP_cells_0&Endothelial_cells_0&Fibroblasts_0&Malignant_cells_0&Myeloid_0'=0,'Adipocytes_0&BP_cells_0&Endothelial_cells_0&Fibroblasts_0&Malignant_cells_0&Myeloid_0&T_cells_0'=0,'Adipocytes_0&BP_cells_0&Endothelial_cells_0&Fibroblasts_0&Malignant_cells_0&Myeloid_0&Normal_epithelial_0'=0,'Adipocytes_0&BP_cells_0&Endothelial_cells_0&Fibroblasts_0&Malignant_cells_0&Myeloid_0&Normal_epithelial_0&T_cells_0'=0)
upset(fromExpression(expressionInput), nsets=8, nintersects=30, show.numbers="yes", main.bar.color="#ea5d4e", sets.bar.color="#317eab", empty.intersections=NULL, order.by = "freq", number.angles = 0, mainbar.y.label ="No. of Intersections", sets.x.label ="Set size")
invisible(dev.off())

library (UpSetR)
pdf("SuppFigure11E Intervene_upset_proximal.pdf", width=14, height=8, onefile=FALSE, useDingbats=FALSE)
expressionInput <- c('T_cells_0'=1348,'Normal_epithelial_0'=135,'Normal_epithelial_0&T_cells_0'=18,'Myeloid_0'=1567,'Myeloid_0&T_cells_0'=337,'Myeloid_0&Normal_epithelial_0'=77,'Myeloid_0&Normal_epithelial_0&T_cells_0'=4,'Malignant_cells_0'=4338,'Malignant_cells_0&T_cells_0'=347,'Malignant_cells_0&Normal_epithelial_0'=776,'Malignant_cells_0&Normal_epithelial_0&T_cells_0'=27,'Malignant_cells_0&Myeloid_0'=410,'Malignant_cells_0&Myeloid_0&T_cells_0'=116,'Malignant_cells_0&Myeloid_0&Normal_epithelial_0'=48,'Malignant_cells_0&Myeloid_0&Normal_epithelial_0&T_cells_0'=3,'Fibroblasts_0'=3745,'Fibroblasts_0&T_cells_0'=423,'Fibroblasts_0&Normal_epithelial_0'=124,'Fibroblasts_0&Normal_epithelial_0&T_cells_0'=5,'Fibroblasts_0&Myeloid_0'=581,'Fibroblasts_0&Myeloid_0&T_cells_0'=175,'Fibroblasts_0&Myeloid_0&Normal_epithelial_0'=24,'Fibroblasts_0&Myeloid_0&Normal_epithelial_0&T_cells_0'=1,'Fibroblasts_0&Malignant_cells_0'=2476,'Fibroblasts_0&Malignant_cells_0&T_cells_0'=445,'Fibroblasts_0&Malignant_cells_0&Normal_epithelial_0'=120,'Fibroblasts_0&Malignant_cells_0&Normal_epithelial_0&T_cells_0'=14,'Fibroblasts_0&Malignant_cells_0&Myeloid_0'=487,'Fibroblasts_0&Malignant_cells_0&Myeloid_0&T_cells_0'=136,'Fibroblasts_0&Malignant_cells_0&Myeloid_0&Normal_epithelial_0'=11,'Fibroblasts_0&Malignant_cells_0&Myeloid_0&Normal_epithelial_0&T_cells_0'=1,'Endothelial_cells_0'=695,'Endothelial_cells_0&T_cells_0'=80,'Endothelial_cells_0&Normal_epithelial_0'=18,'Endothelial_cells_0&Normal_epithelial_0&T_cells_0'=1,'Endothelial_cells_0&Myeloid_0'=147,'Endothelial_cells_0&Myeloid_0&T_cells_0'=28,'Endothelial_cells_0&Myeloid_0&Normal_epithelial_0'=9,'Endothelial_cells_0&Myeloid_0&Normal_epithelial_0&T_cells_0'=0,'Endothelial_cells_0&Malignant_cells_0'=229,'Endothelial_cells_0&Malignant_cells_0&T_cells_0'=16,'Endothelial_cells_0&Malignant_cells_0&Normal_epithelial_0'=29,'Endothelial_cells_0&Malignant_cells_0&Normal_epithelial_0&T_cells_0'=0,'Endothelial_cells_0&Malignant_cells_0&Myeloid_0'=30,'Endothelial_cells_0&Malignant_cells_0&Myeloid_0&T_cells_0'=5,'Endothelial_cells_0&Malignant_cells_0&Myeloid_0&Normal_epithelial_0'=5,'Endothelial_cells_0&Malignant_cells_0&Myeloid_0&Normal_epithelial_0&T_cells_0'=0,'Endothelial_cells_0&Fibroblasts_0'=625,'Endothelial_cells_0&Fibroblasts_0&T_cells_0'=91,'Endothelial_cells_0&Fibroblasts_0&Normal_epithelial_0'=16,'Endothelial_cells_0&Fibroblasts_0&Normal_epithelial_0&T_cells_0'=1,'Endothelial_cells_0&Fibroblasts_0&Myeloid_0'=165,'Endothelial_cells_0&Fibroblasts_0&Myeloid_0&T_cells_0'=36,'Endothelial_cells_0&Fibroblasts_0&Myeloid_0&Normal_epithelial_0'=2,'Endothelial_cells_0&Fibroblasts_0&Myeloid_0&Normal_epithelial_0&T_cells_0'=0,'Endothelial_cells_0&Fibroblasts_0&Malignant_cells_0'=266,'Endothelial_cells_0&Fibroblasts_0&Malignant_cells_0&T_cells_0'=26,'Endothelial_cells_0&Fibroblasts_0&Malignant_cells_0&Normal_epithelial_0'=12,'Endothelial_cells_0&Fibroblasts_0&Malignant_cells_0&Normal_epithelial_0&T_cells_0'=0,'Endothelial_cells_0&Fibroblasts_0&Malignant_cells_0&Myeloid_0'=39,'Endothelial_cells_0&Fibroblasts_0&Malignant_cells_0&Myeloid_0&T_cells_0'=7,'Endothelial_cells_0&Fibroblasts_0&Malignant_cells_0&Myeloid_0&Normal_epithelial_0'=1,'Endothelial_cells_0&Fibroblasts_0&Malignant_cells_0&Myeloid_0&Normal_epithelial_0&T_cells_0'=0,'BP_cells_0'=507,'BP_cells_0&T_cells_0'=627,'BP_cells_0&Normal_epithelial_0'=6,'BP_cells_0&Normal_epithelial_0&T_cells_0'=4,'BP_cells_0&Myeloid_0'=368,'BP_cells_0&Myeloid_0&T_cells_0'=628,'BP_cells_0&Myeloid_0&Normal_epithelial_0'=5,'BP_cells_0&Myeloid_0&Normal_epithelial_0&T_cells_0'=11,'BP_cells_0&Malignant_cells_0'=475,'BP_cells_0&Malignant_cells_0&T_cells_0'=466,'BP_cells_0&Malignant_cells_0&Normal_epithelial_0'=36,'BP_cells_0&Malignant_cells_0&Normal_epithelial_0&T_cells_0'=10,'BP_cells_0&Malignant_cells_0&Myeloid_0'=132,'BP_cells_0&Malignant_cells_0&Myeloid_0&T_cells_0'=224,'BP_cells_0&Malignant_cells_0&Myeloid_0&Normal_epithelial_0'=14,'BP_cells_0&Malignant_cells_0&Myeloid_0&Normal_epithelial_0&T_cells_0'=17,'BP_cells_0&Fibroblasts_0'=207,'BP_cells_0&Fibroblasts_0&T_cells_0'=616,'BP_cells_0&Fibroblasts_0&Normal_epithelial_0'=3,'BP_cells_0&Fibroblasts_0&Normal_epithelial_0&T_cells_0'=1,'BP_cells_0&Fibroblasts_0&Myeloid_0'=140,'BP_cells_0&Fibroblasts_0&Myeloid_0&T_cells_0'=347,'BP_cells_0&Fibroblasts_0&Myeloid_0&Normal_epithelial_0'=0,'BP_cells_0&Fibroblasts_0&Myeloid_0&Normal_epithelial_0&T_cells_0'=2,'BP_cells_0&Fibroblasts_0&Malignant_cells_0'=450,'BP_cells_0&Fibroblasts_0&Malignant_cells_0&T_cells_0'=669,'BP_cells_0&Fibroblasts_0&Malignant_cells_0&Normal_epithelial_0'=11,'BP_cells_0&Fibroblasts_0&Malignant_cells_0&Normal_epithelial_0&T_cells_0'=4,'BP_cells_0&Fibroblasts_0&Malignant_cells_0&Myeloid_0'=164,'BP_cells_0&Fibroblasts_0&Malignant_cells_0&Myeloid_0&T_cells_0'=244,'BP_cells_0&Fibroblasts_0&Malignant_cells_0&Myeloid_0&Normal_epithelial_0'=7,'BP_cells_0&Fibroblasts_0&Malignant_cells_0&Myeloid_0&Normal_epithelial_0&T_cells_0'=10,'BP_cells_0&Endothelial_cells_0'=29,'BP_cells_0&Endothelial_cells_0&T_cells_0'=49,'BP_cells_0&Endothelial_cells_0&Normal_epithelial_0'=1,'BP_cells_0&Endothelial_cells_0&Normal_epithelial_0&T_cells_0'=0,'BP_cells_0&Endothelial_cells_0&Myeloid_0'=29,'BP_cells_0&Endothelial_cells_0&Myeloid_0&T_cells_0'=73,'BP_cells_0&Endothelial_cells_0&Myeloid_0&Normal_epithelial_0'=2,'BP_cells_0&Endothelial_cells_0&Myeloid_0&Normal_epithelial_0&T_cells_0'=0,'BP_cells_0&Endothelial_cells_0&Malignant_cells_0'=22,'BP_cells_0&Endothelial_cells_0&Malignant_cells_0&T_cells_0'=19,'BP_cells_0&Endothelial_cells_0&Malignant_cells_0&Normal_epithelial_0'=1,'BP_cells_0&Endothelial_cells_0&Malignant_cells_0&Normal_epithelial_0&T_cells_0'=0,'BP_cells_0&Endothelial_cells_0&Malignant_cells_0&Myeloid_0'=8,'BP_cells_0&Endothelial_cells_0&Malignant_cells_0&Myeloid_0&T_cells_0'=7,'BP_cells_0&Endothelial_cells_0&Malignant_cells_0&Myeloid_0&Normal_epithelial_0'=0,'BP_cells_0&Endothelial_cells_0&Malignant_cells_0&Myeloid_0&Normal_epithelial_0&T_cells_0'=2,'BP_cells_0&Endothelial_cells_0&Fibroblasts_0'=53,'BP_cells_0&Endothelial_cells_0&Fibroblasts_0&T_cells_0'=58,'BP_cells_0&Endothelial_cells_0&Fibroblasts_0&Normal_epithelial_0'=0,'BP_cells_0&Endothelial_cells_0&Fibroblasts_0&Normal_epithelial_0&T_cells_0'=0,'BP_cells_0&Endothelial_cells_0&Fibroblasts_0&Myeloid_0'=46,'BP_cells_0&Endothelial_cells_0&Fibroblasts_0&Myeloid_0&T_cells_0'=36,'BP_cells_0&Endothelial_cells_0&Fibroblasts_0&Myeloid_0&Normal_epithelial_0'=0,'BP_cells_0&Endothelial_cells_0&Fibroblasts_0&Myeloid_0&Normal_epithelial_0&T_cells_0'=0,'BP_cells_0&Endothelial_cells_0&Fibroblasts_0&Malignant_cells_0'=23,'BP_cells_0&Endothelial_cells_0&Fibroblasts_0&Malignant_cells_0&T_cells_0'=25,'BP_cells_0&Endothelial_cells_0&Fibroblasts_0&Malignant_cells_0&Normal_epithelial_0'=0,'BP_cells_0&Endothelial_cells_0&Fibroblasts_0&Malignant_cells_0&Normal_epithelial_0&T_cells_0'=0,'BP_cells_0&Endothelial_cells_0&Fibroblasts_0&Malignant_cells_0&Myeloid_0'=11,'BP_cells_0&Endothelial_cells_0&Fibroblasts_0&Malignant_cells_0&Myeloid_0&T_cells_0'=4,'BP_cells_0&Endothelial_cells_0&Fibroblasts_0&Malignant_cells_0&Myeloid_0&Normal_epithelial_0'=0,'BP_cells_0&Endothelial_cells_0&Fibroblasts_0&Malignant_cells_0&Myeloid_0&Normal_epithelial_0&T_cells_0'=0,'Adipocytes_0'=82,'Adipocytes_0&T_cells_0'=11,'Adipocytes_0&Normal_epithelial_0'=11,'Adipocytes_0&Normal_epithelial_0&T_cells_0'=0,'Adipocytes_0&Myeloid_0'=37,'Adipocytes_0&Myeloid_0&T_cells_0'=8,'Adipocytes_0&Myeloid_0&Normal_epithelial_0'=3,'Adipocytes_0&Myeloid_0&Normal_epithelial_0&T_cells_0'=2,'Adipocytes_0&Malignant_cells_0'=95,'Adipocytes_0&Malignant_cells_0&T_cells_0'=8,'Adipocytes_0&Malignant_cells_0&Normal_epithelial_0'=32,'Adipocytes_0&Malignant_cells_0&Normal_epithelial_0&T_cells_0'=0,'Adipocytes_0&Malignant_cells_0&Myeloid_0'=16,'Adipocytes_0&Malignant_cells_0&Myeloid_0&T_cells_0'=1,'Adipocytes_0&Malignant_cells_0&Myeloid_0&Normal_epithelial_0'=3,'Adipocytes_0&Malignant_cells_0&Myeloid_0&Normal_epithelial_0&T_cells_0'=0,'Adipocytes_0&Fibroblasts_0'=393,'Adipocytes_0&Fibroblasts_0&T_cells_0'=19,'Adipocytes_0&Fibroblasts_0&Normal_epithelial_0'=25,'Adipocytes_0&Fibroblasts_0&Normal_epithelial_0&T_cells_0'=0,'Adipocytes_0&Fibroblasts_0&Myeloid_0'=68,'Adipocytes_0&Fibroblasts_0&Myeloid_0&T_cells_0'=11,'Adipocytes_0&Fibroblasts_0&Myeloid_0&Normal_epithelial_0'=5,'Adipocytes_0&Fibroblasts_0&Myeloid_0&Normal_epithelial_0&T_cells_0'=1,'Adipocytes_0&Fibroblasts_0&Malignant_cells_0'=108,'Adipocytes_0&Fibroblasts_0&Malignant_cells_0&T_cells_0'=9,'Adipocytes_0&Fibroblasts_0&Malignant_cells_0&Normal_epithelial_0'=19,'Adipocytes_0&Fibroblasts_0&Malignant_cells_0&Normal_epithelial_0&T_cells_0'=0,'Adipocytes_0&Fibroblasts_0&Malignant_cells_0&Myeloid_0'=32,'Adipocytes_0&Fibroblasts_0&Malignant_cells_0&Myeloid_0&T_cells_0'=3,'Adipocytes_0&Fibroblasts_0&Malignant_cells_0&Myeloid_0&Normal_epithelial_0'=2,'Adipocytes_0&Fibroblasts_0&Malignant_cells_0&Myeloid_0&Normal_epithelial_0&T_cells_0'=0,'Adipocytes_0&Endothelial_cells_0'=247,'Adipocytes_0&Endothelial_cells_0&T_cells_0'=16,'Adipocytes_0&Endothelial_cells_0&Normal_epithelial_0'=9,'Adipocytes_0&Endothelial_cells_0&Normal_epithelial_0&T_cells_0'=0,'Adipocytes_0&Endothelial_cells_0&Myeloid_0'=59,'Adipocytes_0&Endothelial_cells_0&Myeloid_0&T_cells_0'=10,'Adipocytes_0&Endothelial_cells_0&Myeloid_0&Normal_epithelial_0'=0,'Adipocytes_0&Endothelial_cells_0&Myeloid_0&Normal_epithelial_0&T_cells_0'=0,'Adipocytes_0&Endothelial_cells_0&Malignant_cells_0'=62,'Adipocytes_0&Endothelial_cells_0&Malignant_cells_0&T_cells_0'=2,'Adipocytes_0&Endothelial_cells_0&Malignant_cells_0&Normal_epithelial_0'=9,'Adipocytes_0&Endothelial_cells_0&Malignant_cells_0&Normal_epithelial_0&T_cells_0'=1,'Adipocytes_0&Endothelial_cells_0&Malignant_cells_0&Myeloid_0'=8,'Adipocytes_0&Endothelial_cells_0&Malignant_cells_0&Myeloid_0&T_cells_0'=0,'Adipocytes_0&Endothelial_cells_0&Malignant_cells_0&Myeloid_0&Normal_epithelial_0'=0,'Adipocytes_0&Endothelial_cells_0&Malignant_cells_0&Myeloid_0&Normal_epithelial_0&T_cells_0'=0,'Adipocytes_0&Endothelial_cells_0&Fibroblasts_0'=402,'Adipocytes_0&Endothelial_cells_0&Fibroblasts_0&T_cells_0'=14,'Adipocytes_0&Endothelial_cells_0&Fibroblasts_0&Normal_epithelial_0'=26,'Adipocytes_0&Endothelial_cells_0&Fibroblasts_0&Normal_epithelial_0&T_cells_0'=0,'Adipocytes_0&Endothelial_cells_0&Fibroblasts_0&Myeloid_0'=102,'Adipocytes_0&Endothelial_cells_0&Fibroblasts_0&Myeloid_0&T_cells_0'=10,'Adipocytes_0&Endothelial_cells_0&Fibroblasts_0&Myeloid_0&Normal_epithelial_0'=7,'Adipocytes_0&Endothelial_cells_0&Fibroblasts_0&Myeloid_0&Normal_epithelial_0&T_cells_0'=0,'Adipocytes_0&Endothelial_cells_0&Fibroblasts_0&Malignant_cells_0'=75,'Adipocytes_0&Endothelial_cells_0&Fibroblasts_0&Malignant_cells_0&T_cells_0'=5,'Adipocytes_0&Endothelial_cells_0&Fibroblasts_0&Malignant_cells_0&Normal_epithelial_0'=10,'Adipocytes_0&Endothelial_cells_0&Fibroblasts_0&Malignant_cells_0&Normal_epithelial_0&T_cells_0'=0,'Adipocytes_0&Endothelial_cells_0&Fibroblasts_0&Malignant_cells_0&Myeloid_0'=16,'Adipocytes_0&Endothelial_cells_0&Fibroblasts_0&Malignant_cells_0&Myeloid_0&T_cells_0'=1,'Adipocytes_0&Endothelial_cells_0&Fibroblasts_0&Malignant_cells_0&Myeloid_0&Normal_epithelial_0'=2,'Adipocytes_0&Endothelial_cells_0&Fibroblasts_0&Malignant_cells_0&Myeloid_0&Normal_epithelial_0&T_cells_0'=0,'Adipocytes_0&BP_cells_0'=8,'Adipocytes_0&BP_cells_0&T_cells_0'=10,'Adipocytes_0&BP_cells_0&Normal_epithelial_0'=1,'Adipocytes_0&BP_cells_0&Normal_epithelial_0&T_cells_0'=0,'Adipocytes_0&BP_cells_0&Myeloid_0'=9,'Adipocytes_0&BP_cells_0&Myeloid_0&T_cells_0'=6,'Adipocytes_0&BP_cells_0&Myeloid_0&Normal_epithelial_0'=1,'Adipocytes_0&BP_cells_0&Myeloid_0&Normal_epithelial_0&T_cells_0'=1,'Adipocytes_0&BP_cells_0&Malignant_cells_0'=9,'Adipocytes_0&BP_cells_0&Malignant_cells_0&T_cells_0'=4,'Adipocytes_0&BP_cells_0&Malignant_cells_0&Normal_epithelial_0'=2,'Adipocytes_0&BP_cells_0&Malignant_cells_0&Normal_epithelial_0&T_cells_0'=0,'Adipocytes_0&BP_cells_0&Malignant_cells_0&Myeloid_0'=6,'Adipocytes_0&BP_cells_0&Malignant_cells_0&Myeloid_0&T_cells_0'=4,'Adipocytes_0&BP_cells_0&Malignant_cells_0&Myeloid_0&Normal_epithelial_0'=0,'Adipocytes_0&BP_cells_0&Malignant_cells_0&Myeloid_0&Normal_epithelial_0&T_cells_0'=0,'Adipocytes_0&BP_cells_0&Fibroblasts_0'=21,'Adipocytes_0&BP_cells_0&Fibroblasts_0&T_cells_0'=5,'Adipocytes_0&BP_cells_0&Fibroblasts_0&Normal_epithelial_0'=0,'Adipocytes_0&BP_cells_0&Fibroblasts_0&Normal_epithelial_0&T_cells_0'=0,'Adipocytes_0&BP_cells_0&Fibroblasts_0&Myeloid_0'=10,'Adipocytes_0&BP_cells_0&Fibroblasts_0&Myeloid_0&T_cells_0'=6,'Adipocytes_0&BP_cells_0&Fibroblasts_0&Myeloid_0&Normal_epithelial_0'=0,'Adipocytes_0&BP_cells_0&Fibroblasts_0&Myeloid_0&Normal_epithelial_0&T_cells_0'=0,'Adipocytes_0&BP_cells_0&Fibroblasts_0&Malignant_cells_0'=6,'Adipocytes_0&BP_cells_0&Fibroblasts_0&Malignant_cells_0&T_cells_0'=8,'Adipocytes_0&BP_cells_0&Fibroblasts_0&Malignant_cells_0&Normal_epithelial_0'=1,'Adipocytes_0&BP_cells_0&Fibroblasts_0&Malignant_cells_0&Normal_epithelial_0&T_cells_0'=0,'Adipocytes_0&BP_cells_0&Fibroblasts_0&Malignant_cells_0&Myeloid_0'=4,'Adipocytes_0&BP_cells_0&Fibroblasts_0&Malignant_cells_0&Myeloid_0&T_cells_0'=2,'Adipocytes_0&BP_cells_0&Fibroblasts_0&Malignant_cells_0&Myeloid_0&Normal_epithelial_0'=1,'Adipocytes_0&BP_cells_0&Fibroblasts_0&Malignant_cells_0&Myeloid_0&Normal_epithelial_0&T_cells_0'=0,'Adipocytes_0&BP_cells_0&Endothelial_cells_0'=8,'Adipocytes_0&BP_cells_0&Endothelial_cells_0&T_cells_0'=9,'Adipocytes_0&BP_cells_0&Endothelial_cells_0&Normal_epithelial_0'=0,'Adipocytes_0&BP_cells_0&Endothelial_cells_0&Normal_epithelial_0&T_cells_0'=1,'Adipocytes_0&BP_cells_0&Endothelial_cells_0&Myeloid_0'=10,'Adipocytes_0&BP_cells_0&Endothelial_cells_0&Myeloid_0&T_cells_0'=10,'Adipocytes_0&BP_cells_0&Endothelial_cells_0&Myeloid_0&Normal_epithelial_0'=0,'Adipocytes_0&BP_cells_0&Endothelial_cells_0&Myeloid_0&Normal_epithelial_0&T_cells_0'=0,'Adipocytes_0&BP_cells_0&Endothelial_cells_0&Malignant_cells_0'=4,'Adipocytes_0&BP_cells_0&Endothelial_cells_0&Malignant_cells_0&T_cells_0'=1,'Adipocytes_0&BP_cells_0&Endothelial_cells_0&Malignant_cells_0&Normal_epithelial_0'=0,'Adipocytes_0&BP_cells_0&Endothelial_cells_0&Malignant_cells_0&Normal_epithelial_0&T_cells_0'=0,'Adipocytes_0&BP_cells_0&Endothelial_cells_0&Malignant_cells_0&Myeloid_0'=2,'Adipocytes_0&BP_cells_0&Endothelial_cells_0&Malignant_cells_0&Myeloid_0&T_cells_0'=2,'Adipocytes_0&BP_cells_0&Endothelial_cells_0&Malignant_cells_0&Myeloid_0&Normal_epithelial_0'=0,'Adipocytes_0&BP_cells_0&Endothelial_cells_0&Malignant_cells_0&Myeloid_0&Normal_epithelial_0&T_cells_0'=0,'Adipocytes_0&BP_cells_0&Endothelial_cells_0&Fibroblasts_0'=28,'Adipocytes_0&BP_cells_0&Endothelial_cells_0&Fibroblasts_0&T_cells_0'=12,'Adipocytes_0&BP_cells_0&Endothelial_cells_0&Fibroblasts_0&Normal_epithelial_0'=1,'Adipocytes_0&BP_cells_0&Endothelial_cells_0&Fibroblasts_0&Normal_epithelial_0&T_cells_0'=0,'Adipocytes_0&BP_cells_0&Endothelial_cells_0&Fibroblasts_0&Myeloid_0'=5,'Adipocytes_0&BP_cells_0&Endothelial_cells_0&Fibroblasts_0&Myeloid_0&T_cells_0'=8,'Adipocytes_0&BP_cells_0&Endothelial_cells_0&Fibroblasts_0&Myeloid_0&Normal_epithelial_0'=0,'Adipocytes_0&BP_cells_0&Endothelial_cells_0&Fibroblasts_0&Myeloid_0&Normal_epithelial_0&T_cells_0'=1,'Adipocytes_0&BP_cells_0&Endothelial_cells_0&Fibroblasts_0&Malignant_cells_0'=7,'Adipocytes_0&BP_cells_0&Endothelial_cells_0&Fibroblasts_0&Malignant_cells_0&T_cells_0'=8,'Adipocytes_0&BP_cells_0&Endothelial_cells_0&Fibroblasts_0&Malignant_cells_0&Normal_epithelial_0'=0,'Adipocytes_0&BP_cells_0&Endothelial_cells_0&Fibroblasts_0&Malignant_cells_0&Normal_epithelial_0&T_cells_0'=0,'Adipocytes_0&BP_cells_0&Endothelial_cells_0&Fibroblasts_0&Malignant_cells_0&Myeloid_0'=1,'Adipocytes_0&BP_cells_0&Endothelial_cells_0&Fibroblasts_0&Malignant_cells_0&Myeloid_0&T_cells_0'=1,'Adipocytes_0&BP_cells_0&Endothelial_cells_0&Fibroblasts_0&Malignant_cells_0&Myeloid_0&Normal_epithelial_0'=0,'Adipocytes_0&BP_cells_0&Endothelial_cells_0&Fibroblasts_0&Malignant_cells_0&Myeloid_0&Normal_epithelial_0&T_cells_0'=0)
upset(fromExpression(expressionInput), nsets=8, nintersects=30, show.numbers="yes", main.bar.color="#ea5d4e", sets.bar.color="#317eab", empty.intersections=NULL, order.by = "freq", number.angles = 0, mainbar.y.label ="No. of Intersections", sets.x.label ="Set size")
invisible(dev.off())

########################################
###### Supplementary Figure 14 #########
########################################
###############################################################################################
############# SUPPLEMENTARY FIGURE 14A - SWARMPLOT PROXIMAL DISTAL ENDOBLAST       #############
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

d1 <- read.table("SCAFE/enichment_proximal_ids/improved_main_ids__Endothelial_cells__proximal.tsv", sep = '\t', header = TRUE, fill = TRUE)
d2 <- read.table("SCAFE/enrichment_distals_ids/improved_main_ids__Endothelial_cells.tsv", sep = '\t', header = TRUE, fill = TRUE)
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
      high = "darkgreen",  
      low = "green3",  
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
      high = "darkgreen",  
      low = "green3",  
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
      high = "darkgreen",  
      low = "green3",  
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

pdf('SuppFigure14A Endothelial_proximal_vs_distal_double.pdf')
print(p)
print(p1)
print(p2)
dev.off()



####################################################################################################
############# SUPPLEMENTARY FIGURE 14B - SWARMPLOT PROXIMAL DISTAL B cells plasmablasts #############
####################################################################################################

library(EnhancedVolcano)
library(ggrepel)
library(ggplot2)
library(miloR)
library(plyr)
library(dplyr)
library(future)
library(stringi)
library(tidyverse)


d1 <- read.table("SCAFE/enichment_proximal_ids/improved_main_ids__BP_cells__proximal.tsv", sep = '\t', header = TRUE, fill = TRUE)
d2 <- read.table("SCAFE/enrichment_distals_ids/improved_main_ids__BP_cells.tsv", sep = '\t', header = TRUE, fill = TRUE)
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
      high = "darkorange3",  
      low = "orange", ,   
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
      high = "darkorange3",  
      low = "orange", ,   
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

pdf('SuppFigure14B BP cells_proximal_vs_distal_double.pdf')
print(p)
print(p1)
print(p2)
dev.off()


########################################
###### Supplementary Figure 15 #########
########################################
###############################################################################################################################
############# Supplementary Figure 15 AC - SWARMPLOT - Baseline - Crossover  - Sensitive - Resistant   ########################
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
## files for resistant clones Supp fig 15 A 
#d2 <- read.table("clones_res_A_vs_B__distal_cre.tsv", sep = '\t', header = TRUE, fill = TRUE)
#d1 <- read.table("clones_res_B_vs_A__distal_cre.tsv", sep = '\t', header = TRUE, fill = TRUE)

## files for sensitive clones Supp fig 15 C 
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

#pdf ("SuppFigure15A swarmplot resistant_clones A_v_B.pdf")
pdf ("SuppFigure15C swarmplot sensitive_clones A_v_B.pdf.pdf")
p
p1
p2  
dev.off()


###############################################################################################################################
############# Supplementary Figure 15 BD - SWARMPLOT - Baseline - EOT  - Sensitive - Resistant   ########################
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
## files for resistant clones Supp fig 15 B 
#d2 <- read.table("clones_res_A_vs_C__distal_cre.tsv", sep = '\t', header = TRUE, fill = TRUE)
#d1 <- read.table("clones_res_C_vs_A__distal_cre.tsv", sep = '\t', header = TRUE, fill = TRUE)
## files for resistant clones Supp fig 15 D 
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

#pdf ("SuppFigure15B swarmplot resistant_clones A_v_C.pdf")
pdf ("SuppFigure15D swarmplot sensitive_clones A_v_C.pdf.pdf")
p
p1
p2  
dev.off()


########################################
###### Supplementary Figure 16 #########
########################################
#######################################################################################
############# SUPP FIG 16A-D - Fraction of enhancer active with ESR1 TFBS across time ##
#######################################################################################
library(Seurat)
library(plyr)
library(dplyr)
library(future)
library (ggplot2)
library (ggrepel)
library (stringi)
library(tidyverse)
library(ggsignif)

### available from PP_snakemake_scafe

dat <- readRDS("response__main.ids.cor__Malignant_cells__ESR1__all__distal_cre.rds") ## file for resistant clones Supp fig 16 A-B
#dat <- readRDS("response__main.ids.cor__Malignant_cells__AR__all__distal_cre.rds") ## file for resistant clones Supp fig 16 C-D
						  
						  
dat$ESR1_nCells <- dat$active_set / dat$active_CRE
#dat_var <- dat[dat$group=="GR",] ### good responders
dat_var <- dat[dat$group=="BR",] ### poor responders

pAB <- wilcox.test(dat_var[dat_var$time_point=="A", colnames (dat_var)=="ESR1_nCells"], dat_var[dat_var$time_point=="B", colnames (dat_var)=="ESR1_nCells"], paired = F)[3]
pBC <- wilcox.test(dat_var[dat_var$time_point=="B", colnames (dat_var)=="ESR1_nCells"], dat_var[dat_var$time_point=="C", colnames (dat_var)=="ESR1_nCells"], paired = F)[3]
pAC <- wilcox.test(dat_var[dat_var$time_point=="A", colnames (dat_var)=="ESR1_nCells"], dat_var[dat_var$time_point=="C", colnames (dat_var)=="ESR1_nCells"], paired = F)[3]
pvalues <- c(pAB, pBC, pAC)

outline_colors <- c("#AA3830", '#9D2E99', '#350D7F') # Box outline colors
box_fill_colors <- c('#FFE2E0', '#FFD1FF', '#DAC3FF') # Box fill colors
point_fill_colors <- c('#FF8077', '#E66CE1', '#7E33FF') # Point fill colors

#pdf('SuppFigure16A Boxplot_ESR1_GR.pdf')
pdf('SuppFigure16B Boxplot_ESR1_BR.pdf')
#pdf('SuppFigure16C Boxplot_AR_GR.pdf')
#pdf('SuppFigure16D Boxplot_AR_BR.pdf')
a <- ggplot(dat_var, aes(time_point, ESR1_nCells, label = patient_number)) + 
 geom_boxplot(color = outline_colors, 
               fill = box_fill_colors, 
               alpha = 0.5, size = 1.5) + 
  geom_point(aes(fill = factor(time_point), colour = factor(time_point)), 
             shape = 21, size = 5, stroke = 1.2) + 
  scale_fill_manual(values = point_fill_colors) +
  scale_colour_manual(values = outline_colors) + 
  ylab('Fraction of active enhancers with ESR1 TFBS')+ 
  #ylab('Fraction of active enhancers with AR TFBS')+ 
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
              map_signif_level = FALSE,
              test = "wilcox.test", 
              step_increase = 0.07)
print(a)
dev.off()




########################################
###### Supplementary Figure 17 #########
########################################
#######################################################################################
############# SUPP FIG 17 A-B - HEATMAP RESISTANT/SENSITIVE CLONES VS ARM #############
#######################################################################################

library(Seurat)
library(plyr)
library(dplyr)
library(future)
library (mixtools)
library (pheatmap)
library (amap)
library (gwasrapidd)
library (ggplot2)
library(ggplot2)
library(reshape2)
library(Seurat)
library(ggplot2)
library(ggrepel)
library (msigdbr)
library (dplyr)
library (clusterProfiler)
library (dplyr)

scores <- as.data.frame(data.frame(read.table("scores_gene_cluster_clones.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)))
scores$time_point <- substr(rownames (scores), 4, 4)
#scores generated  in PP_gene_cluster_correlation
scores_B <- scores[scores$time_point=="B",]

load("metadata_git.RData")
 library (dplyr)
 df_averages <- data.frame (scores_B %>% 
  group_by(clone_id) %>%
  summarise(across(starts_with("cluster"), mean, na.rm = TRUE))
)

rownames (df_averages) <- df_averages[,1]
df_averages <- t(df_averages[, -1])
df_averages <- t(scale (t(df_averages)))


clones_ord <-c( unique (metadata$clone_id[metadata$clone_cat2=="res" &!is.na(metadata$clone_cat2)]), unique (metadata$clone_id[metadata$clone_cat2=="sens" &!is.na(metadata$clone_cat2)]))
clones_cat <-c( rep ('res' , length (unique (metadata$clone_id[metadata$clone_cat2=="res" &!is.na(metadata$clone_cat2)]))), rep ('sens', length (unique (metadata$clone_id[metadata$clone_cat2=="sens" &!is.na(metadata$clone_cat2)]))))

colnames (df_averages) <- gsub ("-", "_", colnames (df_averages))
clones_cat <- clones_cat [clones_ord%in%colnames (df_averages)]
clones_ord <- clones_ord [clones_ord%in%colnames (df_averages)]
df_averages <- df_averages[, match(clones_ord, colnames(df_averages))]

mean_mat <- as.matrix (df_averages)

breaksList = seq(-3, 3, by = 0.05)
col.pal <- colorRampPalette(c("dodgerblue4","white","brown3"))(length(breaksList))

CL <- data.frame(cbind (clones_ord, clones_cat))
annotation_col <- data.frame(CL[,-1])
rownames (annotation_col) <- CL$clones_ord
colnames (annotation_col) <- 'clones_cat'

annotation_col$arm <- ifelse (gsub("^.|_.*", "", rownames (annotation_col)) %in% c('02', '05', '07', '09', '11', '13',  '15', '17', '20', '22', '25'), "EXE_LET",gsub("^.|_.*", "", rownames (annotation_col)))
annotation_col$arm <- ifelse (annotation_col$arm%in% c('01', '03', '04', '06', '08', '10',  '12', '14', '16', '18', '19','21','23','24'), "LET_EXE", annotation_col$arm)
m <- ifelse (annotation_col$arm%in% c('01', '03', '04', '06', '08', '10',  '12', '14', '16', '18', '19','21','23','24'), "LET_EXE", annotation_col$arm)

anno_colors = list(  arm = c("EXE_LET" = "#2A6FB0", "LET_EXE" = "#E86F00"),
clones_cat = c("sens" = "green", "res" = "red"))

############# SUPP FIG 17A - HEATMAP RESISTANT CLONES VS ARM ##########################

mean_mat_res <- mean_mat[, annotation_col$clones_cat=="res"]
annotation_col_res <- annotation_col[annotation_col$clones_cat=="res",]

dcols <- hcluster(t(mean_mat_res), method = 'pearson', link = 'complete')
drows <- hcluster(mean_mat_res, method = 'euclidean', link = 'ward')

resistant_heatmap <- pheatmap(mean_mat_res, color = col.pal, show_rownames = T, show_colnames = T, main = "",  cluster_rows = drows, cluster_cols = dcols, 
breaks = breaksList, annotation_col=annotation_col_res, annotation_colors = anno_colors, fontsize = 10, border_color = NA)

pdf('SuppFigure17A Heatmap_resistant_clones_scores_arm.pdf')
print (resistant_heatmap)
dev.off()

############# SUPP FIG 17B - HEATMAP SENSITIVE CLONES VS ARM ##########################


mean_mat_sens <- mean_mat[, annotation_col$clones_cat=="sens"]
annotation_col_sens <- annotation_col[annotation_col$clones_cat=="sens",]

dcols <- hcluster(t(mean_mat_sens), method = 'pearson', link = 'complete')
drows <- hcluster(mean_mat_sens, method = 'euclidean', link = 'ward')

sensitive_heatmap <- pheatmap(mean_mat_sens, color = col.pal, show_rownames = T, show_colnames = T, main = "",  cluster_rows = drows, cluster_cols = dcols, 
breaks = breaksList, annotation_col=annotation_col_sens, annotation_colors = anno_colors, fontsize = 10, border_color = NA)

pdf('SuppFigure17B Heatmap_sensitive_clones_scores_arm.pdf')
print (sensitive_heatmap)
dev.off()


#######################################################################################
############# SUPP FIG 17C-E - HEATMAP TME VS ARM THREE TIME POINTS ###################
#######################################################################################

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
met_clin <- data.frame(read.table("clinical_data.txt", header=T, sep="\t", stringsAsFactors=F))
##### needs clinical data in supp Table 1

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
dta$arm <- ifelse (dta$paired  %in% c('02', '05', '07', '09', '11', '13',  '15', '17', '20', '22', '25'), "EXE_LET",dta$paired )
dta$arm <- ifelse (dta$arm%in% c('01', '03', '04', '06', '08', '10',  '12', '14', '16', '18', '19','21','23','24'), "LET_EXE", dta$arm)

dta$category <- stri_sub(rownames (dta),-1)
annot <- metadata[, colnames (metadata) %in% c("orig.ident", "perc_malig_surg10")]
annot <- annot[!duplicated (metadata$orig.ident), ]
table (annot$orig.ident == rownames (dta))
dta <- cbind (dta , annot)

############# SUPP FIG 17C - HEATMAP TME VS ARM BASELINE  ##########################
dta1_A <- dta1[dta$category == "A",]
dta_A <- dta[dta$category == "A",]

annotation_row <- data.frame(cbind(dta_A$perc_malig_surg10, dta_A$arm))
rownames(annotation_row) <- rownames(dta_A)
colnames(annotation_row) <- c('Response', 'arm')

max_val <- 20
breaks_vals <- c(seq(0, max_val, length.out = 101))
breaksList = seq(0, 35, by = 0.05)
colors <- colorRampPalette(c("white", "indianred1", "brown3", "brown4"))(length(breaksList))
col.pal <- colorRampPalette(c("dodgerblue4", "white", "white", "brown3", "brown3", "brown3"))(length(breaksList))

anno_colors = list(
  Response = c("GR" = "green", "BR" = "red"),
  arm = c("EXE_LET" = "#2A6FB0", "LET_EXE" = "#E86F00")
)

dcols <- hcluster(t(dta1_A), method = 'spearman', link = 'complete')
drows <- hcluster(dta1_A, method = 'spearman', link = 'complete')

pdf ('SuppFigure17C Heatmap_TME_ARM_Baseline.pdf')
res <- pheatmap(
  t(dta1_A), 
  show_rownames = TRUE, show_colnames = TRUE,  
  cluster_rows = dcols, cluster_cols = drows, fontsize = 9,
  annotation_col = annotation_row, annotation_colors = anno_colors,
  col = colors, 
   width = 10, height = 8 , border_color = NA, main = 'Baseline'
)
print (res)
dev.off()

############# SUPP FIG 17D - HEATMAP TME VS ARM CO  ##########################

dta1_B <- dta1[dta$category == "B",]
dta_B <- dta[dta$category == "B",]

annotation_row <- data.frame(cbind(dta_B$perc_malig_surg10, dta_B$arm))
rownames(annotation_row) <- rownames(dta_B)
colnames(annotation_row) <- c('Response', 'arm')

dcols <- hcluster(t(dta1_B), method = 'spearman', link = 'complete')
drows <- hcluster(dta1_B, method = 'spearman', link = 'complete')

pdf ('SuppFigure17D Heatmap_TME_ARM_CO.pdf')
res <- pheatmap(
  t(dta1_B), 
  show_rownames = TRUE, show_colnames = TRUE,  
  cluster_rows = dcols, cluster_cols = drows, fontsize = 9,
  annotation_col = annotation_row, annotation_colors = anno_colors,
  col = colors, 
   width = 10, height = 8 , border_color = NA, main = 'CO'
)
print (res)
dev.off()


############# SUPP FIG 17E - HEATMAP TME VS ARM EOT  ##########################

dta1_C <- dta1[dta$category == "C",]
dta_C <- dta[dta$category == "C",]

annotation_row <- data.frame(cbind(dta_C$perc_malig_surg10, dta_C$arm))
rownames(annotation_row) <- rownames(dta_C)
colnames(annotation_row) <- c('Response', 'arm')

dcols <- hcluster(t(dta1_C), method = 'spearman', link = 'complete')
drows <- hcluster(dta1_C, method = 'spearman', link = 'complete')

pdf ('SuppFigure17E Heatmap_TME_ARM_EOT.pdf')

res <- pheatmap(
  t(dta1_C), 
  show_rownames = TRUE, show_colnames = TRUE,  
  cluster_rows = dcols, cluster_cols = drows, fontsize = 9,
  annotation_col = annotation_row, annotation_colors = anno_colors,
  col = colors, 
   width = 10, height = 8 , border_color = NA , main = 'EOT'
)
print (res)
dev.off()

#######################################################################################
############# SUPP FIG 17F-H - BOXPLOT EXPANSION CONTRACTION ARM ######################
#######################################################################################


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
load("metadata_git.RData")

merA <- data.frame(read.table("clinical_data.txt", header=T, sep='\t', row.names=))
merA <- merA [!merA$patient_number=="4",]

norm_EC_AB <- cbind (statusAB1_n_EC, merA$response, merA$patient_numb, merA$AI.sequence)
colnames (norm_EC_AB) <- c('contraction','lit_expansion', 'str_expansion',  'response', 'patient_number', 'arm')
norm_EC_AB$category <- 'AB'


combined_data <- rbind(norm_EC_AB)
outline_colors    <- c("#1F4E79", "#A05100")  # darker outlines: blue, orange
box_fill_colors   <- c("#BFD7F6", "#FFDAB3")  # softer/lighter fills
point_fill_colors <- c("#2A6FB0", "#E86F00")  # medium/intense points
		 
############# SUPP FIG 17F - BOXPLOT CONTRACTION VS ARM BASELINE-CO  ##########################
		 
pdf	('SuppFigure17F contraction_AB_arm.pdf')	   
combined_data <- combined_data %>% filter(!is.na(arm) & !is.na(contraction))
pval_AB <- wilcox.test(norm_EC_AB$contraction[norm_EC_AB$arm == 'EXE-LET'],
                       norm_EC_AB$contraction[norm_EC_AB$arm == 'LET-EXE'], 
                       paired = FALSE)$p.value
ggplot(combined_data, aes(x = arm, y = contraction)) + 
  geom_boxplot(aes(fill = arm, colour = arm), 
               alpha = 0.2, outlier.colour = NA, notch = FALSE, size = 1.5) + 
  geom_point(aes(fill = arm, colour = arm), shape = 21, size = 5, , stroke = 1.2) + 
  scale_fill_manual(values = setNames(c(box_fill_colors, point_fill_colors), 
                                      levels(combined_data$arm))) +
  scale_colour_manual(values = setNames(c(outline_colors), 
                                         levels(combined_data$arm))) +  
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
  labs(x = 'arm') +  
  geom_signif(data = subset(combined_data, category == 'AB'), comparisons = list(c('EXE-LET', 'LET-EXE')),
               annotations = sprintf("p = %.3f", pval_AB),
               map_signif_level = FALSE,
               test = "wilcox.test", 
               step_increase = 0.1)
dev.off()


############# SUPP FIG 17G - BOXPLOT EXPANSION VS ARM BASELINE-CO  ##########################
pdf	('SuppFigure17G lit_expansion_AB_arm.pdf')	
combined_data <- combined_data %>% filter(!is.na(arm) & !is.na(lit_expansion))
pval_AB <- wilcox.test(norm_EC_AB$lit_expansion[norm_EC_AB$arm == 'EXE-LET'],
                       norm_EC_AB$lit_expansion[norm_EC_AB$arm == 'LET-EXE'], 
                       paired = FALSE)$p.value

ggplot(combined_data, aes(x = arm, y = lit_expansion)) + 
  geom_boxplot(aes(fill = arm, colour = arm), 
               alpha = 0.2, outlier.colour = NA, notch = FALSE, size = 1.5) + 
  geom_point(aes(fill = arm, colour = arm), shape = 21, size = 5, , stroke = 1.2) + 
  scale_fill_manual(values = setNames(c(box_fill_colors, point_fill_colors), 
                                      levels(combined_data$arm))) +
  scale_colour_manual(values = setNames(c(outline_colors), 
                                         levels(combined_data$arm))) + 
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
  labs(x = 'arm') +  
  geom_signif(data = subset(combined_data, category == 'AB'), comparisons = list(c('EXE-LET', 'LET-EXE')),
               annotations = sprintf("p = %.3f", pval_AB),
               map_signif_level = FALSE,
               test = "wilcox.test", 
               step_increase = 0.1)
dev.off()

############# SUPP FIG 17H - BOXPLOT STRONG EXPANSION VS ARM BASELINE-CO  ##########################
pdf	('SuppFigure17H str_expansion_AB_arm.pdf')	 
combined_data <- combined_data %>% filter(!is.na(arm) & !is.na(str_expansion))
pval_AB <- wilcox.test(norm_EC_AB$str_expansion[norm_EC_AB$arm == 'EXE-LET'],
                       norm_EC_AB$str_expansion[norm_EC_AB$arm == 'LET-EXE'], 
                       paired = FALSE)$p.value
ggplot(combined_data, aes(x = arm, y = str_expansion)) + 
  geom_boxplot(aes(fill = arm, colour = arm), 
               alpha = 0.2, outlier.colour = NA, notch = FALSE, size = 1.5) + 
  geom_point(aes(fill = arm, colour = arm), shape = 21, size = 5, , stroke = 1.2) + 
  scale_fill_manual(values = setNames(c(box_fill_colors, point_fill_colors), 
                                      levels(combined_data$arm))) +
  scale_colour_manual(values = setNames(c(outline_colors), 
                                         levels(combined_data$arm))) + 
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
  geom_signif(data = subset(combined_data, category == 'AB'), comparisons = list(c('EXE-LET', 'LET-EXE')),
               annotations = sprintf("p = %.3f", pval_AB),
               map_signif_level = FALSE,
               test = "wilcox.test", 
               step_increase = 0.1) 
dev.off()



########################################
###### Supplementary Figure 18 #########
########################################
###############################################
######## SUPP FIG 18A EXE-LET - A to B ########
###############################################

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
d2 <- read.table("SCAFE/arm_EXE-LET_A_vs_B__distal_cre.tsv", sep = '\t', header = TRUE, fill = TRUE)
d1 <- read.table("SCAFE/arm_EXE-LET_B_vs_A__distal_cre.tsv", sep = '\t', header = TRUE, fill = TRUE)

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

selected_collections <- c('ESR1', 'TP53', 'GATA3', 'FOXA1', 'AR')

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
      annotate("text", x = x, y = max(d1$logFC) + 1, label = plus_label, hjust = 1, size = 3),
      annotate("text", x = x, y = min(d1$logFC) - 1, label = minus_label, hjust = 1, size = 3)
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

pdf ("SuppFigure18A SCAFE_EXE_LET_A_to_B.pdf")
p
p1
p2  
dev.off()


###############################################
######## SUPP FIG 18B LET-EXE - A to B #########
###############################################

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
d2 <- read.table("SCAFE/arm_LET-EXE_A_vs_B__distal_cre.tsv", sep = '\t', header = TRUE, fill = TRUE)
d1 <- read.table("SCAFE/arm_LET-EXE_B_vs_A__distal_cre.tsv", sep = '\t', header = TRUE, fill = TRUE)

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

selected_collections <- c('ESR1', 'TP53', 'GATA3', 'FOXA1', 'AR')

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
      annotate("text", x = x, y = max(d1$logFC) + 1, label = plus_label, hjust = 1, size = 3),
      annotate("text", x = x, y = min(d1$logFC) - 1, label = minus_label, hjust = 1, size = 3)
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


pdf ("SuppFigure18B SCAFE_LET_EXE_A_to_B.pdf")
p
p1
p2  
dev.off()



###############################################
######## SUPP FIG 18D LET-EXE - A to C #########
###############################################

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
d2 <- read.table("SCAFE/arm_LET-EXE_A_vs_C__distal_cre.tsv", sep = '\t', header = TRUE, fill = TRUE)
d1 <- read.table("SCAFE/arm_LET-EXE_C_vs_A__distal_cre.tsv", sep = '\t', header = TRUE, fill = TRUE)

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

selected_collections <- c('ESR1', 'AR')

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
  p$data$logFC_color <- ifelse((p$data$logFC_color == -2 | p$data$logFC_color == 2), p$data$logFC_color, 0)
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
  p1$data$logFC_color <- ifelse((p1$data$logFC_color == -2 | p1$data$logFC_color == 2), p1$data$logFC_color, 0)
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
  p2$data$logFC_color <- ifelse((p2$data$logFC_color == -2 | p2$data$logFC_color == 2), p2$data$logFC_color, 0)
  p2$data <- p2$data[order(!is.na(p2$data$logFC_color), p2$data$logFC_color), ]

p2 <- p2 + mapply(function(x, plus_label, minus_label) { 
    list(
      annotate("text", x = x, y = max(d1$logFC) + 1, label = plus_label, hjust = 1, size = 3),
      annotate("text", x = x, y = min(d1$logFC) - 1, label = minus_label, hjust = 1, size = 3)
    ) 
  }, seq_along(mer1_filtered$Var1), paste ("(", rev(mer1_filtered$plus), ")"), paste ("(",rev(mer1_filtered$minus), ")"))


pdf ("SuppFigure18D SCAFE_LET_EXE_A_to_C.pdf")
p
p1
p2  
dev.off()




################################################
######## SUPP FIG 18C EXE-LET - A to C #########
################################################
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
d2 <- read.table("SCAFE/arm_EXE-LET_A_vs_C__distal_cre.tsv", sep = '\t', header = TRUE, fill = TRUE)
d1 <- read.table("SCAFE/arm_EXE-LET_C_vs_A__distal_cre.tsv", sep = '\t', header = TRUE, fill = TRUE)

d2$oddsRatio <- (-d2$oddsRatio)

  d <- rbind(d1[,colnames(d1) %in% c('collection', 'pValueLog', 'oddsRatio')], 
             d2[,colnames(d2) %in% c('collection', 'pValueLog', 'oddsRatio')])
  
  max <- max(d$oddsRatio[!d$oddsRatio %in% c('Inf')])
  min <- min(d$oddsRatio[!d$oddsRatio %in% c('-Inf')])
  
  d$pValue <- 10^(-d$pValueLog)

filtered_data <- subset(d, (pValue < 0.05 & oddsRatio > 1.5) | pValue < 0.05 & oddsRatio < -1.5 )
collection_count <- table(filtered_data$collection)

print(head (colnames(d)))
collection_count_plus <- as.data.frame(table (d$collection[d$oddsRatio > 1.75 & d$pValue < 0.05])) ; colnames (collection_count_plus)[2] <- 'plus'
collection_count_minus <- as.data.frame(table (d$collection[!d$oddsRatio > -1.5 & d$pValue < 0.05])) ;colnames (collection_count_minus)[2] <- 'minus'

selected_collections <- c('ESR1', 'AR')

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

pdf ("SuppFigure18C SCAFE_EXE_LET_A_to_C.pdf")
p
p1
p2  
dev.off()

########################################
###### Supplementary Figure 19 #########
########################################
#######################################################################################
############# SUPP FIG 19 - ALLUVIAL PLOT CHANGE IN SUBTYPE between baseline and EOT ##
#######################################################################################
library(ggplot2)
library(ggalluvial)

met <- read.table("clinical_data.txt", header=T, sep='\t')
##### load clinical data in Supp Table 1
met1 <- met[!is.na(met$BAS_subtype) & !is.na(met$EOT_subtype),]

tab <- as.data.frame (table (met1$BAS_subtype, met1$EOT_subtype))
colnames(tab) <- c("BAS_subtype", "EOT_subtype", "Count")

pdf('SuppFigure19 alluvial_plot_subtype_change_PAM50.pdf')

ggplot(tab, aes(axis1 = BAS_subtype, axis2 = EOT_subtype, y = Count)) +
  geom_alluvium(aes(fill = BAS_subtype)) + 
  geom_stratum(width = 0.2) + 
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3, 
            vjust = -0.5, hjust = 0.5) +  
  scale_x_discrete(limits = c("Baseline Subtype", "End of Treatment Subtype"), expand = c(0.05, 0.05)) +
  scale_fill_manual(values = c("LUM-A" = "blue",  "LUM-B" = "light blue", "HER-2-enriched" = "deep pink")) + 
  theme_minimal() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 14, face = "bold", angle = 0, hjust = 0.5, margin = margin(t = 0, r = 5)), 
        axis.text.y = element_text(size = 10, face = "bold"),  
        axis.title.x = element_text(size = 16, face = "bold", margin = margin(t = 10)), 
        axis.title.y = element_text(size = 16, face = "bold"),  
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5)) +  
  labs(title = "Changes in Subtypes from Baseline to End of Treatment", 
       x = "Subtype Status", 
       y = "Counts") +
  geom_alluvium(aes(group = interaction(BAS_subtype, EOT_subtype)), 
                fill = NA, color = "black", size = 0.5)  

dev.off()

########################################
###### Supplementary Figure 20 #########
########################################
########################################################################
############# SUPP FIG 20A - UMAP INFERCNV with normal cells ##
########################################################################

load("annot_umap_git.RData")
#### Need InfercNV UMAP and annotation object obtained from PP_InferCNV

annot_umap$patient_col <- ifelse(annot_umap$patient == "red", "lavender", annot_umap$patient_col)
pdf ('SuppFigure20A_umap_CNV_with_normal.pdf')
ordered_patients <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "red")
ordered_colors <- annot_umap$patient_col[match(ordered_patients, annot_umap$patient)]
plot(mat.umap$layout, col = annot_umap$patient_col, cex = .5, xlim = c(-18, 18), 
     xlab = expression(bold("UMAP 1")), ylab = expression(bold("UMAP 2")))
legend("topright", legend = ordered_patients, pch = 16, col = ordered_colors, cex = 0.9, bty = "n")
dev.off()


############################################################################
############# SUPP FIG 20B - CNA value in malignant and normal epithelial ###
############################################################################

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

pdf('SuppFigure20B_CNA_value_normal_malignant.pdf')
metadata_m <- metadata[metadata$main.ids%in%c('Epithelial_cells'),]
par(mar = c(10, 4, 4, 2) + 0.1)  

colors <- c('purple', '#D8BFD8',  
            'green', '#90EE90',  
            'pink', '#FFDDDD')   


boxplot(CNA_value ~ interaction(main.ids.cor, time_point), data = metadata_m,
        main = "Boxplot of CNA_value by Patient and Group",  
        ylab = "CNA_value",  
        xlab = "",  
        col = colors,  
        border = c('black', 'purple', 'black', 'darkgreen', 'black', 'deeppink'),  # Box border color
        outcol = colors,  
        las = 2,  
        ylim = c(0, 0.5),
		xaxt = 'n',  
        cex.axis = 1.2,  
        cex.lab = 1.4,  
        cex.main = 1.5,  
        font.main = 2,  
        font.lab = 2)  
dev.off()

interaction_labels <- with(metadata_m, interaction(main.ids.cor, time_point))
axis(1, at = seq_along(levels(interaction_labels)), labels = levels(interaction_labels), las = 2, 
     cex.axis = 1.1, font = 2)  

pvalA <- t.test(metadata$CNA_value[metadata$main.ids.cor == "Malignant_cells" & metadata$time_point=="A"], 
metadata$CNA_value[metadata$main.ids.cor == "Normal_epithelial"& metadata$time_point=="A"])[3]

pvalB <- t.test(metadata$CNA_value[metadata$main.ids.cor == "Malignant_cells" & metadata$time_point=="B"], 
metadata$CNA_value[metadata$main.ids.cor == "Normal_epithelial"& metadata$time_point=="B"])[3]

pvalC <- t.test(metadata$CNA_value[metadata$main.ids.cor == "Malignant_cells" & metadata$time_point=="C"], 
metadata$CNA_value[metadata$main.ids.cor == "Normal_epithelial"& metadata$time_point=="C"])[3]


############################################################################
############# SUPP FIG 21 - ELBOW and GAP plots ###
############################################################################
library(cluster)
library(factoextra)

b <- data.frame(read.table("correl_matrix_sensitiv_resistant_genes.txt", header=T, sep="\t", row.names=))

pdf ('SuppFigure21A elbow_plot_correlation_matrix.pdf')
fviz_nbclust(b, kmeans, method = "wss",  k.max = 30) + geom_vline(xintercept = 10, linetype = 2)+ labs(subtitle = "Elbow method")
dev.off()

pdf ('SuppFigure21B gap_plot_correlation_matrix.pdf')
plot <- fviz_nbclust(b, kmeans, nstart = 25,  method = "gap_stat", nboot = 20,  k.max = 30)+ geom_vline(xintercept = 10, linetype = 2) + labs(subtitle = "Gap statistic method") 
plot$layers <- plot$layers[-which(sapply(plot$layers, function(x) "GeomVline" %in% class(x$geom)))]
print(plot)
dev.off()















