########################################
###### Supplementary Figure 11 #########
########################################

################################################################################################################
############# SUPPLEMENTARY FIGURE 11A- DINUCLEOTIDES    #######################################################
################################################################################################################

#### code from Anthony

################################################################################################################
############# SUPPLEMENTARY FIGURE 11B- DINUCLEOTIDES     ######################################################
################################################################################################################

#### code from Anthony


################################################################################################################
############# SUPPLEMENTARY FIGURE 11C- BARPLOT INTERSECT ATAC chromHMM     ####################################
################################################################################################################

#### code from Villlads

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

pdf("SuppFigure11C Intersection with scATAC-seq.pdf", width=14, height=8, onefile=FALSE, useDingbats=FALSE)

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
############# SUPPLEMENTARY FIGURE 11D-E - UPSET PLOTS - proximal distal     ################################
################################################################################################################

#### code from Villlads

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
###### Supplementary Figure 12 #########
########################################
###############################################################################################
############# SUPPLEMENTARY FIGURE 12A - SWARMPLOT PROXIMAL DISTAL ENDOBLAST       #############
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

pdf('SuppFigure12A Endothelial_proximal_vs_distal_double.pdf')
print(p)
print(p1)
print(p2)
dev.off()



###############################################################################################
############# SUPPLEMENTARY FIGURE 12B - SWARMPLOT PROXIMAL DISTAL ENDOBLAST       #############
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

pdf('SuppFigure12B BP cells_proximal_vs_distal_double.pdf')
print(p)
print(p1)
print(p2)
dev.off()

###############################################################################################
############# SUPPLEMENTARY FIGURE 12C - UPSET PLOT of TFBS in hormone sensitive enhancers ####
###############################################################################################

### code from Villada

########################################
###### Supplementary Figure 13 #########
########################################
#######################################################################################
############# SUPP FIG 13 - ALLUVIAL PLOT CHANGE IN SUBTYPE between baseline and EOT ##
#######################################################################################
library(ggplot2)
library(ggalluvial)

met <- read.table("clinical_data.txt", header=T, sep='\t')
##### load clinical data in Supp Table 1
met1 <- met[!is.na(met$BAS_subtype) & !is.na(met$EOT_subtype),]

tab <- as.data.frame (table (met1$BAS_subtype, met1$EOT_subtype))
colnames(tab) <- c("BAS_subtype", "EOT_subtype", "Count")

pdf('SuppFigure13 alluvial_plot_subtype_change_PAM50.pdf')

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
###### Supplementary Figure 14 #########
########################################
########################################################################
############# SUPP FIG 14A - UMAP INFERCNV with normal cells ##
########################################################################

load("annot_umap_git.RData")
#### Need InfercNV UMAP and annotation object obtained from PP_InferCNV

annot_umap$patient_col <- ifelse(annot_umap$patient == "red", "lavender", annot_umap$patient_col)
pdf ('SuppFigure14A_umap_CNV_with_normal.pdf')
ordered_patients <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "red")
ordered_colors <- annot_umap$patient_col[match(ordered_patients, annot_umap$patient)]
plot(mat.umap$layout, col = annot_umap$patient_col, cex = .5, xlim = c(-18, 18), 
     xlab = expression(bold("UMAP 1")), ylab = expression(bold("UMAP 2")))
legend("topright", legend = ordered_patients, pch = 16, col = ordered_colors, cex = 0.9, bty = "n")
dev.off()


############################################################################
############# SUPP FIG 14B - CNA value in malignant and normal epithelial ###
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

pdf('SuppFigure14B_CNA_value_normal_malignant.pdf')
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
############# SUPP FIG 15 - ELBOW and GAP plots ###
############################################################################
library(cluster)
library(factoextra)

b <- data.frame(read.table("correl_matrix_sensitiv_resistant_genes.txt", header=T, sep="\t", row.names=))

pdf ('SuppFigure15A elbow_plot_correlation_matrix.pdf')
fviz_nbclust(b, kmeans, method = "wss",  k.max = 30) + geom_vline(xintercept = 10, linetype = 2)+ labs(subtitle = "Elbow method")
dev.off()

pdf ('SuppFigure15B gap_plot_correlation_matrix.pdf')
plot <- fviz_nbclust(b, kmeans, nstart = 25,  method = "gap_stat", nboot = 20,  k.max = 30)+ geom_vline(xintercept = 10, linetype = 2) + labs(subtitle = "Gap statistic method") 
plot$layers <- plot$layers[-which(sapply(plot$layers, function(x) "GeomVline" %in% class(x$geom)))]
print(plot)
dev.off()















