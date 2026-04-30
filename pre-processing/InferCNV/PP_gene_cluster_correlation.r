####################################################################################################################################################
## Calculates the score for each gene cluster and publically available gene module #################################################################
## Correlates the scores and output the correlation matrix #########################################################################################
####################################################################################################################################################
library(Seurat)
library(plyr)
library(dplyr)

load ('Scaled_NeoLetExe_final_git.RData')
#### Load Seurat object obtained from PP_Seurat_object
cluster_genes_annot <- data.frame(read.table("gene_sensitive_resistant_annotation.txt", header=TRUE, sep="\t", stringsAsFactors=F))
#gene_status found in supplementary Table 7

Scaled_m <- Scaled[,Scaled@meta.data$main.ids.cor=="Malignant_cells"]
unique_clusters <- unique(cluster_genes_annot$cutree1)
for (cluster_name in unique_clusters) {
  gen_list <- list (rownames(cluster_genes_annot)[cluster_genes_annot$cutree1 == cluster_name])
  Scaled_m <- AddModuleScore(object = Scaled_m, features = gen_list, ctrl = length(gen_list), name = paste0(cluster_name, "_score"))

  score_col_name_old <- paste0(cluster_name, "_score1")
  score_col_name_new <- paste0(cluster_name, "_score")
  colnames(Scaled_m@meta.data)[colnames(Scaled_m@meta.data) == score_col_name_old] <- score_col_name_new
 
print (cluster_name)
}


## modules from PMID: 34493872
modules <- read.table("gene_module.txt",sep='\t',header=F,row.names=1,fill=T)
for (col in c(1:7)) {
  gene_list <- unlist(c(modules[col,]))
  Scaled_m <- AddModuleScore(object = Scaled_m, features = list(gene_list), ctrl = length(gene_list), name = paste0('GM',col, "_score"))
  
  score_col_name_old <- paste0('GM',col, "_score1")
  score_col_name_new <- paste0('GM',col, "_score")
  colnames(Scaled_m@meta.data)[colnames(Scaled_m@meta.data) == score_col_name_old] <- score_col_name_new
 
  print (col)
}

## modules from PMID: 37258682
gene_lists <- data.frame(read.table("MP_pancancer_tirosh.txt",sep='\t',header=T,fill=T))
for (col in colnames(gene_lists)) {
  gene_list <- na.omit(gene_lists[[col]])
  Scaled_m <- AddModuleScore(object = Scaled_m, features = list(gene_list), ctrl = length(gene_list), name = paste0(col, "_score"))
  
  score_col_name_old <- paste0(col, "_score1")
  score_col_name_new <- paste0(col, "_score")
  colnames(Scaled_m@meta.data)[colnames(Scaled_m@meta.data) == score_col_name_old] <- score_col_name_new
  
  print (col)
}

scores <- Scaled_m@meta.data[, grepl("^cluster", colnames(Scaled_m@meta.data)) | grepl("^GM", colnames(Scaled_m@meta.data)) | grepl("^MP", colnames(Scaled_m@meta.data)) ]

cor_sco <- cor (scores, method='pearson')
write.table (cor_sco, "correlation_matrices_scores.txt",  sep = "\t")


####################################################################################################################################################
## Calculates the score for each gene cluster  #####################################################################################################
## Averages the score in resistant and sensitive clones ############################################################################################
####################################################################################################################################################
library(Seurat)
library(plyr)
library(dplyr)

load ('Scaled_NeoLetExe_final_git.RData')
#### Load Seurat object obtained from PP_Seurat_object
cluster_genes_annot <- data.frame(read.table("gene_sensitive_resistant_annotation.txt", header=TRUE, sep="\t", stringsAsFactors=F))
#gene_status found in supplementary Table 7

Scaled_m <- Scaled[,Scaled@meta.data$main.ids.cor=="Malignant_cells"]

unique_clusters <- unique(cluster_genes_annot$cutree1)
for (cluster_name in unique_clusters) {
  gen_list <- list (rownames(cluster_genes_annot)[cluster_genes_annot$cutree1 == cluster_name])
  Scaled_m <- AddModuleScore(object = Scaled_m, features = gen_list, ctrl = length(gen_list), name = paste0(cluster_name, "_score"))

  score_col_name_old <- paste0(cluster_name, "_score1")
  score_col_name_new <- paste0(cluster_name, "_score")
  colnames(Scaled_m@meta.data)[colnames(Scaled_m@meta.data) == score_col_name_old] <- score_col_name_new
 
print (cluster_name)
}


clones <- unique (Scaled@meta.data$clone_id[Scaled@meta.data$clone_cat2%in% c("res", "sens") &!is.na(Scaled@meta.data$clone_cat2)])
Subset_Cells <- Scaled_m[,Scaled_m@meta.data$clone_id%in%clones]

scores <- Subset_Cells@meta.data[, grepl("^cluster", colnames(Subset_Cells@meta.data))]

average <- aggregate(get(cluster_name) ~ clone_id, data = Subset_Cells@meta.data, FUN = mean)
scores$clone_id <- Subset_Cells@meta.data$clone_id
 library (dplyr)
 df_averages <- data.frame (scores %>%
  group_by(clone_id) %>%
  summarise(across(starts_with("cluster"), mean, na.rm = TRUE))
)

write.table (scores, "maligmant_cels_metadata_with_modules.txt",  sep = "\t")
write.table (Subset_Cells@meta.data, "all_maligmant_cells_metadata_with_modules.txt",  sep = "\t")
