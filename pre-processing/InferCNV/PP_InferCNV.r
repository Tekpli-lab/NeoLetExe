######################################################################################################################################################
## calling inferred copy number on all epithelial cells ##############################################################################################
## clustering epitelial cells with a reference of known epithelial cells to distinguish malignant and normal epithelial cells#########################
######################################################################################################################################################

library (Seurat)
library (infercnv)
library (Cairo)
options(bitmapType='cairo')

###################### EMBO REFERENCE FOR NORMAL #########################
#### Load reference dataset, PMID: 33950524
load ('Scaled_EMBO_normal.RData')
Scaled_normal <- Scaled
Scaled_normal@meta.data$infer_group <- paste ("norm", Scaled_normal@meta.data$RNA_snn_res.1.2, sep="_")
Scaled_normal <- Scaled_normal[,Scaled_normal@meta.data$infer_group %in% c( 'norm_1', 'norm_19', 'norm_26', 'norm_6', 'norm_28', 'norm_5', 'norm_12', 'norm_0', 'norm_22', 'norm_20', 'norm_23', 'norm_11', 'norm_3', 'norm_29')]

###################### load test dataset #########################
load ("Scaled_NeoLetExe_final_git")
#### Load Seurat object from PP_Seurat_object
Scaled_NL <- Scaled
common <- intersect (rownames (Scaled_NL), rownames (Scaled_normal))

Scaled_NL@meta.data$infer_group <- paste ("NL", Scaled_NL@meta.data$RNA_snn_res.1.2, sep="_")
Scaled_NL <- Scaled_NL[rownames (Scaled_NL) %in% common,Scaled_NL@meta.data$main.ids %in% c("Epithelial_cells")]
Scaled_normal <- Scaled_normal[rownames (Scaled_normal) %in% common,]


as_matrix <- function(mat){

  tmp <- matrix(data=0L, nrow = mat@Dim[1], ncol = mat@Dim[2])
  
  row_pos <- mat@i+1
  col_pos <- findInterval(seq(mat@x)-1,mat@p[-1])+1
  val <- mat@x
    
  for (i in seq_along(val)){
      tmp[row_pos[i],col_pos[i]] <- val[i]
  }
    
  row.names(tmp) <- mat@Dimnames[[1]]
  colnames(tmp) <- mat@Dimnames[[2]]
  return(tmp)
}

Merged <- merge(Scaled_NL, y = c(Scaled_normal), project = "epithelial_infercnv")

a<- as.matrix (as_matrix(Merged@assays$RNA@counts))
annot <- Merged@meta.data[,c(42, 1)]   ### ajust the numbers 
table (annot[,1])

infercnv_obj = CreateInfercnvObject(raw_counts_matrix=a, 
									annotations_file=annot, 
									delim="\t", 
									gene_order_file="/",
									ref_group_names=c('norm_1', 'norm_19', 'norm_26', 'norm_6', 'norm_28', 'norm_5', 'norm_12', 'norm_0', 'norm_22', 'norm_20', 'norm_23', 'norm_11', 'norm_3', 'norm_29'))

infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1,
                             out_dir='/', 
                             cluster_by_groups=T, 
							 plot_steps=F,
                             denoise=T,
                             num_threads=20,
							 window_length=101,
							 max_centered_threshold=3,
							 sd_amplifier=1.3,
							 analysis_mode = "samples",
							 cluster_references = TRUE,
							 output_format = NA,
							 HMM=F)
							 

##################################
### MAKE UMAP ####################
##################################

library(Seurat)
library (umap)
library (mixtools)

target<- data.frame(read.table("infercnv.observations.txt", header=T, sep=" ", row.names=))
normal<- data.frame(read.table("infercnv.references.txt", header=T, sep=" ", row.names=))

mat <- cbind (target, normal)
mat <- as.matrix (mat)
mat.umap = umap(t(mat)) 

library (gdata)
keep (mat.umap,   sure=TRUE)
save.image("annot_umap_git.RData")

#################################################
### Annot_Umap_update Seurat ####################
#################################################
library (Seurat)
library (dbscan)
library (fpc)
load("annot_umap_git.RData")
cl <- sNNclust(mat.umap$layout, k = 1000, eps = 1, minPts = 675)
table (cl$cluster)
annot_umap <- data.frame (cl$cluster)
rownames (annot_umap) <- rownames (mat.umap$data)
rownames (annot_umap) <- gsub ("\\.", "-", rownames(annot_umap))

## color - types ##
annot_umap$type <- substr(rownames (annot_umap),1,1)
annot_umap$type <- gsub ("P", "blue", annot_umap$type)
annot_umap$type <- gsub ("G", "red", annot_umap$type)


## color - patients ##
patient <- substr (rownames (annot_umap),1,3)
patient <- gsub ("G[0-9]", "red", patient)
patient <- gsub ("red[^å]", "red", patient)

patient <- gsub ("P", "", patient) 
annot_umap$patient<- patient

a <- table (annot_umap$cl.cluster, annot_umap$patient)


################################################################################################
## Data without the normal reference and normal epithelial clustering with reference ###########
################################################################################################
mat2 <- mat.umap$data
mat3 <- mat2[substr (rownames (mat2), 1,1)=="P",]
annot_umap2 <- annot_umap[substr (rownames (mat2), 1,1)=="P",]
table (gsub("\\.", "-", rownames (mat3))== rownames (annot_umap2))

mat4 <- mat3[!annot_umap2$cl.cluster==6,]
annot_umap2 <- annot_umap2[!annot_umap2$cl.cluster==6,]
table (gsub("\\.", "-", rownames (mat4))== rownames (annot_umap2))

library(Seurat)
library (umap)
library (mixtools)
mat.umap2 = umap(mat4)

library (gdata)
keep ( mat.umap, mat.umap2, annot_umap, annot_umap2, sure =T)
save.image ('annot_umap_git.RData')



