######################################################################################################################################################
## calling inferred copy number on all cells of each patient #########################################################################################
## creates the outputs to create phylogentic trees and clones ########################################################################################
######################################################################################################################################################

library (Seurat)
library (infercnv)
library (Cairo)
options(bitmapType='cairo')

load ('Scaled_EMBO_normal.RData')
#### Load reference dataset, PMID: 33950524
Scaled_normal <- Scaled
Scaled_normal@meta.data$infer_group <- paste ("norm", Scaled_normal@meta.data$RNA_snn_res.1.2, sep="_")
Scaled_normal <- Scaled_normal[,Scaled_normal@meta.data$infer_group %in% c( 'norm_1', 'norm_19', 'norm_26', 'norm_6', 'norm_28', 'norm_5', 'norm_12', 'norm_0', 'norm_22', 'norm_20', 'norm_23', 'norm_11', 'norm_3', 'norm_29')]



load ("Malignant.RData")
## load malignant cells, obtained from PP_InferCNV_all

Scaled_NL <- Subset_Cells
common <- intersect (rownames (Scaled_NL), rownames (Scaled_normal))
Scaled_NL@meta.data$infer_group <-Scaled_NL@meta.data$orig.ident
Scaled_NL <- Scaled_NL[rownames (Scaled_NL) %in% common,]
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


sample1 <-Scaled_NL[, Scaled_NL@meta.data$patient_number=="08"]
 
Merged <- merge(sample1, y = c(Scaled_normal), project = "epithelial_infercnv")

dim (Merged@meta.data)
head (Merged@meta.data)

a<- as.matrix (Merged@assays$RNA@counts)
annot <- Merged@meta.data[,c(48, 1)]  ### ajust the numbers

table (annot$infer_group)

infercnv_obj = CreateInfercnvObject(raw_counts_matrix=a, 
									annotations_file=annot, 
									delim="\t", 
									gene_order_file="/home/xaviert/genes_location_nodup.txt",
									ref_group_names=c('norm_1', 'norm_19', 'norm_26', 'norm_6', 'norm_28', 'norm_5', 'norm_12', 'norm_0', 'norm_22', 'norm_20', 'norm_23', 'norm_11', 'norm_3', 'norm_29'))

infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1,
                             out_dir='/storage/xaviert/NL_final2/inferCNV/71A', 
                             cluster_by_groups=F, 
							 plot_steps=F,
                             denoise=T,
							 window_length=101,
							 tumor_subcluster_partition_method = "random_trees",
							 max_centered_threshold=3,
							 sd_amplifier=1.3,
							 analysis_mode = "subclusters",
							 cluster_references = F,
							 HMM = T,
							 HMM_type='i6',
							 num_threads=100)