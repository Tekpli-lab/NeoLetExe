library(Seurat)
library(plyr)
library(dplyr)
library(future)
library (mixtools)
library (pheatmap)
library (amap)
library (gwasrapidd)
library (ggplot2)

load("Scaled_NeoLetExe_final_git.RData")
#### Load Seurat object obtained from PP_Seurat_object

### all resistant cells versus all sensitive cells
Subset_Cells <- Scaled[, Scaled@meta.data$clone_cat2 %in% c('res', 'sens')]
Idents(Subset_Cells) <- "clone_cat2"
all_markers_res_vs_sens <- FindAllMarkers(object = Subset_Cells, min.pct=0.2)
analyze_gene_sets(all_markers_res_vs_sens, "markergenes_all_sensi_vs_all_resis_clone_cat2.txt")


### resistant clone each  versus all sensitive
heads <- c("p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj", "cluster", "gene", "ID")
for (clone in unique (Scaled@meta.data$clone_id[Scaled@meta.data$clone_cat2=="res" &!is.na(Scaled@meta.data$clone_cat2)])){
    Subset_Cells <- Scaled[,Scaled@meta.data$clone_id %in% clone | Scaled@meta.data$clone_cat2 =="sens"]
    Idents(Subset_Cells) <- "clone_cat2"
    markers <- FindAllMarkers(object = Subset_Cells, min.pct=0.2)
    markers$ID <- clone
	heads <- rbind (heads, markers)
    write.table(markers, paste("resistant_markergenes", clone, ".txt", sep="_"), sep="\t")
  }
heads  <- heads [ -1, ]
write.table (heads, "markergenes_resistant_each_clone.txt",  sep = "\t")

heads_a <- heads [ heads$cluster=='res' & heads$avg_log2FC > 0.5 & as.numeric(heads$p_val_adj) < 0.05,]
agress_count <- aggregate(data=heads_a[, c(7,8)], ID~.,FUN=paste,collapse=",")
counts <- data.frame (table (heads_a$gene))
table(counts$Var1==agress_count$gene)
agress_count$Freq <- counts$Freq
agress_count <- agress_count[order(mer_ag$Freq , decreasing = T),]
agress_count$clust <- "res"
write.table (agress_count, "counts_markergenes_resistant_each_clone_vs_all_sensitive.txt",  sep = "\t")


### sensitive clone each versus all resistamt
heads <- c("p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj", "cluster", "gene", "ID")
for (clone in unique (Scaled@meta.data$clone_id[Scaled@meta.data$clone_cat2=="sens" &!is.na(Scaled@meta.data$clone_cat2)])){
    Subset_Cells <- Scaled[,Scaled@meta.data$clone_id %in% clone | Scaled@meta.data$clone_cat2 =="res"]
    Idents(Subset_Cells) <- "clone_cat2"
    markers <- FindAllMarkers(object = Subset_Cells, min.pct=0.2)
    markers$ID <- clone
	heads <- rbind (heads, markers)
    write.table(markers, paste("sensitive_markergenes", clone, ".txt", sep="_"), sep="\t") 
  }
heads  <- heads [ -1, ]
write.table (heads, "markergenes_sensitive_each_clone.txt",  sep = "\t")

heads_a <- heads [ heads$cluster=='sens' & heads$avg_log2FC > 0.5 & as.numeric(heads$p_val_adj) < 0.05,]
agress_count <- aggregate(data=heads_a[, c(7,8)], ID~.,FUN=paste,collapse=",")
counts <- data.frame (table (heads_a$gene))
table(counts$Var1==agress_count$gene)
agress_count$Freq <- counts$Freq
agress_count <- agress_count[order(mer_ag$Freq , decreasing = T),]
agress_count$clust <- "res"
write.table (agress_count, "counts_markergenes_sensitive_each_clone_vs_all resistance.txt",  sep = "\t")


#### Correlation between genes upregulated in agressive and sensitive clones 
gene_status <- data.frame(read.table("gene_sensitive_resistant_annotation.txt", header=TRUE, sep="\t", stringsAsFactors=F))
Subset_Cells <- Scaled[rownames (Scaled@assays$RNA)%in%rownames(gene_status), Scaled@meta.data$main.ids.cor=="Malignant_cells"]
a <- as.matrix(Subset_Cells@assays$RNA@data)
b<- cor (t(a), method='pearson')
mean_mat <- as.matrix (b)
write.table (mean_mat, "correlation_matrix.tx",  sep = "\t")
