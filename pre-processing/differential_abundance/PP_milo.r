######################################################################################################################################################
## Use of milo to test for dfifferential abundance between two parameters  ###########################################################################
## creates the outputs to create milo swarmplots #####################################################################################################
######################################################################################################################################################

library (Seurat)
library (miloR)
library(SingleCellExperiment)
library(scater)
library(scran)
library(dplyr)
library(patchwork)

############################################################
# Differential abundabnce oer cell type and time point ####
############################################################
#### example on T / NK cells

load("Scaled_NL_Tcell.RData")
#### Load Seurat object from PP_subclustering
dir.create("milo_time")
setwd("milo_time")

for (clinic in c('time_point')){

Scaled_A_C <- Scaled[,!Scaled@meta.data$time_point=="B"]
Scaled_A_C@project.name <- paste ('Scaled_A_C', clinic, sep="_")

Scaled_A_B <- Scaled[,!Scaled@meta.data$time_point=="C"]
Scaled_A_B@project.name <- paste ('Scaled_A_B', clinic, sep="_")

Scaled_B_C <- Scaled[,!Scaled@meta.data$time_point=="A"]
Scaled_B_C@project.name <- paste ('Scaled_B_C', clinic, sep="_")

scaled_objects <- list(Scaled_A_B,  Scaled_A_C, Scaled_B_C)
data_r <- data.frame(object="_",Nhood=NA, sig_Nhood=NA, ngenes=NA, ncells=NA, stringsAsFactors=F)

lapply(seq_along(scaled_objects), function(i) {
  obj <- as.SingleCellExperiment(scaled_objects[[i]])

  
  milo_obj <- Milo(obj)
  milo_obj <- buildGraph(milo_obj, k = 10, d = 55, reduced.dim = "PCA")
  milo_obj <- makeNhoods(milo_obj, prop = 0.2, k = 10, d = 55, refined = TRUE, reduced_dims = "PCA")
  milo_obj <- countCells(milo_obj, meta.data = as.data.frame(colData(milo_obj)), sample="orig.ident") 
  data_r <- rbind(data_r, c(scaled_objects[[i]]@project.name, dim(nhoodCounts(milo_obj)), dim(milo_obj)))
  milo_obj <- calcNhoodDistance(milo_obj, d = 55, reduced.dim = "PCA")
  milo_obj <- buildNhoodGraph(milo_obj)
  assign(paste('milo', scaled_objects[[i]]@project.name, sep="_"), milo_obj, envir = .GlobalEnv)
  print (scaled_objects[[i]]@project.name)
  print (clinic)  
})


for (ids in c( "minor.ids", "final.ids") {
    for (vari in 'time_point') {
    pdf(paste('swarm_plots_', ids, '_', vari, '.pdf', sep=''))
  for (cate in c( 'Scaled_A_B',  'Scaled_A_C', 'Scaled_B_C')) {
       milo_obj <- get(paste('milo', cate, clinic, sep="_"))
       design_obj <- data.frame(colData(milo_obj))[,c("orig.ident", "time_point", "patient_number")]
       design_obj <- distinct(design_obj)
       rownames(design_obj) <- design_obj$orig.ident
	   da_res <- testNhoods(milo_obj, design = ~ get(vari), design.df = design_obj)
       da_res <- annotateNhoods(milo_obj, da_res, coldata_col = ids)
       assign(paste('da_res', cate, vari, ids, sep="_"), da_res)
			if (any(da_res$FDR < 0.05)) {
			plot(plotDAbeeswarm(da_res, group.by = ids) + ggtitle(paste(gsub('Scaled_', '', cate), vari, sep='_')))
}	  
	   da_res_sig_1 <- da_res[da_res$FDR <0.05 & da_res$logFC <0,]
	   da_res_sig_2 <- da_res[da_res$FDR <0.05 & da_res$logFC >0,]
	   tot <- data.frame (table (data.frame(da_res[, colnames (da_res)%in%ids])))
			if (nrow(da_res_sig_1) > 0) {
			sig1 <- data.frame(table(data.frame(da_res_sig_1[, colnames(da_res_sig_1) %in% ids])))
			} else {
			sig1 <- data.frame(Var1 = character(), Freq = numeric())
}
			if (nrow(da_res_sig_2) > 0) {
			sig2 <- data.frame(table(data.frame(da_res_sig_2[, colnames(da_res_sig_2) %in% ids])))
			} else {
			sig2 <- data.frame(Var1 = character(), Freq = numeric())
}
	   combined_data <- merge(tot, merge(sig1, sig2, by="Var1", all=TRUE), by="Var1", all=TRUE)
	   combined_data[is.na(combined_data)] <- 0
	   combined_data$HG1 <- 1- phyper (combined_data$Freq.x, sum(combined_data$Freq.x)+ sum (combined_data$Freq.y), sum (combined_data$Freq) - (sum(combined_data$Freq.x)+ sum (combined_data$Freq.y)), combined_data$Freq)
	   combined_data$HG2 <- 1- phyper (combined_data$Freq.y, sum(combined_data$Freq.x)+ sum (combined_data$Freq.y), sum (combined_data$Freq) - (sum(combined_data$Freq.x)+ sum (combined_data$Freq.y)), combined_data$Freq)
	   combined_data$FDR1 <-p.adjust(combined_data$HG1, method = 'fdr', n = length( combined_data$Var1))
	   combined_data$FDR2 <- p.adjust(combined_data$HG2, method = 'fdr', n = length (combined_data$Var1))
	   colnames (combined_data) <- c("type", 'Nhood', 'Var1', 'Var2', 'HG1', 'HG2', 'FDR1', 'FDR2')
	   assign(paste('count', cate, sep="_"), combined_data)
    }
    dev.off()
		temp <- data.frame(unique (Scaled@meta.data[,colnames(Scaled@meta.data)%in%ids]))
		temp <- data.frame(temp[order(temp[,1]),])
		colnames (temp) <- 'temp' 

		combi <- merge (temp, count_Scaled_A_B, by.x="temp", by.y="type", all.x=T)
		combi1 <- merge (combi, count_Scaled_A_C, by.x="temp", by.y="type", all.x=T)
		combi <- merge (combi1, count_Scaled_B_C, by.x="temp", by.y="type", all.x=T)
		combi$Row.names <- vari

	assign(paste('count', vari, ids, sep="_"), combi) }
 
}

    write.table (apply(count_time_point_minor.ids,2,as.character), paste('counts', '_', clinic, '_', 'minor.txt', sep=""),  sep = "\t")
	write.table (apply(count_time_point_final.ids,2,as.character), paste('counts', '_', clinic, '_', 'final.txt', sep=""),  sep = "\t")
	write.table (apply(count_time_point_detail.ids,2,as.character), paste('counts', '_', clinic, '_', 'detail.txt', sep=""),  sep = "\t")
    objects <- ls()
	filtered_objects <- objects[grep("^milo_|^da_res_", objects)]
	save(list = filtered_objects, file = paste ("milo_", clinic, '.RData', sep=""))
	
}


############################################################
# Differenttial abundabnce oer cell type and response ######
############################################################
#### example on Endothelial cells

load("Scaled_NL_Endoth.RData")
#### Load Seurat object from PP_subclustering

dir.create("milo_response")
setwd("milo_response")

for (clinic in c('perc_malig_surg10')){


Scaled_C <- Scaled[,Scaled@meta.data$time_point=="C" & !is.na(Scaled@meta.data[, colnames (Scaled@meta.data)%in%clinic])]
Scaled_C@project.name <- paste ('Scaled_C', clinic, sep="_")

Scaled_B <- Scaled[,Scaled@meta.data$time_point=="B" & !is.na(Scaled@meta.data[, colnames (Scaled@meta.data)%in%clinic])]
Scaled_B@project.name <- paste ('Scaled_B', clinic, sep="_")

Scaled_A <- Scaled[,Scaled@meta.data$time_point=="A" & !is.na(Scaled@meta.data[, colnames (Scaled@meta.data)%in%clinic])]
Scaled_A@project.name <- paste ('Scaled_A', clinic, sep="_")

scaled_objects <- list(Scaled_A,  Scaled_B, Scaled_C)
data_r <- data.frame(object="_",Nhood=NA, sig_Nhood=NA, ngenes=NA, ncells=NA, stringsAsFactors=F)

lapply(seq_along(scaled_objects), function(i) {
  obj <- as.SingleCellExperiment(scaled_objects[[i]])
  milo_obj <- Milo(obj)
  milo_obj <- buildGraph(milo_obj, k = 10, d = 55, reduced.dim = "PCA")
  milo_obj <- makeNhoods(milo_obj, prop = 0.2, k = 10, d = 55, refined = TRUE, reduced_dims = "PCA")
  milo_obj <- countCells(milo_obj, meta.data = as.data.frame(colData(milo_obj)), sample="orig.ident") 
  data_r <- rbind(data_r, c(scaled_objects[[i]]@project.name, dim(nhoodCounts(milo_obj)), dim(milo_obj)))
  milo_obj <- calcNhoodDistance(milo_obj, d = 55, reduced.dim = "PCA")
  milo_obj <- buildNhoodGraph(milo_obj)
  assign(paste('milo', scaled_objects[[i]]@project.name, sep="_"), milo_obj, envir = .GlobalEnv)
  print (scaled_objects[[i]]@project.name)
  print (clinic)  
})

    if (clinic == "perc_malig_surg10") {
    varia <- c("perc_malig_surg10")
   } 
   
for (ids in c("minor.ids", "final.ids")) {

    for (vari in varia) {
    pdf(paste('swarm_plots_', ids, '_', vari, '.pdf', sep=''))
  for (cate in c( "Scaled_A",  "Scaled_B", "Scaled_C")) {
       milo_obj <- get(paste('milo', cate, clinic, sep="_"))
       design_obj <- data.frame(colData(milo_obj))[,c('perc_malig_surg10', 'orig.ident')]
       design_obj <- distinct(design_obj)
       rownames(design_obj) <- design_obj$orig.ident
	   da_res <- testNhoods(milo_obj, design = ~ get(vari), design.df = design_obj)
       da_res <- annotateNhoods(milo_obj, da_res, coldata_col = ids)
       assign(paste('da_res', cate, vari, ids, sep="_"), da_res)
if (any(da_res$FDR < 0.05)) {
       plot(plotDAbeeswarm(da_res, group.by = ids) + ggtitle(paste(gsub('Scaled_', '', cate), vari, sep='_')))
}	  
	   da_res_sig_1 <- da_res[da_res$FDR <0.05 & da_res$logFC <0,]
	   da_res_sig_2 <- da_res[da_res$FDR <0.05 & da_res$logFC >0,]
	   tot <- data.frame (table (data.frame(da_res[, colnames (da_res)%in%ids])))
if (nrow(da_res_sig_1) > 0) {
  sig1 <- data.frame(table(data.frame(da_res_sig_1[, colnames(da_res_sig_1) %in% ids])))
} else {
  sig1 <- data.frame(Var1 = character(), Freq = numeric())
}
if (nrow(da_res_sig_2) > 0) {
  sig2 <- data.frame(table(data.frame(da_res_sig_2[, colnames(da_res_sig_2) %in% ids])))
} else {
  sig2 <- data.frame(Var1 = character(), Freq = numeric())
}
	   combined_data <- merge(tot, merge(sig1, sig2, by="Var1", all=TRUE), by="Var1", all=TRUE)
	   combined_data[is.na(combined_data)] <- 0
	   combined_data$HG1 <- 1- phyper (combined_data$Freq.x, sum(combined_data$Freq.x)+ sum (combined_data$Freq.y), sum (combined_data$Freq) - (sum(combined_data$Freq.x)+ sum (combined_data$Freq.y)), combined_data$Freq)
	   combined_data$HG2 <- 1- phyper (combined_data$Freq.y, sum(combined_data$Freq.x)+ sum (combined_data$Freq.y), sum (combined_data$Freq) - (sum(combined_data$Freq.x)+ sum (combined_data$Freq.y)), combined_data$Freq)
	   combined_data$FDR1 <-p.adjust(combined_data$HG1, method = 'fdr', n = length( combined_data$Var1))
	   combined_data$FDR2 <- p.adjust(combined_data$HG2, method = 'fdr', n = length (combined_data$Var1))
	   colnames (combined_data) <- c("type", 'Nhood', 'Var1', 'Var2', 'HG1', 'HG2', 'FDR1', 'FDR2')
	   assign(paste('count', cate, sep="_"), combined_data)
    }
    dev.off()
		temp <- data.frame(unique (Scaled@meta.data[,colnames(Scaled@meta.data)%in%ids]))
		temp <- data.frame(temp[order(temp[,1]),])
		colnames (temp) <- 'temp' 

		combi <- merge (temp, count_Scaled_A, by.x="temp", by.y="type", all.x=T)
		combi1 <- merge (combi, count_Scaled_B, by.x="temp", by.y="type", all.x=T)
		combi <- merge (combi1, count_Scaled_C, by.x="temp", by.y="type", all.x=T)
		combi$Row.names <- vari

	assign(paste('count', vari, ids, sep="_"), combi)
  
  }
 
}
    objects <- ls()
	filtered_objects <- objects[grep("^milo_|^da_res_", objects)]
	save(list = filtered_objects, file = paste ("milo_", clinic, '.RData', sep=""))
	rm(list = filtered_objects)
	if (clinic == 'perc_malig_surg10') {
  	write.table (apply(rbind(count_perc_malig_surg10_minor.ids),2,as.character), paste('counts', '_', clinic, '_', 'minor.txt', sep=""),  sep = "\t")
	write.table (apply(rbind(count_perc_malig_surg10_final.ids),2,as.character), paste('counts', '_', clinic, '_', 'final.txt', sep=""),  sep = "\t")
  } 

}



