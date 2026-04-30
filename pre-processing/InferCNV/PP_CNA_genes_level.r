####################################################################################################################################################
## Calculates the frequency of inferred copy number at the gene levels #############################################################################
## Test using hypergeometric the significant enrichment for amplification and deletion in sensitive and resistant clonses###########################
####################################################################################################################################################

library(Seurat)
library(plyr)
library(dplyr)
library(future)
library (mixtools)
library (tidyverse)
library (stringr)

patient_numbers <- c('01', '02', '03', '05', '06', '07', '08', '09', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23', '24', '25')

load("metadata_git.RData")
#### load the metadata available at EGA

### _cnv_genes.dat are the result of the output PP_InferCNV_per_sample.r
gene_data_list <- list()
for (num in patient_numbers) {
  file_name <- paste0("P", num, "_cnv_genes.dat")  
  gene_data <- read.table(file_name, header=TRUE, sep="\t", stringsAsFactors=FALSE)
  gene_data$pat_num <- paste0("P", num)  
  gene_data_list[[paste0("gene", num)]] <- gene_data  
}

## preparing data
target1 <- do.call(rbind, gene_data_list)
target1$cell_group_name  <- gsub ("all_observations.all_observations.", "", target1$cell_group_name )
target1$clone_id <- paste (target1$pat_num, target1$cell_group_name, sep="_")
target1$GS <- paste (target1$gene, target1$state, sep="_")
target1$GS2 <- gsub ("_5", "_4", target1$GS)
target1$GS2 <- gsub ("_6", "_4", target1$GS2)
target1$GS2 <- gsub ("_1", "_2", target1$GS2)


##############################################################################
## FREQUENCY OF AMPLIFIED OR DELETED GENES IN CLONES ######
##############################################################################

sensit <- metadata$clone_id [metadata$clone_cat=="sens" & !is.na(metadata$clone_cat)]
agress <- metadata$clone_id [metadata$clone_cat=="res" & !is.na(metadata$clone_cat)]

count_patterns_and_create_columns <- function(dataframe, prefix) {
  patterns_to_count <- c('amp', 'del', 'neutr')
  
  for (pattern in patterns_to_count) {
    new_col_name <- paste0(prefix, "_", pattern)
    dataframe[, new_col_name] <- sapply(dataframe$state3, function(row) {
      str_count(row, pattern)
    })
  }
  
  return(dataframe)
}

####### resistant clones #####################
heads <- c("Cytoband","gene", "state", "GS", "GS2","state2", 'state3') 
heads1 <- c("Cytoband","gene", "state", "GS", "GS2","state2", 'state3') 
for (numb in unique(agress)){
tar <- target1[target1$clone_id%in%numb,]
tar_a <- tar[tar$clone_id%in%agress,]
tar_a$state2 <- gsub (".*_", "", tar_a$GS2)
tar_cyto <- tar_a
tar_cyto$state3 <- ifelse (tar_cyto$state2 ==2 , 'del', tar_cyto$state2)
tar_cyto$state3 <- ifelse (tar_cyto$state2 ==4 , 'amp', tar_cyto$state3)
tar_cyto$state3 <- ifelse (tar_cyto$state3 =='amp' | tar_cyto$state3 =='del', tar_cyto$state3, 'neutr')
tar_cyto <- tar_cyto[, colnames (tar_cyto)%in% heads1]

heads <- rbind (heads, tar_cyto)
}
heads <-  heads[-1,]

write.table (heads, "gene_status_agress_clones.txt",  sep = "\t")

result_df <- data.frame(heads %>% group_by(gene) %>%
  summarize(n_gene = n(),          
			state3 = paste(state3, collapse = ", "), .groups = 'drop')) 

result_agress <- count_patterns_and_create_columns(result_df, "agress")

####### sensitive clones #####################
heads <- c("Cytoband","gene", "state", "GS", "GS2","state2", 'state3') 
heads1 <- c("Cytoband","gene", "state", "GS", "GS2","state2", 'state3') 
for (numb in unique(sensit)){
tar <- target1[target1$clone_id%in%numb,]
tar_a <- tar[tar$clone_id%in%sensit,]
tar_a$state2 <- gsub (".*_", "", tar_a$GS2)
tar_cyto <- tar_a
tar_cyto$state3 <- ifelse (tar_cyto$state2 ==2 , 'del', tar_cyto$state2)
tar_cyto$state3 <- ifelse (tar_cyto$state2 ==4 , 'amp', tar_cyto$state3)
tar_cyto$state3 <- ifelse (tar_cyto$state3 =='amp' | tar_cyto$state3 =='del', tar_cyto$state3, 'neutr')
tar_cyto <- tar_cyto[, colnames (tar_cyto)%in% heads1]

heads <- rbind (heads, tar_cyto)
}
heads <-  heads[-1,]

write.table (heads, "gene_status_sensit_clones.txt",  sep = "\t")

result_df <- data.frame(heads %>% group_by(gene) %>%
  summarize(n_gene = n(),         
			state3 = paste(state3, collapse = ", "), .groups = 'drop')) 

result_sensit <- count_patterns_and_create_columns(result_df, "sensit")

####### all clones #####################
heads <- c("Cytoband","gene", "state",  "chr", "start","end", "GS", "GS2","state2", 'state3') 
heads1 <- c("Cytoband","gene", "state",  "chr", "start","end", "GS", "GS2","state2", 'state3') 
for (numb in  unique(metadata$clone_id)){
tar_a <- target1[target1$clone_id%in%numb,]
tar_a$state2 <- gsub (".*_", "", tar_a$GS2)
tar_cyto <- tar_a
tar_cyto$state3 <- ifelse (tar_cyto$state2 ==2 , 'del', tar_cyto$state2)
tar_cyto$state3 <- ifelse (tar_cyto$state2 ==4 , 'amp', tar_cyto$state3)
tar_cyto$state3 <- ifelse (tar_cyto$state3 =='amp' | tar_cyto$state3 =='del', tar_cyto$state3, 'neutr')
tar_cyto <- tar_cyto[, colnames (tar_cyto)%in% heads1]

heads <- rbind(heads, tar_cyto)
}
heads <-  heads[-1,]
write.table (heads, "gene_status_all_clones.txt",  sep = "\t")

result_df <- data.frame(heads %>% group_by(gene) %>%
  summarize(n_gene = n(),         
			state3 = paste(state3, collapse = ", "),
			chr = unique (chr), 
			start= unique (start),
			end= unique (end),
			.groups = 'drop')) 

result_all <- count_patterns_and_create_columns(result_df, "all")

################# enrichment and hypergeo ###############
total <- merge (result_all, result_agress, by.x='gene', by.y='gene', all.x=T)
total2 <- merge (total, result_sensit, by.x='gene', by.y='gene', all.x=T)

total2$enrich_agress_amp <- (total2$agress_amp*total2$n_gene.x)/(total2$all_amp*total2$n_gene.y)
total2$enrich_agress_del <- (total2$agress_del*total2$n_gene.x)/(total2$all_del*total2$n_gene.y)
total2$enrich_sensit_amp <- (total2$sensit_amp*total2$n_gene.x)/(total2$all_amp*total2$n_gene)
total2$enrich_sensit_del <- (total2$sensit_del*total2$n_gene.x)/(total2$all_del*total2$n_gene)

total2$hg_agress_amp <- 1- phyper(total2$agress_amp, total2$all_amp, total2$n_gene.x - total2$all_amp, total2$n_gene.y, lower.tail = TRUE, log.p = FALSE)
total2$hg_agress_del <- 1-phyper(total2$agress_del, total2$all_del, total2$n_gene.x  - total2$all_del, total2$n_gene.y, lower.tail = TRUE, log.p = FALSE)
total2$hg_sensit_amp <- 1-phyper(total2$sensit_amp, total2$all_amp, total2$n_gene.x  - total2$all_amp, total2$n_gene, lower.tail = TRUE, log.p = FALSE)
total2$hg_sensit_del <- 1-phyper(total2$sensit_del, total2$all_del, total2$n_gene.x  - total2$all_del, total2$n_gene, lower.tail = TRUE, log.p = FALSE)

total2$fdr_agress_amp <- p.adjust(total2$hg_agress_amp, method = 'fdr', n = length(total2$hg_agress_amp))
total2$fdr_agress_del <- p.adjust(total2$hg_agress_del, method = 'fdr', n = length(total2$hg_agress_del))
total2$fdr_sensit_amp <- p.adjust(total2$hg_sensit_amp, method = 'fdr', n = length(total2$hg_sensit_amp))
total2$fdr_sensit_del <- p.adjust(total2$hg_sensit_del, method = 'fdr', n = length(total2$hg_sensit_del))

write.table (total2, "results_gene_status_stats.txt",  sep = "\t")