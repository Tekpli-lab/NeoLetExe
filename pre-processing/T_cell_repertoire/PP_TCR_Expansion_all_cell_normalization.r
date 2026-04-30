library(Seurat)
library(plyr)
library(dplyr)
library(mixtools)
library(DropletUtils)
library(singleCellTK)
library(tidyr)
library(scRepertoire)
library(stringi)
library(stringr)


load("metadata_git.RData")
#### load the metadata available at EGA
TCR<- readRDS("combined_metadata_repertoire_tcr_mix1.rds")
#### load the scRepertoire output from  TCR snakemake

count_entries <- function(data_frame_name, count_data_frame_name) {
 data_frame <- as.data.frame(get(data_frame_name))
 
 count_df <- data_frame %>%
  mutate(RowNum = row_number()) %>%
  pivot_longer(cols = -RowNum, names_to = "Column", values_to = "Entry") %>%
  group_by(RowNum, Entry) %>%
  summarize(n = n(), .groups = 'drop') %>%
  pivot_wider(names_from = Entry, values_from = n, values_fill = list(n = 0))
	
	count_df<- as.data.frame (count_df)
	rownames(count_df) <- rownames(data_frame)
	count_df <- count_df[,-1]
	assign(count_data_frame_name, count_df, envir = .GlobalEnv)
}

 cell_n <- as.data.frame(table(metadata$orig.ident))
 rownames(cell_n) <- cell_n$Var1
 colnames(cell_n) <- c('orig.ident', 'cell_n')

 TCRa <- as.data.frame(table(TCR$orig.ident[!is.na(TCR$CTnt)]))
 rownames(TCRa) <- TCRa$Var1
 colnames(TCRa) <- c('orig.ident', 'TCRa')

 Tcell_n <- as.data.frame(table(TCR$orig.ident[!TCR$detail.ids %in% c('NK_cyto', 'CD4_low', 'GDT', 'UNA', 'NK_rest', 'Treg_min', 'CRTAM_GDT')]))
 rownames(Tcell_n) <- Tcell_n$Var1
 colnames(Tcell_n) <- c('orig.ident', 'Tcell_n')

 ME_n <- as.data.frame(table(metadata$orig.ident[!metadata$main.ids.cor %in% c('Malignant_cells', 'Normal_epithelial')]))
 rownames(ME_n) <- ME_n$Var1
 colnames(ME_n) <- c('orig.ident', 'ME_n')

 tab <- as.data.frame.matrix(table(TCR$orig.ident, TCR$CTnt))
 tab[rownames(tab) %in% c('P24A', 'P06A', 'P11B', 'P01C', 'P02C', 'P01B', 'P02B'), ] <- NA

 mer <- merge(cell_n, TCRa, by.x='orig.ident', by.y='orig.ident', all.x=TRUE)
 mer <- merge(mer, Tcell_n, by.x='orig.ident', by.y='orig.ident', all.x=TRUE)
 mer <- merge(mer, ME_n, by.x='orig.ident', by.y='orig.ident', all.x=TRUE)
 
mer$patient_number <- substr (mer$orig.ident,3,5)
mer$patient_number <- gsub ("A", "",mer$patient_number )
mer$patient_number <- gsub ("B", "",mer$patient_number )
mer$patient_number <- gsub ("C", "",mer$patient_number )
mer$time_point <- stri_sub(mer$orig.ident,-1)
rownames (mer) <- mer$orig.ident

TCRtab <- merge (mer, tab, by.x='orig.ident', by.y=0, all.x=T)
rownames (TCRtab) <- TCRtab$orig.ident
TCRtab <- TCRtab[,-(1: dim (mer)[2])]

all(rownames(TCRtab) == rownames(mer))

 nA <- TCRtab[mer$time_point == "A" & !rownames(TCRtab) == "P04A", ]
 nB <- TCRtab[mer$time_point == "B", ]
 nC <- TCRtab[mer$time_point == "C", ]

 merA <- mer[mer$time_point == "A" & !rownames(mer) == "P04A", ]
 merB <- mer[mer$time_point == "B", ]
 merC <- mer[mer$time_point == "C", ]

 ################################# ###########################
 ### Normalization and expansion contraction assignement  ####
 ################################# ###########################
 TCRtab_norm_Tcells <- TCRtab / (mer$cell_n)
 A <- TCRtab_norm_Tcells[mer$time_point == "A" & !rownames(TCRtab_norm_Tcells) == "P04A", ]
 B <- TCRtab_norm_Tcells[mer$time_point == "B", ]
 C <- TCRtab_norm_Tcells[mer$time_point == "C", ]

 B_A <- B - A
 C_A <- C - A
 C_B <- C - B

 statusAB1 <- ifelse(is.na(B_A), NA,
        ifelse(B_A > 0 & nB >= 5, 'str_expansion',
        ifelse(nB - nA > 0 & nB >= 5, 'str_expansion',
        ifelse(B_A > 0 & nB >= 2, 'lit_expansion',
        ifelse(nB - nA > 0 & nB >= 2, 'lit_expansion',
        ifelse(B_A < 0 & nA >= 2, 'contraction',
        ifelse(nA == 1 & nB == 0, 'singletonA',
        ifelse(nB == 1 & nA == 0, 'singletonB',
        ifelse(nA == 0 & nB == 0, NA,
        ifelse(nA == 1 & nB == 1, 'singletonAB',
        ifelse(B_A < 0 & nA == 1 & nB >= 2, 'amb_contraction',
        ifelse(B_A > 0 & nA >= 2 & nB == 1, 'amb_expansion',
        ifelse(B_A == 0, NA, 'stable')))))))))))))
count_entries("statusAB1", "statusAB1_n")
rownames (statusAB1_n) <- str_sub(rownames (statusAB1_n), 1, -2)
statusAB1_n <- statusAB1_n[ , !colnames (statusAB1_n) %in% c('NA')]
statusAB1_n[rownames(statusAB1_n)%in%c('P24', 'P06', 'P11', 'P02', 'P01'), ] <- NA

statusBC1 <-ifelse (is.na (C_B), NA, 
			ifelse (C_B > 0 & nC >= 5, 'str_expansion',
			ifelse (nC-nB > 0 & nC >= 5, 'str_expansion',
			ifelse (C_B > 0 & nC >= 2, 'lit_expansion',
			ifelse (nC-nB > 0 & nC >= 2, 'lit_expansion',				
			ifelse (C_B < 0 & nB >= 2, 'contraction',
			ifelse (nB ==1 & nC==0, 'singletonB', 
			ifelse (nC ==1 & nB==0, 'singletonC', 
			ifelse (nB ==0 & nC==0, NA,
			ifelse (nB ==1 & nC==1, 'singletonBC', 
			ifelse (C_B < 0 & nB == 1 & nC >= 2, 'amb_contraction',
			ifelse (C_B > 0 & nB >= 2 & nC == 1, 'amb_expansion',
			ifelse (C_B == 0, NA, 'stable' )))))))))))))
count_entries("statusBC1", "statusBC1_n")
rownames (statusBC1_n) <- str_sub(rownames (statusBC1_n), 1, -2)
statusBC1_n <- statusBC1_n[ , !colnames (statusBC1_n) %in% c('NA')]
statusBC1_n[rownames(statusBC1_n)%in%c('P24', 'P06', 'P11', 'P02', 'P01'), ] <- NA

statusAC1 <- ifelse (is.na (C_A), NA, 
			ifelse (C_A > 0 & nC >= 5, 'str_expansion',
			ifelse (nC-nA > 0 & nC >= 5, 'str_expansion',
			ifelse (C_A > 0 & nC >= 2, 'lit_expansion',
			ifelse (nC-nA > 0 & nC >= 2, 'lit_expansion',				
			ifelse (C_A < 0 & nA >= 2, 'contraction',
			ifelse (nA ==1 & nC==0, 'singletonA', 
			ifelse (nC ==1 & nA==0, 'singletonC', 
			ifelse (nA ==0 & nC==0, NA,
			ifelse (nA ==1 & nC==1, 'singletonAC', 
			ifelse (C_A < 0 & nA == 1 & nC >= 2, 'amb_contraction',
			ifelse (C_A > 0 & nA >= 2 & nC == 1, 'amb_expansion',
			ifelse (C_A == 0, NA, 'stable' )))))))))))))
count_entries("statusAC1", "statusAC1_n")
rownames (statusAC1_n) <- str_sub(rownames (statusAC1_n), 1, -2)
statusAC1_n <- statusAC1_n[ , !colnames (statusAC1_n) %in% c('NA')]
statusAC1_n[rownames(statusAC1_n)%in%c('P24', 'P06', 'P11', 'P02', 'P01'), ] <- NA

statusABC1 <- data.frame(matrix (paste(statusAB1, statusBC1, sep = "_"), nrow = dim (statusAB1)[1], ncol =dim (statusAB1)[2]))
rownames (statusABC1) <- rownames (statusAB1)
colnames (statusABC1) <- colnames (statusAB1)
substitutions <- c(
"NA_amb_expansion"				= "NA",
"NA_contraction"				= "late_contraction",
"NA_lit_expansion"				= "new_late_expansion",				
"NA_NA"							= "NA",
"NA_singletonB"					= "NA",    						
"NA_singletonBC"				= "NA",    						
"NA_singletonC"					= "singletonC",
"NA_str_expansion"				= "new_late_expansion",				
"lit_expansion_contraction"		="lit_expansion_contraction",
"lit_expansion_lit_expansion"	="lit_expansion_lit_expansion",
"singletonA_NA"					="singletonA",
"singletonA_singletonC"			= "singletonA_singletonC",
"singletonAB_lit_expansion"   = "singleton_lit_expansion", 		
"singletonAB_singletonB"		= "singletonAB",
"singletonB_lit_expansion"		= "singleton_lit_expansion",		
"singletonB_singletonB"			= "singletonB",
"singletonB_singletonBC"		= "singletonBC",
"str_expansion_contraction"		= "str_expansion_contraction",
"str_expansion_lit_expansion"	= "str_expansion_lit_expansion",	
"amb_expansion_singletonB"		= "NA",									
"amb_expansion_singletonBC"		= "singletonBC",
"contraction_contraction"		= "contraction_contraction", 
"contraction_lit_expansion"		= "contraction_lit_expansion", 
"contraction_NA"				= "contraction_singleton", 			
"contraction_singletonB"		= "contraction_singleton",
"contraction_singletonBC"		= "contraction_singleton",
"contraction_singletonC"		= "contraction_singleton",
"lit_expansion_amb_expansion"	= "lit_expansion",					
"singletonA_lit_expansion"		= "singleton_lit_expansion",
"singletonAB_singletonBC"		= "singletonABC",
"str_expansion_str_expansion"	= "str_expansion_str_expansion",
"contraction_amb_expansion"		= "contraction_singleton",
"contraction_str_expansion"		= "contraction_str_expansion",
"str_expansion_amb_expansion"	= "NA", 							
"lit_expansion_str_expansion" 	= "lit_expansion_str_expansion",
"singletonB_str_expansion" 	= "singleton_str_expansion",
"singletonAB_str_expansion" 	= "singleton_str_expansion",
"singletonA_str_expansion"	 	= "singleton_str_expansion",
"amb_expansion_lit_expansion"	= "late_lit_expansion",
"amb_expansion_str_expansion"	= "new_late_expansion",      
"lit_expansion_NA" 				= "lit_expansion_contraction",
"singletonAB_NA" 				="singletonAB",
"singletonB_NA"					="singletonB"
)


replace_entries <- function(x) {
 if (x %in% names(substitutions)) {
  return(substitutions[x])
 } else {
  return(x) 
 }
}
statusABC1 <- apply(statusABC1, c(1, 2), replace_entries)

count_entries("statusABC1", "statusABC1_n")
rownames (statusABC1_n) <- str_sub(rownames (statusABC1_n), 1, -2)
statusABC1_n <- statusABC1_n[ , !colnames (statusABC1_n) %in% c('NA')]
statusABC1_n[rownames(statusABC1_n)%in%c('P24', 'P06', 'P11', 'P02', 'P01'), ] <- NA
print (colnames (statusABC1_n))

statusAB1_n_all <- statusAB1_n /rowSums (statusAB1_n) *100 
statusAB1_n_f <- statusAB1_n[ , colnames (statusAB1_n) %in% c('str_expansion', 'contraction', 'lit_expansion')]
statusAB1_n_EC <- statusAB1_n_f /rowSums (statusAB1_n_f) *100 

statusBC1_n_all <- statusBC1_n /rowSums (statusBC1_n) *100 
statusBC1_n_f <- statusBC1_n[ , colnames (statusBC1_n) %in% c('str_expansion', 'contraction', 'lit_expansion')]
statusBC1_n_EC <- statusBC1_n_f /rowSums (statusBC1_n_f) *100

statusAC1_n_all <- statusAC1_n /rowSums (statusAC1_n) *100 
statusAC1_n_f <- statusAC1_n[ , colnames (statusAC1_n) %in% c('str_expansion', 'contraction', 'lit_expansion')]
statusAC1_n_EC <- statusAC1_n_f /rowSums (statusAC1_n_f) *100

statusABC1_n_all <- statusABC1_n /rowSums (statusABC1_n) *100 
statusABC1_n_f <- statusABC1_n[ , !colnames (statusABC1_n) %in% c('late_contraction', 'lit_expansion', 'singletonA', 'singletonB', 'singletonC','singletonAB','singletonBC','singletonA_singletonC', 'singletonABC')]
statusABC1_n_EC <- statusABC1_n_f /rowSums (statusABC1_n_f) *100

counts <- apply(TCRtab, 1, function(row) sum(row != 0))
counts_exp <- apply(TCRtab, 1, function(row) sum(row >= 2))
counts_exp5 <- apply(TCRtab, 1, function(row) sum(row >= 5))

av_size <- apply(TCRtab, 1, function(row) mean(row[row != 0]))
av_size_exp <- apply(TCRtab, 1, function(row) { filtered_data <- row[row >= 2] 
if (length(filtered_data) == 0) { return(0) } else { return(mean(filtered_data)) } })
av_size_exp5 <- apply(TCRtab, 1, function(row) { filtered_data <- row[row >= 5] 
if (length(filtered_data) == 0) { return(0) } else { return(mean(filtered_data)) } })

med_size_exp <- apply(TCRtab, 1, function(row) { filtered_data <- row[row >= 2] 
if (length(filtered_data) == 0) { return(0) } else { return(median(filtered_data)) } })
med_size_exp5 <- apply(TCRtab, 1, function(row) { filtered_data <- row[row >= 5] 
if (length(filtered_data) == 0) { return(0) } else { return(median(filtered_data)) } })

counts_df <- cbind(counts = counts, 
				  counts_exp = counts_exp,
				  counts_exp5 = counts_exp5, 				  
				  av_size = av_size, 
				  av_size_exp = av_size_exp, 
				  av_size_exp5 = av_size_exp5, 	
				  med_size_exp = med_size_exp,
				  med_size_exp5 = med_size_exp5)
a <- colnames (counts_df) 

counts_df <- as.data.frame(t(apply(counts_df, 1, function(row) {
 if (row[1] == 0 | is.na(row[1])) {return(rep(NA, length(row)))} else {return(row) } } ))) 
colnames (counts_df) <- a

counts_df1 <- counts_df / mer$cell_n *100

counts_A <- counts_df1[mer$time_point=="A" & !rownames (counts_df1)=="P04A",]
counts_B <- counts_df1[mer$time_point=="B" ,]
counts_C <- counts_df1[mer$time_point=="C" ,]
rownames (counts_A) <- str_sub(rownames (counts_A), 1, -2)
rownames (counts_B) <- str_sub(rownames (counts_B), 1, -2)
rownames (counts_C) <- str_sub(rownames (counts_C), 1, -2)

colnames (counts_A) <- paste (colnames (counts_A), 'A', sep="_")
colnames (counts_B) <- paste (colnames (counts_B), 'B', sep="_")
colnames (counts_C) <- paste (colnames (counts_C), 'C', sep="_")

colnames (statusAB1_n) <- paste ('n', 'AB', colnames (statusAB1_n), sep="_")
colnames (statusBC1_n) <- paste ('n','BC', colnames (statusBC1_n), sep="_")
colnames (statusAC1_n) <- paste ('n','AC', colnames (statusAC1_n), sep="_")
colnames (statusABC1_n) <- paste ('n','ABC', colnames (statusABC1_n), sep="_")

colnames (statusAB1_n_all) <- paste ('all','AB', colnames (statusAB1_n_all), sep="_")
colnames (statusBC1_n_all) <- paste ('all','BC', colnames (statusBC1_n_all), sep="_")
colnames (statusAC1_n_all) <- paste ('all','AC', colnames (statusAC1_n_all), sep="_")
colnames (statusABC1_n_all) <- paste ('all','ABC', colnames (statusABC1_n_all), sep="_")

colnames (statusAB1_n_EC) <- paste ('WOs', 'AB', colnames (statusAB1_n_EC), sep="_")
colnames (statusBC1_n_EC) <- paste ('WOs', 'BC',colnames (statusBC1_n_EC), sep="_")
colnames (statusAC1_n_EC) <- paste ('WOs', 'AC',colnames (statusAC1_n_EC), sep="_")
colnames (statusABC1_n_EC) <- paste ('WOs', 'ABC',colnames (statusABC1_n_EC), sep="_")

library(gdata)

TCR_norm_cell_bind <- cbind(counts_A, counts_B, counts_C, 
		statusAB1_n, statusAB1_n_all, statusAB1_n_EC,
		statusAC1_n, statusAC1_n_all, statusAC1_n_EC, 
		statusBC1_n, statusBC1_n_all, statusBC1_n_EC, 
		statusABC1_n, statusABC1_n_all, statusABC1_n_EC)

keep (counts_A, counts_B, counts_C, 
		statusAB1_n, statusAB1_n_all, statusAB1_n_EC, statusAB1, nA,
		statusAC1_n, statusAC1_n_all, statusAC1_n_EC, statusBC1, nB,
		statusBC1_n, statusBC1_n_all, statusBC1_n_EC, statusAC1, nC,
		statusABC1_n, statusABC1_n_all, statusABC1_n_EC, statusABC1,
		TCR_norm_cell_bind, sure=T)
save.image("TCR_norm_cell_bind.RData")
