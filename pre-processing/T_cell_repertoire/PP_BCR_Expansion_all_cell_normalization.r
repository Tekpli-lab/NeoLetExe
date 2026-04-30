
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
BCR<- readRDS("combined_metadata_repertoire_bcr_mix2.rds")
#### load the scRepertoire output from  BCR snakemake


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

 BCRa <- as.data.frame(table(BCR$orig.ident[!is.na(BCR$CTnt)]))
 rownames(BCRa) <- BCRa$Var1
 colnames(BCRa) <- c('orig.ident', 'BCRa')

 Bcell_n <- as.data.frame (table (BCR$orig.ident[!BCR$final.ids%in%c('UNA')]))
 rownames (Bcell_n) <- Bcell_n$Var1
 colnames (Bcell_n) <- c('orig.ident', 'Bcell_n')

 ME_n <- as.data.frame(table(metadata$orig.ident[!metadata$main.ids.cor %in% c('Malignant_cells', 'Normal_epithelial')]))
 rownames(ME_n) <- ME_n$Var1
 colnames(ME_n) <- c('orig.ident', 'ME_n')

 tab <- as.data.frame.matrix(table(BCR$orig.ident, BCR$CTnt))
 tab[rownames(tab) %in% c('P07A', 'P07B', 'P07C'), ] <- NA

 mer <- merge(cell_n, BCRa, by.x='orig.ident', by.y='orig.ident', all.x=TRUE)
 mer <- merge(mer, Bcell_n, by.x='orig.ident', by.y='orig.ident', all.x=TRUE)
 mer <- merge(mer, ME_n, by.x='orig.ident', by.y='orig.ident', all.x=TRUE)
 
mer$patient_number <- substr (mer$orig.ident,3,5)
mer$patient_number <- gsub ("A", "",mer$patient_number )
mer$patient_number <- gsub ("B", "",mer$patient_number )
mer$patient_number <- gsub ("C", "",mer$patient_number )
mer$time_point <- stri_sub(mer$orig.ident,-1)
rownames (mer) <- mer$orig.ident

BCRtab <- merge (mer, tab, by.x='orig.ident', by.y=0, all.x=T)
rownames (BCRtab) <- BCRtab$orig.ident
BCRtab <- BCRtab[,-(1: dim (mer)[2])]

all(rownames(BCRtab) == rownames(mer))

nA <- BCRtab[mer$time_point=="A" & !rownames (BCRtab)=="P04A",]
nB <- BCRtab[mer$time_point=="B" ,]
nC <- BCRtab[mer$time_point=="C" ,]

merA <- mer[mer$time_point=="A" & !rownames (mer)=="P04A",]
merB <- mer[mer$time_point=="B" ,]
merC <- mer[mer$time_point=="C" ,]

 ################################# ###########################
 ### Normalization and expansion contraction assignement  ####
 ################################# ###########################
 BCRtab_norm_Tcells <- BCRtab / (mer$cell_n)
 A <- BCRtab_norm_Tcells[mer$time_point == "A" & !rownames(BCRtab_norm_Tcells) == "P04A", ]
 B <- BCRtab_norm_Tcells[mer$time_point == "B", ]
 C <- BCRtab_norm_Tcells[mer$time_point == "C", ]

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
statusAB1_n[rownames(statusAB1_n)%in%c('NL74'), ] <- NA

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
statusBC1_n[rownames(statusBC1_n)%in%c('NL74'), ] <- NA

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
statusAC1_n[rownames(statusAC1_n)%in%c('NL74'), ] <- NA

statusAB1_n_all <- statusAB1_n /rowSums (statusAB1_n) *100 
statusAB1_n_f <- statusAB1_n[ , colnames (statusAB1_n) %in% c('str_expansion', 'contraction', 'lit_expansion')]
statusAB1_n_EC <- statusAB1_n_f /rowSums (statusAB1_n_f) *100 

statusBC1_n_all <- statusBC1_n /rowSums (statusBC1_n) *100 
statusBC1_n_f <- statusBC1_n[ , colnames (statusBC1_n) %in% c('str_expansion', 'contraction', 'lit_expansion')]
statusBC1_n_EC <- statusBC1_n_f /rowSums (statusBC1_n_f) *100

statusAC1_n_all <- statusAC1_n /rowSums (statusAC1_n) *100 
statusAC1_n_f <- statusAC1_n[ , colnames (statusAC1_n) %in% c('str_expansion', 'contraction', 'lit_expansion')]
statusAC1_n_EC <- statusAC1_n_f /rowSums (statusAC1_n_f) *100


colnames (statusAB1_n) <- paste ('n', 'AB', colnames (statusAB1_n), sep="_")
colnames (statusBC1_n) <- paste ('n','BC', colnames (statusBC1_n), sep="_")
colnames (statusAC1_n) <- paste ('n','AC', colnames (statusAC1_n), sep="_")

colnames (statusAB1_n_all) <- paste ('all','AB', colnames (statusAB1_n_all), sep="_")
colnames (statusBC1_n_all) <- paste ('all','BC', colnames (statusBC1_n_all), sep="_")
colnames (statusAC1_n_all) <- paste ('all','AC', colnames (statusAC1_n_all), sep="_")

colnames (statusAB1_n_EC) <- paste ('WOs', 'AB', colnames (statusAB1_n_EC), sep="_")
colnames (statusBC1_n_EC) <- paste ('WOs', 'BC',colnames (statusBC1_n_EC), sep="_")
colnames (statusAC1_n_EC) <- paste ('WOs', 'AC',colnames (statusAC1_n_EC), sep="_")


library(gdata)

BCR_norm_cell_bind <- cbind(
		statusAB1_n, statusAB1_n_all, statusAB1_n_EC,
		statusAC1_n, statusAC1_n_all, statusAC1_n_EC, 
		statusBC1_n, statusBC1_n_all, statusBC1_n_EC)


for (old_value in names(replacements)) {
  new_value <- replacements[old_value]
  rownames(BCR_norm_cell_bind) <- gsub(old_value, new_value, rownames(BCR_norm_cell_bind))
  rownames(statusAB1_n_EC) <- gsub(old_value, new_value, rownames(statusAB1_n_EC))
  rownames(statusAC1_n_EC) <- gsub(old_value, new_value, rownames(statusAC1_n_EC))
}

keep (BCR_norm_cell_bind,
		statusAB1_n, statusAB1_n_all, statusAB1_n_EC, statusAB1, nA,
		statusAC1_n, statusAC1_n_all, statusAC1_n_EC, statusBC1, nB,
		statusBC1_n, statusBC1_n_all, statusBC1_n_EC, statusAC1, nC, sure=T)
save.image("BCR_norm_cell_bind.RData")
