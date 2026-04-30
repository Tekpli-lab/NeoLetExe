

make_cell_contrast_groups <- function(metadata, contrast, metadata_colum_list){
	### Add conditional statements for contrasts, that will return a column vector with cell ids. The values should be foreground/background/Not_incl. Not_incl for cells which are excluded from the contrast.


	if(contrast == "Malignant_OnTreat_vs_baseline"){
                var = "main.ids.cor"
                ct = "Malignant_cells"

                metadata[[contrast]] <- "Not_incl"
                metadata[metadata[,var] == ct, contrast] <- "foreground"
		metadata[metadata[,var] == ct & metadata[,"Time_point"] == "Baseline", contrast] <- "background" 
                return(metadata[,contrast])
        }

	if(contrast == "Malignant_4mnd_vs_baseline"){
                var = "main.ids.cor"
                ct = "Malignant_cells"

                metadata[[contrast]] <- "Not_incl"
                metadata[metadata[,var] == ct & metadata[,"Time_point"] == "4mnd", contrast] <- "foreground"
                metadata[metadata[,var] == ct & metadata[,"Time_point"] == "Baseline", contrast] <- "background"
                return(metadata[,contrast])
        }

        if(contrast == "Malignant_C_vs_AB"){
                var = "main.ids.cor"
                ct = "Malignant_cells"

                metadata[[contrast]] <- "Not_incl"
                metadata[metadata[,var] == ct, contrast] <- "background"
                metadata[metadata[,var] == ct & metadata[,"Time_point"] == "4mnd", contrast] <- "foreground"
                return(metadata[,contrast])
        }

        if(contrast == "Malignant_AB_vs_C"){
                var = "main.ids.cor"
                ct = "Malignant_cells"

                metadata[[contrast]] <- "Not_incl"
                metadata[metadata[,var] == ct, contrast] <- "foreground"
                metadata[metadata[,var] == ct & metadata[,"Time_point"] == "4mnd", contrast] <- "background"
                return(metadata[,contrast])
        }

        if(contrast == "Malignant_2mnd_vs_baseline"){
                var = "main.ids.cor"
                ct = "Malignant_cells"

                metadata[[contrast]] <- "Not_incl"
                metadata[metadata[,var] == ct & metadata[,"Time_point"] == "2mnd", contrast] <- "foreground"
                metadata[metadata[,var] == ct & metadata[,"Time_point"] == "Baseline", contrast] <- "background"
                return(metadata[,contrast])
        }

        if(contrast == "Malignant_baseline_vs_2mnd"){
                var = "main.ids.cor"
                ct = "Malignant_cells"

                metadata[[contrast]] <- "Not_incl"
                metadata[metadata[,var] == ct & metadata[,"Time_point"] == "2mnd", contrast] <- "background"
                metadata[metadata[,var] == ct & metadata[,"Time_point"] == "Baseline", contrast] <- "foreground"
                return(metadata[,contrast])
        }

        if(contrast == "Malignant_4mnd_vs_2mnd"){
                var = "main.ids.cor"
                ct = "Malignant_cells"

                metadata[[contrast]] <- "Not_incl"
                metadata[metadata[,var] == ct & metadata[,"Time_point"] == "4mnd", contrast] <- "foreground"
                metadata[metadata[,var] == ct & metadata[,"Time_point"] == "2mnd", contrast] <- "background"
                return(metadata[,contrast])
        }

        if(contrast == "Malignant_2mnd_vs_4mnd"){
                var = "main.ids.cor"
                ct = "Malignant_cells"

                metadata[[contrast]] <- "Not_incl"
                metadata[metadata[,var] == ct & metadata[,"Time_point"] == "2mnd", contrast] <- "foreground"
                metadata[metadata[,var] == ct & metadata[,"Time_point"] == "4mnd", contrast] <- "background"
                return(metadata[,contrast])
        }

        if(contrast == "Malignant_baseline_vs_OnTreat"){
                var = "main.ids.cor"
                ct = "Malignant_cells"

                metadata[[contrast]] <- "Not_incl"
                metadata[metadata[,var] == ct, contrast] <- "background"
                metadata[metadata[,var] == ct & metadata[,"Time_point"] == "Baseline", contrast] <- "foreground"
                return(metadata[,contrast])
        }


        if(contrast == "Malignant_baseline_vs_4mnd"){
                var = "main.ids.cor"
                ct = "Malignant_cells"

                metadata[[contrast]] <- "Not_incl"
                metadata[metadata[,var] == ct & metadata[,"Time_point"] == "4mnd", contrast] <- "background"
                metadata[metadata[,var] == ct & metadata[,"Time_point"] == "Baseline", contrast] <- "foreground"
                return(metadata[,contrast])
        }



	
	### Response longitudinal analysis

	#GR A vs B
        if(contrast == "response_GR_A_vs_B"){ 
                var = "response"

                metadata[[contrast]] <- "Not_incl"
                metadata[metadata[,var] == "GR" & metadata[,"main.ids.cor"] == "Malignant_cells" & metadata[,"time_point"] == "B", contrast] <- "background"
                metadata[metadata[,var] == "GR" & metadata[,"main.ids.cor"] == "Malignant_cells" & metadata[,"time_point"] == "A", contrast] <- "foreground"
                return(metadata[,contrast])
        }

        if(contrast == "response_GR_B_vs_A"){ 
                var = "response"

                metadata[[contrast]] <- "Not_incl"
                metadata[metadata[,var] == "GR" & metadata[,"main.ids.cor"] == "Malignant_cells" & metadata[,"time_point"] == "A", contrast] <- "background"
                metadata[metadata[,var] == "GR" & metadata[,"main.ids.cor"] == "Malignant_cells" & metadata[,"time_point"] == "B", contrast] <- "foreground"
                return(metadata[,contrast])
        }

	#GR B vs C
	if(contrast == "response_GR_C_vs_B"){ 
                var = "response"

                metadata[[contrast]] <- "Not_incl"
                metadata[metadata[,var] == "GR" & metadata[,"main.ids.cor"] == "Malignant_cells" & metadata[,"time_point"] == "B", contrast] <- "background"
                metadata[metadata[,var] == "GR" & metadata[,"main.ids.cor"] == "Malignant_cells" & metadata[,"time_point"] == "C", contrast] <- "foreground"
                return(metadata[,contrast])
        }

        if(contrast == "response_GR_B_vs_C"){ 
                var = "response"

                metadata[[contrast]] <- "Not_incl"
                metadata[metadata[,var] == "GR" & metadata[,"main.ids.cor"] == "Malignant_cells" & metadata[,"time_point"] == "C", contrast] <- "background"
                metadata[metadata[,var] == "GR" & metadata[,"main.ids.cor"] == "Malignant_cells" & metadata[,"time_point"] == "B", contrast] <- "foreground"
                return(metadata[,contrast])
        }

	#GR A vs C
        if(contrast == "response_GR_C_vs_A"){ 
                var = "response"

                metadata[[contrast]] <- "Not_incl"
                metadata[metadata[,var] == "GR" & metadata[,"main.ids.cor"] == "Malignant_cells" & metadata[,"time_point"] == "A", contrast] <- "background"
                metadata[metadata[,var] == "GR" & metadata[,"main.ids.cor"] == "Malignant_cells" & metadata[,"time_point"] == "C", contrast] <- "foreground"
                return(metadata[,contrast])
        }

        if(contrast == "response_GR_A_vs_C"){ 
                var = "response"

                metadata[[contrast]] <- "Not_incl"
                metadata[metadata[,var] == "GR" & metadata[,"main.ids.cor"] == "Malignant_cells" & metadata[,"time_point"] == "C", contrast] <- "background"
                metadata[metadata[,var] == "GR" & metadata[,"main.ids.cor"] == "Malignant_cells" & metadata[,"time_point"] == "A", contrast] <- "foreground"
                return(metadata[,contrast])
        }


        #BR A vs B
        if(contrast == "response_BR_A_vs_B"){ 
                var = "response"

                metadata[[contrast]] <- "Not_incl"
                metadata[metadata[,var] == "BR" & metadata[,"main.ids.cor"] == "Malignant_cells" & metadata[,"time_point"] == "B", contrast] <- "background"
                metadata[metadata[,var] == "BR" & metadata[,"main.ids.cor"] == "Malignant_cells" & metadata[,"time_point"] == "A", contrast] <- "foreground"
                return(metadata[,contrast])
        }

        if(contrast == "response_BR_B_vs_A"){ 
                var = "response"

                metadata[[contrast]] <- "Not_incl"
                metadata[metadata[,var] == "BR" & metadata[,"main.ids.cor"] == "Malignant_cells" & metadata[,"time_point"] == "A", contrast] <- "background"
                metadata[metadata[,var] == "BR" & metadata[,"main.ids.cor"] == "Malignant_cells" & metadata[,"time_point"] == "B", contrast] <- "foreground"
                return(metadata[,contrast])
        }

        #BR B vs C
        if(contrast == "response_BR_C_vs_B"){ 
                var = "response"

                metadata[[contrast]] <- "Not_incl"
                metadata[metadata[,var] == "BR" & metadata[,"main.ids.cor"] == "Malignant_cells" & metadata[,"time_point"] == "B", contrast] <- "background"
                metadata[metadata[,var] == "BR" & metadata[,"main.ids.cor"] == "Malignant_cells" & metadata[,"time_point"] == "C", contrast] <- "foreground"
                return(metadata[,contrast])
        }

        if(contrast == "response_BR_B_vs_C"){ 
                var = "response"

                metadata[[contrast]] <- "Not_incl"
                metadata[metadata[,var] == "BR" & metadata[,"main.ids.cor"] == "Malignant_cells" & metadata[,"time_point"] == "C", contrast] <- "background"
                metadata[metadata[,var] == "BR" & metadata[,"main.ids.cor"] == "Malignant_cells" & metadata[,"time_point"] == "B", contrast] <- "foreground"
                return(metadata[,contrast])
        }

        #BR A vs C
        if(contrast == "response_BR_C_vs_A"){ 
                var = "response"

                metadata[[contrast]] <- "Not_incl"
                metadata[metadata[,var] == "BR" & metadata[,"main.ids.cor"] == "Malignant_cells" & metadata[,"time_point"] == "A", contrast] <- "background"
                metadata[metadata[,var] == "BR" & metadata[,"main.ids.cor"] == "Malignant_cells" & metadata[,"time_point"] == "C", contrast] <- "foreground"
                return(metadata[,contrast])
        }

        if(contrast == "response_BR_A_vs_C"){ 
                var = "response"

                metadata[[contrast]] <- "Not_incl"
                metadata[metadata[,var] == "BR" & metadata[,"main.ids.cor"] == "Malignant_cells" & metadata[,"time_point"] == "C", contrast] <- "background"
                metadata[metadata[,var] == "BR" & metadata[,"main.ids.cor"] == "Malignant_cells" & metadata[,"time_point"] == "A", contrast] <- "foreground"
                return(metadata[,contrast])
        }


	### Clones longitudinal analysis

        if(contrast == "clones_res_A_vs_B"){
                var = "clone_cat2"

                metadata[[contrast]] <- "Not_incl"
                metadata[metadata[,var] == "res" & metadata[,"main.ids.cor"] == "Malignant_cells" & metadata[,"time_point"] == "B", contrast] <- "background"
                metadata[metadata[,var] == "res" & metadata[,"main.ids.cor"] == "Malignant_cells" & metadata[,"time_point"] == "A", contrast] <- "foreground"
                return(metadata[,contrast])
        }

        if(contrast == "clones_res_A_vs_C"){
                var = "clone_cat2"
                
		metadata[[contrast]] <- "Not_incl"
                metadata[metadata[,var] == "res" & metadata[,"main.ids.cor"] == "Malignant_cells" & metadata[,"time_point"] == "C", contrast] <- "background"
                metadata[metadata[,var] == "res" & metadata[,"main.ids.cor"] == "Malignant_cells" & metadata[,"time_point"] == "A", contrast] <- "foreground"
                return(metadata[,contrast])
        }

        if(contrast == "clones_res_B_vs_A"){
                var = "clone_cat2"

                metadata[[contrast]] <- "Not_incl"
                metadata[metadata[,var] == "res" & metadata[,"main.ids.cor"] == "Malignant_cells" & metadata[,"time_point"] == "A", contrast] <- "background"
                metadata[metadata[,var] == "res" & metadata[,"main.ids.cor"] == "Malignant_cells" & metadata[,"time_point"] == "B", contrast] <- "foreground"
                return(metadata[,contrast])
        }

        if(contrast == "clones_res_C_vs_A"){
                var = "clone_cat2"

                metadata[[contrast]] <- "Not_incl"
                metadata[metadata[,var] == "res" & metadata[,"main.ids.cor"] == "Malignant_cells" & metadata[,"time_point"] == "A", contrast] <- "background"
                metadata[metadata[,var] == "res" & metadata[,"main.ids.cor"] == "Malignant_cells" & metadata[,"time_point"] == "C", contrast] <- "foreground"
                return(metadata[,contrast])
        }

        if(contrast == "clones_sens_A_vs_B"){
                var = "clone_cat2"

                metadata[[contrast]] <- "Not_incl"
                metadata[metadata[,var] == "sens" & metadata[,"main.ids.cor"] == "Malignant_cells" & metadata[,"time_point"] == "B", contrast] <- "background"
                metadata[metadata[,var] == "sens" & metadata[,"main.ids.cor"] == "Malignant_cells" & metadata[,"time_point"] == "A", contrast] <- "foreground"
                return(metadata[,contrast])
        }

        if(contrast == "clones_sens_A_vs_C"){
                var = "clone_cat2"

                metadata[[contrast]] <- "Not_incl"
                metadata[metadata[,var] == "sens" & metadata[,"main.ids.cor"] == "Malignant_cells" & metadata[,"time_point"] == "C", contrast] <- "background"
                metadata[metadata[,var] == "sens" & metadata[,"main.ids.cor"] == "Malignant_cells" & metadata[,"time_point"] == "A", contrast] <- "foreground"
                return(metadata[,contrast])
        }

        if(contrast == "clones_sens_B_vs_A"){
                var = "clone_cat2"

                metadata[[contrast]] <- "Not_incl"
                metadata[metadata[,var] == "sens" & metadata[,"main.ids.cor"] == "Malignant_cells" & metadata[,"time_point"] == "A", contrast] <- "background"
                metadata[metadata[,var] == "sens" & metadata[,"main.ids.cor"] == "Malignant_cells" & metadata[,"time_point"] == "B", contrast] <- "foreground"
                return(metadata[,contrast])
        }

        if(contrast == "clones_sens_C_vs_A"){
                var = "clone_cat2"

                metadata[[contrast]] <- "Not_incl"
                metadata[metadata[,var] == "sens" & metadata[,"main.ids.cor"] == "Malignant_cells" & metadata[,"time_point"] == "A", contrast] <- "background"
                metadata[metadata[,var] == "sens" & metadata[,"main.ids.cor"] == "Malignant_cells" & metadata[,"time_point"] == "C", contrast] <- "foreground"
                return(metadata[,contrast])
        }

	### ARM ###
	#EXE-LET A vs B
        if(contrast == "arm_EXE-LET_A_vs_B"){ 
                var = "arm"

                metadata[[contrast]] <- "Not_incl"
                metadata[metadata[,var] == "EXE-LET" & metadata[,"improved_main_ids"] == "Malignant_cells" & metadata[,"time_point"] == "B", contrast] <- "background"
                metadata[metadata[,var] == "EXE-LET" & metadata[,"improved_main_ids"] == "Malignant_cells" & metadata[,"time_point"] == "A", contrast] <- "foreground"
                return(metadata[,contrast])
        }

        if(contrast == "arm_EXE-LET_B_vs_A"){ 
                var = "arm"

                metadata[[contrast]] <- "Not_incl"
                metadata[metadata[,var] == "EXE-LET" & metadata[,"improved_main_ids"] == "Malignant_cells" & metadata[,"time_point"] == "A", contrast] <- "background"
                metadata[metadata[,var] == "EXE-LET" & metadata[,"improved_main_ids"] == "Malignant_cells" & metadata[,"time_point"] == "B", contrast] <- "foreground"
                return(metadata[,contrast])
        }

        #EXE-LET B vs C
        if(contrast == "arm_EXE-LET_C_vs_B"){ 
                var = "arm"

                metadata[[contrast]] <- "Not_incl"
                metadata[metadata[,var] == "EXE-LET" & metadata[,"improved_main_ids"] == "Malignant_cells" & metadata[,"time_point"] == "B", contrast] <- "background"
                metadata[metadata[,var] == "EXE-LET" & metadata[,"improved_main_ids"] == "Malignant_cells" & metadata[,"time_point"] == "C", contrast] <- "foreground"
                return(metadata[,contrast])
        }

        if(contrast == "arm_EXE-LET_B_vs_C"){ 
                var = "arm"

                metadata[[contrast]] <- "Not_incl"
                metadata[metadata[,var] == "EXE-LET" & metadata[,"improved_main_ids"] == "Malignant_cells" & metadata[,"time_point"] == "C", contrast] <- "background"
                metadata[metadata[,var] == "EXE-LET" & metadata[,"improved_main_ids"] == "Malignant_cells" & metadata[,"time_point"] == "B", contrast] <- "foreground"
                return(metadata[,contrast])
        }

        #EXE-LET A vs C
        if(contrast == "arm_EXE-LET_C_vs_A"){ 
                var = "arm"

                metadata[[contrast]] <- "Not_incl"
                metadata[metadata[,var] == "EXE-LET" & metadata[,"improved_main_ids"] == "Malignant_cells" & metadata[,"time_point"] == "A", contrast] <- "background"
                metadata[metadata[,var] == "EXE-LET" & metadata[,"improved_main_ids"] == "Malignant_cells" & metadata[,"time_point"] == "C", contrast] <- "foreground"
                return(metadata[,contrast])
        }

        if(contrast == "arm_EXE-LET_A_vs_C"){ 
                var = "arm"

                metadata[[contrast]] <- "Not_incl"
                metadata[metadata[,var] == "EXE-LET" & metadata[,"improved_main_ids"] == "Malignant_cells" & metadata[,"time_point"] == "C", contrast] <- "background"
                metadata[metadata[,var] == "EXE-LET" & metadata[,"improved_main_ids"] == "Malignant_cells" & metadata[,"time_point"] == "A", contrast] <- "foreground"
                return(metadata[,contrast])
        }

	#LET-EXE A vs B
        if(contrast == "arm_LET-EXE_A_vs_B"){ 
                var = "arm"

                metadata[[contrast]] <- "Not_incl"
                metadata[metadata[,var] == "LET-EXE" & metadata[,"improved_main_ids"] == "Malignant_cells" & metadata[,"time_point"] == "B", contrast] <- "background"
                metadata[metadata[,var] == "LET-EXE" & metadata[,"improved_main_ids"] == "Malignant_cells" & metadata[,"time_point"] == "A", contrast] <- "foreground"
                return(metadata[,contrast])
        }

        if(contrast == "arm_LET-EXE_B_vs_A"){ 
                var = "arm"

                metadata[[contrast]] <- "Not_incl"
                metadata[metadata[,var] == "LET-EXE" & metadata[,"improved_main_ids"] == "Malignant_cells" & metadata[,"time_point"] == "A", contrast] <- "background"
                metadata[metadata[,var] == "LET-EXE" & metadata[,"improved_main_ids"] == "Malignant_cells" & metadata[,"time_point"] == "B", contrast] <- "foreground"
                return(metadata[,contrast])
        }

        #LET-EXE B vs C
        if(contrast == "arm_LET-EXE_C_vs_B"){ 
                var = "arm"

                metadata[[contrast]] <- "Not_incl"
                metadata[metadata[,var] == "LET-EXE" & metadata[,"improved_main_ids"] == "Malignant_cells" & metadata[,"time_point"] == "B", contrast] <- "background"
                metadata[metadata[,var] == "LET-EXE" & metadata[,"improved_main_ids"] == "Malignant_cells" & metadata[,"time_point"] == "C", contrast] <- "foreground"
                return(metadata[,contrast])
        }

        if(contrast == "arm_LET-EXE_B_vs_C"){ 
                var = "arm"

                metadata[[contrast]] <- "Not_incl"
                metadata[metadata[,var] == "LET-EXE" & metadata[,"improved_main_ids"] == "Malignant_cells" & metadata[,"time_point"] == "C", contrast] <- "background"
                metadata[metadata[,var] == "LET-EXE" & metadata[,"improved_main_ids"] == "Malignant_cells" & metadata[,"time_point"] == "B", contrast] <- "foreground"
                return(metadata[,contrast])
        }

        #LET-EXE A vs C
        if(contrast == "arm_LET-EXE_C_vs_A"){ 
                var = "arm"

                metadata[[contrast]] <- "Not_incl"
                metadata[metadata[,var] == "LET-EXE" & metadata[,"improved_main_ids"] == "Malignant_cells" & metadata[,"time_point"] == "A", contrast] <- "background"
                metadata[metadata[,var] == "LET-EXE" & metadata[,"improved_main_ids"] == "Malignant_cells" & metadata[,"time_point"] == "C", contrast] <- "foreground"
                return(metadata[,contrast])
        }

        if(contrast == "arm_LET-EXE_A_vs_C"){ 
                var = "arm"

                metadata[[contrast]] <- "Not_incl"
                metadata[metadata[,var] == "LET-EXE" & metadata[,"improved_main_ids"] == "Malignant_cells" & metadata[,"time_point"] == "C", contrast] <- "background"
                metadata[metadata[,var] == "LET-EXE" & metadata[,"improved_main_ids"] == "Malignant_cells" & metadata[,"time_point"] == "A", contrast] <- "foreground"
                return(metadata[,contrast])
        }


	### MAIN IDS VS REST ###
	for (var in metadata_column_list){
		if(grepl(var, contrast)){
			cell_type_list <- unique(metadata[,var])
			for (ct in cell_type_list){
				if(grepl(ct, contrast)){
					metadata[[contrast]] <- "background"
					metadata[metadata[,var] == ct, contrast] <- "foreground"
					return(metadata[,contrast])
				}
			}
		}
	}

}



if(!exists("let_me_source")){ ### Ignore the 'if' wrapper. With snakemake, this block should always be executed. Just wraps the executeable code, such that the functions can be sourced to another script without running the code inside the block.
	library(optparse)

        ###########################
        #### THE EXECUTED CODE ####
        ###########################

        ### COMMAND LINE OPTIONS ###
        option_list = list(make_option(c("-o", "--output"), type = "character", default = NULL, metavar = "character"),
			   make_option(c("-m", "--sc_metadata"), type = "character", default = NULL, metavar = "character"),
			   make_option(c("-g", "--contrast"), type = "character", default = NULL, help = "Name of the contrast", metavar = "character"),
			   make_option(c("-v", "--metadata_column_list"), type = "character", default = NULL, metavar = "character"))

	opt = parse_args(OptionParser(option_list = option_list))
	contrast <- opt$contrast
	output <- opt$output
	sc_metadata_path <- opt$sc_metadata
	metadata_column_list <- strsplit(opt$metadata_column_list, " ")[[1]]
	
	metadata <- read.csv(sc_metadata_path, header = T, row.names = 1)


	####

	contrast_vector <- as.data.frame(make_cell_contrast_groups(metadata, contrast, metadata_column_list))
	rownames(contrast_vector) <- rownames(metadata)
	colnames(contrast_vector) <- contrast
	write.csv(contrast_vector, file = output, row.names=T, col.names=T)
}	

