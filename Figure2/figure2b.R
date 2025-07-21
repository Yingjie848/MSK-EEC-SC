# Make heatmap for genotype at any of Wild-type (0), Muated (1), and Missing values (3)

rm(list=ls())
library(ggplot2)
library(gtools)
library(ggpubr)
library(ComplexHeatmap)

# This function make heatmap for genotype at any of 0:Wild-type,1:Mutated,3:Missing
'make_genotype_heatmap' <- function(
		genotypes, # row is variant, column is cell, value is genotype
		sample_name = NULL,
		show_row_names = FALSE,
		show_column_names = TRUE,
		cluster_rows = FALSE,
		cluster_columns = FALSE,
		width = 5, # width of the heatmap in cm
		use_raster = FALSE){

	# order cells from top variants to bottom variants
    genotypes <- genotypes[,rownames(psychTools::dfOrder(t(genotypes)))]
	
	# get genotypes and define colors
	unique_genotypes <- as.character(sort(unique(as.vector(genotypes))))
	
	genotype_name <- c("0"="Wild-type", "1"="Mutated","3"="Missing")
	legend_params <- list(title = "Genotype", 
							at = unique_genotypes, border = "black", labels = genotype_name[unique_genotypes], 
							title_gp = grid::gpar(fontsize = 8, fontface = "bold"), 
							labels_gp = grid::gpar(fontsize = 8), 
							legend_height = 4, legend_width = 6, 
							grid_height = grid::unit(4, "mm"), grid_width = grid::unit(8, "mm"))

	col_fun <- c("0"="#787877","1"="#CD0002","3"="gray")
	col_fun <- col_fun[unique_genotypes]

	# count mutated cells for each variant
	count_mutated_cells <- function(genotypes){

			mutated_cells <- data.frame(
											variant=rownames(genotypes),
											n_of_mutated_cells = apply(genotypes,1,function(x) sum(x==1 | x==2,na.rm=T)),
											n_of_wt_cells      = apply(genotypes,1,function(x) sum(x==0,na.rm=T)),
											n_of_cells_with_missing_values = apply(genotypes,1,function(x) sum(x==3,na.rm=T))
										) %>% 
								mutate(     n_of_genotyped_cells=n_of_mutated_cells+n_of_wt_cells,
											perc_of_genotyped_cells=round(n_of_genotyped_cells/ncol(genotypes)*100,1),
											perc_of_mutated_cells=round(n_of_mutated_cells/ncol(genotypes)*100,1)
								)
			
	}
	mutated_cells <- count_mutated_cells(genotypes)

	# Create heatmap annotation for mutated cells
	ra <- HeatmapAnnotation(
        mutated_cells = anno_text(paste0(mutated_cells$perc_of_mutated_cells,'%'), gp=gpar(fontsize = 8)), which='row')

	# Create heatmap
	p <- Heatmap(genotypes,name="Genotype",
				cluster_rows = cluster_rows,
				cluster_columns = cluster_columns,
				show_row_names = show_row_names,row_names_side = "left",show_column_names = show_column_names,
				col = col_fun,
				show_row_dend = F,show_column_dend = F,border=T,
				row_names_gp = gpar(fontsize = 8,fontface="italic"),
				column_names_gp = gpar(fontsize = 8),
				column_title = paste0(sample_name,"\n",ncol(genotypes)," cells"),
				column_title_gp = gpar(fontsize = 12),
				heatmap_legend_param = legend_params,
				right_annotation=ra,
				use_raster = use_raster,raster_device="CairoPNG",width=unit(width, "cm"))

	p
  
}

variants_HEC6_ISHI <- c(
            # both
			"PDGFRA:p.V824=",
            "PIK3R1:p.F392=",
            "PTEN:p.V290*", # Frame_Shift_Del
            "CTCF:p.T204Nfs*26", # Frame_Shift_Ins
            # HEC6, sorted by % mutated cells in HEC6-ISHI
            "PIK3R1:p.Y73=",
            "KDR:p.V297I",
            "CTNNB1:p.D32V",
            "MET:p.N375S",
            "MET:p.S178=",
            "PIK3CA:p.R108H",
            "TP53:p.R273H",
            "RB1:p.R661W",
            "KDR:p.T293I",
            "PTEN:p.X85_splice",
            "EGFR:p.A289V",
            "PIK3CA:p.E545=",
            "ESR1:p.A239D",
            "MLH1:p.V384D",
            "EGFR:p.V786M",
            "CSF1R:p.S937G",
            # ISHI, sorted by % mutated cells in ISHI-HEC6
            # ISHI indel
            "PTEN:p.T319*", # Frame_Shift_Del
            # "ATM:p.E1313Dfs*7", # Frame_Shift_Del, failed in HEC6-ISHI
            # "NOTCH1:p.S2486Rfs*103", # Frame_Shift_Del, failed in ISHI-HEC6
            "PIK3R1:p.L570P",
            "RB1:p.R320*",
            "TP53:p.M246V",
            "ARID1A:p.Q548=",
            "TP53:p.D49H",
            "HNF1A:p.A298=",
            "CSF1R:p.N968S",
            "JAK1:p.P733=",
            "MTOR:p.T1903I",
            "IDH1:p.N184D"

		)

main <- function(){
	
	outdir="../output_final/figures/"
	if(!dir.exists(outdir)){
		dir.create(outdir, recursive = TRUE)
	}

	# HEC6-ISHI
	hec6_ishi_genotypes <- readRDS("../output_final/mutation_data/HEC6-ISHI.heatmap_somatic_tumor.filtered_doublets.rds")

	hec6_ishi_genotypes <- hec6_ishi_genotypes[variants_HEC6_ISHI,]

	hec6_ishi_ht <- make_genotype_heatmap(
			hec6_ishi_genotypes,
			sample_name="HEC6-ISHI",
			show_row_names=T,
			show_column_names=F,
			cluster_rows=F,
			cluster_columns=F,
			width=5,
			use_raster=TRUE)

	# ISHI-HEC6
	ishi_hec6_genotypes <- readRDS("../output_final/mutation_data/ISHI-HEC6.heatmap_somatic_tumor.filtered_doublets.rds")

	ishi_hec6_genotypes <- ishi_hec6_genotypes[variants_HEC6_ISHI,]

	ishi_hec6_ht <- make_genotype_heatmap(
			ishi_hec6_genotypes,
			sample_name="ISHI-HEC6",
			show_row_names=T,
			show_column_names=F,
			cluster_rows=F,
			cluster_columns=F,
			width=5,
			use_raster=TRUE)

	# Save heatmaps
	pdf(file.path(outdir,"/figure2b_HEC6-ISHI.pdf"), width=8, height=6)
	draw(hec6_ishi_ht)
	dev.off()

	pdf(file.path(outdir,"/figure2b_ISHI-HEC6.pdf"), width=8, height=6)
	draw(ishi_hec6_ht)
	dev.off()

}

main()