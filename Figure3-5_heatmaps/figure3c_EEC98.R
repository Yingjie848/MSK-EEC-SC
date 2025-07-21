# Make heatmap for genotype at any of Wild-type (0), Muated (1), and Missing values (3)

rm(list=ls())
library(dplyr)
library(ggplot2)
library(gtools)
library(ggpubr)
library(ComplexHeatmap)

# sort genotype matrix before making heatmap
'sort_genotypes' <- function(genotypes,noSortVariants=FALSE){

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

        mutated_cells <- count_mutated_cells(genotypes) %>% arrange(desc(perc_of_mutated_cells))
        mutated_cells_reverse <- mutated_cells %>% arrange(perc_of_mutated_cells)
        if(!isTRUE(noSortVariants)){
            gt <- t(t(genotypes)[,mutated_cells_reverse$variant])
            genotypes.2 <- t(t(genotypes[,rownames(psychTools::dfOrder(t(gt)))])[,mutated_cells$variant])
        }else{
            gt <- t(t(genotypes)[,rev(rownames(genotypes))])
            genotypes.2 <- t(t(genotypes[,rownames(psychTools::dfOrder(t(gt)))]))
        }
        
        #make_genotype_heatmap(genotypes.2,sample_name="tumor",show_row_names=T,show_column_names=F,cluster_columns=F,width=5,use_raster=FALSE)
        list(variants_order=match(rownames(genotypes.2),rownames(genotypes)),
             cells_order=match(colnames(genotypes.2),colnames(genotypes))
        )
}

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
    #genotypes <- genotypes[,rownames(psychTools::dfOrder(t(genotypes)))]

    genotypes_orders <- sort_genotypes(genotypes)
    genotypes <- genotypes[genotypes_orders$variants_order,genotypes_orders$cells_order]

	
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

main <- function(){
	
	outdir="../output_final/figures/"
	if(!dir.exists(outdir)){
		dir.create(outdir, recursive = TRUE)
	}

	# load final genotypes
	genotypes <- readRDS("../output_final/mutation_data/EEC98.filtered_small_clusters.rds")$genotypes_tumor

	ht <- make_genotype_heatmap(
			genotypes,
			sample_name="EEC98",
			show_row_names=T,
			show_column_names=F,
			cluster_rows=F,
			cluster_columns=F,
			width=5,
			use_raster=TRUE)

	# Save heatmaps
	pdf(file.path(outdir,"/figure3c_EEC98.pdf"), width=8, height=6)
	draw(ht)
	dev.off()

}

main()