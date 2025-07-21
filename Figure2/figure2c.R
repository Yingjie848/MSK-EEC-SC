# Combine HEC6-only, ISHI-only, HEC6-ISHI, ISHI-HEC6 cells
# visualize in PCA/UMAP plot
# make cell component plot
# make scatter plot for obs. vs. exp.
# however, this is not working, pca cannot run for cells with NA

# Identify HEC6 and ISHI cells
"identify_HEC6_ISHI" <- function(genotypes){

    HEC6_variants <- c('PIK3R1:p.Y73=',
                    'KDR:p.V297I',
                    'RB1:p.R661W',
                    'MET:p.S178=',
                    'CTNNB1:p.D32V',
                    'TP53:p.R273H',
                    'MET:p.N375S',
                    'PIK3CA:p.R108H',
                    'MLH1:p.V384D',
                    'KDR:p.T293I',
                    'PTEN:p.X85_splice',
                    'EGFR:p.A289V',
                    'ESR1:p.A239D',
                    'PIK3CA:p.E545=',
                    'EGFR:p.V786M',
                    'CSF1R:p.S937G')
    ISHI_variants <- c(
                'RB1:p.R320*',
                'PIK3R1:p.L570P',
                'TP53:p.D49H',
                'JAK1:p.P733=',
                #'PIK3CA:p.T1025=',
                'HNF1A:p.A298=',
                'TP53:p.M246V',
                'CSF1R:p.N968S',
                'ARID1A:p.Q548=',
                'MTOR:p.T1903I',
                'IDH1:p.N184D',
                "PTEN:p.T319*" # Frame_Shift_Del
                )

    if(sum(HEC6_variants %in% rownames(genotypes))>0 & sum(ISHI_variants %in% rownames(genotypes))>0){
        cell_identity_HEC6_ISHI <- data.frame(cell = colnames(genotypes),
                                                HEC6_mutations = apply(genotypes[HEC6_variants,],2,function(x) sum(x==1|x==2)),
                                                ISHI_mutations = apply(genotypes[ISHI_variants,],2,function(x) sum(x==1|x==2))) %>% 
                                    dplyr::mutate(group=ifelse(HEC6_mutations>0 & ISHI_mutations==0,'HEC6',
                                                        ifelse(HEC6_mutations==0 & ISHI_mutations>0,'ISHI',
                                                        ifelse(HEC6_mutations>0 & ISHI_mutations>0, 'HEC6 & ISHI','Others')))
                                    )
    }else if(sum(HEC6_variants %in% rownames(genotypes))>0){
        cell_identity_HEC6_ISHI <- data.frame(cell = colnames(genotypes),
                                                HEC6_mutations = apply(genotypes[HEC6_variants,],2,function(x) sum(x==1|x==2)),
                                                ISHI_mutations = 0) %>% 
                                    dplyr::mutate(group=ifelse(HEC6_mutations>0 & ISHI_mutations==0,'HEC6','Others')
                                    )
    }else if(sum(ISHI_variants %in% rownames(genotypes))>0){
        cell_identity_HEC6_ISHI <- data.frame(cell = colnames(genotypes),
                                                HEC6_mutations = 0,
                                                ISHI_mutations = apply(genotypes[ISHI_variants,],2,function(x) sum(x==1|x==2))) %>% 
                                    dplyr::mutate(group=ifelse(HEC6_mutations==0 & ISHI_mutations>0,'ISHI','Others')
                                    )
    }

    cell_identity_HEC6_ISHI <- cell_identity_HEC6_ISHI %>% dplyr::filter(group != "HEC6 & ISHI")

}


# Perform dimensional reduction using PCA, then visualize by UMAP
"create_umap" <- function(genotypes, group = NULL, outDir = "outdir") {

        Do_PCA <- function(genotypes) {
                pca_res <- prcomp(genotypes)

                save(pca_res, file = paste0(outDir, "/pca.Rdata"))

                pca_res$x
        }

        dir.create(outDir)

        pca_data <- Do_PCA(t(genotypes))

        umap_res <- umap(pca_data, n_components = 2)
        umap_layout <- data.frame(umap_res$layout) %>% dplyr::mutate(PC1 = pca_data[, 1], PC2 = pca_data[, 2])
        colnames(umap_layout) <- c("UMAP_1", "UMAP_2", "PC1", "PC2")
        umap_layout$cell <- rownames(umap_layout)

        umap_layout <- umap_layout %>%
                left_join(group, by = "cell") %>%
                sample_frac(1)

        # HEC6 and ISHI samples
        p <- ggplot(umap_layout %>% dplyr::filter(sample %in% c("HEC6", "ISHI")), aes(x = UMAP_1, y = UMAP_2)) +
                geom_point(aes(color = group), position = position_jitter(h = 1, w = 1), shape = 16, size = 1.5) +
                theme_classic() +
                scale_color_manual("Cell Identity", values = c("HEC6" = "#E77745", "ISHI" = "#2A2D7C", "HEC6 & ISHI" = "gray"))

        pdf(paste0(outDir, "/umap_colored_by_group.HEC6_and_ISHI.pdf"), 6, 5)
        print(p)
        dev.off()

        # HEC6-ISHI sample
        p <- ggplot(umap_layout %>% dplyr::filter(sample %in% c("HEC6-ISHI")), aes(x = UMAP_1, y = UMAP_2)) +
                geom_point(aes(color = group), position = position_jitter(h = 1, w = 1), shape = 16, size = 1.5) +
                theme_classic() +
                scale_color_manual("Cell Identity", values = c("HEC6" = "#E77745", "ISHI" = "#2A2D7C", "HEC6 & ISHI" = "gray"))

        pdf(paste0(outDir, "/umap_colored_by_group.HEC6-ISHI.pdf"), 6, 5)
        print(p)
        dev.off()

        # ISHI-HEC6 sample
        p <- ggplot(umap_layout %>% dplyr::filter(sample %in% c("ISHI-HEC6")), aes(x = UMAP_1, y = UMAP_2)) +
                geom_point(aes(color = group), position = position_jitter(h = 1, w = 1), shape = 16, size = 1.5) +
                theme_classic() +
                scale_color_manual("Cell Identity", values = c("HEC6" = "#E77745", "ISHI" = "#2A2D7C", "HEC6 & ISHI" = "gray"))

        pdf(paste0(outDir, "/umap_colored_by_group.ISHI-HEC6.pdf"), 6, 5)
        print(p)
        dev.off()

}


main <- function() {
        
        outdir <- "make_cell_line_mixture_umap"
        dir.create(outdir)

        # Get mutations data for HEC6-ISHI and ISHI-HEC6
        hec6_ishi <- readRDS("HEC6-ISHI.heatmaps_for_filters.ordered_from_top2bottom.rds")$genotypes_no_filter
        ishi_hec6 <- readRDS("ISHI-HEC6.heatmaps_for_filters.ordered_from_top2bottom.rds")$genotypes_no_filter

        # Get mutations data for HEC6 and ISHI
        rdata_hec6 <- readRDS("HEC6.genotypes.rds")
        rdata_ishi <- readRDS("ISHI.genotypes.rds")

        # Get variant name and id mapping
        variant_name_map <- readRDS("HEC6-ISHI.genotypes.rds")$somatic_variant_name_map

        # remain same order for hec6 and ishi as in HEC6-ISHI
        variant_ids <- variant_name_map[match(rownames(hec6_ishi), variant_name_map$variant_name), ]$variant_id
        # chr10.89720798.GTACT.G in ISHI, but not called in HEC6 by GATK
        variant_ids_common <- variant_ids[!variant_ids %in% "chr10.89720798.GTACT.G"]

        hec6 <- t(rdata_hec6$genotypes[, variant_ids_common])
        tmp <- matrix(rep(0, ncol(hec6)), 1, ncol(hec6), byrow = TRUE)
        rownames(tmp) <- "chr10.89720798.GTACT.G"
        hec6 <- rbind(hec6, tmp)

        ishi <- t(rdata_ishi$genotypes[, c(variant_ids_common, "chr10.89720798.GTACT.G")])

        rownames(hec6) <- variant_name_map[match(rownames(hec6), variant_name_map$variant_id), ]$variant_name
        rownames(ishi) <- variant_name_map[match(rownames(ishi), variant_name_map$variant_id), ]$variant_name

        # create cell id
        colnames(hec6) <- paste0("HEC6_", colnames(hec6))
        colnames(ishi) <- paste0("ISHI_", colnames(ishi))
        colnames(hec6_ishi) <- paste0("HEC6-ISHI_", colnames(hec6_ishi))
        colnames(ishi_hec6) <- paste0("ISHI-HEC6_", colnames(ishi_hec6))

        # Get cell identity
        hec6_cell_identity <- identify_HEC6_ISHI(hec6) %>% dplyr::mutate(sample = "HEC6")
        ishi_cell_identity <- identify_HEC6_ISHI(ishi) %>% dplyr::mutate(sample = "ISHI")
        hec6_ishi_cell_identity <- identify_HEC6_ISHI(hec6_ishi) %>% dplyr::mutate(sample = "HEC6-ISHI")
        ishi_hec6_cell_identity <- identify_HEC6_ISHI(ishi_hec6) %>% dplyr::mutate(sample = "ISHI-HEC6")

        cell_identity <- rbind(hec6_cell_identity, ishi_cell_identity, hec6_ishi_cell_identity, ishi_hec6_cell_identity)

        # combine genotypes
        tidy_hec6 <- reshape2::melt(hec6)
        tidy_ishi <- reshape2::melt(ishi)
        tidy_hec6_ishi <- reshape2::melt(hec6_ishi)
        tidy_ishi_hec6 <- reshape2::melt(ishi_hec6)

        tidy_genotypes <- rbind(tidy_hec6, tidy_ishi, tidy_hec6_ishi, tidy_ishi_hec6)

        genotypes <- reshape2::acast(tidy_genotypes, Var1 ~ Var2, value.var = "value")
        genotypes[is.na(genotypes)] <- 0

        # Get cells with cell identity
        genotypes <- genotypes[, colnames(genotypes) %in% cell_identity$cell]
        
        # create umap
        create_umap(genotypes, group = cell_identity, outDir = outdir)

}