# Compare VAF for bulk data with single cell pseudo-bulk data

compare_vaf_sc_and_bulk <- function(SAMPLE_NAME,genotypes,rdata,outDir="sc_vs_bulk"){

        dir.create(outDir,recursive=TRUE)

        #genotypes=rdata$genotype.somatic.tumor
        variant_id <- rdata$somatic_variant_name_map[match(rownames(genotypes),rdata$somatic_variant_name_map$variant_name),]$variant_id

        genotypes_ <- t(genotypes)
        ad_ <- t(t(rdata$ad[,variant_id])[,colnames(genotypes)])
        dp_ <- t(t(rdata$dp[,variant_id])[,colnames(genotypes)])
        vaf_ <- ad_/dp_

        # calculate mean vaf for each variant
        mean_vaf <- data.frame(variant_id=variant_id,
                                cells=nrow(genotypes_),
                                mean_vaf_pseudo_bulk=sapply(1:ncol(ad_),function(i){
                                        sum(ad_[,i])/sum(dp_[,i])
                                }),
                                mean_vaf_across_cells=sapply(1:ncol(vaf_),function(i){
                                        mean(vaf_[,i],na.rm=TRUE)
                                }),
                                mutated_cells=sapply(1:ncol(ad_),function(i){
                                        sum(genotypes_[,i]==1 | genotypes_[,i]==2)
                                })
                                ) %>%
                    dplyr::mutate(pct_mutated_cells=mutated_cells/cells*100)

        df_vaf <- merge(mean_vaf,rdata$maf,by='variant_id') # %>% dplyr::select(-index)
        write.table(df_vaf,paste0(outDir,"/",SAMPLE_NAME,".sc_vs_bulk.txt"),quote=F,sep="\t",row.names=F)

        df_vaf <- df_vaf %>% dplyr::filter(!is.na(GBC.vaf))

        p1 <- ggplot(df_vaf,aes(x=GBC.vaf*100,y=mean_vaf_pseudo_bulk*100,color=variant_name)) + 
                                geom_point(size=3) + geom_abline() + 
                                xlab("VAF (%) in bulk data ") + ylab("VAF (%) in single cell data") + 
                                ggtitle(paste0("R = ",round(cor(df_vaf$GBC.vaf,df_vaf$mean_vaf_pseudo_bulk),3))) + theme_minimal() + coord_cartesian(xlim=c(0,100),ylim=c(0,100))
        p2 <- ggplot(df_vaf,aes(x=GBC.vaf*100,y=mean_vaf_pseudo_bulk*100,color=variant_name)) + 
                                geom_point(size=3) + geom_abline() + 
                                xlab("VAF (%) in bulk data ") + ylab("VAF (%) in single cell data") + 
                                ggtitle(paste0("R = ",round(cor(df_vaf$GBC.vaf,df_vaf$mean_vaf_pseudo_bulk),3))) + theme_minimal() + coord_cartesian(xlim=c(0,100),ylim=c(0,100)) + theme(legend.position="none")
                                
        pdf(paste0(outDir,"/",SAMPLE_NAME,".sc_vs_bulk.pdf"),8,5); print(p1);dev.off()
        pdf(paste0(outDir,"/",SAMPLE_NAME,".sc_vs_bulk.no_legend.pdf"),5,5); print(p2);dev.off()
        

}

# Compare VAF between single cell and bulk in HEC6-ISHI sample
compare_vaf_btw_sc_and_bulk(SAMPLE_NAME = "HEC6-ISHI",
                            genotypes = readRDS("../output_final/mutation_data/HEC6-ISHI.genotypes.rds"), # genotypes for final mutations and cells
                            rdata = readRDS("../output_final/mutation_data_all/HEC6-ISHI.genotypes.rds"), # rdata include ad, dp, maf
                            outDir = "../output_final/figures/figure2f")

# Compare VAF between single cell and bulk in ISHI-HEC6 sample
compare_vaf_btw_sc_and_bulk(SAMPLE_NAME = "ISHI-HEC6",
                            genotypes = readRDS("../output_final/mutation_data/ISHI-HEC6.genotypes.rds"), # genotypes for final mutations and cells
                            rdata = readRDS("../output_final/mutation_data_all/ISHI-HEC6.genotypes.rds"), # rdata include ad, dp, maf
                            outDir = "../output_final/figures/figure2f")