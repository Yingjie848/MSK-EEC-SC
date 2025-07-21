# take genotypes.csv file, extract x% of cells, and save into multiple files required by infSCITE
Args <- commandArgs()
genotypes_file <- Args[6]
outdir <- Args[7]
pct_cells <- as.numeric(Args[8])

if(is.na(pct_cells)){
    pct_cells=100
}

input <- read.csv(genotypes_file,check.names=F)

dir.create(outdir)

output <- t(input[,2:ncol(input)])
colnames(output) <- input[,1]
dim(output)

# extract cells randomly
if(pct_cells<100){
   ncells_tobe_extracted=round(ncol(output)*pct_cells/100)
   cell_idx <- sample(1:ncol(output),ncells_tobe_extracted,replace=FALSE)
}else{
   cell_idx <- 1:ncol(output)
}
output <- output[,cell_idx]
cat(dim(output),file=paste0(outdir,"/summary.txt"))

write.table(output,paste0(outdir,"/data.txt"),sep=" ",row.names=F,col.names = F,quote=F)
write.table(rownames(output),paste0(outdir,"/data.geneNames"),sep=" ",row.names=F,col.names = F,quote=F)
write.table(colnames(output),paste0(outdir,"/data.sample"),sep=" ",row.names=F,col.names = F,quote=F)
write.csv(output,paste0(outdir,"/genotypes.csv"))