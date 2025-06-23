#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

if( length(args) != 3 ){
  stop("Usage: counts metadata output")
}

# define useful variable names from input arguments
counts=args[1]
metadata=args[2]
out=args[3]

library(edgeR)
# read in counts
data <- read.table(counts,header=T,row.names=1)
# subset counts to genes with >20 total counts
data_subset <- data[rowSums(data)>20,]
# define sample groups; reorder data_subset to match metadata order
# NOTE: this will remove samples in data_subset that aren't in metadata!
# also save metadata as a plain vector after matching the names
metadata <- read.table(metadata,header=T,row.names=1)
data_subset = data_subset[,rownames(metadata)]
metadata <- metadata[,1]
# create edgeR object
genes <- DGEList(counts=data_subset, group=metadata)
# calculate TMM factors
genes <- calcNormFactors(genes)
# estimate dispersion
genes <- estimateDisp(genes)
# run differential
stats <- exactTest(genes)
# run FDR correction
QValue <- p.adjust(stats$table$PValue,method="BH")
# add to the test$table data frame
stats$table <- cbind(stats$table,QValue)
# calculate normalized expression
norm <- cpm(genes)
# write tables to files
# col.names=NA adds an extra tab to the first row so that the columns line up
outdiff = paste(out,".diff.txt",sep="")
outnorm = paste(out,".norm.txt",sep="")
write.table(stats$table,outdiff,sep="\t",quote=F,col.names=NA)
write.table(norm,outnorm,sep="\t",quote=F,col.names=NA)
