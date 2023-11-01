#note this data has already been log+1 transformed
ex <- read.csv('TCGA/ex.csv')
ex <- ex[,2:ncol(ex)]#get rid of first column of gene names
fi <- read.csv('TCGA/fi.csv')
si <- read.csv('TCGA/si.csv')

# restrict to top 500 highly expressed and variable genes
means <- apply(ex,1,mean) # average expression for each genes
ex_temp <- ex[means > quantile(means,.5),] # take top 75% of genes based on mean expression
stds <- apply(ex_temp,1,sd) # Take std of expression for each genes 
ex_filtered <- ex_temp[stds>quantile(stds,1-(500/length(stds))),] # Keep top 5000 most variable genes

fi_filtered <- fi[rownames(fi) %in% rownames(ex_filtered),]

#now filter out subjects without survival info
si_filtered <- si[!is.na(si$status) & si$survival != 0,]
samps <- paste0('X',gsub('-','.',si_filtered$sampleID))
ex_filtered <- ex_filtered[,colnames(ex_filtered) %in% samps]

write.csv(ex_filtered,file='TCGA/ex_filtered.csv',row.names = FALSE)
write.csv(fi_filtered,file='TCGA/fi_filtered.csv',row.names = FALSE)
write.csv(si_filtered,file='TCGA/si_filtered.csv',row.names = FALSE)

library(survNMF)
#generate cv folds
X <- list()
X[[1]] <- ex_filtered
folds <- generate_folds(123,X,5)
save(folds,file='TCGA/folds.RData')
