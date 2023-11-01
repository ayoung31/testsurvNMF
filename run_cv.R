args <- commandArgs(trailingOnly = TRUE)
f <- as.numeric(args[1])
k <- as.numeric(args[2])

library(survNMF)

alpha=c(0,1)
lambda=c(0,.001,.005,seq(from=.01,to=.1,by=.01))

ex <- read.csv('TCGA/ex_filtered.csv')
si <- read.csv('TCGA/si_filtered.csv')

#format as lists
X <- list()
X[[1]] <- as.matrix(ex)
y <- list()
y[[1]] <- si$survival
delta <- list()
delta[[1]] <- si$status

load('TCGA/folds.RData')

if(!file.exists(paste0('TCGA/cv_res_fold_',f,'_k_',k,'.RData'))){
test <- cv(X,y,delta,theta=1,nfold=5,alpha=alpha,lambda=lambda,seed=123,folds,f,k)
save(test,file=paste0('TCGA/cv_res_fold_',f,'_k_',k,'.RData'))
}
