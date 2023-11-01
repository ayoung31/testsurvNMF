library(survNMF)
library(ggplot2)
library(NMF)

ex <- read.csv('TCGA/ex_filtered.csv')
si <- read.csv('TCGA/si_filtered.csv')
fi <- read.csv('TCGA/fi_filtered.csv')


# The primary inputs need to be lists for survNMF
# I designed it this way for multi-study handling
X <- list()
X[[1]] <- as.matrix(ex)
y <- list()
y[[1]] <- si$survival
delta <- list()
delta[[1]] <- si$status

# Run standard NMF on TCGA
fit0 <- nmf(X[[1]],5)

# Initialize H with NMF runs
#H0 <- init_H(X,5)

#initialize random H
H0 <- list()
H0[[1]] <- matrix(rexp(length(y[[1]])*5),nrow=5,ncol=144)

#Run our method survNMF
fit <- optimize_loss(X=X,H0=H0,k=5,y=y,delta=delta,theta=1,alpha=1,lambda=0,tol=.02,maxit = 15000)

#Here is a nice way to compare the fitted H to the initialization
plot(fit$H[[1]][4,],H0[[1]][4,]) # compares factor 4, you can change this
abline(coef = c(0, 1))


# ID factors for survNMF fit
load('cmbSubtypes.RData')
names(subtypeGeneList) <- schemaList
top_genes <- subtypeGeneList[["DECODER"]]
source('helper_functions.R')
W <- fit$W
rownames(W) <- fi$x
# get top genes in each factor
tops <- get_top_genes(j,W)
# create top genes table
create_table(tops,top_genes,genes_from,j)

# can also ID factors for regular NMF fit
Wnmf <- fit$W
rownames(Wnmf) <- fi$x
tops_nmf <- get_top_genes(j,Wnmf)
create_table(tops_nmf,top_genes,genes_from,j)

# read the full sample info dataset
si_full <- read.csv('TCGA/si.csv')
# these are the samples that have adequate survival info
keep <- !is.na(si_full$status) & si_full$survival != 0


# This is the result for Richard's NMF
H_2step <- read.delim('5-200-Full-sample-Table.tsv')
H_2step <- H_2step[,5:ncol(H_2step)] 
H_2step <- H_2step[keep,] # remove samples without adequate survival info

H_joint <- fit$H[[1]]


# compare H from joint model and H from Richard's NMF
# Look at top genes table to figure out which factors correspond to eachother
# Factor 3 of Richard's NMF run is basal
comp <- data.frame(h2step=H_2step[,3],hjoint=fit$H[[1]][4,]) 
library(ggplot2)
ggplot(comp,aes(x=h2step,y=hjoint))+
  geom_point()+
  geom_abline(intercept=0)


# here is some code to run additional survival model on fitted H
fitsurv <- cv.glmnet(t(fit$H[[1]]),Surv(y[[1]],delta[[1]]),
                     family='cox',type.measure = "C",alpha=.0005)
plot(fitsurv)
fitsurv$cvm[fitsurv$lambda==fitsurv$lambda.1se]
coef(fitsurv,s='lambda.min')

dat <- data.frame(survival=y[[1]],status=delta[[1]],h=scale(fit$H[[1]][4,]))
fits <- coxph(Surv(survival,status)~h,dat)
summary(fits)





