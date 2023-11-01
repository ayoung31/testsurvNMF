library(gt)
library(gtExtras)
# This function returns the top genes in each factor,
# which factor is mostly basal,
# which factor is mostly classical,
# genes from these two factors to use in clustering step
get_top_genes <- function(j,avgW){
  dat <- data.frame(row=1:25)
  # avgW <- factors[,4:ncol(factors)]
  # rownames(avgW) <- factors$symbol
  
  for(i in 1:ncol(avgW)){
    m <- apply(avgW[,setdiff(1:ncol(avgW),i),drop=FALSE],1,max)
    diff <- avgW[,i] - m
    #maxdiff <- apply(diff,1,max)
    names <- names(diff[order(diff,decreasing = TRUE)][1:25])
    #weights <- as.numeric(avgW[names,i])
    #top_genes[[colnames(avgW[,i,drop=FALSE])]] <- names
    dat[[paste0('factor',i)]] <- names
  }
  dat$row <- NULL
  # b <- numeric(j)
  # c <- numeric(j)
  # merged_genes_for_cluster <- list()
  # for(i in 1:j){
  #   test <- mutate_all(factors[factors$factor==i,4:ncol(factors),drop=FALSE],as.numeric)
  #   temp <- apply(test,1,function(x){
  #     n <- 1:length(x)
  #     x[i]-max(x[setdiff(n,i)])})
  #   temp2 <- factors$symbol[factors$factor==i]
  #   temp3 <- temp2[order(temp,decreasing=TRUE)]
  #   dat[paste0('factor',i)] <- temp3[1:25]
  #   c[i] <- length(temp3[temp3 %in% clas])
  #   b[i] <- length(temp3[temp3 %in% bas])
  #   
  # }
  # dat$row <- NULL
  # if(sum(c)==0){
  #   genes_for_cluster <- dat[,which.max(b)]
  # }else if(sum(b)==0){
  #   genes_for_cluster <- dat[,which.max(c)]
  # }else if((sum(c) + sum(b)) ==0){
  #   print('WARNING!!!!\n')
  #   genes_for_cluster <- NULL
  # }else{
  #   genes_for_cluster <- unique(c(dat[,which.max(c)],dat[,which.max(b)]))
  # }
  # print(genes_for_cluster)
  
  # return(list(top_genes_per_factor=dat,classical_factor=which.max(c),
  #             basal_factor=which.max(b),genes_for_cluster=genes_for_cluster))
  return(top_genes_per_factor=dat)
  
}

## This function creates a table of top genes in each factor with 
## basal and classical genes highlighted
create_table <- function(dat,top_genes,name,j){
  builder <- function(x, genes){cells_body(columns = !!rlang::sym(x), rows = !!rlang::sym(x) %in% genes)}
  
  tab1 <- dat %>% gt() %>% 
    tab_style(cell_fill(color='orange'),
              locations=lapply(colnames(dat), builder, genes = top_genes[['BasalTumor']])) %>%
    tab_style(style=list(cell_fill(color='blue'),
                         cell_text(color='white')),
              locations=lapply(colnames(dat), builder, genes = top_genes[['ClassicalTumor']])) %>%
    tab_style(cell_fill(color='pink'),
              locations=lapply(colnames(dat), builder, genes = top_genes[['iCAF']])) %>%
    tab_style(cell_fill(color='green'),
              locations=lapply(colnames(dat), builder, genes = top_genes[['myCAF']])) %>%
    tab_style(cell_fill(color='yellow'),
              locations=lapply(colnames(dat), builder, genes = top_genes[['NormalStroma']])) %>%
    tab_style(style=list(cell_fill(color='purple'),
                         cell_text(color='white')),
              locations=lapply(colnames(dat), builder, genes = top_genes[['Exocrine']])) %>%
    tab_style(style=list(cell_fill(color='red'),
                         cell_text(color='white')),
              locations=lapply(colnames(dat), builder, genes = top_genes[['Endocrine']]))
  
  return(tab1)
}

## This function widens data for heatmap preprocessing
widen <- function(h,from,fac){
  h_temp <- h %>% filter((genes_from == from & nfactor == fac)) %>%
    select(-genes_from,-nfactor) %>%
    arrange(nclus) 
  h_temp2 <- h_temp %>% 
    pivot_wider(id_cols = sampleID, 
                names_from = nclus, 
                values_from = cluster, 
                names_prefix = from,
                unused_fn = unique)
  return(h_temp2)
}


postProcess <- function(genes_from,K,top_genes){
  # start loop over total number of factors used in NMF
  for(j in K){
    
    # input weight matrix for top 200 differentially expressed genes
    factors <- read.delim(
      paste0(genes_from,'/K=',j,'/',j,'-200-Differential-gene-table.tsv'),header = TRUE,sep='\t')

    # load(paste0(genes_from,'/K=',j,'/standard_NMF.RData'))
    # symbols <- read.csv(paste0(genes_from,'/full_symbols.csv'))
    # 
    # 
    # factors <- nmf_fit@fit@W
    # rownames(factors) <- symbols$x
    
    imp_genes <- factors %>% dplyr::select(symbol,factor )
    rows <- factors$symbol
    factors <- factors[,4:ncol(factors)]
    rownames(factors) <- rows
    
    # get top genes in each factor
    dat <- get_top_genes(j,factors)
    
    # create top genes table
    create_table(dat,top_genes,genes_from,j)
    
    # save results
    save(dat,file=paste0(genes_from,'/K=',j,'/top_genes_for_clustering.RData'))
    
  } # end loop over total number of factors
}

cluster <- function(ex_sub,genes_from,clus_out,j,C){
  #clustering output folder
  
  if(!dir.exists(clus_out)){
    dir.create(clus_out)
  }
  title <- clus_out
  
  
  # perform clustering on this expression matrix and output class assignments
  clus_res <- ConsensusClusterPlus(ex_sub,maxK=C,reps=50,pItem=0.8,pFeature=1,
                                   title=title,clusterAlg="hc",distance="pearson",seed=126,plot="png")
  
  clus <- data.frame(sampleID=names(clus_res[[2]][["consensusClass"]]))
  for(i in 2:C){
    clus[,i] <- clus_res[[i]][["consensusClass"]]
  }

  colnames(clus) <- c('sampleID',paste0('K',j,'_C',2:(ncol(clus))))
  return(clus)
}


bas_or_clas <- function(small_si){
  tab <- as.matrix(table(small_si$tsub, small_si$curcol))
  print(tab)
  labs <- numeric(ncol(tab))
  for(i in 1:ncol(tab)){
    cur <- tab[,i,drop=FALSE]

    if(which.max(cur)==1){
      labs[i] <- 'Basal-like'
    }else{
      labs[i] <- 'Classical'
    }
    small_si$curcol[small_si$curcol==i] <- labs[i]
  }
  
  return(small_si$curcol)
}

row_normalize <- function(ex){
  ex <- t(ex)
  means <- colMeans(ex,na.rm = TRUE)
  sds <- apply(ex,2,function(x){sd(x,na.rm = TRUE)})
  test <- scale(ex,center=means,scale=sds)
  ex_norm <- t(test)
  ex_norm[is.nan(ex_norm)] <- 0
  return(ex_norm)
}

