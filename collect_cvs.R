cv <- list()
i <- 1
for(k in 2:12){
  for(f in 1:5){
    if(file.exists(paste0("TCGA/cv_res_fold_",f,'_k_',k,'.RData'))){
      load(paste0("TCGA/cv_res_fold_",f,'_k_',k,'.RData'))
      cv[[i]] <- test
    }else{
      cv[[i]] <- NULL
    }
    i <- i + 1
  }
}
cv_all <- do.call('rbind',cv)

library(dplyr)

cv_all$surv_perc <- 100*(cv_all$surv_loss/cv_all$loss)

cv_avg <- cv_all %>% #filter(alpha==.25) %>%
  group_by(k,lambda,alpha) %>% summarise(aloss=mean(loss),
                                         mean_beta=mean(nbeta),
                                         med_beta=median(nbeta),
                                         anmf=mean(nmf_loss),
                                         asurv=mean(surv_loss),
                                         apen=mean(pen_loss),
                                         asperc=mean(surv_perc),
                                         avg_cindex=mean(cindex)) %>%
  arrange(avg_cindex)

temp <- cv_avg %>% group_by(alpha,lambda) %>% 
  summarise(k=k[which.max(avg_cindex)],
            c=avg_cindex[which.max(avg_cindex)])%>%
  arrange(alpha)
temp %>% filter(alpha==.2)

library(ggplot2)

ggplot(cv_avg, aes(x=k,y=anmf,color=as.factor(alpha)))+
  geom_point()+
  geom_line()+
  facet_wrap(vars(lambda))+
  labs(color='alpha',y='Average Reconstruction Error')

ggplot(cv_avg, aes(x=k,y=avg_cindex,color=as.factor(alpha)))+
  geom_point()+
  geom_line()+
  facet_wrap(vars(lambda))+
  labs(color='alpha',y='Average c-index')

ggplot(cv_avg, aes(x=k,y=med_beta,color=as.factor(alpha)))+
  geom_point()+
  geom_line()+
  facet_wrap(vars(lambda))+
  labs(color='alpha',y='median # beta kept')


ggplot(cv_avg, aes(x=k,y=aloss,color=as.factor(alpha)))+
  geom_point()+
  geom_line()+
  facet_wrap(vars(lambda))+
  labs(color='alpha',y='Average loss')


ggplot(cv_avg, aes(x=k,y=asperc,color=as.factor(alpha)))+
  geom_point()+
  geom_line()+
  facet_wrap(vars(lambda))+
  labs(color='alpha',y='average survival contribution to loss (%)')

ggplot(cv_avg, aes(x=k,y=asurv,color=as.factor(alpha)))+
  geom_point()+
  geom_line()+
  facet_wrap(vars(lambda))+
  labs(color='alpha',y='average partial likelihood')


cv_avg_red <- cv_avg %>% filter(k %in% 2:10)


ggplot(cv_avg_red, aes(x=k,y=avg_cindex,color=as.factor(alpha)))+
  geom_point()+
  geom_line()+
  facet_wrap(vars(lambda))+
  labs(color='alpha')

cv_k6 <- cv_avg %>% filter(k==6)

ggplot(cv_k6, aes(x=lambda,y=avg_cindex,color=as.factor(alpha)))+
  geom_point()+
  geom_line()+
  labs(color='alpha')

cv_all_k6 <- cv_all %>% filter(k==6)

ggplot(cv_all_k6, aes(x=as.factor(lambda),y=cindex,group=interaction(lambda,alpha),fill=as.factor(alpha)))+
  geom_boxplot()+
  labs(fill='alpha',x='lambda')

### average difference
cv_diff <- cv_all %>% group_by(fold, k, lambda) %>%
  summarise(cdiff = cindex[alpha==1]-cindex[alpha==0],
            rediff = nmf_loss[alpha==1] - nmf_loss[alpha==0])

ggplot(cv_diff, aes(x=as.factor(k),y=cdiff))+
  geom_boxplot()+
  facet_wrap(vars(lambda))+
  geom_hline(yintercept=0,color='red')+
  labs(x='k',y='Difference in c-index between alpha=1 and alpha=0')

ggplot(cv_diff, aes(x=as.factor(k),y=rediff))+
  geom_boxplot()+
  facet_wrap(vars(lambda))+
  geom_hline(yintercept=0,color='red')+
  labs(x='k',y='Difference in recon error between alpha=1 and alpha=0')


cv_diff_avg <- cv_diff %>% group_by(k,lambda) %>%
  summarise(acindex=mean(cdiff),arecon=mean(rediff))

ggplot(cv_diff_avg,aes(x=k,y=acindex,color=as.factor(lambda)))+
  geom_point()+
  geom_line()


ggplot(cv_diff_avg,aes(x=k,y=arecon,color=as.factor(lambda)))+
  geom_point()+
  geom_line()
