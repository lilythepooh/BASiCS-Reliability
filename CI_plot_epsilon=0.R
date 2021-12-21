library(matrixStats)
library(bayestestR)
library(dplyr)
library(ggplot2)
setwd()#set data path to where you stored the posterior results.
data.path<-getwd()
n_run<-200 #set the number of run to the number of replications you ran
pdf("0_CI.pdf",paper="USr",width=17)
for (index_run in (1:n_run)){
  posterior_delta<-read.table(paste0(data.path,"/sim_delta_",index_run,".txt"))
  posterior_delta<-as.matrix(posterior_delta)
  colnames(posterior_delta)<-gsub("Gene.","",colnames(posterior_delta))
  delta1<-read.table(paste0(data.path,"/simulate_delta1.txt"))
  ci_hdi_delta_89<-rep(0,ncol(posterior_delta))
  ci_hdi_delta_low_89<-rep(0,ncol(posterior_delta))
  ci_hdi_delta_high_89<-rep(0,ncol(posterior_delta))
  ci_hdi_delta_50<-rep(0,ncol(posterior_delta))
  ci_hdi_delta_low_50<-rep(0,ncol(posterior_delta))
  ci_hdi_delta_high_50<-rep(0,ncol(posterior_delta))
  for (i in 1:ncol(posterior_delta)){
    cihdi_89<-ci(posterior_delta[,i],ci=0.89,method="HDI")
    ci_hdi_delta_89[i]=cihdi_89$CI
    ci_hdi_delta_low_89[i]=cihdi_89$CI_low
    ci_hdi_delta_high_89[i]=cihdi_89$CI_high
    cihdi_50<-ci(posterior_delta[,i],ci=0.5,method="HDI")
    ci_hdi_delta_50[i]=cihdi_50$CI
    ci_hdi_delta_low_50[i]=cihdi_50$CI_low
    ci_hdi_delta_high_50[i]=cihdi_50$CI_high
  }
  length(which(delta1<ci_hdi_delta_low_89))
  length(which(delta1>ci_hdi_delta_high_89))
  length(which(delta1>ci_hdi_delta_high_50))
  length(which(delta1<ci_hdi_delta_low_50))
  delta1$true.delta<-delta1$x
  delta1$gene.name<-colnames(posterior_delta)
  hdi_89<-data.frame(gene.name=colnames(posterior_delta),upper=ci_hdi_delta_high_89,lower=ci_hdi_delta_low_89,Median=colMedians(posterior_delta),Mean=colMeans(posterior_delta))
  hdi_50<-data.frame(gene.name=colnames(posterior_delta),upper=ci_hdi_delta_high_50,lower=ci_hdi_delta_low_50,Median=colMedians(posterior_delta),Mean=colMeans(posterior_delta))
  #cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  p <- ggplot() + # fill=name allow to automatically dedicate a color for each group
    geom_errorbar(data=hdi_89,aes(x=gene.name, ymin=lower,ymax=upper),size=1,width=0.5,color="royalblue")+
    geom_errorbar(data=hdi_50,aes(x=gene.name, ymin=lower,ymax=upper),size=1,width=0.5,color="orange")+
    geom_point(data=delta1,aes(x=gene.name,y=true.delta))+
    geom_point(data=hdi_89,aes(x=gene.name,y=Median,shape=17,colour="green4"))+
    geom_point(data=hdi_89,aes(x=gene.name,y=Mean,shape=15,colour="#009E73"))+
    scale_shape_identity()+
    scale_colour_discrete(name  ="Estimation",
                          breaks=c("green4", "#009E73"),
                          labels=c("Median", "Mean")) +
    #scale_shape_discrete(name  ="Estimation",
    #                     breaks=c(17, 15),
    #                    labels=c("Median", "Mean"))+
    labs(title = "posterior Credible Interval of delta",
         subtitle = paste0("run ",index_run),
         y="delta")+
    theme(axis.text.x = element_text(size=5)) 
  print(p)
  
  posterior_mu<-read.table(paste0(data.path,"/sim_mu_",index_run,".txt"))
  posterior_mu<-as.matrix(posterior_mu)
  colnames(posterior_mu)<-gsub("Gene.","",colnames(posterior_mu))
  mu1<-read.table(paste0(data.path,"/simulate_mu1.txt"))
  mu1<-as.matrix(mu1)
  mu1<-mu1[1:ncol(posterior_delta)]
  mu1<-data.frame(gene.name=colnames(posterior_mu),true.mu=mu1)
  ci_hdi_mu_89<-rep(0,ncol(posterior_mu))
  ci_hdi_mu_low_89<-rep(0,ncol(posterior_mu))
  ci_hdi_mu_high_89<-rep(0,ncol(posterior_mu))
  ci_hdi_mu_50<-rep(0,ncol(posterior_mu))
  ci_hdi_mu_low_50<-rep(0,ncol(posterior_mu))
  ci_hdi_mu_high_50<-rep(0,ncol(posterior_mu))
  for (i in 1:ncol(posterior_mu)){
    cihdi_89<-ci(posterior_mu[,i],method="HDI")
    ci_hdi_mu_89[i]=cihdi_89$CI
    ci_hdi_mu_low_89[i]=cihdi_89$CI_low
    ci_hdi_mu_high_89[i]=cihdi_89$CI_high
    cihdi_50<-ci(posterior_mu[,i],ci=0.5,method="HDI")
    ci_hdi_mu_50[i]=cihdi_50$CI
    ci_hdi_mu_low_50[i]=cihdi_50$CI_low
    ci_hdi_mu_high_50[i]=cihdi_50$CI_high
  }
  length(which(mu1$true.mu<ci_hdi_mu_low_89))
  length(which(mu1$true.mu>ci_hdi_mu_high_89))
  length(which(mu1$true.mu>ci_hdi_mu_high_50))
  length(which(mu1$true.mu<ci_hdi_mu_low_50))
  #mu1$true.mu<-mu1$x
  #mu1$gene.name<-colnames(posterior_mu)
  hdi_89<-data.frame(gene.name=colnames(posterior_mu),upper=ci_hdi_mu_high_89,lower=ci_hdi_mu_low_89,Median=colMedians(posterior_mu),Mean=colMeans(posterior_mu))
  hdi_50<-data.frame(gene.name=colnames(posterior_mu),upper=ci_hdi_mu_high_50,lower=ci_hdi_mu_low_50,Median=colMedians(posterior_mu),Mean=colMeans(posterior_mu))
  #cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  p <- ggplot() + # fill=name allow to automatically dedicate a color for each group
    geom_errorbar(data=hdi_89,aes(x=gene.name, ymin=lower,ymax=upper),size=1,width=0.5,color="royalblue")+
    geom_errorbar(data=hdi_50,aes(x=gene.name, ymin=lower,ymax=upper),size=1,width=0.5,color="orange")+
    geom_point(data=mu1,aes(x=gene.name,y=true.mu))+
    geom_point(data=hdi_89,aes(x=gene.name,y=Median,shape=17,colour="green4"))+
    geom_point(data=hdi_89,aes(x=gene.name,y=Mean,shape=15,colour="#009E73"))+
    scale_shape_identity()+
    scale_colour_discrete(name  ="Estimation",
                          breaks=c("green4", "#009E73"),
                          labels=c("Median", "Mean")) +
    #scale_shape_discrete(name  ="Estimation",
    #                     breaks=c(17, 15),
    #                    labels=c("Median", "Mean"))+
    labs(title = "posterior Credible Interval of mu",
         subtitle = paste0("run ",index_run),
         y="mu")+
    theme(axis.text.x = element_text(size=5)) 
  print(p)
  
  posterior_nu<-read.table(paste0(data.path,"/sim_nu_",index_run,".txt"))
  posterior_nu<-as.matrix(posterior_nu)
  colnames(posterior_nu)<-gsub("Cell.","",colnames(posterior_nu))
  colnames(posterior_nu)<-gsub("_Batch1","",colnames(posterior_nu))
  nu1<-read.table(paste0(data.path,"/simulate_nu1.txt"))
  ci_hdi_nu_89<-rep(0,ncol(posterior_nu))
  ci_hdi_nu_low_89<-rep(0,ncol(posterior_nu))
  ci_hdi_nu_high_89<-rep(0,ncol(posterior_nu))
  ci_hdi_nu_50<-rep(0,ncol(posterior_nu))
  ci_hdi_nu_low_50<-rep(0,ncol(posterior_nu))
  ci_hdi_nu_high_50<-rep(0,ncol(posterior_nu))
  for (i in 1:ncol(posterior_nu)){
    cihdi_89<-ci(posterior_nu[,i],method="HDI")
    ci_hdi_nu_89[i]=cihdi_89$CI
    ci_hdi_nu_low_89[i]=cihdi_89$CI_low
    ci_hdi_nu_high_89[i]=cihdi_89$CI_high
    cihdi_50<-ci(posterior_nu[,i],ci=0.5,method="HDI")
    ci_hdi_nu_50[i]=cihdi_50$CI
    ci_hdi_nu_low_50[i]=cihdi_50$CI_low
    ci_hdi_nu_high_50[i]=cihdi_50$CI_high
  }
  nu1$true.nu<-nu1$x
  nu1$cell.name<-colnames(posterior_nu)
  hdi_89<-data.frame(cell.name=colnames(posterior_nu),upper=ci_hdi_nu_high_89,lower=ci_hdi_nu_low_89,Median=colMedians(posterior_nu),Mean=colMeans(posterior_nu))
  hdi_50<-data.frame(cell.name=colnames(posterior_nu),upper=ci_hdi_nu_high_50,lower=ci_hdi_nu_low_50,Median=colMedians(posterior_nu),Mean=colMeans(posterior_nu))
  #cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  p <- ggplot() + # fill=name allow to automatically dedicate a color for each group
    geom_errorbar(data=hdi_89,aes(x=cell.name, ymin=lower,ymax=upper),size=1,width=0.5,color="royalblue")+
    geom_errorbar(data=hdi_50,aes(x=cell.name, ymin=lower,ymax=upper),size=1,width=0.5,color="orange")+
    geom_point(data=nu1,aes(x=cell.name,y=true.nu))+
    geom_point(data=hdi_89,aes(x=cell.name,y=Median,shape=17,colour="green4"))+
    geom_point(data=hdi_89,aes(x=cell.name,y=Mean,shape=15,colour="#009E73"))+
    scale_shape_identity()+
    scale_colour_discrete(name  ="Estimation",
                          breaks=c("green4", "#009E73"),
                          labels=c("Median", "Mean")) +
    #scale_shape_discrete(name  ="Estimation",
    #                     breaks=c(17, 15),
    #                    labels=c("Median", "Mean"))+
    labs(title = "posterior Credible Interval of nu",
         subtitle = paste0("run ",index_run),
         y="nu")+
    theme(axis.text.x = element_text(size=5)) 
  print(p)
  
  posterior_phi<-read.table(paste0(data.path,"/sim_phi_",index_run,".txt"))
  posterior_phi<-as.matrix(posterior_phi)
  colnames(posterior_phi)<-gsub("Cell.","",colnames(posterior_phi))
  colnames(posterior_phi)<-gsub("_Batch1","",colnames(posterior_phi))
  phi1<-read.table(paste0(data.path,"/simulate_phi1.txt"))
  ci_hdi_phi_89<-rep(0,ncol(posterior_phi))
  ci_hdi_phi_low_89<-rep(0,ncol(posterior_phi))
  ci_hdi_phi_high_89<-rep(0,ncol(posterior_phi))
  ci_hdi_phi_50<-rep(0,ncol(posterior_phi))
  ci_hdi_phi_low_50<-rep(0,ncol(posterior_phi))
  ci_hdi_phi_high_50<-rep(0,ncol(posterior_phi))
  for (i in 1:ncol(posterior_phi)){
    cihdi_89<-ci(posterior_phi[,i],method="HDI")
    ci_hdi_phi_89[i]=cihdi_89$CI
    ci_hdi_phi_low_89[i]=cihdi_89$CI_low
    ci_hdi_phi_high_89[i]=cihdi_89$CI_high
    cihdi_50<-ci(posterior_phi[,i],ci=0.5,method="HDI")
    ci_hdi_phi_50[i]=cihdi_50$CI
    ci_hdi_phi_low_50[i]=cihdi_50$CI_low
    ci_hdi_phi_high_50[i]=cihdi_50$CI_high
  }
  phi1$true.phi<-phi1$x
  phi1$cell.name<-colnames(posterior_phi)
  hdi_89<-data.frame(cell.name=colnames(posterior_phi),upper=ci_hdi_phi_high_89,lower=ci_hdi_phi_low_89,Median=colMedians(posterior_phi),Mean=colMeans(posterior_phi))
  hdi_50<-data.frame(cell.name=colnames(posterior_phi),upper=ci_hdi_phi_high_50,lower=ci_hdi_phi_low_50,Median=colMedians(posterior_phi),Mean=colMeans(posterior_phi))
  #cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  p <- ggplot() + # fill=name allow to automatically dedicate a color for each group
    geom_errorbar(data=hdi_89,aes(x=cell.name, ymin=lower,ymax=upper),size=1,width=0.5,color="royalblue")+
    geom_errorbar(data=hdi_50,aes(x=cell.name, ymin=lower,ymax=upper),size=1,width=0.5,color="orange")+
    geom_point(data=phi1,aes(x=cell.name,y=true.phi))+
    geom_point(data=hdi_89,aes(x=cell.name,y=Median,shape=17,colour="green4"))+
    geom_point(data=hdi_89,aes(x=cell.name,y=Mean,shape=15,colour="#009E73"))+
    scale_shape_identity()+
    scale_colour_discrete(name  ="Estimation",
                          breaks=c("green4", "#009E73"),
                          labels=c("Median", "Mean")) +
    #scale_shape_discrete(name  ="Estimation",
    #                     breaks=c(17, 15),
    #                    labels=c("Median", "Mean"))+
    labs(title = "posterior Credible Interval of phi",
         subtitle = paste0("run ",index_run),
         y="phi")+
    theme(axis.text.x = element_text(size=5)) 
  print(p)
  
  posterior_s<-read.table(paste0(data.path,"/sim_s_",index_run,".txt"))
  posterior_s<-as.matrix(posterior_s)
  colnames(posterior_s)<-gsub("Cell.","",colnames(posterior_s))
  colnames(posterior_s)<-gsub("_Batch1","",colnames(posterior_s))
  s1<-read.table(paste0(data.path,"/simulate_s1.txt"))
  ci_hdi_s_89<-rep(0,ncol(posterior_s))
  ci_hdi_s_low_89<-rep(0,ncol(posterior_s))
  ci_hdi_s_high_89<-rep(0,ncol(posterior_s))
  ci_hdi_s_50<-rep(0,ncol(posterior_s))
  ci_hdi_s_low_50<-rep(0,ncol(posterior_s))
  ci_hdi_s_high_50<-rep(0,ncol(posterior_s))
  for (i in 1:ncol(posterior_s)){
    cihdi_89<-ci(posterior_s[,i],method="HDI")
    ci_hdi_s_89[i]=cihdi_89$CI
    ci_hdi_s_low_89[i]=cihdi_89$CI_low
    ci_hdi_s_high_89[i]=cihdi_89$CI_high
    cihdi_50<-ci(posterior_s[,i],ci=0.5,method="HDI")
    ci_hdi_s_50[i]=cihdi_50$CI
    ci_hdi_s_low_50[i]=cihdi_50$CI_low
    ci_hdi_s_high_50[i]=cihdi_50$CI_high
  }
  s1$true.s<-s1$x
  s1$cell.name<-colnames(posterior_s)
  hdi_89<-data.frame(cell.name=colnames(posterior_s),upper=ci_hdi_s_high_89,lower=ci_hdi_s_low_89,Median=colMedians(posterior_s),Mean=colMeans(posterior_s))
  hdi_50<-data.frame(cell.name=colnames(posterior_s),upper=ci_hdi_s_high_50,lower=ci_hdi_s_low_50,Median=colMedians(posterior_s),Mean=colMeans(posterior_s))
  #cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  p <- ggplot() + # fill=name allow to automatically dedicate a color for each group
    geom_errorbar(data=hdi_89,aes(x=cell.name, ymin=lower,ymax=upper),size=1,width=0.5,color="royalblue")+
    geom_errorbar(data=hdi_50,aes(x=cell.name, ymin=lower,ymax=upper),size=1,width=0.5,color="orange")+
    geom_point(data=s1,aes(x=cell.name,y=true.s))+
    geom_point(data=hdi_89,aes(x=cell.name,y=Median,shape=17,colour="green4"))+
    geom_point(data=hdi_89,aes(x=cell.name,y=Mean,shape=15,colour="#009E73"))+
    scale_shape_identity()+
    scale_colour_discrete(name  ="Estimation",
                          breaks=c("green4", "#009E73"),
                          labels=c("Median", "Mean")) +
    #scale_shape_discrete(name  ="Estimation",
    #                     breaks=c(17, 15),
    #                    labels=c("Median", "Mean"))+
    labs(title = "posterior Credible Interval of s",
         subtitle = paste0("run ",index_run),
         y="s")+
    theme(axis.text.x = element_text(size=5)) 
  print(p)
  
  posterior_theta<-read.table(paste0(data.path,"/sim_theta_",index_run,".txt"))
  posterior_theta<-as.matrix(posterior_theta)
  colnames(posterior_theta)<-gsub("Cell.","",colnames(posterior_theta))
  colnames(posterior_theta)<-gsub("_Batch1","",colnames(posterior_theta))
  theta1<-read.table(paste0(data.path,"/simulate_theta1.txt"))
  ci_hdi_theta_89<-rep(0,ncol(posterior_theta))
  ci_hdi_theta_low_89<-rep(0,ncol(posterior_theta))
  ci_hdi_theta_high_89<-rep(0,ncol(posterior_theta))
  ci_hdi_theta_50<-rep(0,ncol(posterior_theta))
  ci_hdi_theta_low_50<-rep(0,ncol(posterior_theta))
  ci_hdi_theta_high_50<-rep(0,ncol(posterior_theta))
  for (i in 1:ncol(posterior_theta)){
    cihdi_89<-ci(posterior_theta[,i],method="HDI")
    ci_hdi_theta_89[i]=cihdi_89$CI
    ci_hdi_theta_low_89[i]=cihdi_89$CI_low
    ci_hdi_theta_high_89[i]=cihdi_89$CI_high
    cihdi_50<-ci(posterior_theta[,i],ci=0.5,method="HDI")
    ci_hdi_theta_50[i]=cihdi_50$CI
    ci_hdi_theta_low_50[i]=cihdi_50$CI_low
    ci_hdi_theta_high_50[i]=cihdi_50$CI_high
  }
  theta1$true.theta<-theta1$x
  theta1$cell.name<-colnames(posterior_theta)
  hdi_89<-data.frame(cell.name=colnames(posterior_theta),upper=ci_hdi_theta_high_89,lower=ci_hdi_theta_low_89,Median=colMedians(posterior_theta),Mean=colMeans(posterior_theta))
  hdi_50<-data.frame(cell.name=colnames(posterior_theta),upper=ci_hdi_theta_high_50,lower=ci_hdi_theta_low_50,Median=colMedians(posterior_theta),Mean=colMeans(posterior_theta))
  #cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  p <- ggplot() + # fill=name allow to automatically dedicate a color for each group
    geom_errorbar(data=hdi_89,aes(x=cell.name, ymin=lower,ymax=upper),size=1,width=0.5,color="royalblue")+
    geom_errorbar(data=hdi_50,aes(x=cell.name, ymin=lower,ymax=upper),size=1,width=0.5,color="orange")+
    geom_point(data=theta1,aes(x=cell.name,y=true.theta))+
    geom_point(data=hdi_89,aes(x=cell.name,y=Median,shape=17,colour="green4"))+
    geom_point(data=hdi_89,aes(x=cell.name,y=Mean,shape=15,colour="#009E73"))+
    scale_shape_identity()+
    scale_colour_discrete(name  ="Estimation",
                          breaks=c("green4", "#009E73"),
                          labels=c("Median", "Mean")) +
    #scale_thetahape_discrete(name  ="Estimation",
    #                     breaks=c(17, 15),
    #                    labels=c("Median", "Mean"))+
    labs(title = "posterior Credible Interval of theta",
         subtitle = paste0("run ",index_run),
         y="theta",x="batch")+
    theme(axis.text.x = element_text(size=5)) 
  print(p)
  
}
dev.off()




########

posterior_delta<-read.table(paste0(data.path,"/sim_delta_",index_run,".txt"))
posterior_delta<-as.matrix(posterior_delta)
colnames(posterior_delta)<-gsub("Gene.","",colnames(posterior_delta))
delta1<-read.table(paste0(data.path,"/simulate_delta1.txt"))
ci_hdi_delta_89<-rep(0,ncol(posterior_delta))
ci_hdi_delta_low_89<-rep(0,ncol(posterior_delta))
ci_hdi_delta_high_89<-rep(0,ncol(posterior_delta))
ci_hdi_delta_50<-rep(0,ncol(posterior_delta))
ci_hdi_delta_low_50<-rep(0,ncol(posterior_delta))
ci_hdi_delta_high_50<-rep(0,ncol(posterior_delta))
for (i in 1:ncol(posterior_delta)){
  cihdi_89<-ci(posterior_delta[,i],method="HDI")
  ci_hdi_delta_89[i]=cihdi_89$CI
  ci_hdi_delta_low_89[i]=cihdi_89$CI_low
  ci_hdi_delta_high_89[i]=cihdi_89$CI_high
  cihdi_50<-ci(posterior_delta[,i],ci=0.5,method="HDI")
  ci_hdi_delta_50[i]=cihdi_50$CI
  ci_hdi_delta_low_50[i]=cihdi_50$CI_low
  ci_hdi_delta_high_50[i]=cihdi_50$CI_high
}
delta1$true.delta<-delta1$x
delta1$gene.name<-colnames(posterior_delta)
hdi_89<-data.frame(gene.name=colnames(posterior_delta),upper=ci_hdi_delta_high_89,lower=ci_hdi_delta_low_89,Median=colMedians(posterior_delta),Mean=colMeans(posterior_delta))
hdi_50<-data.frame(gene.name=colnames(posterior_delta),upper=ci_hdi_delta_high_50,lower=ci_hdi_delta_low_50,Median=colMedians(posterior_delta),Mean=colMeans(posterior_delta))
#cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
p <- ggplot() + # fill=name allow to automatically dedicate a color for each group
  geom_errorbar(data=hdi_89,aes(x=gene.name, ymin=lower,ymax=upper),size=1,width=0.5,color="royalblue")+
  geom_errorbar(data=hdi_50,aes(x=gene.name, ymin=lower,ymax=upper),size=1,width=0.5,color="orange")+
  geom_point(data=delta1,aes(x=gene.name,y=true.delta))+
  geom_point(data=hdi_89,aes(x=gene.name,y=Median,shape=17,colour="green4"))+
  geom_point(data=hdi_89,aes(x=gene.name,y=Mean,shape=15,colour="#009E73"))+
  scale_shape_identity()+
  scale_colour_discrete(name  ="Estimation",
                          breaks=c("green4", "#009E73"),
                          labels=c("Median", "Mean")) +
  #scale_shape_discrete(name  ="Estimation",
  #                     breaks=c(17, 15),
   #                    labels=c("Median", "Mean"))+
  labs(title = "posterior Credible Interval of delta",
       subtitle = paste0("run ",index_run),
       y="delta")+
  theme(axis.text.x = element_text(size=5)) 
p

posterior_mu<-read.table(paste0(data.path,"/sim_mu_",index_run,".txt"))
posterior_mu<-as.matrix(posterior_mu)
colnames(posterior_mu)<-gsub("Gene.","",colnames(posterior_mu))
mu1<-read.table(paste0(data.path,"/simulate_mu1.txt"))
ci_hdi_mu_89<-rep(0,ncol(posterior_mu))
ci_hdi_mu_low_89<-rep(0,ncol(posterior_mu))
ci_hdi_mu_high_89<-rep(0,ncol(posterior_mu))
ci_hdi_mu_50<-rep(0,ncol(posterior_mu))
ci_hdi_mu_low_50<-rep(0,ncol(posterior_mu))
ci_hdi_mu_high_50<-rep(0,ncol(posterior_mu))
for (i in 1:ncol(posterior_mu)){
  cihdi_89<-ci(posterior_mu[,i],method="HDI")
  ci_hdi_mu_89[i]=cihdi_89$CI
  ci_hdi_mu_low_89[i]=cihdi_89$CI_low
  ci_hdi_mu_high_89[i]=cihdi_89$CI_high
  cihdi_50<-ci(posterior_mu[,i],ci=0.5,method="HDI")
  ci_hdi_mu_50[i]=cihdi_50$CI
  ci_hdi_mu_low_50[i]=cihdi_50$CI_low
  ci_hdi_mu_high_50[i]=cihdi_50$CI_high
}
mu1$true.mu<-mu1$x
mu1$gene.name<-colnames(posterior_mu)
hdi_89<-data.frame(gene.name=colnames(posterior_mu),upper=ci_hdi_mu_high_89,lower=ci_hdi_mu_low_89,Median=colMedians(posterior_mu),Mean=colMeans(posterior_mu))
hdi_50<-data.frame(gene.name=colnames(posterior_mu),upper=ci_hdi_mu_high_50,lower=ci_hdi_mu_low_50,Median=colMedians(posterior_mu),Mean=colMeans(posterior_mu))
#cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
p <- ggplot() + # fill=name allow to automatically dedicate a color for each group
  geom_errorbar(data=hdi_89,aes(x=gene.name, ymin=lower,ymax=upper),size=1,width=0.5,color="royalblue")+
  geom_errorbar(data=hdi_50,aes(x=gene.name, ymin=lower,ymax=upper),size=1,width=0.5,color="orange")+
  geom_point(data=mu1,aes(x=gene.name,y=true.mu))+
  geom_point(data=hdi_89,aes(x=gene.name,y=Median,shape=17,colour="green4"))+
  geom_point(data=hdi_89,aes(x=gene.name,y=Mean,shape=15,colour="#009E73"))+
  scale_shape_identity()+
  scale_colour_discrete(name  ="Estimation",
                        breaks=c("green4", "#009E73"),
                        labels=c("Median", "Mean")) +
  #scale_shape_discrete(name  ="Estimation",
  #                     breaks=c(17, 15),
  #                    labels=c("Median", "Mean"))+
  labs(title = "posterior Credible Interval of mu",
       subtitle = paste0("run ",index_run),
       y="mu")+
  theme(axis.text.x = element_text(size=5)) 
p

posterior_nu<-read.table(paste0(data.path,"/sim_nu_",index_run,".txt"))
posterior_nu<-as.matrix(posterior_nu)
colnames(posterior_nu)<-gsub("Cell.","",colnames(posterior_nu))
colnames(posterior_nu)<-gsub("_Batch1","",colnames(posterior_nu))
nu1<-read.table(paste0(data.path,"/simulate_",index_run,"_nu1.txt"))
ci_hdi_nu_89<-rep(0,ncol(posterior_nu))
ci_hdi_nu_low_89<-rep(0,ncol(posterior_nu))
ci_hdi_nu_high_89<-rep(0,ncol(posterior_nu))
ci_hdi_nu_50<-rep(0,ncol(posterior_nu))
ci_hdi_nu_low_50<-rep(0,ncol(posterior_nu))
ci_hdi_nu_high_50<-rep(0,ncol(posterior_nu))
for (i in 1:ncol(posterior_nu)){
  cihdi_89<-ci(posterior_nu[,i],method="HDI")
  ci_hdi_nu_89[i]=cihdi_89$CI
  ci_hdi_nu_low_89[i]=cihdi_89$CI_low
  ci_hdi_nu_high_89[i]=cihdi_89$CI_high
  cihdi_50<-ci(posterior_nu[,i],ci=0.5,method="HDI")
  ci_hdi_nu_50[i]=cihdi_50$CI
  ci_hdi_nu_low_50[i]=cihdi_50$CI_low
  ci_hdi_nu_high_50[i]=cihdi_50$CI_high
}
nu1$true.nu<-nu1$x
nu1$cell.name<-colnames(posterior_nu)
hdi_89<-data.frame(cell.name=colnames(posterior_nu),upper=ci_hdi_nu_high_89,lower=ci_hdi_nu_low_89,Median=colMedians(posterior_nu),Mean=colMeans(posterior_nu))
hdi_50<-data.frame(cell.name=colnames(posterior_nu),upper=ci_hdi_nu_high_50,lower=ci_hdi_nu_low_50,Median=colMedians(posterior_nu),Mean=colMeans(posterior_nu))
#cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
p <- ggplot() + # fill=name allow to automatically dedicate a color for each group
  geom_errorbar(data=hdi_89,aes(x=cell.name, ymin=lower,ymax=upper),size=1,width=0.5,color="royalblue")+
  geom_errorbar(data=hdi_50,aes(x=cell.name, ymin=lower,ymax=upper),size=1,width=0.5,color="orange")+
  geom_point(data=nu1,aes(x=cell.name,y=true.nu))+
  geom_point(data=hdi_89,aes(x=cell.name,y=Median,shape=17,colour="green4"))+
  geom_point(data=hdi_89,aes(x=cell.name,y=Mean,shape=15,colour="#009E73"))+
  scale_shape_identity()+
  scale_colour_discrete(name  ="Estimation",
                        breaks=c("green4", "#009E73"),
                        labels=c("Median", "Mean")) +
  #scale_shape_discrete(name  ="Estimation",
  #                     breaks=c(17, 15),
  #                    labels=c("Median", "Mean"))+
  labs(title = "posterior Credible Interval of nu",
       subtitle = paste0("run ",index_run),
       y="nu")+
  theme(axis.text.x = element_text(size=5)) 
p

posterior_phi<-read.table(paste0(data.path,"/sim_phi_",index_run,".txt"))
posterior_phi<-as.matrix(posterior_phi)
colnames(posterior_phi)<-gsub("Cell.","",colnames(posterior_phi))
colnames(posterior_phi)<-gsub("_Batch1","",colnames(posterior_phi))
phi1<-read.table(paste0(data.path,"/simulate_",index_run,"_phi1.txt"))
ci_hdi_phi_89<-rep(0,ncol(posterior_phi))
ci_hdi_phi_low_89<-rep(0,ncol(posterior_phi))
ci_hdi_phi_high_89<-rep(0,ncol(posterior_phi))
ci_hdi_phi_50<-rep(0,ncol(posterior_phi))
ci_hdi_phi_low_50<-rep(0,ncol(posterior_phi))
ci_hdi_phi_high_50<-rep(0,ncol(posterior_phi))
for (i in 1:ncol(posterior_phi)){
  cihdi_89<-ci(posterior_phi[,i],method="HDI")
  ci_hdi_phi_89[i]=cihdi_89$CI
  ci_hdi_phi_low_89[i]=cihdi_89$CI_low
  ci_hdi_phi_high_89[i]=cihdi_89$CI_high
  cihdi_50<-ci(posterior_phi[,i],ci=0.5,method="HDI")
  ci_hdi_phi_50[i]=cihdi_50$CI
  ci_hdi_phi_low_50[i]=cihdi_50$CI_low
  ci_hdi_phi_high_50[i]=cihdi_50$CI_high
}
phi1$true.phi<-phi1$x
phi1$cell.name<-colnames(posterior_phi)
hdi_89<-data.frame(cell.name=colnames(posterior_phi),upper=ci_hdi_phi_high_89,lower=ci_hdi_phi_low_89,Median=colMedians(posterior_phi),Mean=colMeans(posterior_phi))
hdi_50<-data.frame(cell.name=colnames(posterior_phi),upper=ci_hdi_phi_high_50,lower=ci_hdi_phi_low_50,Median=colMedians(posterior_phi),Mean=colMeans(posterior_phi))
#cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
p <- ggplot() + # fill=name allow to automatically dedicate a color for each group
  geom_errorbar(data=hdi_89,aes(x=cell.name, ymin=lower,ymax=upper),size=1,width=0.5,color="royalblue")+
  geom_errorbar(data=hdi_50,aes(x=cell.name, ymin=lower,ymax=upper),size=1,width=0.5,color="orange")+
  geom_point(data=phi1,aes(x=cell.name,y=true.phi))+
  geom_point(data=hdi_89,aes(x=cell.name,y=Median,shape=17,colour="green4"))+
  geom_point(data=hdi_89,aes(x=cell.name,y=Mean,shape=15,colour="#009E73"))+
  scale_shape_identity()+
  scale_colour_discrete(name  ="Estimation",
                        breaks=c("green4", "#009E73"),
                        labels=c("Median", "Mean")) +
  #scale_shape_discrete(name  ="Estimation",
  #                     breaks=c(17, 15),
  #                    labels=c("Median", "Mean"))+
  labs(title = "posterior Credible Interval of phi",
       subtitle = paste0("run ",index_run),
       y="phi")+
  theme(axis.text.x = element_text(size=5)) 
p

posterior_s<-read.table(paste0(data.path,"/sim_s_",index_run,".txt"))
posterior_s<-as.matrix(posterior_s)
colnames(posterior_s)<-gsub("Cell.","",colnames(posterior_s))
colnames(posterior_s)<-gsub("_Batch1","",colnames(posterior_s))
s1<-read.table(paste0(data.path,"/simulate_",index_run,"_s1.txt"))
ci_hdi_s_89<-rep(0,ncol(posterior_s))
ci_hdi_s_low_89<-rep(0,ncol(posterior_s))
ci_hdi_s_high_89<-rep(0,ncol(posterior_s))
ci_hdi_s_50<-rep(0,ncol(posterior_s))
ci_hdi_s_low_50<-rep(0,ncol(posterior_s))
ci_hdi_s_high_50<-rep(0,ncol(posterior_s))
for (i in 1:ncol(posterior_s)){
  cihdi_89<-ci(posterior_s[,i],method="HDI")
  ci_hdi_s_89[i]=cihdi_89$CI
  ci_hdi_s_low_89[i]=cihdi_89$CI_low
  ci_hdi_s_high_89[i]=cihdi_89$CI_high
  cihdi_50<-ci(posterior_s[,i],ci=0.5,method="HDI")
  ci_hdi_s_50[i]=cihdi_50$CI
  ci_hdi_s_low_50[i]=cihdi_50$CI_low
  ci_hdi_s_high_50[i]=cihdi_50$CI_high
}
s1$true.s<-s1$x
s1$cell.name<-colnames(posterior_s)
hdi_89<-data.frame(cell.name=colnames(posterior_s),upper=ci_hdi_s_high_89,lower=ci_hdi_s_low_89,Median=colMedians(posterior_s),Mean=colMeans(posterior_s))
hdi_50<-data.frame(cell.name=colnames(posterior_s),upper=ci_hdi_s_high_50,lower=ci_hdi_s_low_50,Median=colMedians(posterior_s),Mean=colMeans(posterior_s))
#cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
p <- ggplot() + # fill=name allow to automatically dedicate a color for each group
  geom_errorbar(data=hdi_89,aes(x=cell.name, ymin=lower,ymax=upper),size=1,width=0.5,color="royalblue")+
  geom_errorbar(data=hdi_50,aes(x=cell.name, ymin=lower,ymax=upper),size=1,width=0.5,color="orange")+
  geom_point(data=s1,aes(x=cell.name,y=true.s))+
  geom_point(data=hdi_89,aes(x=cell.name,y=Median,shape=17,colour="green4"))+
  geom_point(data=hdi_89,aes(x=cell.name,y=Mean,shape=15,colour="#009E73"))+
  scale_shape_identity()+
  scale_colour_discrete(name  ="Estimation",
                        breaks=c("green4", "#009E73"),
                        labels=c("Median", "Mean")) +
  #scale_shape_discrete(name  ="Estimation",
  #                     breaks=c(17, 15),
  #                    labels=c("Median", "Mean"))+
  labs(title = "posterior Credible Interval of s",
       subtitle = paste0("run ",index_run),
       y="s")+
  theme(axis.text.x = element_text(size=5)) 
p

posterior_theta<-read.table(paste0(data.path,"/sim_theta_",index_run,".txt"))
posterior_theta<-as.matrix(posterior_theta)
colnames(posterior_theta)<-gsub("Cell.","",colnames(posterior_theta))
colnames(posterior_theta)<-gsub("_Batch1","",colnames(posterior_theta))
theta1<-read.table(paste0(data.path,"/simulate_",index_run,"_theta1.txt"))
ci_hdi_theta_89<-rep(0,ncol(posterior_theta))
ci_hdi_theta_low_89<-rep(0,ncol(posterior_theta))
ci_hdi_theta_high_89<-rep(0,ncol(posterior_theta))
ci_hdi_theta_50<-rep(0,ncol(posterior_theta))
ci_hdi_theta_low_50<-rep(0,ncol(posterior_theta))
ci_hdi_theta_high_50<-rep(0,ncol(posterior_theta))
for (i in 1:ncol(posterior_theta)){
  cihdi_89<-ci(posterior_theta[,i],method="HDI")
  ci_hdi_theta_89[i]=cihdi_89$CI
  ci_hdi_theta_low_89[i]=cihdi_89$CI_low
  ci_hdi_theta_high_89[i]=cihdi_89$CI_high
  cihdi_50<-ci(posterior_theta[,i],ci=0.5,method="HDI")
  ci_hdi_theta_50[i]=cihdi_50$CI
  ci_hdi_theta_low_50[i]=cihdi_50$CI_low
  ci_hdi_theta_high_50[i]=cihdi_50$CI_high
}
theta1$true.theta<-theta1$x
theta1$cell.name<-colnames(posterior_theta)
hdi_89<-data.frame(cell.name=colnames(posterior_theta),upper=ci_hdi_theta_high_89,lower=ci_hdi_theta_low_89,Median=colMedians(posterior_theta),Mean=colMeans(posterior_theta))
hdi_50<-data.frame(cell.name=colnames(posterior_theta),upper=ci_hdi_theta_high_50,lower=ci_hdi_theta_low_50,Median=colMedians(posterior_theta),Mean=colMeans(posterior_theta))
#cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
p <- ggplot() + # fill=name allow to automatically dedicate a color for each group
  geom_errorbar(data=hdi_89,aes(x=cell.name, ymin=lower,ymax=upper),size=1,width=0.5,color="royalblue")+
  geom_errorbar(data=hdi_50,aes(x=cell.name, ymin=lower,ymax=upper),size=1,width=0.5,color="orange")+
  geom_point(data=theta1,aes(x=cell.name,y=true.theta))+
  geom_point(data=hdi_89,aes(x=cell.name,y=Median,shape=17,colour="green4"))+
  geom_point(data=hdi_89,aes(x=cell.name,y=Mean,shape=15,colour="#009E73"))+
  scale_shape_identity()+
  scale_colour_discrete(name  ="Estimation",
                        breaks=c("green4", "#009E73"),
                        labels=c("Median", "Mean")) +
  #scale_thetahape_discrete(name  ="Estimation",
  #                     breaks=c(17, 15),
  #                    labels=c("Median", "Mean"))+
  labs(title = "posterior Credible Interval of theta",
       subtitle = paste0("run ",index_run),
       y="theta",x="batch")+
  theme(axis.text.x = element_text(size=5)) 
p
