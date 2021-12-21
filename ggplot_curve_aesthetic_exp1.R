library(lemon)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(reshape2)
library(ggridges)
library(viridis)
library(hrbrthemes)
library(grDevices)
library(gridExtra)
library(gtable)
library(grid)
#library(lattice)
setwd("./non_regression_BASiCS/fixed_dataset/epsilon=0")
data.path.epsilon0<-getwd()
setwd("./non_regression_BASiCS/fixed_dataset/epsilon=0.25")
data.path.epsilon1quarter<-getwd()
setwd("./non_regression_BASiCS/fixed_dataset/epsilon=0.5")
data.path.epsilon2quarter<-getwd()
setwd("./non_regression_BASiCS/fixed_dataset/epsilon=0.75")
data.path.epsilon3quarter<-getwd()
setwd("./non_regression_BAsiCS/fixed_dataset/epsilon=1")
data.path.epsilon4quarter<-getwd()
setwd("./non_regression_BASiCS/fixed_dataset/epsilons")
data.path.epsilons<-getwd()
delta1<-read.table(paste0(data.path.epsilon0,"/simulate_delta1.txt"))
delta1<-as.matrix(delta1)
mu1<-read.table(paste0(data.path.epsilon0,"/simulate_mu1.txt"))
mu1<-as.matrix(mu1)
n_bio_gene=nrow(delta1)
mu1<-mu1[1:n_bio_gene,]
#pdf("0_bio_gene_curve_stochsticity.pdf",paper="USr",width=17)
cairo_pdf("0_bio_gene_curve_stochsticity.pdf",width=13,height=8,onefile=T)
for (i in (1:n_bio_gene)){
  i_posterior_delta.epsilon0<-read.table(paste0(data.path.epsilons,"/",i,"_posterior_delta.epsilon0.txt"))
  i_df_posterior_delta.epsilon0<-as.data.frame(t(as.matrix(i_posterior_delta.epsilon0)))
i_df_posterior_delta.epsilon0$name<-sprintf("run.%s",seq(1:200))
i_df_posterior_delta.epsilon0$epsilon<-rep(0,200)
i_df_posterior_delta.epsilon0.long = melt(i_df_posterior_delta.epsilon0, id.vars= "name")
i_df_posterior_delta.epsilon0.long$epsilon<-rep(0,length(i_df_posterior_delta.epsilon0.long$value))

i_posterior_delta.epsilon1q<-read.table(paste0(data.path.epsilons,"/",i,"_posterior_delta.epsilon1q.txt"))
i_df_posterior_delta.epsilon1q<-as.data.frame(t(as.matrix(i_posterior_delta.epsilon1q)))
i_df_posterior_delta.epsilon1q$name<-sprintf("run.%s",seq(1:200))
#i_df_posterior_delta.epsilon1q$epsilon<-rep(0.25,200)
i_df_posterior_delta.epsilon1q.long = melt(i_df_posterior_delta.epsilon1q, id.vars= "name")
i_df_posterior_delta.epsilon1q.long$epsilon=rep(0.25,length(i_df_posterior_delta.epsilon1q.long$value))

i_posterior_delta.epsilon2q<-read.table(paste0(data.path.epsilons,"/",i,"_posterior_delta.epsilon2q.txt"))
i_df_posterior_delta.epsilon2q<-as.data.frame(t(as.matrix(i_posterior_delta.epsilon2q)))
i_df_posterior_delta.epsilon2q$name<-sprintf("run.%s",seq(1:200))
i_df_posterior_delta.epsilon2q.long = melt(i_df_posterior_delta.epsilon2q, id.vars= "name")
i_df_posterior_delta.epsilon2q.long$epsilon=rep(0.5,length(i_df_posterior_delta.epsilon2q.long$value))

i_posterior_delta.epsilon3q<-read.table(paste0(data.path.epsilons,"/",i,"_posterior_delta.epsilon3q.txt"))
i_df_posterior_delta.epsilon3q<-as.data.frame(t(as.matrix(i_posterior_delta.epsilon3q)))
i_df_posterior_delta.epsilon3q$name<-sprintf("run.%s",seq(1:200))
i_df_posterior_delta.epsilon3q.long = melt(i_df_posterior_delta.epsilon3q, id.vars= "name")
i_df_posterior_delta.epsilon3q.long$epsilon=rep(0.75,length(i_df_posterior_delta.epsilon3q.long$value))

i_posterior_delta.epsilon4q<-read.table(paste0(data.path.epsilons,"/",i,"_posterior_delta.epsilon4q.txt"))
i_df_posterior_delta.epsilon4q<-as.data.frame(t(as.matrix(i_posterior_delta.epsilon4q)))
i_df_posterior_delta.epsilon4q$name<-sprintf("run.%s",seq(1:200))
i_df_posterior_delta.epsilon4q.long = melt(i_df_posterior_delta.epsilon4q, id.vars= "name")
i_df_posterior_delta.epsilon4q.long$epsilon=rep(1,length(i_df_posterior_delta.epsilon4q.long$value))


i_df_posterior_delta.long=rbind(i_df_posterior_delta.epsilon0.long,
                                i_df_posterior_delta.epsilon1q.long,
                                i_df_posterior_delta.epsilon2q.long,
                                i_df_posterior_delta.epsilon3q.long,
                                i_df_posterior_delta.epsilon4q.long)
i_df_posterior_delta.long<- data.frame(i_df_posterior_delta.long,
                                       Df = rep(c("\u03B5=0","\u03B5=0.25",
                                                  "\u03B5=0.5","\u03B5=0.75",
                                                  "\u03B5=1"),
                              times=c(nrow(i_df_posterior_delta.epsilon0.long),
                                      nrow(i_df_posterior_delta.epsilon1q.long),
                                      nrow(i_df_posterior_delta.epsilon2q.long),
                                      nrow(i_df_posterior_delta.epsilon3q.long),
                                      nrow(i_df_posterior_delta.epsilon4q.long))))
p_delta<-ggplot(aes(x=value, group=name,color=as.factor(epsilon),shape=as.factor(epsilon)), data=i_df_posterior_delta.long)+
  geom_density()+
 scale_color_manual(name=expression(epsilon),
                  values =c("0"="indianred2", "0.25"="Orange", 
                            "0.5"="seagreen3", "0.75"="royalblue",
                            "1"="deeppink4"), 
                  labels=c("\u03B5=0","\u03B5=0.25",
                          "\u03B5=0.5","\u03B5=0.75",
                          "\u03B5=1"))+
  guides(color = guide_legend(override.aes = list(shape=21,fill=c("0"="indianred2", "0.25"="Orange", 
                                                                  "0.5"="seagreen3", "0.75"="royalblue",
                                                                  "1"="deeppink4")))) +
  geom_vline(xintercept=delta1[i],color="black", size=1)+
  theme(axis.text.y = element_text(),
        plot.margin=unit(c(0.3,0,0,0.3), "cm"))+
  theme(axis.text.x = element_text(),axis.title.x = element_text())+
  labs(x=expression(delta[i]))+
  facet_wrap(~Df,nrow=5,strip.position='right')
  

i_posterior_mu.epsilon0<-read.table(paste0(data.path.epsilons,"/",i,"_posterior_mu.epsilon0.txt"))
i_df_posterior_mu.epsilon0<-as.data.frame(t(as.matrix(i_posterior_mu.epsilon0)))
i_df_posterior_mu.epsilon0$name<-sprintf("run.%s",seq(1:200))
i_df_posterior_mu.epsilon0.long = melt(i_df_posterior_mu.epsilon0, id.vars= "name")
i_df_posterior_mu.epsilon0.long$epsilon<-rep(0,length(i_df_posterior_mu.epsilon0.long$value))

i_posterior_mu.epsilon1q<-read.table(paste0(data.path.epsilons,"/",i,"_posterior_mu.epsilon1q.txt"))
i_df_posterior_mu.epsilon1q<-as.data.frame(t(as.matrix(i_posterior_mu.epsilon1q)))
i_df_posterior_mu.epsilon1q$name<-sprintf("run.%s",seq(1:200))
i_df_posterior_mu.epsilon1q.long = melt(i_df_posterior_mu.epsilon1q, id.vars= "name")
i_df_posterior_mu.epsilon1q.long$epsilon<-rep(0.25,length(i_df_posterior_mu.epsilon1q.long$value))

i_posterior_mu.epsilon2q<-read.table(paste0(data.path.epsilons,"/",i,"_posterior_mu.epsilon2q.txt"))
i_df_posterior_mu.epsilon2q<-as.data.frame(t(as.matrix(i_posterior_mu.epsilon2q)))
i_df_posterior_mu.epsilon2q$name<-sprintf("run.%s",seq(1:200))
i_df_posterior_mu.epsilon2q.long = melt(i_df_posterior_mu.epsilon2q, id.vars= "name")
i_df_posterior_mu.epsilon2q.long$epsilon<-rep(0.5,length(i_df_posterior_mu.epsilon2q.long$value))

i_posterior_mu.epsilon3q<-read.table(paste0(data.path.epsilons,"/",i,"_posterior_mu.epsilon3q.txt"))
i_df_posterior_mu.epsilon3q<-as.data.frame(t(as.matrix(i_posterior_mu.epsilon3q)))
i_df_posterior_mu.epsilon3q$name<-sprintf("run.%s",seq(1:200))
i_df_posterior_mu.epsilon3q.long = melt(i_df_posterior_mu.epsilon3q, id.vars= "name")
i_df_posterior_mu.epsilon3q.long$epsilon<-rep(0.75,length(i_df_posterior_mu.epsilon3q.long$value))

i_posterior_mu.epsilon4q<-read.table(paste0(data.path.epsilons,"/",i,"_posterior_mu.epsilon4q.txt"))
i_df_posterior_mu.epsilon4q<-as.data.frame(t(as.matrix(i_posterior_mu.epsilon4q)))
i_df_posterior_mu.epsilon4q$name<-sprintf("run.%s",seq(1:200))
i_df_posterior_mu.epsilon4q.long = melt(i_df_posterior_mu.epsilon4q, id.vars= "name")
i_df_posterior_mu.epsilon4q.long$epsilon=rep(1,length(i_df_posterior_mu.epsilon4q.long$value))

i_df_posterior_mu.long=rbind(i_df_posterior_mu.epsilon0.long,
                             i_df_posterior_mu.epsilon1q.long,
                             i_df_posterior_mu.epsilon2q.long,
                             i_df_posterior_mu.epsilon3q.long,
                             i_df_posterior_mu.epsilon4q.long)
i_df_posterior_mu.long<- data.frame(i_df_posterior_mu.long,
                                       Df = rep(c("\u03B5=0","\u03B5=0.25",
                                                  "\u03B5=0.5","\u03B5=0.75",
                                                  "\u03B5=1"),
                                                times=c(nrow(i_df_posterior_mu.epsilon0.long),
                                                        nrow(i_df_posterior_mu.epsilon1q.long),
                                                        nrow(i_df_posterior_mu.epsilon2q.long),
                                                        nrow(i_df_posterior_mu.epsilon3q.long),
                                                        nrow(i_df_posterior_mu.epsilon4q.long))))
p_mu<-ggplot(aes(x=value, group=name,color=as.factor(epsilon),shape=as.factor(epsilon)), data=i_df_posterior_mu.long)+
  geom_density()+
  scale_color_manual(name=expression(epsilon),
                     values =c("0"="indianred2", "0.25"="Orange", 
                               "0.5"="seagreen3", "0.75"="royalblue",
                               "1"="deeppink4"), 
                     labels=c("\u03B5=0","\u03B5=0.25",
                              "\u03B5=0.5","\u03B5=0.75",
                              "\u03B5=1"))+
  guides(color = guide_legend(override.aes = list(shape=21,fill=c("0"="indianred2", "0.25"="Orange", 
                                                                  "0.5"="seagreen3", "0.75"="royalblue",
                                                                  "1"="deeppink4")))) +
  geom_vline(xintercept=mu1[i],color="black", size=1)+
  theme(axis.text.y = element_text(),
        plot.margin=unit(c(0.3,0,0,0.3), "cm"))+
  theme(axis.text.x = element_text(),axis.title.x = element_text())+
  labs(x=expression(mu[i]))+
  facet_wrap(~Df,nrow=5,strip.position='right')

#grid.arrange(arrangeGrob(p_delta+theme(legend.position="none"), p_mu+theme(legend.position="bottom"),ncol=2),
                   #      top="Posteriors of gene-specific parameters from 200 replications")
grid_arrange_shared_legend(p_delta+theme(strip.background = element_blank(),strip.text.y = element_blank()), 
                           gridExtra::arrangeGrob(p_mu+theme(legend.position="none",axis.title.y = element_blank()), ncol=1), ncol=2, nrow=1, 
                           top=paste("Posteriors of gene-specific parameters for gene i=",i,", from 200 replications"))
}
dev.off()


#pdf
cairo_pdf("0_cell_curve_stochsticity.pdf",width=13,height=8,onefile=T)
nu1<-read.table(paste0(data.path.epsilon0,"/simulate_nu1.txt"))
nu1<-as.matrix(nu1)
phi1<-read.table(paste0(data.path.epsilon0,"/simulate_phi1.txt"))
phi1<-as.matrix(phi1)
s1<-read.table(paste0(data.path.epsilon0,"/simulate_s1.txt"))
s1<-as.matrix(s1)
n_cell<-nrow(nu1)
for (j in (1:n_cell)){
  j_posterior_nu.epsilon0<-read.table(paste0(data.path.epsilons,"/",j,"_posterior_nu.epsilon0.txt"))
  j_df_posterior_nu.epsilon0<-as.data.frame(t(as.matrix(j_posterior_nu.epsilon0)))
  j_df_posterior_nu.epsilon0$name<-sprintf("run.%s",seq(1:200))
  j_df_posterior_nu.epsilon0.long = melt(j_df_posterior_nu.epsilon0, id.vars= "name")
  j_df_posterior_nu.epsilon0.long$epsilon=rep(0,length(j_df_posterior_nu.epsilon0.long$value))
  
  j_posterior_nu.epsilon1q<-read.table(paste0(data.path.epsilons,"/",j,"_posterior_nu.epsilon1q.txt"))
  j_df_posterior_nu.epsilon1q<-as.data.frame(t(as.matrix(j_posterior_nu.epsilon1q)))
  j_df_posterior_nu.epsilon1q$name<-sprintf("run.%s",seq(1:200))
  j_df_posterior_nu.epsilon1q.long = melt(j_df_posterior_nu.epsilon1q, id.vars= "name")
  j_df_posterior_nu.epsilon1q.long$epsilon=rep(0.25,length(j_df_posterior_nu.epsilon1q.long$value))
 
  j_posterior_nu.epsilon2q<-read.table(paste0(data.path.epsilons,"/",j,"_posterior_nu.epsilon2q.txt"))
  j_df_posterior_nu.epsilon2q<-as.data.frame(t(as.matrix(j_posterior_nu.epsilon2q)))
  j_df_posterior_nu.epsilon2q$name<-sprintf("run.%s",seq(1:200))
  j_df_posterior_nu.epsilon2q.long = melt(j_df_posterior_nu.epsilon2q, id.vars= "name")
  j_df_posterior_nu.epsilon2q.long$epsilon=rep(0.5,length(j_df_posterior_nu.epsilon2q.long$value))
  
  j_posterior_nu.epsilon3q<-read.table(paste0(data.path.epsilons,"/",j,"_posterior_nu.epsilon3q.txt"))
  j_df_posterior_nu.epsilon3q<-as.data.frame(t(as.matrix(j_posterior_nu.epsilon3q)))
  j_df_posterior_nu.epsilon3q$name<-sprintf("run.%s",seq(1:200))
  j_df_posterior_nu.epsilon3q.long = melt(j_df_posterior_nu.epsilon3q, id.vars= "name")
  j_df_posterior_nu.epsilon3q.long$epsilon=rep(0.75,length(j_df_posterior_nu.epsilon3q.long$value))
  
  j_posterior_nu.epsilon4q<-read.table(paste0(data.path.epsilons,"/",j,"_posterior_nu.epsilon4q.txt"))
  j_df_posterior_nu.epsilon4q<-as.data.frame(t(as.matrix(j_posterior_nu.epsilon4q)))
  j_df_posterior_nu.epsilon4q$name<-sprintf("run.%s",seq(1:200))
  j_df_posterior_nu.epsilon4q.long = melt(j_df_posterior_nu.epsilon4q, id.vars= "name")
  j_df_posterior_nu.epsilon4q.long$epsilon=rep(1,length(j_df_posterior_nu.epsilon4q.long$value))
  
  j_df_posterior_nu.long=rbind(j_df_posterior_nu.epsilon0.long,
                               j_df_posterior_nu.epsilon1q.long,
                               j_df_posterior_nu.epsilon2q.long,
                               j_df_posterior_nu.epsilon3q.long,
                               j_df_posterior_nu.epsilon4q.long)
  j_df_posterior_nu.long<- data.frame(j_df_posterior_nu.long,
                                      Df = rep(c("\u03B5=0","\u03B5=0.25",
                                                 "\u03B5=0.5","\u03B5=0.75",
                                                 "\u03B5=1"),
                                               times=c(nrow(j_df_posterior_nu.epsilon0.long),
                                                       nrow(j_df_posterior_nu.epsilon1q.long),
                                                       nrow(j_df_posterior_nu.epsilon2q.long),
                                                       nrow(j_df_posterior_nu.epsilon3q.long),
                                                       nrow(j_df_posterior_nu.epsilon4q.long))))
                                                       
  p_nu<-ggplot(aes(x=value, group=name,color=as.factor(epsilon),shape=as.factor(epsilon)), data=j_df_posterior_nu.long)+
    geom_density()+
    scale_color_manual(name=expression(epsilon),
                       values =c("0"="indianred2", "0.25"="Orange", 
                                 "0.5"="seagreen3", "0.75"="royalblue",
                                 "1"="deeppink4"), 
                       labels=c("\u03B5=0","\u03B5=0.25",
                                "\u03B5=0.5","\u03B5=0.75",
                                "\u03B5=1"))+
    guides(color = guide_legend(override.aes = list(shape=21,fill=c("0"="indianred2", "0.25"="Orange", 
                                                                    "0.5"="seagreen3", "0.75"="royalblue",
                                                                    "1"="deeppink4")))) +
    geom_vline(xintercept=nu1[j],color="black", size=1)+
    theme(axis.text.y = element_text(),
          plot.margin=unit(c(0.3,0,0,0.3), "cm"))+
    theme(axis.text.x = element_text(),axis.title.x = element_text())+
    labs(x=expression(nu[j]))+
    facet_wrap(~Df,nrow=5,strip.position='right')

  
  j_posterior_phi.epsilon0<-read.table(paste0(data.path.epsilons,"/",j,"_posterior_phi.epsilon0.txt"))
  j_df_posterior_phi.epsilon0<-as.data.frame(t(as.matrix(j_posterior_phi.epsilon0)))
  j_df_posterior_phi.epsilon0$name<-sprintf("run.%s",seq(1:200))
  j_df_posterior_phi.epsilon0.long = melt(j_df_posterior_phi.epsilon0, id.vars= "name")
  j_df_posterior_phi.epsilon0.long$epsilon<-rep(0,length(j_df_posterior_phi.epsilon0.long$value))
  
  j_posterior_phi.epsilon1q<-read.table(paste0(data.path.epsilons,"/",j,"_posterior_phi.epsilon1q.txt"))
  j_df_posterior_phi.epsilon1q<-as.data.frame(t(as.matrix(j_posterior_phi.epsilon1q)))
  j_df_posterior_phi.epsilon1q$name<-sprintf("run.%s",seq(1:200))
  j_df_posterior_phi.epsilon1q.long = melt(j_df_posterior_phi.epsilon1q, id.vars= "name")
  j_df_posterior_phi.epsilon1q.long$epsilon<-rep(0.25,length(j_df_posterior_phi.epsilon1q.long$value))
  
  j_posterior_phi.epsilon2q<-read.table(paste0(data.path.epsilons,"/",j,"_posterior_phi.epsilon2q.txt"))
  j_df_posterior_phi.epsilon2q<-as.data.frame(t(as.matrix(j_posterior_phi.epsilon2q)))
  j_df_posterior_phi.epsilon2q$name<-sprintf("run.%s",seq(1:200))
  j_df_posterior_phi.epsilon2q.long = melt(j_df_posterior_phi.epsilon2q, id.vars= "name")
  j_df_posterior_phi.epsilon2q.long$epsilon<-rep(0.5,length(j_df_posterior_phi.epsilon2q.long$value))
  
  j_posterior_phi.epsilon3q<-read.table(paste0(data.path.epsilons,"/",j,"_posterior_phi.epsilon3q.txt"))
  j_df_posterior_phi.epsilon3q<-as.data.frame(t(as.matrix(j_posterior_phi.epsilon3q)))
  j_df_posterior_phi.epsilon3q$name<-sprintf("run.%s",seq(1:200))
  j_df_posterior_phi.epsilon3q.long = melt(j_df_posterior_phi.epsilon3q, id.vars= "name")
  j_df_posterior_phi.epsilon3q.long$epsilon<-rep(0.75,length(j_df_posterior_phi.epsilon3q.long$value))
  
  j_posterior_phi.epsilon4q<-read.table(paste0(data.path.epsilons,"/",j,"_posterior_phi.epsilon4q.txt"))
  j_df_posterior_phi.epsilon4q<-as.data.frame(t(as.matrix(j_posterior_phi.epsilon4q)))
  j_df_posterior_phi.epsilon4q$name<-sprintf("run.%s",seq(1:200))
  j_df_posterior_phi.epsilon4q.long = melt(j_df_posterior_phi.epsilon4q, id.vars= "name")
  j_df_posterior_phi.epsilon4q.long$epsilon<-rep(1,length(j_df_posterior_phi.epsilon4q.long$value))
  
  j_df_posterior_phi.long=rbind(j_df_posterior_phi.epsilon0.long,
                               j_df_posterior_phi.epsilon1q.long,
                               j_df_posterior_phi.epsilon2q.long,
                               j_df_posterior_phi.epsilon3q.long,
                               j_df_posterior_phi.epsilon4q.long)
  j_df_posterior_phi.long<- data.frame(j_df_posterior_phi.long,
                                      Df = rep(c("\u03B5=0","\u03B5=0.25",
                                                 "\u03B5=0.5","\u03B5=0.75",
                                                 "\u03B5=1"),
                                               times=c(nrow(j_df_posterior_phi.epsilon0.long),
                                                       nrow(j_df_posterior_phi.epsilon1q.long),
                                                       nrow(j_df_posterior_phi.epsilon2q.long),
                                                       nrow(j_df_posterior_phi.epsilon3q.long),
                                                       nrow(j_df_posterior_phi.epsilon4q.long))))
  
  p_phi<-ggplot(aes(x=value, group=name,color=as.factor(epsilon),shape=as.factor(epsilon)), data=j_df_posterior_phi.long)+
    geom_density()+
    scale_color_manual(name=expression(epsilon),
                       values =c("0"="indianred2", "0.25"="Orange", 
                                 "0.5"="seagreen3", "0.75"="royalblue",
                                 "1"="deeppink4"), 
                       labels=c("\u03B5=0","\u03B5=0.25",
                                "\u03B5=0.5","\u03B5=0.75",
                                "\u03B5=1"))+
    guides(color = guide_legend(override.aes = list(shape=21,fill=c("0"="indianred2", "0.25"="Orange", 
                                                                    "0.5"="seagreen3", "0.75"="royalblue",
                                                                    "1"="deeppink4")))) +
    geom_vline(xintercept=phi1[j],color="black", size=1)+
    theme(axis.text.y = element_text(),
          plot.margin=unit(c(0.3,0,0,0.3), "cm"))+
    theme(axis.text.x = element_text(),axis.title.x = element_text())+
    labs(x=expression(Phi[j]))+
    facet_wrap(~Df,nrow=5,strip.position='right')
  

  j_posterior_s.epsilon0<-read.table(paste0(data.path.epsilons,"/",j,"_posterior_s.epsilon0.txt"))
  j_df_posterior_s.epsilon0<-as.data.frame(t(as.matrix(j_posterior_s.epsilon0)))
  j_df_posterior_s.epsilon0$name<-sprintf("run.%s",seq(1:200))
  j_df_posterior_s.epsilon0.long = melt(j_df_posterior_s.epsilon0, id.vars= "name")
  j_df_posterior_s.epsilon0.long$epsilon<-rep(0,length(j_df_posterior_s.epsilon0.long$value))
  
  
  j_posterior_s.epsilon1q<-read.table(paste0(data.path.epsilons,"/",j,"_posterior_s.epsilon1q.txt"))
  j_df_posterior_s.epsilon1q<-as.data.frame(t(as.matrix(j_posterior_s.epsilon1q)))
  j_df_posterior_s.epsilon1q$name<-sprintf("run.%s",seq(1:200))
  j_df_posterior_s.epsilon1q.long = melt(j_df_posterior_s.epsilon1q, id.vars= "name")
  j_df_posterior_s.epsilon1q.long$epsilon<-rep(0.25,length(j_df_posterior_s.epsilon1q.long$value))
  
  
  j_posterior_s.epsilon2q<-read.table(paste0(data.path.epsilons,"/",j,"_posterior_s.epsilon2q.txt"))
  j_df_posterior_s.epsilon2q<-as.data.frame(t(as.matrix(j_posterior_s.epsilon2q)))
  j_df_posterior_s.epsilon2q$name<-sprintf("run.%s",seq(1:200))
  j_df_posterior_s.epsilon2q.long = melt(j_df_posterior_s.epsilon2q, id.vars= "name")
  j_df_posterior_s.epsilon2q.long$epsilon<-rep(0.5,length(j_df_posterior_s.epsilon2q.long$value))
  
  j_posterior_s.epsilon3q<-read.table(paste0(data.path.epsilons,"/",j,"_posterior_s.epsilon3q.txt"))
  j_df_posterior_s.epsilon3q<-as.data.frame(t(as.matrix(j_posterior_s.epsilon3q)))
  j_df_posterior_s.epsilon3q$name<-sprintf("run.%s",seq(1:200))
  j_df_posterior_s.epsilon3q.long = melt(j_df_posterior_s.epsilon3q, id.vars= "name")
  j_df_posterior_s.epsilon3q.long$epsilon<-rep(0.75,length(j_df_posterior_s.epsilon3q.long$value))
  
  
  j_posterior_s.epsilon4q<-read.table(paste0(data.path.epsilons,"/",j,"_posterior_s.epsilon4q.txt"))
  j_df_posterior_s.epsilon4q<-as.data.frame(t(as.matrix(j_posterior_s.epsilon4q)))
  j_df_posterior_s.epsilon4q$name<-sprintf("run.%s",seq(1:200))
  j_df_posterior_s.epsilon4q.long = melt(j_df_posterior_s.epsilon4q, id.vars= "name")
  j_df_posterior_s.epsilon4q.long$epsilon<-rep(1,length(j_df_posterior_s.epsilon4q.long$value))
  
  j_df_posterior_s.long=rbind(j_df_posterior_s.epsilon0.long,
                               j_df_posterior_s.epsilon1q.long,
                               j_df_posterior_s.epsilon2q.long,
                               j_df_posterior_s.epsilon3q.long,
                               j_df_posterior_s.epsilon4q.long)
  j_df_posterior_s.long<- data.frame(j_df_posterior_s.long,
                                       Df = rep(c("\u03B5=0","\u03B5=0.25",
                                                  "\u03B5=0.5","\u03B5=0.75",
                                                  "\u03B5=1"),
                                                times=c(nrow(j_df_posterior_s.epsilon0.long),
                                                        nrow(j_df_posterior_s.epsilon1q.long),
                                                        nrow(j_df_posterior_s.epsilon2q.long),
                                                        nrow(j_df_posterior_s.epsilon3q.long),
                                                        nrow(j_df_posterior_s.epsilon4q.long))))
  
  p_s<-ggplot(aes(x=value, group=name,color=as.factor(epsilon),shape=as.factor(epsilon)), data=j_df_posterior_s.long)+
    geom_density()+
    scale_color_manual(name=expression(epsilon),
                       values =c("0"="indianred2", "0.25"="Orange", 
                                 "0.5"="seagreen3", "0.75"="royalblue",
                                 "1"="deeppink4"), 
                       labels=c("\u03B5=0","\u03B5=0.25",
                                "\u03B5=0.5","\u03B5=0.75",
                                "\u03B5=1"))+
    guides(color = guide_legend(override.aes = list(shape=21,fill=c("0"="indianred2", "0.25"="Orange", 
                                                                    "0.5"="seagreen3", "0.75"="royalblue",
                                                                    "1"="deeppink4")))) +
    geom_vline(xintercept=s1[j],color="black", size=1)+
    theme(axis.text.y = element_text(),
          plot.margin=unit(c(0.3,0,0,0.3), "cm"))+
    theme(axis.text.x = element_text(),axis.title.x = element_text())+
    labs(x=expression(s[j]))+
    facet_wrap(~Df,nrow=5,strip.position='right')
  
  grid_arrange_shared_legend(p_nu+theme(strip.background = element_blank(),strip.text.y = element_blank()), 
                             p_phi+theme(strip.background = element_blank(),strip.text.y = element_blank()),
                             gridExtra::arrangeGrob(p_s+theme(legend.position="none",axis.title.y = element_blank()), ncol=1), ncol=3, nrow=1, 
                             top=paste("Posteriors of cell-specific parameters for cell j=",j,", from 200 replications"))
  
}
dev.off()
