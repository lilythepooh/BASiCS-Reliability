library(distr)
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
delta_0<-rlnorm(n=1000,meanlog=0,sdlog=sqrt(0.5))

epsilon1<-0.25
pdf1<- function(x) (epsilon1*dgamma(x,shape=1,rate=1)+(1-epsilon1)*dlnorm(x,meanlog=0,sdlog=sqrt(0.5))) # probability density function
dist1 <-AbscontDistribution(d=pdf1)  # signature for a dist with pdf ~ p
rdist1 <- r(dist1)                 # function to create random variates from p
delta_1<-rdist1(1000)


epsilon2<-0.5
pdf2<- function(x) (epsilon2*dgamma(x,shape=1,rate=1)+(1-epsilon2)*dlnorm(x,meanlog=0,sdlog=sqrt(0.5))) # probability density function
dist2 <-AbscontDistribution(d=pdf2)  # signature for a dist with pdf ~ p
rdist2 <- r(dist2)                 # function to create random variates from p
delta_2<-rdist2(1000)

epsilon3<-0.75
pdf3<- function(x) (epsilon3*dgamma(x,shape=1,rate=1)+(1-epsilon3)*dlnorm(x,meanlog=0,sdlog=sqrt(0.5))) # probability density function
dist3 <-AbscontDistribution(d=pdf3)  # signature for a dist with pdf ~ p
rdist3 <- r(dist3)                 # function to create random variates from p
delta_3<-rdist3(1000)

delta_4<-rgamma(n=1000,shape=1,rate=1)

df<-data.frame(epsilon=factor(rep(c("0","0.25","0.5","0.75","1"),each=1000)),
               delta=c(delta_0,delta_1,delta_2,delta_3,delta_4))
df<- data.frame(df,Df = rep(c("\u03B5=0","\u03B5=0.25","\u03B5=0.5",
                              "\u03B5=0.75","\u03B5=1"),
                            times=c(length(delta_0),
                                    length(delta_1),
                                    length(delta_2),
                                    length(delta_3),
                                    length(delta_4))))
den0<-estimate_density(delta_0)
den2<-estimate_density(delta_2)
den4<-estimate_density(delta_4)
#df<-data.frame(epsilon=factor(rep(c("0","0.5","1"),each=length(den0$x))),
               #delta=c(den0$x,den2$x,den4$x),density=c(den0$y,den2$y,den4$y))
df<-data.frame(epsilon=factor(rep(c("0","0.5","1"),each=length(delta_0))),
               delta=c(delta_0,delta_2,delta_4))
df<- data.frame(df,Df = rep(c("\u03B5=0","\u03B5=0.5","\u03B5=1"),
                          times=c(length(delta_0),
                                length(delta_2),
                               length(delta_4))))
#df$density<-c((estimate_density(delta_0))$y,estimate_density(delta_2),estimate_density(delta_4))
prior_delta<-ggplot(df,aes(x=delta,color=epsilon,linetype=epsilon)) +
    theme_classic() +
  geom_density(size=1,show.legend = F)+
  stat_density(geom="line",position="identity", size = 0) + 
    # The posterior1
    #geom_area(color = "indianred2", fill = NA, size =1, linetype="dotted",
             # data = estimate_density(delta_0, extend = TRUE)) +
    # The posterior2
    #geom_area(color = "seagreen3", fill = NA, size = 1,
             # data = estimate_density(delta_2, extend = TRUE))+
    # The prior
    #geom_area(color = "deeppink4", fill = NA, size = 1,linetype="dashed",
   #           data = estimate_density(delta_4, extend = TRUE))+
  scale_x_continuous(limits = c(0.000001, 3.1))+
  scale_linetype_manual(name=expression(epsilon),
                        values = c("dotted", "solid","dashed"),
                        labels=c("\u03B5=0",
                         "\u03B5=0.5",
                      "\u03B5=1")) +
  scale_color_manual(name=expression(epsilon),
                 values =c("indianred2",
                         "seagreen3",
                         "deeppink4"), 
                labels=c("\u03B5=0",
                         "\u03B5=0.5",
                        "\u03B5=1"))+
 # guides(color = guide_legend(override.aes = list(linetype=c("0"="indianred2",
  #                                                           "0.5"="seagreen3",
   #                                                          "1"="deeppink4"),
    #                                              fill=c("0"="indianred2",
      #                                                         "0.5"="seagreen3",
    #                                                             "1"="deeppink4")))) +
  theme(axis.text.y = element_text(),
        plot.margin=unit(c(0.3,0,0,0.3), "cm"))+
  theme(axis.text.x = element_text(),axis.title.x = element_text())+
  labs(x=expression(delta[i]),title=paste("priors with different \u03B5"),color  = expression(epsilon), linetype = expression(epsilon))
print(prior_delta)

p_delta<-ggplot(aes(x=delta, group=epsilon,color=as.factor(epsilon),shape=as.factor(epsilon)), data=df)+
  geom_density(size=1)+
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
  theme(axis.text.y = element_text(),
        plot.margin=unit(c(0.3,0,0,0.3), "cm"))+
  theme(axis.text.x = element_text(),axis.title.x = element_text())+
  labs(x=expression(delta[i]),title=paste("priors with different \u03B5"))+
  facet_wrap(~Df,nrow=5,strip.position='right')
print(p_delta)


p_delta<-ggplot(aes(x=delta, group=epsilon,color=as.factor(epsilon),shape=as.factor(epsilon)), data=df)+
  geom_density(size=1)+
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
  theme(axis.text.y = element_text(),
        plot.margin=unit(c(0.3,0,0,0.3), "cm"))+
  theme(axis.text.x = element_text(),axis.title.x = element_text())+
  labs(x=expression(delta[i]),title=paste("priors with different \u03B5"))+
  facet_wrap(~Df,nrow=5,strip.position='right')
print(p_delta)
