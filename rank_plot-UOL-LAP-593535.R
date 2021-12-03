plot(freq_true_delta_in_89hdi,type="p",pch=3,col="blue",ylim=c(0.5,1))
points(freq_true_mu_in_89hdi,type="p",pch=4,col="red")
abline(h=0.89)
library(bayesplot)
#setwd("C:/Users/hksca/OneDrive/Bayesian Big Data/BASiCS and Beyond/BASiCS Robustness/HPC BASiCS/2021/Simulation Experiment/Too many Simulations/Local Running")
#setwd("C:/Users/hksca/OneDrive/Bayesian Big Data/BASiCS and Beyond/BASiCS Robustness/HPC BASiCS/2021/Simulation Experiment/NewBASiCS/resimulate")
setwd("C:/Users/mmsli/OneDrive/Bayesian Big Data/BASiCS and Beyond/BASiCS Robustness/HPC BASiCS/2021/Simulation Experiment/NewBASiCS/resimulate")
data.path<-getwd()
ppc_ecdf_overlay_2 <- function (y, yrep, ..., pad = TRUE, size = 0.25, alpha = 0.7) 
{
        
        y <- bayesplot:::validate_y(y)
        yrep <- bayesplot:::validate_yrep(yrep, y)
        ggplot(bayesplot:::melt_yrep(yrep), aes_(x = ~value)) + hline_at(c(0, 
                                                                           0.5, 1), size = c(0.2, 0.1, 0.2), linetype = 2, color = bayesplot:::get_color("dh")) + 
                stat_ecdf(mapping = aes_(group = ~rep_id, color = "yrep"), 
                          geom = "step", size = size, alpha = alpha, pad = pad) + 
                stat_ecdf(data = data.frame(value = y), mapping = aes_(color = "y"), 
                          geom = c("step"), size = 1, pad = pad) + bayesplot:::scale_color_ppc_dist() + 
                xlab(bayesplot:::y_label()) + 
                scale_x_continuous(limits=c(0,50),expand=c(0,0)) + scale_y_continuous(limits=c(0,1),expand=c(0,0),breaks = c(0, 0.5, 
                                                                                                                              1)) + yaxis_title(FALSE) + xaxis_title(FALSE) + yaxis_ticks(FALSE)
}
i=1
j=1
#CI = qbinom(c(0.005,0.5,0.995), size=1000,prob  =  1/101)

rank_delta<-read.table(paste0(data.path,"/rank_delta.txt"))
rank_delta<-rank_delta[1:47,]
hist(rank_delta[!is.na(rank_delta[,i]),i],xlab=expression(rank(delta[i])))
abline(h=10,col="red",lwd=2,lty="dashed")
samps = matrix(sample(c(1:50),size = length(rank_delta[!is.na(rank_delta[,i]),i])*500, replace=T),500,length(rank_delta[!is.na(rank_delta[,i]),i]))
ppc_ecdf_overlay_2(rank_delta[!is.na(rank_delta[,i]),i],samps) + geom_abline(slope=1/50,intercept=0,colour="grey45",linetype="dashed")

rank_mu<-read.table(paste0(data.path,"/rank_mu.txt"))
hist(rank_mu[!is.na(rank_mu[,i]),i],xlab=expression(rank(mu[i])))
abline(h=10,col="red",lwd=2,lty="dashed")
#legend(x = "topright",          # Position
 #      legend = "uniform[0,10]",  # Legend texts
  #     lty = "dashed",           # Line types
   #    col = "red",           # Line colors
    #   lwd = 2)                 # Line width
samps = matrix(sample(c(1:50),size = length(rank_mu[!is.na(rank_mu[,i]),i])*500, replace=T),500,length(rank_mu[!is.na(rank_mu[,i]),i]))
ppc_ecdf_overlay_2(rank_mu[!is.na(rank_mu[,i]),i],samps) + geom_abline(slope=1/50,intercept=0,colour="grey45",linetype="dashed")

rank_nu<-read.table(paste0(data.path,"/rank_nu.txt"))
hist(rank_nu[!is.na(rank_nu[,j]),j],xlab=expression(rank(nu[j])))
abline(h=10,col="red",lwd=2,lty="dashed")
samps = matrix(sample(c(1:50),size = length(rank_nu[!is.na(rank_nu[,j]),i])*500, replace=T),500,length(rank_nu[!is.na(rank_nu[,j]),j]))
ppc_ecdf_overlay_2(rank_nu[!is.na(rank_nu[,j]),j],samps) + geom_abline(slope=1/50,intercept=0,colour="grey45",linetype="dashed")

rank_phi<-read.table(paste0(data.path,"/rank_phi.txt"))
hist(rank_phi[!is.na(rank_phi[,j]),j],xlab=expression(rank(Phi[j])))
abline(h=10,col="red",lwd=2,lty="dashed")
samps = matrix(sample(c(1:50),size = length(rank_phi[!is.na(rank_phi[,j]),i])*500, replace=T),500,length(rank_phi[!is.na(rank_phi[,j]),j]))
ppc_ecdf_overlay_2(rank_phi[!is.na(rank_phi[,j]),j],samps) + geom_abline(slope=1/50,intercept=0,colour="grey45",linetype="dashed")

rank_s<-read.table(paste0(data.path,"/rank_s.txt"))
rank_s<-rank_s[1:47,]
hist(rank_s[!is.na(rank_s[,j]),j],xlab=expression(rank(s[j])))
abline(h=10,col="red",lwd=2,lty="dashed")
samps = matrix(sample(c(0:50),size = length(rank_s[!is.na(rank_s[,j]),i])*500, replace=T),500,length(rank_s[!is.na(rank_s[,j]),j]))
ppc_ecdf_overlay_2(rank_s[!is.na(rank_s[,j]),j],samps) + geom_abline(slope=1/50,intercept=0,colour="grey45",linetype="dashed")

rank_theta<-read.table(paste0(data.path,"/rank_theta.txt"))
rank_theta<-as.matrix(rank_theta)
rank_theta<-rank_theta[1:47,]
hist(rank_theta[!is.na(rank_theta)],xlab=expression(rank(theta)))
abline(h=10,col="red",lwd=2,lty="dashed")
samps = matrix(sample(c(1:50),size = length(rank_theta)*500, replace=T),500,length(rank_theta))
ppc_ecdf_overlay_2(rank_theta[,1],samps) + geom_abline(slope=1/50,intercept=0,colour="grey45",linetype="dashed")
