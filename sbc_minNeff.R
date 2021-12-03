#g1d<-function(delta,mu,nu,phi,s,theta){
 # g_value<-sum(delta)+sum(mu)+theta+sum(nu*phi)+sum(s)
  #return(g_value)
#}
#Compute g for the prior value
#g0<-g1d(delta1,mu1,nu1,phi1,s1,theta1)
#Compute g for the posterior samples
#g_samp<-rep(0,500)
library(LaplacesDemon)
library(stats)
library(BASiCS)
library(SingleCellExperiment)
setwd("C:/Users/hksca/OneDrive/Bayesian Big Data/BASiCS and Beyond/BASiCS Robustness/HPC BASiCS/2021/Simulation Experiment/Too many Simulations/Local Running")
data.path<-getwd()
s2.mu<-0.8
s2.delta<-0.5 #1 3 5 #biological variation #
a.theta<-1 #2 3 4 5 #start from here #
b.theta<-1 
#theta global technical noise Gamma(a,b)
#a/b^2
a.s<-1
b.s<-1
n_run<-100
for (index_run in (70:n_run)){
sim_Output_delta<-read.table(paste0(data.path,"/sim_delta_",index_run,".txt"))
#sim_Output_delta<-as.matrix(sim_Output_delta)
sim_Output_mu<-read.table(paste0(data.path,"/sim_mu_",index_run,".txt"))
sim_Output_nu<-read.table(paste0(data.path,"/sim_nu_",index_run,".txt"))
sim_Output_phi<-read.table(paste0(data.path,"/sim_phi_",index_run,".txt"))
sim_Output_s<-read.table(paste0(data.path,"/sim_s_",index_run,".txt"))
sim_Output_theta<-read.table(paste0(data.path,"/sim_theta_",index_run,".txt"))
#for (samp in (1:500)){
 
 # post_delta<-sim_Output_delta[samp,]
  #post_mu<-sim_Output_mu[samp,]
  #post_nu<-sim_Output_nu[samp,]
  #post_phi<-sim_Output_phi[samp,]
  #post_s<-sim_Output_s[samp,]
  #post_theta<-sim_Output_theta[samp,]
 # g_samp[samp]<-g1d(post_delta,post_mu,post_nu,post_phi,post_s,post_theta)
#}

#compute log-m autocorrelation
#autoCorr<-rep(0,499)
#for (lag in (0:499)){
# autoCorr[lag]<-0
#for (k in (1:(500-lag))){
# autoCorr[lag]<-autoCorr[lag]+g_samp[k]*g_samp[k+lag]
#}
#autoCorr[lag]<-autoCorr[lag]/(500-lag)
#}

#Compute the number of effective samples
#N_eff=500/(1+sum(autoCorr))
Ess<-c(ESS(sim_Output_delta),ESS(sim_Output_mu),ESS(sim_Output_nu),ESS(sim_Output_phi),ESS(sim_Output_s),ESS(sim_Output_theta))
write.table(Ess,paste0("ESS_",index_run,".txt"))
N_eff=min(Ess)
#Define the function for rank statistic calculation
#rank_statistic<-function(g0,g_samp){
#  rank<-0
 # for (l in (1:L)){
  #  if (g_samp[l]<g0){
   #   rank=rank+1
    #}
    #else{
     # rank=rank
    #}
#  }
 # return(rank)
#}

#Set L=N_samp/10=50, 
L=50
#then
if (N_eff>=L){
  #Thin the sample to L samples
  #thin_seq<-seq(1,length(g_samp),10)
  #g_samp_new<-g_samp[thin_seq]
  #Compute the rank statistic
  #r<-rank_statistic(g0,g_samp_new) 
  print(index_run)
  #write.table(g_samp,paste0("sim_g_samp_",index_run,".txt"))
  #write.table(g0,paste0("sim_g0_",index_run,".txt"))
  #write.table(r,paste0("sim_rank_",index_run,".txt"))
}else {
  N_new=2*round((50*15000/N_eff))+10000
  n_new=10000
  X1<-read.table(paste0(data.path,"/simulate_",index_run,"_X1.txt"))
  X1<-as.matrix(X1)
  #Produce the vector indicating whether the gene is a spike-in gene
  Tech<-read.table(paste0(data.path,"/simulate_",index_run,"_Tech.txt"))
  Tech<-as.matrix(Tech)
  SpikesInfo<-read.table(paste0(data.path,"/simulate_",index_run,"_SpikesInfo.txt"))
  spike_input1<-SpikesInfo$SpikeInput
  #Produce the SingleCellExperiment object to use in BASiCS
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = X1[!Tech,])
  )
  spike_sce<- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = X1[Tech,])
  )
  altExp(sce, "spike-ins") <- spike_sce
  rowData(altExp(sce, "spike-ins")) <- data.frame(Name = rownames(X1[Tech,]), Molecules = spike_input1, row.names = rownames(X1[Tech,]))
  Data<-sce
  prior.param<-BASiCS_PriorParam(Data,k = 12,
                                 s2.mu = s2.mu,
                                 s2.delta = s2.delta,
                                 a.delta = 1,
                                 b.delta = 1,
                                 a.s = a.s,
                                 b.s = b.s,
                                 a.theta = a.theta,
                                 b.theta = b.theta,
                                 a.sigma2 = 2,
                                 b.sigma2 = 2,
                                 PriorDelta = 0)
  sim_Output<-BASiCS_MCMC(Data=Data,N=N_new,Thin=2,Burn=10000,
                          Regression=FALSE,WithSpikes = TRUE,
                          PriorParam = prior.param)
  sim_Output_mu<-sim_Output@parameters$mu
  write.table(sim_Output_mu,paste0("sim_mu_",index_run, ".txt"))
  sim_Output_delta<-sim_Output@parameters$delta
  write.table(sim_Output_delta,paste0("sim_delta_",index_run, ".txt"))
  sim_Output_theta<-sim_Output@parameters$theta
  write.table(sim_Output_theta,paste0("sim_theta_",index_run, ".txt"))
  sim_Output_s<-sim_Output@parameters$s
  write.table(sim_Output_s,paste0("sim_s_",index_run, ".txt"))
  sim_Output_phi<-sim_Output@parameters$phi
  write.table(sim_Output_phi,paste0("sim_phi_",index_run, ".txt"))
  sim_Output_nu<-sim_Output@parameters$nu
  write.table(sim_Output_nu,paste0("sim_nu_",index_run, ".txt"))
  #Compute g for the prior value
  #g0<-g1d(delta1,mu1,nu1,phi1,s1,theta1)
  #Compute g for the posterior samples
  #g_samp<-rep(0,L)
  #for (samp in (1:L)){
   # post_delta<-sim_Output_delta[samp,]
    #post_mu<-sim_Output_mu[samp,]
    #post_nu<-sim_Output_nu[samp,]
    #post_phi<-sim_Output_phi[samp,]
    #post_s<-sim_Output_s[samp,]
    #post_theta<-sim_Output_theta[samp,]
    #g_samp[samp]<-g1d(post_delta,post_mu,post_nu,post_phi,post_s,post_theta)
  #}
  #Thin the sample to L samples
  #thin_seq<-seq(1,length(g_samp),length(g_samp)/L)
  
  #g_samp_new<-g_samp[thin_seq]
  #r<-rank_statistic(g0,g_samp_new)
 # write.table(r,paste0("sim_rank_",index_run,".txt"))
 # write.table(g_samp,paste0("sim_g_samp_",index_run,".txt"))
 # write.table(g0,paste0("sim_g0_",index_run,".txt"))
 # write.table(N_eff,paste0("N_eff_",index_run,".txt"))
}
}




