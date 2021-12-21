library(LaplacesDemon)
library(stats)
library(BASiCS)
library(SingleCellExperiment)
setwd("./regression_BASiCS/resimulate")
data.path<-getwd()
s2.mu<-0.8
s2.delta<-0.5 #1 3 5 #biological variation #
a.theta<-1 #2 3 4 5 #start from here #
b.theta<-1 
#theta global technical noise Gamma(a,b)
a.s<-1
b.s<-1
n_run<-100
for (index_run in (2:n_run)){
sim_Output_delta<-read.table(paste0(data.path,"/sim_delta_",index_run,".txt"))
#sim_Output_delta<-as.matrix(sim_Output_delta)
sim_Output_mu<-read.table(paste0(data.path,"/sim_mu_",index_run,".txt"))
sim_Output_nu<-read.table(paste0(data.path,"/sim_nu_",index_run,".txt"))
sim_Output_phi<-read.table(paste0(data.path,"/sim_phi_",index_run,".txt"))
sim_Output_s<-read.table(paste0(data.path,"/sim_s_",index_run,".txt"))
sim_Output_theta<-read.table(paste0(data.path,"/sim_theta_",index_run,".txt"))
Ess<-c(ESS(sim_Output_delta[,colSums(is.na(sim_Output_delta)) == 0]),ESS(sim_Output_mu[,colSums(is.na(sim_Output_mu)) == 0]),ESS(sim_Output_nu[,colSums(is.na(sim_Output_nu)) == 0]),ESS(sim_Output_phi[,colSums(is.na(sim_Output_phi)) == 0]),ESS(sim_Output_s[,colSums(is.na(sim_Output_s)) == 0]),ESS(sim_Output_theta[,colSums(is.na(sim_Output_theta)) == 0]))
write.table(Ess,paste0("ESS_",index_run,".txt"))
N_eff=min(Ess)
#Set L=N_samp/10=50, 
L=50
#then
if (N_eff>=L){
  print(index_run)
}else {
  N_new=2*(round(50*15000/N_eff))+10000
  X1<-read.table(paste0(data.path,"/simulate_",index_run,"_X1.txt"))
  X1<-as.matrix(X1)
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
                                 b.sigma2 = 2)
  sim_Output<-BASiCS_MCMC(Data=Data,N=N_new,Thin=2,Burn=10000,
                          Regression=TRUE,WithSpikes = TRUE,
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
}
}

