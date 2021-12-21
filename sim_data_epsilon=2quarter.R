setwd("./non_regression_BASiCS/fixed_dataset/epsilon=0.5")
#set data path t where you store the synthetic data from non-regression BASiCS
#and where you would like to store your posterior 
#from modified BASiCS MCMC with epsilon=0.5
library(LaplacesDemon)
library(stats)
library(BASiCS)
library(SingleCellExperiment)
#args = commandArgs(trailingOnly = TRUE)
#arrayNumber=as.numeric
epsilon=0.5
args<-as.numeric(args)
n_bio_gene<-100
n_spike_gene<-10
n_cell<-50
n_all_gene<-n_bio_gene+n_spike_gene
s2.mu<-0.8
s2.delta<-0.5 #1 3 5 #biological variation #
a.theta<-1 #2 3 4 5 #start from here #
b.theta<-1 
#theta global technical noise Gamma(a,b)
#a/b^2
a.s<-1
b.s<-1

p<-rep(1,n_cell)

X1<-read.table("simulate_X1.txt")
X1<-as.matrix(X1)
#Produce the vector indicating whether the gene is a spike-in gene
Tech<-read.table("simulate_Tech.txt")
Tech<-as.matrix(Tech)
SpikesInfo<-read.table("simulate_SpikesInfo.txt")
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

#Assume the same prior parameter
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
                               PriorDelta = epsilon)
simBASiCS<-function(args){
  sim_Output<-BASiCS_MCMC(Data=Data,N=15000,Thin=5,Burn=10000,
                          Regression=FALSE,WithSpikes = TRUE,
                          PriorParam = prior.param)
  BASiCS_DiagHist(sim_Output)
  BASiCS_DiagPlot(sim_Output)
  sim_Output_mu<-sim_Output@parameters$mu
  write.table(sim_Output_mu,paste0("sim_mu_",args, ".txt"))
  sim_Output_delta<-sim_Output@parameters$delta
  write.table(sim_Output_delta,paste0("sim_delta_",args, ".txt"))
  sim_Output_theta<-sim_Output@parameters$theta
  write.table(sim_Output_theta,paste0("sim_theta_",args, ".txt"))
  sim_Output_s<-sim_Output@parameters$s
  write.table(sim_Output_s,paste0("sim_s_",args, ".txt"))
  sim_Output_phi<-sim_Output@parameters$phi
  write.table(sim_Output_phi,paste0("sim_phi_",args, ".txt"))
  sim_Output_nu<-sim_Output@parameters$nu
  write.table(sim_Output_nu,paste0("sim_nu_",args, ".txt"))
}
for (args in (1:200)){
  simBASiCS(args)
}

