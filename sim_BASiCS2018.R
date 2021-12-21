library(LaplacesDemon)
library(MASS)
library(stats)
setwd()#set data path to where you choose to store the simulated data
k<-12
eta<-5
n_bio_gene<-100
n_spike_gene<-10
n_cell<-50
n_gene<-n_bio_gene+n_spike_gene
s2.delta<-0.5
s2.mu<-0.8
a.s=1
b.s=1
a.theta<-1 
b.theta<-1
p<-rep(1,n_cell)
a.sigma2<-2
b.sigma2<-2
logmu.mean=rep(0,n_gene)
#equation 40
logmu.Sigma=s2.mu*(diag(n_gene)-rep(1,n_gene)%*%t(rep(1,n_gene))/n_gene)
#equation 3
logmu<-mvrnorm(n = 1, mu=logmu.mean, Sigma=logmu.Sigma, tol = 1e-6, empirical = FALSE)
mu<-exp(logmu)
sigma2<-rinvgamma(n=1, shape=a.sigma2, scale=b.sigma2)
m0.beta<-numeric(k)
V0.beta<-diag(k)
#equation 7
beta<-mvrnorm(n=1,mu=m0.beta,Sigma=sigma2*V0.beta)
#L info given by Alan
L=k-2
logmu.bio<-logmu[1:n_bio_gene]
logmu.max<-max(logmu.bio)
logmu.min<-min(logmu.bio)
ran<-logmu.max-logmu.min
d=ran/(k-1)
RBFlocations<-seq((from=logmu.min+ran/(k-3)),to=(logmu.max+ran/(k-3)),by=ran/(k-3))
#from line 74 of BASICS_PriorParam.R, variance=1.2
variance=1.2
#h in equation 6
#from line 19 of updatesReg.h
h<-(ran/(k-3))*variance
#from line 22-25 of updatesReg.h
#XX the design matrix
#use notation XX instead of X to avoid confusion with the count matrix
XX<-matrix(0,nrow=n_bio_gene,ncol=k)

XX[,1]<-logmu.bio
XX[,2]<-exp(-0.5*(logmu.bio-logmu.min)^2)
for (i in (3:k)){
  XX[,i]<-exp(-0.5*(logmu.bio-RBFlocations[i-2])^2/(h^2))
}
#epsilon<-rt(n=n_bio_gene,df=eta)
#logdelta<-XX%*%beta+sqrt(s2.delta)*epsilon
#equation 16
lambda<-rgamma(n=n_bio_gene,shape=eta/2,rate=eta/2)
logdelta<-rep(0,n_bio_gene)
mean.logdelta<-XX%*%beta

v.rnorm<-Vectorize(rnorm)
logdelta<-v.rnorm(1,mean=mean.logdelta,sd=sqrt(s2.delta/lambda))
delta<-exp(logdelta)
s<-rgamma(n_cell,shape=a.s,rate=b.s)
theta<-rgamma(1,shape=a.theta,rate=b.theta)
phi<-as.vector(n_cell*rdirichlet(alpha=p,n=1))

#simulate the gene expression count matrix of all genes
nu<-rep(0,n_cell)
for (j in (1:n_cell)){
  nu[j]<-rgamma(1,shape=1/theta,rate=(1/(theta*s[j])))
}

#Simulate rho_{ij}~Gamma(1/delta_i,1/delta_i)
rho<-matrix(0,ncol=n_cell,nrow=n_bio_gene)
for (j in (1:n_cell)){
  for (i in (1:n_bio_gene)){
    rho[i,j]<-rgamma(1,shape=1/delta[i],rate=1/delta[i])
  }
}

#simulate the gene expression count of biological genes
#X_{ij}~Poisson(phi_j*nu_j*mu_i*rho_{ij})
#the expression count of biological gene i in cell j, i=1,...,n_bio_gene, j=1,
X<-matrix(0,nrow=n_gene,ncol=n_cell)
pois.Mean<-matrix(0,nrow=n_bio_gene,ncol=n_cell)
for (j in (1:n_cell)){
  for (i in (1:n_bio_gene)){
    pois.Mean[i,j]<-mu[i]*nu[j]*phi[j]*rho[i,j]
    X[i,j]<-rpois(1,pois.Mean[i,j])
  }
}

#simulate the gene expression count for spike-in genes
#X_{ij}~Poisson(mu_i*nu_j)
#the expression count of spike-in gene i in cell j, 
#i=n_bio_gene+1,...n_all_gene, j=1,...,n_cell
for (j in (1:n_cell)){
  for (i in ((n_bio_gene+1):n_gene)){
    X[i,j]<-rpois(1,mu[i]*nu[j])
    #snr<-expectation/noise
    #noise<-X/snr
    #X<-X+noise
  }
}

gene_name<-sprintf("Gene[%s]",seq(1:n_bio_gene))
spike_name<-sprintf("ERCC[%s]",seq(1:n_spike_gene))
cell_name<-sprintf("Cell[%s]",seq(1:n_cell))
rownames(X)<-c(gene_name,spike_name)
colnames(X)<-cell_name
which_bio_gene_count_more_than_0.1_per_cell<-which(rowSums(X[1:n_bio_gene,])>5)
which_gene_count_more_than_0.1_per_cell<-which(rowSums(X)>5)
X1<-X[which_gene_count_more_than_0.1_per_cell,]
Tech1<-ifelse(1:nrow(X1) %in% 
                grep("ERCC",rownames(X1)),T,F)
which_cell_count_more_than_0_bio<-which(colSums(X1[!Tech1,])>0)
which_cell_count_more_than_0_spike<-which(colSums(X1[Tech1,])>0)
which_cell_count_good<-intersect(which_cell_count_more_than_0_bio,which_cell_count_more_than_0_spike)

X1<-X1[,which_cell_count_good]

delta1<-delta[which_bio_gene_count_more_than_0.1_per_cell]
mu1<-mu[which_gene_count_more_than_0.1_per_cell]
nu1<-nu[which_cell_count_good]
phi1<-phi[which_cell_count_good]
s1<-s[which_cell_count_good]
theta1<-theta


write.table(delta1,paste0("simulate_delta1.txt"))
write.table(mu1,paste0("simulate_mu1.txt"))
write.table(nu1,paste0("simulate_nu1.txt"))
write.table(phi1,paste0("simulate_phi1.txt"))
write.table(s1,paste0("simulate_s1.txt"))
write.table(theta1,paste0("simulate_theta1.txt"))
write.table(X1,paste0("simulate_X1.txt"))

#Produce the vector indicating whether the gene is a spike-in gene
Tech<-ifelse(1:nrow(X1) %in% 
               grep("ERCC",rownames(X1)),T,F)
spike_name1<-rownames(X1[Tech,])
mu_spike<-mu[(n_bio_gene+1):n_gene]
spike_input1<-mu_spike[which(spike_name %in% spike_name1)]
SpikesInfo1<-data.frame("SpikeID"=spike_name1,"SpikeInput"=spike_input1)
write.table(Tech,paste0("simulate_Tech.txt"))
write.table(SpikesInfo1,paste0("simulate_SpikesInfo.txt"))

