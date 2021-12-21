setwd()#set data path to where you want to store your synthetic data from non-regression BASiCS
#args = commandArgs(trailingOnly = TRUE)
#arrayNumber=as.numeric
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
#simulate gene specific expected count mu across cells
log_mu<-rnorm(n_all_gene,mean=0,sd=sqrt(s2.mu))
mu<-exp(log_mu)

#simulate biological variation factor delta
log_delta<-rnorm(n_bio_gene,mean=0,sd=sqrt(s2.delta))
delta<-exp(log_delta)

#simulate general technical noise factor theta
theta<-rgamma(1,a.theta,b.theta)
#theta<-theta/snr
#noise<-rnorm()
#X<-X+noise?
#snr=1,2,3,4,5,6,7
###
#snr=expectation/noise
#theta is the global technical noise
#mu is expected expression
#theta<-mu/snr

#simulate cell specific technical noise s
s<-rgamma(n_cell,a.s,b.s)

#simulate cell size factor phi
phi<-as.vector(n_cell*rdirichlet(alpha=p,n=1))

#simulate nu_j~Gamma(1/theta,1/(s_j*theta)),j=1,...,n_cell
nu<-rep(0,n_cell)
for (j in (1:n_cell)){
  nu[j]<-rgamma(1,shape=1/theta,rate=(1/(theta*s[j])))
}

#simulate the gene expression count matrix of all genes
#X is a matrix with n_all_gene rows and n_cell columns
X<-matrix(0,ncol=n_cell,nrow=n_all_gene)
mean<-matrix(0,ncol=n_cell,nrow=n_bio_gene)

#simulate the gene expression count of biological genes
#X_{ij}~Neg-binomial(1/delta_i,phi_j*nu_j*mu_i/(phi_j*nu_j*mu_i+1/delta_i))
#the expression count of biological gene i in cell j, i=1,...,n_bio_gene, j=1,
#for (j in (1:n_cell)){
#  for (i in (1:n_bio_gene)){
#    mean[i,j]<-(phi[j,1]*nu[j]*mu[i])/(phi[j,1]*nu[j]*mu[i]+1/delta[i])
#    X[i,j]<-rnbinom(1,size=1/delta[i],mu=mean)
#  }
#}

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
lambda<-matrix(0,nrow=n_bio_gene,ncol=n_cell)
for (j in (1:n_cell)){
  for (i in (1:n_bio_gene)){
    lambda[i,j]<-mu[i]*nu[j]*phi[j]*rho[i,j]
    X[i,j]<-rpois(1,lambda[i,j])
  }
}

#simulate the gene expression count for spike-in genes
#X_{ij}~Poisson(mu_i*nu_j)
#the expression count of spike-in gene i in cell j, 
#i=n_bio_gene+1,...n_all_gene, j=1,...,n_cell
for (j in (1:n_cell)){
  for (i in ((n_bio_gene+1):n_all_gene)){
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
mu_spike<-mu[(n_bio_gene+1):n_all_gene]
spike_input1<-mu_spike[which(spike_name %in% spike_name1)]
SpikesInfo1<-data.frame("SpikeID"=spike_name1,"SpikeInput"=spike_input1)
write.table(Tech,paste0("simulate_Tech.txt"))
write.table(SpikesInfo1,paste0("simulate_SpikesInfo.txt"))

