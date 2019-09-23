\name{update_all2}
\alias{update_all2}
\title{update_all2}
\description{
A model-based approach with expectation-maximization (EM) implementation that is specifically designed for single-cell genomic data and can jointly cluster single-cell chromatin accessibility and single-cell gene expression data.
}

\usage{
update_all2(phi_atac,f1,f0,w_acc,phi_rna,g1,g0,w_exp,pi_exp,qi,ql,phi_1,po)
}

\arguments{
  \item{phi_atac}{the initialization of the probability for accessibility data belonging to cluster k (i.e. z[i,k]=1)}

  \item{f1}{f1[i,r] the known pdf of accessibility data x[i,r] when ut[i,r]=1}

  \item{f0}{f0[i,r] the known pdf of accessibility data x[i,r] when ut[i,r]=0}

  \item{w_acc}{the initialization of the probability of regulatory region r to be open in cluster k}

  \item{phi_rna}{the initialization of the probability for gene expression data belonging to cluster k (i.e. z[l,k]=1)}

  \item{g1}{g1[l,g] the known pdf of gene expression data y[l,g] when vt[l,g]=1}

  \item{g0}{g1[l,g] the known pdf of gene expression data y[l,g] when vt[l,g]=0}


  \item{w_exp}{the initialization of the probability of gene g to be expressed in cluster k}

  \item{pi_exp}{the initialization of the probability of the gene being expressed in cell i}

  \item{qi}{the initialization of the parameter of probability for a Bernoulli distribution which is representing the capture rate of regulatory element opening status signal in the experimental process}

  \item{ql}{the initialization of the parameter of probability for a Bernoulli distribution which is representing the capture rate of gene expressing status signal in the experiment process}

  \item{pi_exp}{the initialization of the probability of the gene being expressed in cell i}

  \item{phi_1}{the precision parameter of beta distribution, which represents how well the two data types are coupled}

  \item{po}{the number of linked features, assuming 1:po is linked}

}


\value{
phi_atac: the probability for accessibility data belonging to cluster k

phi_rna: the probability for gene expression data belonging to cluster k

w_exp: the probability of gene g to be expressed in cluster k

w_acc: the probability of regulatory region r to be open in cluster k

pi_exp: the probability of the gene being expressed in cell i

qi: the parameter of probability for a Bernoulli distribution which is representing the capture rate of regulatory element opening status signal

ql: the parameter of probability for a Bernoulli distribution which is representing the capture rate of gene expressing status signal

phi_1: the precision parameter for beta distribution of w_exp linked to w_acc

postprob: posterior probability of all the parameters
}

\author{
Jingwen Yang, Jiaxuan WANGWU, Jingyu Li
}

\examples{
library(mixtools)
library(scACE)

k0 <- 2
i0 <- nrow(LZX_data$data_acc)
l0 <- nrow(LZX_data$data_exp)
p0 <- ncol(LZX_data$data_exp)
po <- 500

# f1 and f0
f1 <- (LZX_data$data_acc!=0) + 0
f0 <- 1 - f1

# g1 and g0
data_exp_non0 <- data_exp[which(LZX_data$data_exp!=0)]
if (length(data_exp_non0)>50000){
  data_exp_non0 <- sample(data_exp_non0, 50000)
}
data_exp_non0 <- log(data_exp_non0 + 1)
fit <- gammamixEM(x=data_exp_non0, k = 2, epsilon=2, verb=F)
data_exp <- log(LZX_data$data_exp+1)
dens_scRNA_gg <- cal_density_gg_thres(gg_fit=fit, data=data_exp, thres_low=-1, thres_high=20)
g1 <- dens_scRNA_gg$f1
g0 <- dens_scRNA_gg$f0


# Random initialization
set.seed(1)
phi_atac <- c(1/k0,1/k0)
phi_rna <- c(1/k0,1/k0)
phi_1 <- 10

w_acc <- (rep(0.5,k0*p0),k0,p0)
w_exp <- matrix(runif(k0*p0,min=0.2,max=0.8),k0,p0)
pi_exp <- matrix(c(rep(0.6,l0),rep(0.1,l0)),k0,l0,byrow=TRUE)
qi<-rep(0.5,i0)
ql<-rep(0.5,l0)


niter <- 500
temp <- update_all2(phi_atac=phi_atac,f1=f1,f0=f0,w_acc=w_acc,
                    phi_rna=phi_rna,g1=g1,g0=g0,w_exp=w_exp,pi_exp=pi_exp,
                    qi=qi,ql=ql,phi_1=phi_1,po=po)
initialpostprob <- temp$postprob

postprob <- c()
for(i in c(1:niter)){
  temp<-update_all2(phi_atac=temp$phi_atac,f1=f1,f0=f0,w_acc=temp$w_acc,
                    phi_rna=temp$phi_rna,g1=g1,g0=g0,w_exp=temp$w_exp,pi_exp=temp$pi_exp,
                    qi=temp$qi,ql=temp$ql,phi_1=temp$phi_1,po=po)
  postprob <- c(postprob,temp$postprob)
  print(c(i,temp$postprob))
}
}
