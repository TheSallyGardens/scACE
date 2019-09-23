\name{getClusterGibbs}
\alias{getClusterGibbs}
\title{getClusterGibbs}
\description{
A model-based approach with Markov Chain Monte Carlo (MCMC) implementation that is specifically designed for single-cell genomic data and can jointly cluster single-cell chromatin accessibility and single-cell gene expression data.
}
\usage{
getClusterGibbs(data_acc, data_exp, overlap_seq_acc, overlap_seq_exp, mu0 = 0.5, mu1 = 0.5, nu0 = 2, nCluster, cluster_acc_Ini = NULL, cluster_exp_Ini = NULL, u_acc_tilda_Ini = NULL, v_exp_tilda_Ini = NULL, omega_acc_Ini = NULL, omega_exp_Ini = NULL, nu1_Ini = 1, q_acc_Ini = NULL, q_exp_Ini = NULL, pi_exp_Ini = NULL, niter, quiet = F, eval_logPost = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data_acc}{chromatin accessibility data matrix, cells by regions}

  \item{data_exp}{gene expression data matrix, cells by genes}

  \item{overlap_seq_acc}{the linked features in chromatin accessibility data}

  \item{overlap_seq_exp}{the linked features in gene expression data}

  \item{mu0}{the mean value of beta distribution for omega_acc, which is set as 0.5}

  \item{mu1}{the mean value of beta distribution for omega_exp from genes that are not linked to the accessibility data, which is set as 0.5}

  \item{nu0}{the precision parameter of beta distribution for omega_acc, which represents how well the chromatin accessibility data and gene expression data are coupled and is set as 2}

  \item{nCluster}{the number of clustering assignments}

  \item{cluster_acc_Ini}{the initialization of clustering accessibility data, which is optional}

  \item{cluster_exp_Ini}{the initialization of clustering expression data, which is optional}

  \item{u_acc_tilda_Ini}{the initialization of "contaminated" status of a regulatory element in the single-cell experiment due to the loss of biological signal in the experiment process, which is optional}

  \item{v_exp_tilda_Ini}{the initialization of "contaminated" gene status due to the loss of biological signal in the experiment process, which is optional}

  \item{omega_acc_Ini}{the initialization of the probability of regulatory region r to be open in cluster k, which is optional}

  \item{omega_exp_Ini}{the initialization of the probability of gene g which isn't linked to accessibility data to be expressed in cluster k, which is optional}

  \item{nu1_Ini}{the initialization of the precision parameter of beta distribution, which represents how well the two data types are coupled and is set as 1}

  \item{q_acc_Ini}{the initialization of the parameter of probability for a Bernoulli distribution which is representing the capture rate of regulatory element opening status signal in the experimental process, which is optinal}

  \item{q_exp_Ini}{the initialization of the parameter of probability for a Bernoulli distribution which is representing the capture rate of gene expressing status signal in the experiment process, which is optional}

  \item{pi_exp_Ini}{the initialization of the probability of the gene being expressed in cell i, which is optional}

  \item{niter}{the nummber of iteration}

  \item{quiet}{when it is set as 'TRUE', it will report the number of iteration; otherwise it won't}

  \item{eval_logPost}{when it is set as 'TRUE', it will calculate log posterior probability; otherwise it won't}
}


\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
z_acc_prob: a vector of marginal cluster assignment probability for single-cell chromatin accessibility data

z_exp_prob: a vector of marginal cluster assignment probability for single-cell gene expression data

cluster_acc: a vector of cluster assignment (based on z_acc_prob) for single-cell chromatin accessibility data

cluster_exp: a vector of cluster assignment (based on z_exp_prob) for single-cell gene expression data

gamma_mixture_exp: the estimated gamma-gamma mixture parameters for single-cell gene expression data

z_acc_Trace value depends on whether the cell in the accessibility data belongs to cluster or not in every iteration

q_acc_Trace the parameter of probability for a Bernoulli distribution which is representing the capture rate of regulatory element opening status signal in the experimental process in every iteration

omega_acc_Trace the probability of regulatory region r to be open in the cluster in every iteration

z_exp_Trace value depends on whether the cell in the gene expression data belongs to cluster or not in every iteration

q_exp_Trace the parameter of probability for a Bernoulli distribution which is representing the capture rate of gene expressing status signal in the experiment process in every iteration

omega_exp_Trace the probability of gene g which isn't linked to accessibility data to be expressed in cluster k in every iteration

pi_exp_Trace the probability of the gene being expressed in cell i in every iteration

h_Trace auxiliary variable in every iteration, representing how gene expression clusters are matched to the accessibility clusters

logpost_Trace whether or not calculate posterior probability in every iteration

nu1_Trace the precision parameter of beta distribution, which represents how well the two data types are coupled in every iteration
}


\author{
Zhixiang LIN, Jing Huang
}

\examples{
library(mixtools)
library(scACE)
d <- scACE::LZX_data
result_test <- getClusterGibbs(data_acc=d$data_acc ,data_exp=d$data_exp ,overlap_seq_acc=d$overlap_seq_acc ,overlap_seq_exp=d$overlap_seq_exp ,nCluster=2 ,niter=1000)
}
