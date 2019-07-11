\name{simData}
\alias{simData}
\title{simData}
\description{
Simulate single-cell genomic data by model-based approach, including single-cell chromatin accessibility and single-cell gene expression data for 2 clusters.
}
\usage{
simData(n1, n2, p, Xprobs, muG1, sigmaG1, muG2, sigmaG2, mu0, mu1, nu0, nu1, q_acc, q_exp, pi_exp, overlap_prop = 0.5, diff_prop, diff_prop_atac, diff_prop_rna, cutoff = 10^-6, high1 = 0.9, low1 = 0.1, high2 = 0.9, low2 = 0.1)
}

\arguments{
  \item{n1}{sample size of cluster 1}

  \item{n2}{sample size of cluster 2}

  \item{p}{number of features}

  \item{Xprobs}{numeric non-negative vector of length 2, specifying the probability for the 2 classes; is internally normalized to sum 1. Infinite and missing values are not allowed}

  \item{muG1}{numeric vector of length 2, specifying the mean of fitting normal mixture for gene chromatin accessibility data}

  \item{sigmaG1}{numeric vector of length 2, specifying the standard deviation of fitting normal mixture for gene chromatin accessibility data}

  \item{muG2}{numeric vector of length 2, specifying the mean of fitting normal mixture for gene expression data}

  \item{sigmaG2}{numeric vector of length 2, specifying the standard deviation of fitting normal mixture for gene expression data}

  \item{mu0}{the mean value of beta distribution for omega_acc}

  \item{mu1}{the mean value of beta distribution for omega_exp from genes that are not linked to the accessibility data}

  \item{nu0}{the precision parameter of beta distribution for omega_acc and omega_exp from genes that are not linked to the accessibility data}

  \item{nu1}{the precision parameter of beta distribution for omega_exp from genes that are linked to the accessibility data}

  \item{q_acc}{probability for a Bernoulli distribution which is representing the capture rate of regulatory element opening status signal in the experimental process}

  \item{q_exp}{probability for a Bernoulli distribution which is representing the capture rate of gene expressing status signal in the experiment process}

  \item{pi_exp}{probability of the gene being expressed in cell i}

  \item{overlap_prop}{the proportion of the linked features}

  \item{diff_prop}{the proportion of different features among the linked features}

  \item{diff_prop_atac}{the proportion of different features among the linked features in chromatin accessibility data}

  \item{diff_prop_rna}{the proportion of different features among the linked features in gene expression data}

  \item{cutoff}{the cutoff value for omega_acc and omega_exp}

  \item{high1}{value 1 for omega_acc for different features acorss the clusters}

  \item{low1}{value 2 for omega_acc for different features acorss the clusters}

  \item{high2}{value 1 for omega_exp for non-overlap different features acorss the clusters}

  \item{low2}{value 2 for omega_exp for non-overlap different features acorss the clusters}

}


\value{
  \item{f1}{normal mixture component 1 for the chromatin accessibility data}
  \item{f0}{normal mixture component 2 for the chromatin accessibility data}

  \item{g1}{normal mixture component 1 for the gene expression data}
  \item{g0}{normal mixture component 2 for the gene expression data}

  \item{cluster_atac}{a vector of cluster assignment for the chromatin accessibility data}
  \item{cluster_rna}{a vector of cluster assignment for of the gene expression data}

  \item{x_atac}{cluster assignment z_ik, whether the cell in the chromatin accessibility data belongs to cluster or not}
  \item{x_rna}{cluster assignment z_lk, whether the cell in the gene expression data belongs to cluster or not}

  \item{u_atac}{true status of a regulatory element for the chromatin accessibility data}
  \item{u_rna}{true status of a regulatory element for the linked data, and true gene status for the unlinked gene expression data}

  \item{v_atac}{the "contaminated" status of a regulatory element in the single-cell experiment due to the loss of biological signal in the experiment process}
  \item{v_rna}{the "contaminated" gene status due to the loss of biological signal in the experiment process}

  \item{Data_atac}{simulated chromatin accessibility data}
  \item{Data_rna}{simulated gene expression data}

  \item{pi_pk_atac}{probability of regulatory region r to be open in cluster k}
  \item{pi_pk_rna}{probability of gene g to be expressed in cluster k}

  \item{q_acc}{probability for a Bernoulli distribution which is representing the capture rate of regulatory element opening status signal in the experimental process}
  \item{q_exp}{probability for a Bernoulli distribution which is representing the capture rate of gene expressing status signal in the experiment process}

  \item{pi_exp}{probability of the gene being expressed in cell i}

}
\references{
Lin Z, Zamanighomi M, Daley T, Ma S and Wong WH(2019). "Model-based approach to the joint analysis of single-cell data on chromatin accessibility and gene expression." Statistical Science.}
\author{
Zhixiang LIN, Jing Huang
}


\examples{
library(scACE)

n1 <- 100
n2 <- 100
p <- 100
Xprobs <- c(0.5, 0.5)
muG1 <- c(0, 2)
sigmaG1 <- c(0.8, 0.8)
muG2 <- c(0, 2)
sigmaG2 <- c(0.8, 0.8)
overlap_prop <- 0.8
diff_prop <- 0.2
diff_prop_atac <- 0.001
diff_prop_rna <- 0.001

mu0 <- 0.5
mu1 <- 0.5
nu0 <- 2
nu1 <- 10
q_acc <- rep(0.5, n1)
q_exp <- rep(0.5, n2)
pi_exp <- rbind(rep(0.7, n2), 0.1)


# Simulation 1
simdata1 <- simData(n1=n1, n2=n2, p=p, Xprobs=Xprobs,
                    muG1=muG1, sigmaG1=sigmaG1, muG2=muG2, sigmaG2=sigmaG2,
                    mu0=mu0, mu1=mu1, nu0=nu0, nu1=nu1,
                    q_acc=q_acc, q_exp=q_exp, pi_exp=pi_exp,
                    overlap_prop=overlap_prop, diff_prop=diff_prop,
                    diff_prop_atac=diff_prop_atac, diff_prop_rna=diff_prop_rna,
                    cutoff=10^-6,
                    high1=0.9, low1=0.1, high2=0.9, low2=0.1)
