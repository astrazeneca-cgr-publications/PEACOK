\name{run.logistic2}
\alias{run.logistic2}
\title{
Running logistic regression with confounders too.
}
\description{
Giving input phenotype and genotype vectors as well as coresponding covariates this function fits a logistic regression model phenotype ~ genotype + confounders and returns beta, standard error and p value associated with the genotype. This implementaiton is almost identical to run.logistic but with less output.
}
\usage{
run.logistic2(pheno, geno, confounders)
}
\arguments{
   \item{pheno}{a binary vector with 0 for cases, 1 for controls and NA for misisng values.}
  \item{geno}{a numeric vector, typically binary vector in this package.}
  \item{confounders}{A dataframe with rows for samples and columns for confounder covariates. }
}
\details{
This function simply fits a logistic regresion phenotype ~ genotype + confounders and return test statistics associated with the genotype vector. It assumes that the orders of samples from all input are synced.
}
\value{
 \item{nSamples }{number of samples that have non misisng phenotype and genotype values}
  \item{numCases }{number of samples that have non misisng phenotype and genotype values and have non zero genotypes. Note that in the case of input genotype is nonbinary, the definition here might not be in line with typical difinition for cases.}
  \item{numControls }{number of samples that have non misisng phenotype and genotype values and zhave ero genotypes.}
  \item{p }{The p value for testing the coefficient for genotype is not zero.}
  \item{beta }{the beta coefficient for genotype.}
  \item{se }{the standarr error for the beta coefficient for genotype.}
}

