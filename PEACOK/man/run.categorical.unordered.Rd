\name{run.categorical.unordered}
\alias{run.categorical.unordered}

\title{Fitting multinomial log-linear models with confounders.}
\description{
Given input phenotype and genotype vectors as well as coresponding covariates this function fits a multinomial log-linear model, assuming that the input phenotypes have more than two levels with unordered data.
}
\usage{
run.categorical.unordered(pheno, geno, confounders)
}
\arguments{
  \item{pheno}{a vector with more than two distinct values as unordered categorical values.}
  \item{geno}{a numeric vector, typically binary vector in this package.}
  \item{confounders}{A dataframe with rows for samples and columns for confounder covariates. }
}
