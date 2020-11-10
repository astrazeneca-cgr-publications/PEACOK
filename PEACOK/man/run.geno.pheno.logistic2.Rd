\name{run.geno.pheno.logistic2}
\alias{run.geno.pheno.logistic2}
\title{
Run association tests for a set of input genotypes and binary phenotypes using logistic regression with fewer output.
}
\description{
Run logistic association tests for a set of input genotypes and binary phenotypes with matching confounders. If an annotation information is also provided, it will be used for added annotation. User can also provide the range/block of genotypes and phenotypes so association tests will be performed for choosen gentype/phenotype blocks.
}
\usage{
run.geno.pheno.logistic2(genotypes, phenotypes, file.confounder, file.annotaiton = NULL, geno_start = NULL, geno_end = NULL, pheno_start = NULL, pheno_end = NULL, ignoreConfounder = FALSE, verbose = 0)
}
\arguments{
\item{genotypes}{The input genotype file.(add format info)}
  \item{phenotypes}{The input phenotype file.}
  \item{file.confounder}{The companion confounder file for the input phenotype file.(add format info)}
  \item{file.annotaiton}{Optional. If provided, the info in annotaion file will be used to annotate the output.}
  \item{geno_start}{Optional. Default to the 1.See argument geno_end}
  \item{geno_end}{Optional. Default to the maxinum number of genotypes. (geno_start, geno_end), if provided, together will define the range of genotypes that will be used to run the association tests.}
  \item{pheno_start}{Optional. Default to the 1.See argument pheno_end}
  \item{pheno_end}{Optional. Default to the maxinum number of phenotypes. (pheno_start, pheno_end), if provided, together will define the range of phenotypes that will be used to run the association tests.}
  \item{ignoreConfounder}{Optional, default to FALSE. If setting to TRUE, the confounders will not be used in analysis even if it was provided.}
  \item{verbose}{Optional Default to 0. It currently takes values 0, 1 or 2 and the default value is 0. If the value is 0, no progress information will be printed. If the value is 1, it will print progress at the genotype levele. If the value is set to 2, i will print progress at the (genotype, phenotype) level.}
}
