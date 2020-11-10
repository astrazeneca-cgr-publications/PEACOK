\name{test.categorical.unordered}
\alias{test.categorical.unordered}
\title{
Perform unordered categorical phenotype with multinomial regression
}
\description{
Perform unordered categorical phenotype with multinomial regression and 
saves this result in the multinomial logistic results file
}

\usage{
test.categorical.unordered(opt, vl, varName, varType, thisdata, phenoStartIdx)
}
\arguments{
  \item{opt}{The list of input options provided by user.}
  \item{vl}{The list that holds input phenotype variable list and input data coding list.}
  \item{varName}{The phenotype/variable to be tested, in the form of FIELDID. For example 21022.}
  \item{varType}{The value type of the variable. Should always be "Categorical single" here and used for logging purpose only.}
  \item{thisdata}{The data frame object holds all variables including phenotypes, confounders, depeinding variables and trait of interest.}
  \item{phenoStartIdx}{The column index of the first phenotype.}
}
\details{
The variable will be skipped if less than 500 samples have value.
}
\value{
Return nothing and all output and lof information are written to output directory specified by user.
}


