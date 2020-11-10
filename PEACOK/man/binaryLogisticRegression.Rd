\name{binary.logistic.regression}
\alias{binary.logistic.regression}
\title{
Perform binary logistic regression
}
\description{
Performs binary logistic regression and stores result in 'results-logistic-binary' results file.
}
\usage{
binary.logistic.regression(opt, varName, varType, thisdata, isExposure, phenoStartIdx)
}
\arguments{
  \item{opt}{The list of input options provided by user.}
  \item{varName}{The phenotype/variable to be tested, in the form of FIELDID. For example 21022.}
  \item{varType}{The value type of the variable. Should always be "Categorical single" here and used for logging purpose only.}
  \item{thisdata}{The data frame object holds all variables including phenotypes, confounders, depeinding variables and trait of interest.}
  \item{isExposure}{If the variable is labeled as trait of interest.}
  \item{phenoStartIdx}{The column index of the first phenotype.}
}
\details{
The variable will be skipped if less than 500 samples have value or there are more than two levels or at least one level has less than 10 observations.
}
\value{
Return nothing and all output and lof information are written to output directory specified by user.
}


