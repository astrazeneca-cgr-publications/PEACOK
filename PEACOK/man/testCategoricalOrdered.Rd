\name{test.categorical.ordered}
\alias{test.categorical.ordered}
\title{
Performs ordered logistic regression test
}
\description{
Performs ordered logistic regression test and saves results in ordered logistic results file.
}
\usage{
test.categorical.ordered(opt, vl, varName, varType, thisdata, phenoStartIdx, 
                         orderStr = "")
}
\arguments{
  \item{opt}{The list of input options provided by user.}
  \item{vl}{The list that holds input phenotype variable list and input data coding list.}
  \item{varName}{The phenotype/variable to be tested, in the form of FIELDID. For example 21022.}
  \item{varType}{The value type of the variable. Should always be "Categorical single" here and used for logging purpose only.}
  \item{thisdata}{The data frame object holds all variables including phenotypes, confounders, depeinding variables and trait of interest.}
  \item{phenoStartIdx}{The column index of the first phenotype.}
   \item{orderStr}{The Order String for the variable.}
}


\details{
The variable will be skipped if less than 500 samples have value or there are too many levels.
}
\value{
Return nothing and all output and lof information are written to output directory specified by user.
}


