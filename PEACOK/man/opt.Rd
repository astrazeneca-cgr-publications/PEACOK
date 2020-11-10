\name{opt}
\alias{opt}
\docType{data}
\title{
Default options for "test" dataset
}
\description{
Default options to be used when user runs the program with "--test" option in command line. The files used in "--test" mode will be updated from inst/extdata folder of the package.
}


\details{
The fields are:
List of 15
  \enumerate{
     \item{userId             :} { chr "userId"}
     \item{test               :} { logi TRUE}
     \item{sensitivity        :} { logi FALSE}
     \item{genetic            :} { logi TRUE}
     \item{save               :} { logi FALSE}
     \item{confidenceintervals:} { logi TRUE}
     \item{standardise        :} { logi TRUE}
     \item{help               :} { logi FALSE}
     \item{resDir             :} { chr ""}
     \item{phenofile          :} { chr "inst/extdata/phenotypes.csv"}
     \item{variablelistfile   :} { chr "inst/extdata/variable-lists/outcome-info.tsv"}
     \item{datacodingfile     :} { chr "inst/extdata/variable-lists/data-coding-ordinal-info.txt"}
     \item{traitofinterestfile:} { chr "inst/extdata/exposure.csv"}
     \item{traitofinterest    :} { chr "exposure"}
     \item{varTypeArg         :} { chr "all"}
  }
Please note that the resDir field will be replaced by actual working directory.
}
\examples{
data(opt)
}
\keyword{datasets}
