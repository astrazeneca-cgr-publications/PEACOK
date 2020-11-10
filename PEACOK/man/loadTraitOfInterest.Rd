\name{load.trait.of.interest}
\alias{load.trait.of.interest}
\title{
Load trait of interest
}
\description{
Load trait of interest, either from separate trait of interest file, or from phenotype file.
}
\usage{
load.trait.of.interest(opt, phenotypes)
}
\arguments{
  \item{opt}{The list of input options provided by user.}
  \item{phenotypes}{The phenotype data from load.phenotypes function.}
}
\details{
Load trait of interest, either from separate trait of interest file, or from phenotype file.
}
\value{
Returns a data frame of two columns for User ID and trait of interest.
}

\examples{
#toi <- load.trait.of.interest(opt, phenotype)
}
