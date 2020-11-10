\name{choose.reference.category}
\alias{choose.reference.category}
\title{
Choose Reference Category 
}
\description{
Convert an integer vector with possible NA's to a factor and use the category/factor with most number of examples as reference.
}
\usage{
choose.reference.category(pheno)
}
\arguments{
  \item{pheno}{ An integer vector with possible NA's. }
}

\value{
  A factor with the most number of examples as reference.
}
