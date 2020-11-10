\name{get.is.exposure}
\alias{get.is.exposure}
\title{
Check field exposure 
}
\description{
Check whether a field denotes the trait of interest, as specified.
}
\usage{
get.is.exposure(vl, varName)
}
\arguments{
  \item{vl}{variable list object
}
  \item{varName}{ the field ID to be checked
}
}

\value{
Return true if the field has a "YES" value in TRAIT_OF_INTEST column.
}
