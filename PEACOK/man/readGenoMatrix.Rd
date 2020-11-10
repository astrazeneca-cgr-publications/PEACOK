\name{read.geno.matrix}
\alias{read.geno.matrix}
\title{
Read the whole or a trunk of the Collapsing Matrix by rows. 
}
\description{
Read the whole or a trunk of Collapsing Matrix by rows. 
}
\usage{
read.geno.matrix(matrix_file, gene_start = NULL, gene_end = NULL)
}
\arguments{
  \item{matrix_file}{
  The collapsing matrix to be read.
}
  \item{gene_start}{
  The start gene of the matrix if read by trunk.
}
  \item{gene_end}{
  The end gene of the matrix if read by trunk.
}
}

\value{
Returns a data frame that holds the whole or a trunk of the collapsing matrix.
}
