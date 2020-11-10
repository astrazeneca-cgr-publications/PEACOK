\name{annotate.tree.code}
\alias{annotate.tree.code}
\title{
Add annotation columns to data code that has tree structures
}
\description{
Add annotation columns to data code that has tree structures. For a given datacode file with tree structures, it annotate it by add more columns indicating tree level, root node meaning as well as replacing comma in orignal meaning column with pipe (|) for better csv ewxperience.
}
\usage{
annotate.tree.code(data_code_file)
}
\arguments{
  \item{data_code_file}{an original datacoding file downloaded from UKB website which has a tree structure in node. }
}

\value{
Returns a dataframe with extra annotation columns beside the original columns. 
}

\examples{
  #datacode <- annotate.tree.code(data_code_file)
}
