\name{downsample.by.gender}
\alias{downsample.by.gender}
\title{
Downsample Phenotypes by Gender
}
\description{
  Downsample controls or cases based on gender to reduce the effect of gender imbalance. 
}
\usage{
downsample.by.gender(samples, phenotypes, alpha = 0.05, seed = 1234)
}
\arguments{
  \item{samples}{A dataframe with rows for samples and two columns. The first column is "userID", identifying the samples. The second column is 'is.female", of boolean type, indicating if a sample is female or not.}
  \item{phenotypes}{A dataframe with rows for samples and columns for phenotypes. The first column is "userID", identifying the samples.
}
  \item{alpha}{The significant level(p value) when testing if a phenotype is gender imbalanced. The default is 0.05.}
  \item{seed}{The random seed for reproducibility during downsampling. The default value is 1234.
}
}
\details{
For each input phenotype, 1. If all cases are female, only females are retained in controls. 2. If all cases are male, only mailes are retained in controls. 3). If there are both females and males in cases, it checks gender imbalance using Fisher Exact Test at the signifiant level alpha defined by user(default value is 0.05). 4). If gender imbalance is detected, it  downsamples either males or females in controls so the ratio of controls over cases are matched as closely as possible for gender.
}
\value{
Return a dataframe that has the same structure as input phenotypes dataframe but with excluded samples marked as NA.
}


