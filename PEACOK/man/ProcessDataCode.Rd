\name{process.data.code}
\alias{process.data.code}
\title{
Download and parse UKB data codings
}
\description{
Optionally Download the UKB data coding files  and then parse UKB into a format that can be used by the package.
}
\usage{
process.data.code(variable_info_base, outcome_info_file, download)
}
\arguments{
  \item{variable_info_base}{The folder while all coding will be stored.}
  \item{outcome_info_file}{The "outcone-info"" file as ddefined in PHESANT.}
  \item{download}{Indicator if all data coding files should be downloaded or not. If false, variable_info_base should have an copy of previously downloaded data coding files.
}
}
\details{
This function allows user to download and parse the data codings for "Categorical single" and "Categorical multiple" fields based on the infomation provided in "outcome-info" file, which lists all the fields to be recognized by the package. 
}
\value{
A data frame that holds extra info about each data codings, which can be then merged with previously curated data coding file.
}
