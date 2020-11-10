\name{auto.update.variable.info}
\alias{auto.update.variable.info}

\title{
Update and consolidate a previous version of variable info file with current version of UKB Data Dictionary.
}
\description{
Taking a previous version of the variable info file, and a copy of the UKB Data Dictionary, this function will consolidate and update the previous version based on the new Data Dictionary and output a new version that generally needs to be curated. It also outputs  list of fields that are added or removed from the data dictionary.
}
\usage{
auto.update.variable.info(old_variable_info_file, UKB_data_dictionary_file)
}
\arguments{
  \item{old_variable_info_file}{ A previous versionn of variable info file, which was curated to have PHESANT specific columns for all varibles.
}
  \item{UKB_data_dictionary_file}{ A newer version of UKB Data Dictionary, which can be ontained by calling download.data.dictionary function.}
}

\value{
Return a list of two data frames: updated_fields for all added/removed fields and variable_info for consolidated variable info. 
}

\examples{
    #library(PEACOK)
    #args <- commandArgs(T)
    
    #variable_info_base <- "../variable-info"
    #output_path <- file.path(variable_info_base, 'update-outcome-info')
    #download.data.dictionary(output_path)
    
    #old_variable_info <- file.path(variable_info_base,'outcome-info.tsv')
    #UKB_data_dictionary <- file.path(output_path,"Data_Dictionary_Showcase.csv")
    #updated_variable_info <- file.path(output_path,'new-field-list.tsv')
    #new_variable_info <- file.path(output_path,'outcome-info-new.tsv')
    
    #updated_info <- auto.update.variable.info(old_variable_info, UKB_data_dictionary)
    #if (is.null(updated_info$updated_fields)) {
    #    cat("\nInfo: New new fields were added or deleted")
    #} else {
    #    cat(paste0("Info: Updated fields are writen to file: ", updated_variable_info))
    #    fwrite(updated_info$updated_fields, sep='\t', file=updated_variable_info)
    #}
    # Write out, and make sure it's tab separated.
    #fwrite(updated_info$variable_info, sep='\t', file = new_variable_info)
    #cat(paste0("\nInfo: Updated variable info is writen to file: ", new_variable_info))
}

\keyword{ variable info }
\keyword{ UKB Data dictionary }
