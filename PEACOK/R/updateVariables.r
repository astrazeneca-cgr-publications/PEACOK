# The MIT License (MIT)
# Copyright (c) 2020 Quanli Wang, Center for Genomic Research, AstraZeneca
#
# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated 
# documentation files (the "Software"), to deal in the Software without restriction, including without 
# limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of 
# the Software, and to permit persons to whom the Software is furnished to do so, subject to the following 
# conditions:
#
# The above copyright notice and this permission notice shall be included in all copies or substantial portions 
# of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED 
# TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL 
# THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF 
# CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
# DEALINGS IN THE SOFTWARE.

download.data.dictionary <- function(output_path) {
    output.file <- file.path(output_path,'Data_Dictionary_Showcase.csv')
    tryCatch({
        f <- "http://biobank.ctsu.ox.ac.uk/%7Ebbdatan/Data_Dictionary_Showcase.csv"
        download.file(f, destfile = output.file, quiet = TRUE)
        if (file.exists(output.file)) {
              cat("Info: Download Success\n")
              data.dict <- read.csv(output.file, na.strings = "NONASTRING")
              data.dict[is.na(data.dict)] <-''
              data.dict <- data.frame(lapply(data.dict, gsub, pattern = '"', replacement = '', fixed = TRUE))
              write.table(data.dict, file = output.file, 
                          sep = '\t', quote = FALSE, row.names = FALSE, na = '', eol = '\n')
              return(TRUE)
        } else {
              cat(paste0("\nERROR: Downloading ", f, " failed"))
              return(FALSE)
        }
    ## END TRYCATCH
    }, error = function(e) {
        cat(paste("\nERROR:", gsub("[\r\n]", "", e), sep=" "))
        return(FALSE)
    })
}

.getFieldUpdate <- function (df) {
    # Before we do that, find out which phenotypes have been added, and write this subset to disk to determine which of these to include.
    df_added <- df[is.na(df$Field.x),]
    df_removed <- df[is.na(df$Field.y),]
    
    # Now write the FieldID and the Field to disk.
    names(df_added)[names(df_added) == 'Field.y'] <- 'Field'
    names(df_removed)[names(df_removed) == 'Field.x'] <- 'Field'
    
    if (dim(df_added)[1] >0) {
        df_added$action <- 'added'
    }
    if (dim(df_removed)[1] >0) {
        df_removed$action <- 'removed'
    }
    if (dim(df_added)[1] >0 && dim(df_removed)[1] >0) {
        df_to_check <- rbind(df_added[,c('FieldID', 'Field','action')], df_removed[,c('FieldID', 'Field','action')])
    } else if (dim(df_added)[1] >0) {
        df_to_check <- df_added[,c('FieldID', 'Field','action')]
    } else if (dim(df_removed)[1] >0) {
        df_to_check <- df_removed[,c('FieldID', 'Field','action')]
    } else {
        df_to_check = NULL
    }
    return(df_to_check)
}


auto.update.variable.info <- function(old_variable_info_file, UKB_data_dictionary_file) {
    # Merge the two data frames.
    old_df <- fread(old_variable_info_file, header=TRUE, data.table=FALSE)
    new_df <- fread(UKB_data_dictionary_file, header=TRUE, data.table=FALSE)
    
    PHSEANT_value_types <- c("Integer","Categorical single","Continuous","Categorical multiple") 
    new_df <- new_df[new_df$ValueType %in% PHSEANT_value_types,]
    
    df <- merge(x=old_df, y=new_df, by="FieldID", all=TRUE)

    updated_fields <- .getFieldUpdate(df)
    
    #remove fields that no longer exists in new UKB Data Dictionary
    df <- df[!is.na(df$Field.y),]
    
    #add "X" to PEACOK specific columns for newly added fields
    PEACOK_specific <- c('TRAIT_OF_INTEREST','EXCLUDED','CAT_MULT_INDICATOR_FIELDS','CAT_SINGLE_TO_CAT_MULT','DATA_CODING')
    df[is.na(df$Field.x),PEACOK_specific] <- 'X'
    
    #remove old columns and use new columns from new UKBData Dictionary
    df <- df[,-grep('\\.x', names(df))]
    names(df) <- gsub('\\.y', '', names(df))
    
    return(list(updated_fields = updated_fields, variable_info = df))
}




