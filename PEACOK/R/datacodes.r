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

process.data.code <- function(variable_info_base, outcome_info_file, download) {
    # avoid annoying golbal binding from CRAN check
    ValueType <- NULL
    Coding <- NULL
    
    output_path <- file.path(variable_info_base, 'data-codes')
    variable_info <- fread(outcome_info_file, header=TRUE, data.table=TRUE)
    output_info <- variable_info %>% filter(ValueType %in% c("Categorical single","Categorical multiple")) %>% 
                            dplyr::select(ValueType,Coding) %>%
                            unique() %>% arrange(Coding) %>%
                            mutate(ordinal = -1) %>%
                            mutate(ordering = "") %>%
                            mutate(reassignments = "") %>%
                            mutate(default_value = "") %>%
                            mutate(default_related_field = "") %>%
                            mutate(tree = FALSE) %>%
                            mutate(num_fields = 2) %>%
                            mutate(number_values = 0) %>%
                            mutate(values = "") 
    for (i in 1:dim(output_info)[1]) {
        code <- output_info[i,"Coding"]
        var_type <- output_info[i,"ValueType"]
        output.file <- file.path(output_path,paste0("datacode-",code,".tsv"))
        tryCatch({
            if (download) {
                f <- paste0("http://biobank.ctsu.ox.ac.uk/showcase/codown.cgi?id=",code, "&btn_glow=Download")
                download.file(f, destfile = output.file, quiet = TRUE)
            }
            if (file.exists(output.file)) {
                  cat(paste0("\nInfo: Parsing ", paste0("datacode-",code,".tsv")))
                  data.dict <- fread(output.file, header=TRUE, data.table=TRUE) 
                  output_info$tree[i] <- "parent_id" %in%  names(data.dict)
                  output_info$num_fields[i] <- dim(data.dict)[2]
                  output_info$number_values[i] <- dim(data.dict)[1]
                  output_info$values[i] <- paste0(data.dict$coding, collapse = "|")
                  if (output_info[i,"ValueType"] == "Categorical single") {
                      output_info$ordinal[i] <- "0" # assume unordered
                      output_info$ordering[i] <- ""
                  } else {
                      output_info$ordinal[i] <- -1
                  }
            } else {
                  if (download) {
                      cat(paste0("\nERROR: Downloading ", f, " failed"))
                  } else {
                      cat(paste0("\nERROR: Data code File ", f, " was not found. Consider running with download option."))
                  }
            }
        ## END TRYCATCH
        }, error = function(e) {
            cat(paste("\nERROR:", gsub("[\r\n]", "", e), sep=" "))
        })
    }
    cat("\n")
    output_info <- output_info %>% rename(dataCode = Coding)
    return(output_info)
}

is.tree.code <- function(data_code_file) {
    datacode <- fread(data_code_file, sep='\t', header=TRUE, data.table=TRUE, na.strings=c("", "NA"))
    return ("parent_id" %in%  names(datacode));
}

annotate.tree.code <- function(data_code_file) {
    # avoid annoying golbal binding from CRAN check
    parent_id <- NULL
    level <- NULL
    meaning <- NULL
    node_id <- NULL
    root_node_id <- NULL
    
    datacode <- fread(data_code_file, sep='\t', header=TRUE, data.table=TRUE, na.strings=c("", "NA"))
    if (! "parent_id" %in%  names(datacode)) {
        stop(paste0("parent_id column not found in ", data_code_file), call.=FALSE)
    }
    
    #parse tree structure
    datacode <- datacode %>% mutate(meaning = gsub(",","|",meaning)) %>% 
                             mutate(level = ifelse(parent_id ==0,0,NA)) %>% 
                             mutate(root_node_id = ifelse(parent_id ==0,node_id,NA))
    all_done <- FALSE
    current_level <- 0
    while (!all_done) {
        parent_nodes <- datacode$node_id[which(datacode$level==current_level)]
        datacode <- datacode %>% mutate(level = ifelse(parent_id %in% parent_nodes, current_level+1, level))  
        for (pnode in parent_nodes) {
            root_row <- datacode %>% filter(node_id == pnode)
            root_value <- root_row$root_node_id
            datacode <- datacode %>%  mutate(root_node_id = ifelse(parent_id %in% pnode,root_value, root_node_id))
        }
        all_done <- !any(is.na(datacode$level))
        current_level <- current_level + 1
    }
    
    topNodes <- datacode %>% filter(parent_id == 0) %>% dplyr::select(root_node_id,meaning) %>% rename(root = meaning)
    
    datacode <- datacode %>% inner_join(topNodes)
    return(datacode)
  
}

.getDataCode <- function(opt, vl, varName) {
    # avoid annoying golbal binding from CRAN check
    dataPheno <- vl$phenoInfo[which(vl$phenoInfo$FieldID==varName),]
    dataCode <- dataPheno$DATA_CODING
    data_code_base <- file.path(dirname(opt$datacodingfile),'data-codes')
    if (!file.exists(data_code_base)) {
        stop(paste0("Required dota codes folder ", data_code_base, " not found" ), call.=FALSE)
    }
    data_code_file <- file.path(data_code_base, paste0("datacode-", dataCode, ".tsv"))
    if (!file.exists(data_code_file)) {
        stop(paste0("Required dota codes file ", data_code_file, " not found" ), call.=FALSE)
    }
    datacode <- fread(data_code_file, sep='\t', header=TRUE, data.table=TRUE, na.strings=c("", "NA"))
    if (is.tree.code(data_code_file)) {
        datacode <- annotate.tree.code(data_code_file)
    }
    return(datacode)
}