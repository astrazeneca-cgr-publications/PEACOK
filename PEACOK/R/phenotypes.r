# The MIT License (MIT)
# Copyright (c) 2017 Louise AC Millard, MRC Integrative Epidemiology Unit, University of Bristol
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

## load phenotypes from phenotype file
load.phenotypes <- function(opt) {
    ## is not running 'all' then we determine the start and end idxs of phenotypes that we test, so that we can parallelise into multiple jobs
    if (opt$varTypeArg == "all") {
         # reading all data at once
        data <- fread(opt$phenofile, sep=',', header=TRUE, data.table=FALSE, na.strings=c("", "NA"))
        names(data) <- .UKB2PEACOK(names(data), opt$userId)
    } else {
        if (opt$varTypeArg == "user") {
            phenosToTest <- .getPheno2TestForFieldList(opt)
        } else {
            phenosToTest <- .getPheno2TestForPart(opt)
        }
        ## read in the right table columns - a subset of the data file
        data <- fread(opt$phenofile, select=.PEACOK2UKB(phenosToTest,opt$userId), sep=',', header=TRUE, data.table=FALSE, na.strings=c("", "NA"))
        names(data) <- phenosToTest
    } 

  	## this is type conversion as used in the read.table function (that we used to use ((this was changed because read.table cannot read column subsets))
  	#data <- data.frame(lapply(data, function(x) type.convert(as.character(x), as.is = TRUE)))
    data <- data.frame(lapply(data, function(x) as.numeric(x)))
    colnames(data)[1] <- "userID"
    return(data)
}

# Validate the contents of the phenotype file
validate.phenotype.input.header <- function(opt) {
  	print("Validating phenotype data ...")
  	## get just first row so we can check the column names
  	phenoNames <- fread(opt$phenofile, header=FALSE, data.table=FALSE, nrows=1, sep=',')
  	
  	### pheno file validation
  	print(paste("Number of columns in phenotype file: ", length(phenoNames),sep=""))
  	## check user id exists in pheno file
  	idx1 <- which(phenoNames == opt$userId)
  	if (length(idx1)==0) {
        stop(paste("phenotype file doesn't contain userID colunn:", opt$userId), call.=FALSE)
    }
  
  	# we only need the confounders if we are actually running the tests
  	if (opt$save==FALSE && is.null(opt$confounderfile)) {
  	    phenoNames <- .UKB2PEACOK(phenoNames, opt$userId)
      	## confounder variables exist in pheno file
      	if (! "x21022_0_0" %in% phenoNames) {
            stop("phenotype file doesn't contain required age colunn: x21022_0_0", call.=FALSE)
      	}
      	
        if (! "x31_0_0" %in% phenoNames ) {
            stop("phenotype file doesn't contain required sex colunn: x31_0_0", call.=FALSE)
        }
      
      	if (opt$genetic ==TRUE) {
          	if (! "x22000_0_0" %in% phenoNames) {
          	    stop("phenotype file doesn't contain required genetic batch colunn: x22000_0_0", call.=FALSE)
          	}	
      	}
      
      	## if running with sensitivity option then check extra columns exist in pheno file (genetic PCs and assessment centre)
      	if (opt$sensitivity==TRUE) {
        		if (opt$genetic ==TRUE) {
          			## check first 10 genetic PCs exist
          			for (i in 1:10) {
            				if (! paste("x22009_0_", i, sep="") %in% phenoNames) {
               			    stop(paste("phenotype file doesn't contain required genetic principal component colunn: x22009_0_", i, sep=""), call.=FALSE)
              			}
          			}
        		}
        		## assessment centre field
          	if (! "x54_0_0" %in% phenoNames) {
          	    stop("phenotype file doesn't contain required assessment centre colunn: x54_0_0", call.=FALSE)
          	}
        }
  	}
  	print("Phenotype file validated")
}

.getFields <- function(phenoFile, userId) {
    phenoVars <- fread(phenoFile, header=FALSE, data.table=FALSE, nrows=1, sep=',')
    phenoVars<- .UKB2PEACOK(phenoVars, userId)
    phenoVars <- phenoVars[which(phenoVars!=userId)]
    return(phenoVars)
}

.getPheno2TestForFieldList <- function(opt) {
    phenoVars <- .getFields(opt$phenofile, opt$userId)
    fields <- paste0("x",strsplit(opt$fieldlist,",")[[1]])
    ## user ID always included
    phenosToTest <- c(opt$userId)
    for (f in fields) {
        phenosToTest <- c(phenosToTest,phenoVars[startsWith(phenoVars, paste0(f,"_"))])
    }
    return(phenosToTest)
}

.getPheno2TestForPart <- function(opt) {
    phenoVars <- .getFields(opt$phenofile, opt$userId)
    #deal with 90002 fields while instances runs from -1, 0, and 1.
    fields <- unique(gsub("_", "", gsub("_[0-9]+_[0-9]+$", "", phenoVars)))
    
    partSize <- ceiling(length(fields)/opt$numParts)
    partStart <- (opt$partIdx-1)*partSize + 1
    if (opt$partIdx == opt$numParts) {
        partEnd <- length(phenoVars)
    } else {
        partEnd <- partStart + partSize - 1
    }
    
    ## user ID always included
    phenosToTest <- c(opt$userId)
    for (f in fields[partStart:partEnd]) {
        phenosToTest <- c(phenosToTest,phenoVars[startsWith(phenoVars, paste0(f,"_"))])
    }
    return(phenosToTest)
}

.fixOddFieldsToCatMul <- function(vl, data) {
    # examples are variables: 40006, 40011, 40012, 40013
    # get all variables that need their instances changing to arrays
    dataPheno <- vl$phenoInfo[which(vl$phenoInfo$CAT_SINGLE_TO_CAT_MULT=="YES-INSTANCES"),]
    for (i in 1:nrow(dataPheno)) {
        varID <- dataPheno[i,]$FieldID		
        varidString <- paste("x",varID,"_", sep="")			
        
        # get all columns in data dataframe for this variable	
        colIdxs <- which(grepl(varidString,names(data)))
        
        # change format from xvarid_0_0, xvarid_1_0, xvarid_2_0, to xvarid_0_0, xvarid_0_1, xvarid_0_2
        count <- 0
        for (j in colIdxs) {	
            colnames(data)[j] <- paste(varidString, "0_", count, sep="")
            count <- count + 1
        }				
    }
    return(data)
}

