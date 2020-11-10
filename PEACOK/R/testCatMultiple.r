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


test.categorical.multiple <- function(opt, vl, varName, varType, thisdata, phenoStartIdx) {
   # avoid annoying golbal binding from CRAN check
    node_id <- NULL
    coding <- NULL
    level <- NULL
    parent_id <- NULL
    
    datacode <- .getDataCode(opt, vl,varName)
    if (! "tree" %in% names(vl$dataCodeInfo) || opt$forceflat || !"parent_id" %in%  names(datacode)) {
        .test.categorical.multipleFlat(opt, vl, varName, varType, thisdata, phenoStartIdx)
    } else {
        commondata <- thisdata[,1:(phenoStartIdx -1)]
  	    userID <- thisdata[,"userID", drop=FALSE]
  	    pheno <- thisdata[,phenoStartIdx:ncol(thisdata), drop=FALSE]
  	    pheno <- .reassignValue(vl, pheno, varName)
  	    
  	    
  	    max_tree_level <- max(datacode$level)
        for (i in (max_tree_level:0)) {
            cat(paste0("CAT-MULTIPLE-TREE-", i)," || ")
  	        ## get unique values from all columns of this variable
  	        uniqueValues <- sort(unique(na.omit(as.vector(as.matrix(pheno)))))
  	        current_codes <- datacode$coding[datacode$level ==i]
  	        
  	        #only work with current level
  	        uniqueValues <- intersect(current_codes,uniqueValues) 
          	.testImp(opt, vl, varName, varType, pheno, commondata, userID, uniqueValues,datacode)
          	
          	#update pheno to parent levels
          	pheno <- as.matrix(pheno)
          	for (v in uniqueValues) {
          	    codeRow <- datacode %>% filter(coding == v)
          	    if (dim(codeRow)[1] > 1) {
          	        cat("multiple nodes have same coding, will not colapse to parent nodes\n")
          	        #return(NULL)
          	    } else {
          	        if (codeRow$parent_id != 0) {
              	        codeRowParent <- datacode %>% filter(node_id == codeRow$parent_id)
              	        pheno[pheno == v] <- codeRowParent$"coding"
              	    }
          	    }
          	    
          	}
          	pheno <- as.data.frame(pheno)
        }
    }
}
  
# Performs preprocessing of categorical (multiple) fields, namely:
# 1) Reassigning values as specified in data coding file
# 2) Generating binary variable for each category in field, restricting to correct set of participants as specified
# in CAT_MULT_INDICATOR_FIELDS field of variable info file (either NO_NAN, ALL or a field ID)
# 3) Checking derived variable has at least catmultcutoff cases in each group
# 4) Calling binary.logistic.regression function for this derived binary variable
.test.categorical.multipleFlat <- function(opt, vl, varName, varType, thisdata, phenoStartIdx) {
  	cat("CAT-MULTIPLE || ")
  	pheno <- thisdata[,phenoStartIdx:ncol(thisdata), drop=FALSE]
  	pheno <- .reassignValue(vl, pheno, varName)
    
  	commondata <- thisdata[,1:(phenoStartIdx -1)]
  	userID <- thisdata[,"userID", drop=FALSE]
  	
  	## get unique values from all columns of this variable
  	uniqueValues <- sort(unique(na.omit(as.vector(as.matrix(pheno)))))
  	.testImp(opt, vl, varName, varType, pheno, commondata, userID, uniqueValues, datacode = NULL)
}

.testImp <- function(opt, vl, varName, varType, pheno, commondata, userID, uniqueValues, datacode) {
    # avoid annoying golbal binding from CRAN check
    coding <- NULL
   
    
    numCols <- ncol(pheno)
  	numRows <- nrow(pheno)
  	isNumeric <- is.numeric(as.matrix(pheno))
  	phenoStartIdx <- dim(commondata)[2] + 1
  	## for each value create a binary variable and test this
  	for (variableVal in uniqueValues) {
  		  ## numeric negative values we assume are missing - check this
  		  if(is.numeric(variableVal) & variableVal<0) {
  			    cat("SKIP_val:[Fix Me]", variableVal," < 0", sep="")
  			    next
  		  }
  	
  		  cat(" CAT-MUL-BINARY-VAR ", variableVal, " || ", sep="")
  		  .incrementCounter("catMul.binary")
  		
  		  # make zero vector and set 1s for those with this variable value
  		  # make variable for this value
  		  varBinary <- rep.int(0,numRows)
  		  varBinary[which(pheno == variableVal, arr.ind=TRUE)[,"row"]] <- 1
  		  varBinaryFactor <- factor(varBinary)
  
  		  ## data for this new binary variable
  		  newthisdata <- cbind.data.frame(commondata, varBinaryFactor)
  
  		  ## one of 3 ways to decide which examples are negative
        idxsToRemove <- .restrictSample(vl, varName, pheno, variableVal, userID,isNumeric)
  		  if (!is.null(idxsToRemove) & length(idxsToRemove) > 0) {
  			    newthisdata <- newthisdata[-idxsToRemove,]
  		  }
  
  		  facLevels <- levels(newthisdata[,phenoStartIdx])		
  		  idxTrue <- length(which(newthisdata[,phenoStartIdx]==facLevels[1]))
  	    idxFalse <- length(which(newthisdata[,phenoStartIdx]==facLevels[2]))
                  
  	    if (idxTrue < opt$catmultcutoff || idxFalse < opt$catmultcutoff) {
  	        cat("CAT-MULT-SKIP-", opt$catmultcutoff, " (", idxTrue, " vs ", idxFalse, ") || ", sep="")
  			    .incrementCounter(paste("catMul.", opt$catmultcutoff, sep=""))
  	    } else {
  			    isExposure <- .getIsCatMultExposure(vl, varName, variableVal)
  			    .incrementCounter(paste("catMul.over", opt$catmultcutoff, sep=""))
  		     	# binary - so logistic regression
  			    long_name <- paste(varName, gsub(" ", "",variableVal),sep="#")
  			    if (!is.null(datacode)) {
  			        codeRow  <- datacode %>% filter(coding == variableVal)
  			        meaning <- gsub(" ", "",codeRow$meaning)
  			        meaning <- gsub(",", "|",codeRow$meaning)
  			        long_name <- paste(long_name, meaning, sep="#")
  			    } 
  			    binary.logistic.regression(opt, long_name, 
  			                             varType, newthisdata, isExposure, phenoStartIdx)
  		  }
  	}
}
# restricts sample based on value in CAT_MULT_INDICATOR_FIELDS column of variable info file,
# either NO_NAN, ALL or a field ID
# returns idx's that should be removed from the sample
.restrictSample <- function(vl, varName,pheno,variableVal, userID,isNumeric) {
  	# get definition for sample for this variable either NO_NAN, ALL or a variable ID
  	varIndicator <- vl$phenoInfo$CAT_MULT_INDICATOR_FIELDS[which(vl$phenoInfo$FieldID==varName)]
  	return(.restrictSample2(vl, varName,pheno,varIndicator,variableVal, userID,isNumeric))
}

.restrictSample2 <- function(vl, varName,pheno, varIndicator,variableVal, userID, isNumeric) {
  	if (varIndicator=="NO_NAN") { # remove NAs
    		## remove all people with no value for this variable
    		# row indexes with NA in all columns of this cat mult field		
    		ind <- apply(pheno, 1, function(x) all(is.na(x)))
    		naIdxs <- which(ind==TRUE)
    		cat("NO_NAN Remove NA participants ", length(naIdxs), " || ", sep="")
  	} else if (varIndicator=="ALL") {
    		# use all people (no missing assumed) so return empty vector
    		# e.g. hospital data and death registry
    		naIdxs = cbind()
    		cat("ALL || ")
  	} else if (varIndicator!="") {
    		# remove people who have no value for indicator variable
    		indName <- paste("x",varIndicator,"_0_0",sep="");
    		cat("Indicator name ", indName, " || ", sep="");
    		indvarx <- merge(userID, vl$indicatorFields, by="userID", all.x=TRUE, all.y=FALSE, sort=FALSE)		
    		indicatorVar <- indvarx[,indName]
    
    		# remove participants with NA value in this related field
    		indicatorVar <- .replaceNaN(indicatorVar)
    		naIdxs <- which(is.na(indicatorVar))
    		cat("Remove indicator var NAs: ", length(naIdxs), " || ", sep="");
    
    		if (is.numeric(as.matrix(indicatorVar))) {
      			# remove participants with value <0 in this related field - assumed missing indicators
      			lessZero <- which(indicatorVar<0)
      			naIdxs <- union(naIdxs, lessZero)
      			cat("Remove indicator var <0: ", length(lessZero), " || ", sep="")
    		}
  	} else {
  		  stop("Categorical multiples variables need a value for CAT_MULT_INDICATOR_FIELDS", call.=FALSE)
  	}
  
  	## remove people with pheno<0 if they aren't a positive example for this variable indicator
  	## because we can't know if they are a negative example or not
  	if (isNumeric) {
    		idxForVar <- which(pheno == variableVal, arr.ind=TRUE)
    		idxMissing <- which(pheno < 0, arr.ind=TRUE)
    
    		# all people with <0 value and not variableVal
    		naMissing <- setdiff(idxMissing,idxForVar)
    		
    		# add these people with unknowns to set to remove from sample
    		naIdxs <- union(naIdxs, naMissing)
    		cat(paste("Removed ", length(naMissing) ," examples != ", variableVal, " but with missing value (<0) || ", sep=""))
  	} else {
  		  cat("Not numeric || ")
  	}
	  return(naIdxs)
}
