# The MIT License (MIT)
# Copyright (c) 2017 Louise AC Millard, MRC Integrative Epidemiology Unit, University of Bristol
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

# Processing integer fields, namely:
# 1) Reassigning values as specified in the data code information file
# 2) Generate a single value if there are several values (arrays) by taking the mean
# 3) Treating this field as continuous if at least 20 distinct values.
# Otherwise treat as binary or ordered categorical if 2 or more than two values. 
test.integer <- function(opt, vl, varName, varType, thisdata, phenoStartIdx) {
  	cat("INTEGER || ")
  	pheno <- thisdata[,phenoStartIdx:ncol(thisdata)]
  	isExposure <- get.is.exposure(vl, varName)
  	if (is.numeric(as.matrix(pheno))) {
  	    pheno <- .reassignValue(vl, pheno, varName)
  	    
      	## average if multiple columns
      	if (!is.null(dim(pheno))) {
            phenoAvg <- rowMeans(pheno, na.rm=TRUE)
      		  # if participant only has NA values then NaN is generated so we convert back to NA
      		  phenoAvg <- .replaceNaN(phenoAvg)
        } else {
              phenoAvg <- pheno
        }
  
      	uniqVar <- unique(na.omit(phenoAvg))
      	# if >=20 separate values then treat as continuous
      	if (length(uniqVar)>=20) {
      		  thisdatanew <- cbind.data.frame(thisdata[,1:(phenoStartIdx -1)], phenoAvg)
      		  .test.continuous2(opt, vl, varName, varType, thisdatanew, phenoStartIdx)
      		  .incrementCounter( "int.continuous")
      	} else {
      		  ## remove categories if < 10 examples
      	    phenoAvg <- .testNumExamples(opt, phenoAvg)
      		  ## binary if 2 distinct values, else ordered categorical
      		  phenoFactor <- factor(phenoAvg)
      		  numLevels <- length(levels(phenoFactor))
      		  if (numLevels<=1) {
      			    cat("SKIP (number of levels: ",numLevels,")",sep="")
      		      .incrementCounter("int.onevalue")
      		  } else if (numLevels==2) {
      		      .incrementCounter("int.binary")
      			    # binary
      			    thisdatanew <- cbind.data.frame(thisdata[,1:(phenoStartIdx -1)], phenoFactor);
      			    binary.logistic.regression(opt, varName, varType, thisdatanew, isExposure,phenoStartIdx)
      		  } else {
      		      .incrementCounter("int.catord")
      			    cat("3-20 values || ")
      			    # we don't use equal sized bins just the original integers (that have >=10 examples) as categories
      			    thisdatanew <- cbind.data.frame(thisdata[,1:(phenoStartIdx -1)], phenoFactor)
      
      			    # treat as ordinal categorical
      			    test.categorical.ordered(opt, vl, varName, varType, thisdatanew, phenoStartIdx)
      		  }
      	}
  	} else {
  	    cat("SKIP Integer type but not numeric",sep="")
  	}
}
