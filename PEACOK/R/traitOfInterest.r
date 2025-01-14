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


##
## load trait of interest, either from separate trait of interest file, or from phenotype file

load.trait.of.interest <- function(opt, phenotypes) {
    if (opt$save==TRUE) { 
	      # saving not running tests so we don't have a trait of interest
        # add pretend trait of interest so other code doesn't break
        numRows <- nrow(phenotypes)
	      data <- cbind.data.frame(phenotypes$userID, rep(-1, numRows))
    } else {
	      # load the trait of interest specified by the user
        if (is.null(opt$traitofinterestfile)) {
              print("Extracting trait of interest from pheno file ...")
		          data = fread(opt$phenofile, select=c(opt$userId, opt$traitofinterest), sep=',', header=TRUE, data.table=FALSE)
        } else {
              print("Loading trait of interest file ...")
              data = fread(opt$traitofinterestfile, select=c(opt$userId, opt$traitofinterest), sep=',', header=TRUE, data.table=FALSE)
        }
        userId.index <- which(names(data) == opt$userId)
        names(data) <- paste0("x", gsub("[-|.]","_",names(data)))
        names(data)[userId.index] <- opt$userId
	      data <- data.frame(lapply(data,function(x) type.convert(as.character(x))))
    }

    colnames(data)[1] <- "userID"
    colnames(data)[2] <- "geno"
    return(data)
}

# Validate the contents of the trait of interest file
validate.trait.input.header <- function(opt) {
    if (opt$save!=TRUE) {
      	### get header of trait of interest file or pheno file (if no trait of interest file is specified
        print("Validating trait of interest data ...")
      	if (is.null(opt$traitofinterestfile)) {
      		  snpIn  <- fread(opt$phenofile, header=FALSE, data.table=FALSE, nrows=1, sep=',')
      	} else {
      		  snpIn <- fread(opt$traitofinterestfile, header=FALSE, data.table=FALSE, nrows=1, sep=',') 
      	}
      
      	### trait of interest file validation
      	print(paste("Number of columns in trait of interest file:", ncol(snpIn),sep=""))
      	## check user id exists in snp file
      	idx1 <- which(snpIn == opt$userId)
      	if (length(idx1)==0) {
        		if (is.null(opt$traitofinterestfile)) {
        			stop(paste("Phenotype file doesn't contain userID colunn:", opt$userId), call.=FALSE)
        		} else {
        			stop(paste("Trait of interest file doesn't contain userID colunn:", opt$userId), call.=FALSE)
        		}
      	}
      	
      	## check trait of interest exists in trait of interest file
      	idx2 <- which(snpIn == opt$traitofinterest)
      	if (length(idx2)==0) {
        		if (is.null(opt$traitofinterestfile)) {
        			  stop(paste("No trait of interest file specified, and phenotypes file doesn't contain trait of interest variable column:", opt$traitofinterest), call.=FALSE)
        		} else {
        			  stop(paste("Trait of interest file doesn't contain trait of interest variable column:", opt$traitofinterest), call.=FALSE)
        		}
      	}
      	print("Trait of interest file validated")
      }
}

