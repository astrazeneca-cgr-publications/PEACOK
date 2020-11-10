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

.storeNewVar <- function(userIDData, phenoData, varName, type) {
    # add pheno to dataframe
    newdata = data.frame(userID=userIDData, newvar=phenoData)
    names(newdata)[names(newdata)=="newvar"] = varName
    if (type == "bin") {
        pkg.env$derivedBinary <- merge(pkg.env$derivedBinary, newdata, by="userID", all=TRUE)
    } else if (type == "cont") {
        pkg.env$derivedCont <- merge(pkg.env$derivedCont, newdata, by="userID", all=TRUE)
    } else if (type == "catOrd") {
        pkg.env$derivedCatOrd <- merge(pkg.env$derivedCatOrd, newdata, by="userID", all=TRUE)
    } else if (type == "catUnord") {
        pkg.env$derivedCatUnord <- merge(pkg.env$derivedCatUnord, newdata, by="userID", all=TRUE)
    }
}

## makes a smaller data frame containing the data for a particular test
.makeTestDataFrame <- function(datax, confounders, currentVarValues) {
		thisdata <- datax[,c("geno", "userID")]
		thisdata <- merge(thisdata, confounders, by="userID", all.x=TRUE, all.y=FALSE, sort=FALSE)
		currentVarValues <- cbind.data.frame(datax$userID, currentVarValues)
    colnames(currentVarValues)[1] <- "userID"
		thisdata <- merge(thisdata, currentVarValues, by="userID", all.x=TRUE, all.y=FALSE, sort=FALSE)
		return(thisdata)
}

process.options <- function(test = FALSE) {
    opt_parser <- OptionParser(option_list=option_list)
    if (test) {
        #data(test.opt) 
        # opt will be saved in sysdata
    } else {
        opt <- parse_args(opt_parser)
    }
    opt <- .processArgs(opt, opt_parser)
    return(opt)
}

# Parse the arguments input by the user
# if argument 'test' is used then run test phenome scan
.processArgs <- function(opt, opt_parser) {
    if (opt$test==TRUE) {
        resdir <- file.path(getwd(), "results")
        if (!file.exists(resdir)){
            dir.create(resdir)
            print(paste0("results folder created: ",resdir))
        } else {
            print(paste0("Using existing results folder: ", resdir))
        }
        
        datadir <- system.file("extdata", package="PEACOK")
        
      	# set up the test phenome scan settings
      	opt$resDir <- paste0(resdir, "/")
        opt$userId <- 'userId'
      	opt$phenofile <-  file.path(datadir,'phenotypes.csv')
        opt$variablelistfile <- file.path(datadir,'variable-lists/outcome-info.tsv') 
        opt$datacodingfile <- file.path(datadir,'variable-lists/data-coding-ordinal-info.txt')
      	opt$confidenceintervals <- TRUE	
      
      	if(opt$save == FALSE) {
      		opt$traitofinterestfile <- file.path(datadir,'exposure.csv') 
      		opt$traitofinterest <- 'exposure'
      		opt$sensitivity <- FALSE
      		opt$genetic <- TRUE
      	}
    	  opt <- .processParts(opt, opt_parser,opt$partIdx, opt$numParts)
    } else {
    	## check arguments are supplied correctly
    	if (is.null(opt$phenofile)){
      	  print_help(opt_parser)
      	  stop("phenofile argument must be supplied", call.=FALSE)
    	} else if (!file.exists(opt$phenofile)) {
          stop(paste("phenotype data file phenofile=", opt$phenofile, " does not exist", sep=""), call.=FALSE)
      }
    
    	if (opt$save==FALSE && !is.null(opt$traitofinterestfile) && !file.exists(opt$traitofinterestfile)) {
          stop(paste("trait of interest data file traitofinterestfile=", opt$traitofinterestfile, " does not exist", sep=""), call.=FALSE)
      }
    
    	if (opt$save==FALSE && !is.null(opt$confounderfile) && !file.exists(opt$confounderfile)) {
          stop(paste("confounder data file confounderfile=", opt$confounderfile, " does not exist", sep=""), call.=FALSE)
      }
    
    	if (is.null(opt$variablelistfile)){
      	  print_help(opt_parser)
      	  stop("variablelistfile argument must be supplied", call.=FALSE)
    	} else if (!file.exists(opt$variablelistfile)) {
          stop(paste("variable listing file variablelistfile=", opt$variablelistfile, " does not exist", sep=""), call.=FALSE)
      }
    
    	if (is.null(opt$datacodingfile)){
      	  print_help(opt_parser)
      	  stop("datacodingfile argument must be supplied", call.=FALSE)
    	} else if (!file.exists(opt$datacodingfile)) {
          stop(paste("data coding file datacodingfile=", opt$datacodingfile, " does not exist", sep=""), call.=FALSE)
      }
    
    	if (opt$save==FALSE && is.null(opt$traitofinterest)){
      	  print_help(opt_parser)
      	  stop("traitofinterest argument must be supplied", call.=FALSE)
    	}
    
    	if (is.null(opt$resDir)){
      	  print_help(opt_parser)
      	  stop("resDir argument must be supplied", call.=FALSE)
    	}
    	else if (!file.exists(opt$resDir)) {
    		  stop(paste("results directory resDir=", opt$resDir, " does not exist", sep=""), call.=FALSE)
    	}
    	opt <-.processParts(opt, opt_parser,opt$partIdx, opt$numParts)
    }
    if (opt$save==TRUE) {
    	  print("Saving phenotypes to file. Tests of association will not run!")
    }
    return(opt)
}

# Parse the 'part' arguments and check they are valid
.processParts <- function(opt, opt_parser, pIdx, nParts) {
      if (!is.null(opt$fieldlist)) {
          if (file.exists(opt$fieldlist)) { # read from a file
              fields <- as.character(fread(opt$fieldlist, sep=',', header=FALSE, data.table=FALSE)[,])
          } else {
              fields <- strsplit(opt$fieldlist,",")[[1]]
          }
          phenoVarsAll <- as.character(fread(opt$phenofile, header=FALSE, data.table=FALSE, nrows=1, sep=','))
          for (f in fields) {
              if (!any(startsWith(phenoVarsAll, paste0(f,"-")))) {
                  stop(paste0("Field ", f , " was not found in genotype file"), call.=FALSE)
              }
          }
          opt$varTypeArg <- "user"
          opt$fieldlist <- paste(fields, sep = "", collapse = ",")
          print(paste("Running with field(s): ", opt$fieldlist))
      } else {
        	if (is.null(pIdx) && is.null(nParts)) {
              opt$varTypeArg <- "all"
        		  print(paste("Running with all traits in phenotype file:", opt$phenofile));
          } else if (is.null(pIdx)) {
              print_help(opt_parser)
              stop("pIdx argument must be supplied when nParts argument is supplied", call.=FALSE)
          } else if (is.null(nParts)) {
              print_help(opt_parser)
              stop("nParts argument must be supplied when pIdx argument is supplied", call.=FALSE)
          } else if (pIdx<1 || pIdx>nParts) {
        		  print_help(opt_parser)
              stop("pIdx arguments must be between 1 and nParts inclusive", call.=FALSE)
        	} else {
              opt$varTypeArg <- paste(pIdx, "-", nParts, sep="");
              print(paste("Running with part",pIdx,"of",nParts," in phenotype file:", opt$phenofile))
        	}
      }
      return(opt)
}

.PEACOK2UKB <- function(pname, userId) {
    UKB.names <- pname
    userId.index <- which(pname == userId)
    UKB.names <- sub("__","_-",UKB.names) #deal with field 90002, which has instances -1, 0 and 1
	  UKB.names <- sub("_","-",UKB.names) #replace first "_" with "-"
	  UKB.names <- sub("_",".",UKB.names) #replace second "_" with "."
	  UKB.names <- sub("^x", "", UKB.names)
	  if (length(userId.index) > 0) {
	      UKB.names[userId.index] <- userId
	  }
	  return(UKB.names)
}

.UKB2PEACOK <- function(pname, userId) {
    userId.index <- which(pname == userId)
    PEACOK.name <- gsub("[-|.]","_",pname)
    PEACOK.name <- paste0('x',PEACOK.name)
    if (length(userId.index) > 0) {
        PEACOK.name[userId.index] <- userId
    }
	  return(PEACOK.name)
}