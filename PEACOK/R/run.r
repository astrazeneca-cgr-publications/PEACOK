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


run <- function(opt) {
    input <- init.data(opt)
    print("LOADING DONE")
    if(!is.null(opt$extract.only) && opt$extract.only) {
        return(0)
    }
    #initilize package level variables, primarily for passing variables between functions and logging
    .initEnv(opt, input$data)
    
    currentVar <- ""
    currentVarShort <- ""
    first <- TRUE
    
    phenoIdx=0; # zero because then the idx is the position of the previous variable, i.e. the var in currentVar
    for (var in input$phenoVars) { 
        sink()
        sink(pkg.env$resLogFile, append=TRUE)
        
        varx <- gsub("^x", "", var)
        varx <- gsub("_[0-9]+$", "", varx)
        
        varxShort <- gsub("^x", "", var)
        varxShort <- gsub("_[0-9]+_[0-9]+$", "", varxShort)
        
        ## test this variable
        if (currentVar == varx) {
          thisCol <- input$data[,eval(var)]
          thisCol <- .replaceNaN(thisCol)
          currentVarValues <- cbind.data.frame(currentVarValues, thisCol)
        } else if (currentVarShort == varxShort) {
          ## different time point of this var so skip
        } else {
          ## new variable so run test for previous (we have collected all the columns now)
          if (first==FALSE) {
            thisdata <- .makeTestDataFrame(input$data, input$confounders, currentVarValues)
            test.associations(opt, input$vl, currentVar, currentVarShort, thisdata, input$phenoStartIdx)
          }
          first <- FALSE
          
          ## new variable so set values
          currentVar <- varx
          currentVarShort <- varxShort
          
          currentVarValues <- input$data[,eval(var)]
          currentVarValues <- .replaceNaN(currentVarValues)
        }
        phenoIdx <- phenoIdx + 1
    }
    
    if (phenoIdx>0){
        # last variable so test association
        thisdata = .makeTestDataFrame(input$data, input$confounders, currentVarValues)
        test.associations(opt, input$vl, currentVar, currentVarShort, thisdata, input$phenoStartIdx)
    }
    sink()
    
    # save counters of each path in variable flow
    .saveCounts(opt)
    if (opt$save == TRUE) {
        write.table(pkg.env$derivedBinary, file=paste(opt$resDir,"data-binary-",opt$varTypeArg,".csv", sep=""), append=FALSE, quote=FALSE, sep=",", na="", row.names=FALSE, col.names=TRUE);
        write.table(pkg.env$derivedCont, file=paste(opt$resDir,"data-cont-",opt$varTypeArg,".csv", sep=""), append=FALSE, quote=FALSE, sep=",", na="", row.names=FALSE, col.names=TRUE);
        write.table(pkg.env$derivedCatOrd, file=paste(opt$resDir,"data-catord-",opt$varTypeArg,".csv", sep=""), append=FALSE, quote=FALSE, sep=",", na="", row.names=FALSE, col.names=TRUE);
        write.table(pkg.env$derivedCatUnord, file=paste(opt$resDir,"data-catunord-",opt$varTypeArg,".csv", sep=""), append=FALSE, quote=FALSE, sep=",", na="", row.names=FALSE, col.names=TRUE);
    }
  
}