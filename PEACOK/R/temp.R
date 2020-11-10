run.categorical.unordered <- function(pheno, geno, confounders) {
  ###### BEGIN TRYCATCH
  tryCatch({
    numNotNA <- length(which(!is.na(pheno)))
    nref <- NULL
    nalt <- NULL
    
    
    # check there are not too many levels and skip if there are
    numUnique <- length(unique(na.omit(pheno)))
    # num outcome values * (num confounders and trait of interest and bias term)
    if (numUnique > opt$maxunorderedcategories) {
      cat("Too many levels: ", numUnique, " > ", opt$maxunorderedcategories, 
          "(num outcomes values: ", numUnique, ") || SKIP ", sep="")
      .incrementCounter("unordCat.cats")
      return(NULL)
    }
    
    phenoFactor <- choose.reference.category(pheno)
    reference <- levels(phenoFactor)[1]
    
    fit <- multinom(phenoFactor ~ geno + ., data=confounders, maxit=1000)
    ## baseline model with only confounders, to which we compare the model above
    fitB <- multinom(phenoFactor ~ ., data=confounders, maxit=1000)
    
    ## compare model to baseline model
    lres <- lrtest(fit, fitB)
    modelP <- lres[2,"Pr(>Chisq)"]
    
    ## save result to file
    maxFreq <- length(which(phenoFactor==reference))
    numNotNA <- length(which(!is.na(pheno)))
    output <- data.frame(pairs = paste(reference,"#lrtest",sep=""), 
                         nSamples = numNotNA, nref = maxFreq, nalt = "", 
                         beta = -999, lower = -999, upper = -999, p = modelP,stringsAsFactors = FALSE)
    sumx <- summary(fit)
    z <- sumx$coefficients/sumx$standard.errors
    p <- (1 - pnorm(abs(z), 0, 1))*2			
    ci <- confint(fit, "geno", level=0.95)
    ci <- data.frame(ci)
    
    ## get result for each variable category
    uniqVar <- unique(na.omit(pheno))
    for (u in uniqVar) {
      ## no coef for baseline value, and values <0 are assumed to be missing
      if (u == reference || u<0) {
        next
      }
      pvalue <- p[paste(eval(u),sep=""),"geno"]				
      beta <- sumx$coefficients[paste(eval(u),sep=""),"geno"]
      
      if (opt$confidenceintervals == TRUE) { 
        lower <- ci[1, paste("X2.5...", u, sep="")]
        upper <-	ci[1, paste("X97.5...", u, sep="")]
      } else {
        lower <- NA
        upper <- NA
      }
      numThisValue <- length(which(phenoFactor==u))
      ## save result to file
      output <- rbind(output, data.frame(pairs = paste(reference,"#",u,sep=""), 
                                         nSamples = numNotNA, nref = maxFreq, nalt = numThisValue,
                                         beta = beta, lower = lower, upper = upper, p = pvalue,stringsAsFactors = FALSE))
    }
  }, error = function(e) {
  })
  return(output)
}