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

.DoBinary <- function(phenoFactor, genoFactor) {
    tryCatch({
        mylogit <- glm(phenoFactor ~ genoFactor, family="binomial")
        sumx <- summary(mylogit)
        p <- sumx$coefficients['genoFactor1','Pr(>|z|)']
        beta <- sumx$coefficients['genoFactor1','Estimate']
        suppressMessages(cis <- confint(mylogit, 'genoFactor1', level=0.95))
        lower <- cis["2.5 %"]
    	  upper <- cis["97.5 %"]
    ## END TRYCATCH
    }, error = function(e) {
        p <- NA
        beta <- NA
        lower <- NA
        upper <- NA
    }) 
    return(list(p = p, beta = beta, lower = lower, upper = upper))
}

run.logistic2 <- function(pheno, geno, confounders) {
  
  ###### BEGIN TRYCATCH
  tryCatch({
    sumx <- summary(glm( pheno ~ geno + ., data=confounders, family = binomial(link = logit)))
    if ("geno" %in% row.names(sumx$coefficients)) {
      pvalue <- sumx$coefficients[,'Pr(>|z|)']
      beta <- sumx$coefficients[,"Estimate"]
    } else {
      pvalue <- NA
      beta <- NA
    }
    NA_geno_pheno <- !is.na(pheno) & !is.na(geno)
    numNotNA <- length(which(NA_geno_pheno))
    numCases <- length(which(NA_geno_pheno & geno>0))
    numControls <- length(which(NA_geno_pheno & geno==0))
    
    ## END TRYCATCH
  }, error = function(e) {
  })
  return(list(nSamples = numNotNA, nCases = numCases, nControls = numControls, p = pvalue, beta = beta))
}

run.logistic <- function(pheno, geno, confounders) {
  ###### BEGIN TRYCATCH
  tryCatch({
    sumx <- summary(glm( pheno ~ geno + ., data=confounders, family = binomial(link = logit)))
    if ("geno" %in% row.names(sumx$coefficients)) {
      pvalue <- sumx$coefficients['geno','Pr(>|z|)']
      beta <- sumx$coefficients["geno","Estimate"]
      se <- sumx$coefficients["geno","Std. Error"]
      lower <- NA
      upper <- NA
    } else {
      pvalue <- NA
      beta <- NA
      lower <- NA
      upper <- NA
      se <- NA
    }
    NA_geno_pheno <- !is.na(pheno) & !is.na(geno)
    numNotNA <- length(which(NA_geno_pheno))
    numCases <- length(which(NA_geno_pheno & geno>0))
    numControls <- length(which(NA_geno_pheno & geno==0))
    if (numCases > 0) {
      MedCases <- median(pheno[NA_geno_pheno & geno>0])
    } else {
      MedCases <- NA
    }
    if (numControls > 0) {
      MedControls <- median(pheno[NA_geno_pheno & geno==0])
    } else {
      MedControls <- NA
    }
    
    ## END TRYCATCH
  }, error = function(e) {
  })
  return(list(nSamples = numNotNA, nCases = numCases, nControls = numControls,
              p = pvalue, beta = beta, lower = lower, upper = upper, se = se,
              MedCases = MedCases, MedControls = MedControls))
}

run.continuous <- function(pheno, geno, confounders) {
    ###### BEGIN TRYCATCH
    tryCatch({
        if (is.null(confounders)) {
            fit <- lm(pheno ~ geno)
        } else {
            fit <- lm(pheno ~ geno + ., data=confounders)
        }
        sumx <- summary(fit)
        if ("geno" %in% row.names(sumx$coefficients)) {
              pvalue <- sumx$coefficients['geno','Pr(>|t|)']
              beta <- sumx$coefficients["geno","Estimate"]
              se <- sumx$coefficients["geno","Std. Error"]
              cis <- confint(fit, level=0.95)
	            lower <- cis["geno", "2.5 %"]
              upper <- cis["geno", "97.5 %"]
        } else {
              pvalue <- NA
              beta <- NA
	            lower <- NA
              upper <- NA
              se <- NA
        }
        NA_geno_pheno <- !is.na(pheno) & !is.na(geno)
        
        numNotNA <- length(which(NA_geno_pheno))
        numCases <- length(which(NA_geno_pheno & geno>0))
        numControls <- length(which(NA_geno_pheno & geno==0))
        if (numCases > 0) {
            MedCases <- median(pheno[NA_geno_pheno & geno>0])
        } else {
            MedCases <- NA
        }
        if (numControls > 0) {
            MedControls <- median(pheno[NA_geno_pheno & geno==0])
        } else {
            MedControls <- NA
        }

    ## END TRYCATCH
    }, error = function(e) {
    })
    return(list(nSamples = numNotNA, nCases = numCases, nControls = numControls,
                p = pvalue, beta = beta, lower = lower, upper = upper, se = se,
                MedCases = MedCases, MedControls = MedControls))
}

run.categorical.ordered <- function(pheno, geno, confounders) {
    ###### BEGIN TRYCATCH
    tryCatch({
        numNotNA <- length(which(!is.na(pheno)))
        numCases <- length(which(!is.na(pheno) & geno>0))
        #numControls <- numNotNA - numCases
        numControls <- length(which(!is.na(pheno) & geno==0))
        if (numCases > 0) {
            MedCases <- median(pheno[!is.na(pheno) & geno>0])
        } else {
            MedCases <- NA
        }
        if (numControls > 0) {
            MedControls <- median(pheno[!is.na(pheno) & geno==0])
        } else {
            MedControls <- NA
        }
        
        phenoFactor <- factor(pheno)
        fit <- polr(phenoFactor ~ geno + ., data=confounders, Hess=TRUE)
		    ctable <- coef(summary(fit))
		    ct <- coeftest(fit)
		    if ("geno" %in% row.names(ct)) {
  		    pvalue <- ct["geno","Pr(>|t|)"]
  		    beta <- ctable["geno", "Value"]
  		    se <- ctable["geno", "Std. Error"]
  	      lower <- beta - 1.96*se
          upper <- beta + 1.96*se
          if (pvalue == 0) { # ill model
              pvalue <- NA
              beta <- NA
  		        se <- NA
  	          lower <- NA
              upper <- NA
          }
		    } else {
		       pvalue <- NA
		       beta <- NA
		       se <- NA
	         lower <- NA
           upper <- NA
		    }	    
    ## END TRYCATCH
    }, error = function(e) {
    })
    return(list(nSamples = numNotNA, nCases = numCases, nControls = numControls,
                p = pvalue, beta = beta, lower = lower, upper = upper, se = se,
                MedCases = MedCases, MedControls = MedControls))
}

