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

run.geno.pheno.logistic <- function(genotypes, phenotypes, file.confounder, file.annotaiton = NULL, 
                                   geno_start = NULL, geno_end = NULL, pheno_start = NULL, pheno_end = NULL, ignoreConfounder = FALSE, 
                                  verbose = 0) {
    
    if (is(genotypes, "Matrix")) { #If it is  Matrix, we assume it is the varaint matrix from function variant.to.genotype.model
        samples <- colnames(genotypes)
        genos <- rownames(genotypes)
        genotypes <- as.data.frame(t(as.matrix(genotypes))) 
        genotypes <- genotypes %>% mutate(userID = samples) %>% select(userID, genos)
        rownames(genotypes) <- samples
        
    } else {
        # read a trunk of the collapsing matrix
        suppressWarnings(genotypes <- read.geno.matrix(genotypes,gene_start = geno_start,gene_end = geno_end))
    }
    # read a trunk of the phenotype matrix
    if (is(phenotypes, "character")) {
        phenotypes <- read.pheno.matrix(phenotypes, pheno_start = pheno_start, pheno_end = pheno_end)
    }
    # read confounders
    confounders <- fread(file.confounder, header=TRUE)
    
    if (class(genotypes$userID) != class(phenotypes$userID)) {
        suppressWarnings(genotypes$userID <- as(genotypes$userID, class(phenotypes$userID)))
    }
    # harmonize to common samples
    userID <- NULL #trick CRAN check
    suppressMessages(common.sample.IDs <- confounders %>% dplyr::select(userID) %>% inner_join(phenotypes %>% dplyr::select(userID)) %>%
                      inner_join(genotypes %>% dplyr::select(userID)) %>%  arrange(userID))
    
    suppressMessages(confounders <- confounders %>% inner_join(common.sample.IDs) %>% arrange(userID))
    suppressMessages(phenotypes <- phenotypes %>% inner_join(common.sample.IDs) %>% arrange(userID))
    suppressMessages(genotypes <- genotypes %>% inner_join(common.sample.IDs) %>% arrange(userID))
    
    # up to this point, we no longer need userIDs
    confounders <- confounders %>% select(-userID)
    phenotypes <- phenotypes  %>% select(-userID)
    genotypes <- genotypes  %>% select(-userID)
    
    if (!is.null(file.annotaiton)) {
      annotations <- fread(file.annotaiton, header = TRUE)
    } else {
      annotations <- NULL
    }
    
    p <- NULL
    system.time({
    Result = data.frame(geno = character(0), pheno = character(0), Field = character(0), nSamples = numeric(0),
                        nCases = numeric(0), nControls = numeric(0),
                        p = numeric(0), beta = numeric(0), lower = numeric(0), upper = numeric(0), 
                        se = numeric(0), MedCases = numeric(0), MedControls = numeric(0), stringsAsFactors = FALSE)
    
    #df <- cbind.data.frame(common.sample.IDs,genotypes,phenotypes,confounders)
    #write.table(df, file="test.csv",  append=FALSE, quote=FALSE, sep=",", na="", row.names=FALSE, col.names=TRUE)
    for (geno_index in 1: dim(genotypes)[2]) {
        geno <- genotypes[,geno_index]
        if (verbose == 1) {
            print(geno_index)
        }  
        
        for (pheno_index in 1: dim(phenotypes)[2]) {
            if (verbose == 2) {
                print(c(geno_index,pheno_index))
            }  
            pheno <- phenotypes[,pheno_index]
            if (ignoreConfounder) {
                confounders = NULL
            }
            tryCatch({
                result <- run.logistic(pheno, geno, confounders)
                Field <- NULL
                FieldID <- NULL
                if (!is.null(annotations)) {
                  FieldRow <- annotations %>% filter(FieldID == names(phenotypes)[pheno_index])
                  if (!is.null(FieldRow)) {
                    Field <- FieldRow$Field
                    Field <- gsub(",", "|", Field)
                  }
                }
                if (is.null(Field) || length(Field) ==0) {
                  Field <- ""
                }
                
                Result <- Result %>% rbind(data.frame(geno = names(genotypes)[geno_index],
                                                      pheno = names(phenotypes)[pheno_index],
                                                      Field = Field,
                                                      nSamples = result$nSamples,
                                                      nCases = result$nCases,
                                                      nControls = result$nControls,
                                                      p = result$p, 
                                                      beta = result$beta,
                                                      lower = result$lower,
                                                      upper = result$upper, 
                                                      se = result$se,
                                                      MedCases = result$MedCases,
                                                      MedControls = result$MedControls,
                                                      stringsAsFactors = FALSE))
            }, error = function(e) {
            }) 
            
        }
    }
    })
    Result <- Result %>%arrange(p)
    return(Result)
}

run.geno.pheno.logistic2 <- function(genotypes, phenotypes, file.confounder, file.annotaiton = NULL, 
                                    geno_start = NULL, geno_end = NULL, pheno_start = NULL, pheno_end = NULL, ignoreConfounder = FALSE, 
                                    verbose = 0) {
  
  if (is(genotypes, "Matrix")) { #If it is  Matrix, we assume it is the varaint matrix from function variant.to.genotype.model
    samples <- colnames(genotypes)
    genos <- rownames(genotypes)
    genotypes <- as.data.frame(t(as.matrix(genotypes))) 
    genotypes <- genotypes %>% mutate(userID = samples) %>% select(userID, genos)
    rownames(genotypes) <- samples
    
  } else {
    # read a trunk of the collapsing matrix
    suppressWarnings(genotypes <- read.geno.matrix(genotypes,gene_start = geno_start,gene_end = geno_end))
  }
  # read a trunk of the phenotype matrix
  if (is(phenotypes, "character")) {
    phenotypes <- read.pheno.matrix(phenotypes, pheno_start = pheno_start, pheno_end = pheno_end)
  }
  # read confounders
  confounders <- fread(file.confounder, header=TRUE)
  
  if (class(genotypes$userID) != class(phenotypes$userID)) {
    suppressWarnings(genotypes$userID <- as(genotypes$userID, class(phenotypes$userID)))
  }
  # harmonize to common samples
  userID <- NULL #trick CRAN check
  suppressMessages(common.sample.IDs <- confounders %>% dplyr::select(userID) %>% inner_join(phenotypes %>% dplyr::select(userID)) %>%
                     inner_join(genotypes %>% dplyr::select(userID)) %>%  arrange(userID))
  
  suppressMessages(confounders <- confounders %>% inner_join(common.sample.IDs) %>% arrange(userID))
  suppressMessages(phenotypes <- phenotypes %>% inner_join(common.sample.IDs) %>% arrange(userID))
  suppressMessages(genotypes <- genotypes %>% inner_join(common.sample.IDs) %>% arrange(userID))
  
  # up to this point, we no longer need userIDs
  confounders <- confounders %>% select(-userID)
  phenotypes <- phenotypes  %>% select(-userID)
  genotypes <- genotypes  %>% select(-userID)
  
  if (!is.null(file.annotaiton)) {
    annotations <- fread(file.annotaiton, header = TRUE)
  } else {
    annotations <- NULL
  }
  
  p <- NULL
  system.time({
    Result = NULL
    
    #df <- cbind.data.frame(common.sample.IDs,genotypes,phenotypes,confounders)
    #write.table(df, file="test.csv",  append=FALSE, quote=FALSE, sep=",", na="", row.names=FALSE, col.names=TRUE)
    for (geno_index in 1: dim(genotypes)[2]) {
      geno <- genotypes[,geno_index]
      if (verbose == 1) {
        print(geno_index)
      }  
      
      for (pheno_index in 1: dim(phenotypes)[2]) {
        if (verbose == 2) {
          print(c(geno_index,pheno_index))
        }  
        pheno <- phenotypes[,pheno_index]
        if (ignoreConfounder) {
          confounders = NULL
        }
        tryCatch({
          result <- run.logistic2(pheno, geno, confounders)
          Field <- NULL
          FieldID <- NULL
          if (!is.null(annotations)) {
            FieldRow <- annotations %>% filter(FieldID == names(phenotypes)[pheno_index])
            if (!is.null(FieldRow)) {
              Field <- FieldRow$Field
              Field <- gsub(",", "|", Field)
            }
          }
          if (is.null(Field) || length(Field) ==0) {
            Field <- ""
          }
          
          df <- data.frame(geno = names(genotypes)[geno_index],
                           pheno = names(phenotypes)[pheno_index],
                           Field = Field,
                           nSamples = result$nSamples,
                           nCases = result$nCases,
                           nControls = result$nControls,
                           stringsAsFactors = FALSE)
          names(result$p)[1] <- "p_intercept"
          names(result$beta)[1] <- "beta_intercept"
          for (p_index in 2: length(names(result$p))) {
              names(result$p)[p_index] <-  paste0("p_", names(result$p)[p_index])
              names(result$beta)[p_index] <-  paste0("beta_", names(result$beta)[p_index]) 
          }
          df <- df %>% cbind(as.data.frame(t(result$p)), as.data.frame(t(result$beta)))
          if (is.null(Result)) {
              Result <- df
          } else {
              Result <- Result %>% rbind(df)
          }
        }, error = function(e) {
        }) 
        
      }
    }
  })
  return(Result)
}