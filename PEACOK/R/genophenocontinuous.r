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

run.geno.pheno.nonbinary <- function(genotypes, phenotypes, file.confounder, file.annotaiton = NULL, 
                                      ordered = FALSE, unordered = FALSE,
                                      geno_start = NULL, geno_end = NULL, pheno_start = NULL, pheno_end = NULL, 
                                      ignoreConfounder = FALSE, pheno_pheno = FALSE, verbose = 0) {
    if (unordered && ordered) {
        stop("Can only specify ordered or unordered but not both")
    }
  
    if (is(genotypes, "Matrix")) { #If it is  Matrix, we assume it is the varaint matrix from function variant.to.genotype.model
        samples <- colnames(genotypes)
        genos <- rownames(genotypes)
        genotypes <- as.data.frame(t(as.matrix(genotypes))) 
        genotypes <- genotypes %>% mutate(userID = samples) %>% select(userID, genos)
        rownames(genotypes) <- samples
        
    } else {
        # read a trunk of the collapsing matrix
        if (pheno_pheno) {
            genotypes <- read.pheno.matrix(genotypes, pheno_start = geno_start, pheno_end = geno_end)
        } else {
            suppressWarnings(genotypes <- read.geno.matrix(genotypes,gene_start = geno_start,gene_end = geno_end))
        }
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
    if (unordered) {
        Result <-data.frame(geno = character(0), pheno = character(0), Field = character(0),pairs = character(0),
                            nSamples = numeric(0),nref = numeric(0), nalt = numeric(0),
                            p = numeric(0), beta = numeric(0), lower = numeric(0), upper = numeric(0), stringsAsFactors = FALSE)
    } else {
        Result = data.frame(geno = character(0), pheno = character(0), Field = character(0), 
                            nSamples = numeric(0), nCases = numeric(0), nControls = numeric(0),
                            p = numeric(0), beta = numeric(0), lower = numeric(0), upper = numeric(0), se = numeric(0), 
                            MedCases = numeric(0), MedControls = numeric(0), stringsAsFactors = FALSE)
    }
    
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
                if (ordered) {
                    result <- run.categorical.ordered(pheno, geno, confounders)
                } else if(unordered) {
                    result <- run.categorical.unordered(pheno, geno, confounders)
                } else {
                    result <- run.continuous(pheno, geno, confounders)
                }
                Field <- NULL
                FieldID <- NULL
                if (!is.null(annotations)) {
                    FieldRow <- annotations %>% filter(FieldID == names(phenotypes)[pheno_index])
                    if (!is.null(FieldRow)) {
                        Field <- FieldRow$Field
                        Field <- gsub(",", "|", Field)
                    }
                }
                if (is.null(Field)) {
                    Field <- ""
                }
                if (unordered) {
                    result <- result <- mutate(geno = names(genotypes)[geno_index]) %>%
                                        mutate(pheno = names(phenotypes)[pheno_index]) %>%
                                        mutate(Field = Field) %>%
                                        select(names(Result))
                    Result <- Result %>% rbind(result)
                } else {
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
                }
            }, error = function(e) {
            }) 
            
        }
    }
    })
    Result <- Result %>%arrange(p)
    return(Result)
}