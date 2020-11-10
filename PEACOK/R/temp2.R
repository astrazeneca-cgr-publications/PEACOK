run.geno.pheno.nonbinary2 <- function(genotypes, phenotypes, file.confounder, file.annotaiton = NULL, 
                                      geno_start = NULL, geno_end = NULL, pheno_start = NULL, pheno_end = NULL, 
                                      ignoreConfounder = FALSE, verbose = 0) {
  genotypes <- t(genotypes)
  # read a trunk of the phenotype matrix
  if (is(phenotypes, "character")) {
    phenotypes <- read.pheno.matrix(phenotypes, pheno_start = pheno_start, pheno_end = pheno_end)
  }
  # read confounders
  confounders <- fread(file.confounder, header=TRUE)
  userID <- rownames(genotypes)
  if (class(userID) != class(phenotypes$userID)) {
    suppressWarnings(userID <- as(userID, class(phenotypes$userID)))
  }
  
  # harmonize to common samples
  common.sample.IDs <- intersect(intersect(confounders$userID, phenotypes$userID),userID)
 
  suppressMessages(confounders <- confounders %>% filter(userID %in% common.sample.IDs) %>% arrange(userID))
  suppressMessages(phenotypes <- phenotypes %>% filter(userID %in% common.sample.IDs) %>% arrange(userID))
  suppressMessages(genotypes <- genotypes[as.character(phenotypes$userID),])
  
  # up to this point, we no longer need userIDs
  confounders <- confounders %>% select(-userID)
  phenotypes <- phenotypes  %>% select(-userID)

  if (!is.null(file.annotaiton)) {
    annotations <- fread(file.annotaiton, header = TRUE)
  } else {
    annotations <- NULL
  }
  
  p <- NULL
  system.time({
    
    Result = data.frame(geno = character(0), pheno = character(0), Field = character(0), 
                        nSamples = numeric(0), nCases = numeric(0), nControls = numeric(0),
                        p = numeric(0), beta = numeric(0), lower = numeric(0), upper = numeric(0), se = numeric(0), 
                        MedCases = numeric(0), MedControls = numeric(0), stringsAsFactors = FALSE)
    
    
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
          
          result <- run.continuous(pheno, geno, confounders)
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
          
          Result <- Result %>% rbind(data.frame(geno = colnames(genotypes)[geno_index],
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