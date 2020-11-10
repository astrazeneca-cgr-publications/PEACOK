read.variant.annotations <- function(matrix_file, output.file, annotation.columns = c(1,3,4)) {
    annotations <- fread(matrix_file, na.strings = c("NA",""), stringsAsFactors = FALSE, sep = "\t", header = TRUE, select = annotation.columns)
    write.table(annotations,file = output.file, append=FALSE, quote=FALSE, sep=",", na="", row.names=FALSE, col.names=TRUE)
}

partition.variant.matrix <- function(matrix_file, output.folder, trunk.size = 1000) {
  is.gz <- grepl("tsv.gz$",matrix_file)
  if (is.gz) {
      #system(paste0("gunzip -k ", matrix_file)) not all gunzip support -k option
      old_matrix_file <- matrix_file
      matrix_file <- gsub("tsv.gz", "tsv", matrix_file)
      system(paste0("gunzip -c ", old_matrix_file, " > ", matrix_file))
  }
  ## figure out the number of columns/samples
  columns <- fread(matrix_file, na.strings = c("NA",""), stringsAsFactors = FALSE, sep = "\t", header = TRUE, nrows = 0)
  if (length(names(columns)) < 5 || !(names(columns)[4] == "Most Damaging Effect")) { #magic word here 
      stop("Unrecognized varaint matrix format")
  }
  columns <- names(columns[,-c(1,3,4)])
  
  not.done <- TRUE
  skip <- 1
  last_block_size <- 10 #a number greater than zero
  trunk.id <- 1
  while (not.done && last_block_size > 0) {
      print(paste0("Block ID: ", trunk.id))
      last_block_size <- 0
      tryCatch({
        block <- fread(matrix_file, na.strings = c("NA",""), stringsAsFactors = FALSE, sep = "\t", header = FALSE, skip = skip, nrows = trunk.size)
        last_block_size <- dim(block)[1]
        if (last_block_size > 0) { #at least one row read
            variants <- unlist(block[,1])
            var.matrix <- as(as.matrix(block[,-c(1,3,4)]), "sparseMatrix") 
            colnames(var.matrix) <- columns
            rownames(var.matrix) <- variants
            save(var.matrix, file = file.path(output.folder, paste0("trunk_", trunk.id, ".RData")), version = "2")
            trunk.id <- trunk.id + 1
            skip <- skip + last_block_size
            if (last_block_size < trunk.size) {
              not.done <- FALSE
            }
        } 
      }, error = function(e) {
      }) 
  }
  read.variant.annotations(matrix_file, file.path(output.folder, "annotations.csv"))
  if (is.gz) {
      file.remove(matrix_file)
  }
}

## this is the easy to understand but slow implementation
.variant.to.genotype.model.original <- function(var.matrix, model, verbose = FALSE) {
  
    models <- c("dominant", "recessive", "genotypic") # only deal with these three for now
    
    # constants from CAF
    REF <- 0
    HET <- 1
    HOM <- 2
    REF_MALE <- 3
    HOM_MALE <- 4
    
    MISSING <- -1
    MISSING2 <- 5
    
    if (!model %in% models ) {
        stop(paste0("Genetic model " ,model, " is not recognized!"))
    }
    
    var.matrix[var.matrix == MISSING | var.matrix == MISSING2] = NA
    
    transform <- function(variant, model) {
        is.Minor.allele <- variant[1]
        
        v.new <- variant
        if (model == "dominant") {
            if (is.Minor.allele > 0) { 
                v.new[variant == REF | variant == REF_MALE | variant == HET] <- 1
                v.new[variant == HOM | variant == HOM_MALE] <- 0
            } else {
                v.new[variant == HOM | variant == HOM_MALE | variant == HET] <- 1
                v.new[variant == REF | variant == REF_MALE] <- 0
            }
        } else if (model == "recessive") {
            if (is.Minor.allele > 0) {
                v.new[variant == REF | variant == REF_MALE] <- 1
                v.new[variant == HOM | variant == HOM_MALE | variant == HET] <- 0
            } else {
                v.new[variant == HOM | variant == HOM_MALE] <- 1
                v.new[variant == REF | variant == REF_MALE | variant == HET] <- 0
            }
        } else if (model == "genotypic") {
            if (is.Minor.allele > 0) {
                v.new[variant == HOM | variant == HOM_MALE] <- 0
                v.new[ variant == HET] <- 1
                v.new[variant == REF | variant == REF_MALE] <- 2
            } else {
                v.new[variant == REF | variant == REF_MALE] <- 0
                v.new[variant == HET] <- 1
                v.new[variant == HOM | variant == HOM_MALE] <- 2
            }
        }
        v.new[1] <- is.Minor.allele
        return(v.new)
    }
    
    n.variants <- dim(var.matrix)[1]
    for (v in 1:n.variants) {
        if (verbose) {
          print(v)
        }
        var.matrix[v,] <- transform(var.matrix[v,],model)
    }
    return(var.matrix)
}

variant.to.genotype.model <- function(var.matrix, model) {
  
  models <- c("dominant", "recessive", "genotypic") # only deal with these three for now
  
  # constants from CAF
  REF <- 0
  HET <- 1
  HOM <- 2
  REF_MALE <- 3
  HOM_MALE <- 4
  
  MISSING <- -1
  MISSING2 <- 5
  
  if (!model %in% models ) {
    stop(paste0("Genetic model " ,model, " is not recognized!"))
  }
  
  var.matrix[var.matrix == MISSING | var.matrix == MISSING2] = NA
  
  minor <- which(var.matrix[,1] >0)
  if (length(minor)> 0) {
      var.matrix.minor <- var.matrix[minor,]
      var.matrix.minor <- as.matrix(var.matrix.minor)
      if (length(minor) == 1) {
          dim(var.matrix.minor) <- c(1,length(var.matrix.minor))
          rownames(var.matrix.minor)[1] <- rownames(var.matrix)[minor]
      }
      v.minor <- var.matrix.minor
      if (model == "dominant") {
          v.minor[var.matrix.minor == REF | var.matrix.minor == REF_MALE | var.matrix.minor == HET] <- 1
          v.minor[var.matrix.minor == HOM | var.matrix.minor == HOM_MALE] <- 0
      } else if (model == "recessive") {
          v.minor[var.matrix.minor == REF | var.matrix.minor == REF_MALE] <- 1
          v.minor[var.matrix.minor == HOM | var.matrix.minor == HOM_MALE | var.matrix.minor == HET] <- 0
      } else if (model == "genotypic") {
          v.minor[var.matrix.minor == HOM | var.matrix.minor == HOM_MALE] <- 0
          v.minor[ var.matrix.minor == HET] <- 1
          v.minor[var.matrix.minor == REF | var.matrix.minor == REF_MALE] <- 2
      }
      v.minor[,1] <- 1
      rownames(v.minor) <- unname(rownames(v.minor))
      rm(var.matrix.minor)
  } else {
      v.minor <- NULL
  }
  major <- which(var.matrix[,1] ==0)
  if (length(major)> 0) {   #### Do not delete the comments below, they arw hints for optimization
      var.matrix.major <- var.matrix[major,]
      v.major <- var.matrix.major
      if (model == "dominant") {
          #v.major[var.matrix.major == HOM | var.matrix.major == HOM_MALE | var.matrix.major == HET] <- 1
          v.major[var.matrix.major == HOM | var.matrix.major == HOM_MALE] <- 1 #no HET as it is already 1
          #v.major[var.matrix.major == REF | var.matrix.major == REF_MALE] <- 0
          v.major[var.matrix.major == REF_MALE] <- 0 #no REF, as it is already zero
      } else if (model == "recessive") {
          v.major[var.matrix.major == HOM | var.matrix.major == HOM_MALE] <- 1 
          #v.major[var.matrix.major == REF | var.matrix.major == REF_MALE | var.matrix.major == HET] <- 0
          v.major[var.matrix.major == REF_MALE | var.matrix.major == HET] <- 0  #no REF, as it is already zero
      } else if (model == "genotypic") {
          #v.major[var.matrix.major == REF | var.matrix.major == REF_MALE] <- 0 
          v.major[var.matrix.major == REF_MALE] <- 0 #no REF, as it is already zero
          #v.major[var.matrix.major == HET] <- 1     #no HET as it is already 1
          #v.major[var.matrix.major == HOM | var.matrix.major == HOM_MALE] <- 2
          v.major[var.matrix.major == HOM_MALE] <- 2 #no HOM, as it is already 2
      }
      v.major[,1] <- 0
      rownames(v.major) <- unname( rownames(v.major))
      rm(var.matrix.major)
  } else {
      v.major = NULL
  }
  return(rbind(v.minor,v.major))
  #return(list(var.matrix.minor = v.minor, var.matrix.major = v.major))
}

genotype.model.to.geno <- function(genotype.model.data, geno_start = NULL, geno_end = NULL) {
    total_genos <- dim(genotype.model.data)[1]
    if (!is.null(geno_end) && geno_end > total_genos) {
      gene_end <- total_genos
    }
    if ((is.null(geno_start) && !is.null(geno_end)) || (!is.null(geno_start) && is.null(geno_end))) {
      stop("geno block not defined correctly")
    }
    
    if (!is.null(geno_start)) {
        genotype.model.data <- genotype.model.data[geno_start:geno_end,]
    }
    samples <- colnames(genotype.model.data)
    genos <- rownames(genotype.model.data)
    genotype.model.data <- as.data.frame(t(as.matrix(genotype.model.data))) 
    
    userID <- NULL
    genotype.model.data <- genotype.model.data %>% mutate(userID = samples) %>% select(userID, genos)
    rownames(genotype.model.data) <- samples
    return(genotype.model.data)
}