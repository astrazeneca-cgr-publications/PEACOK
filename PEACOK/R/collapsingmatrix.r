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

read.geno.matrix <- function(matrix_file, gene_start = NULL, gene_end = NULL ) {
    if (!file.exists(matrix_file)) {
        stop("geno matrix not found.")
    }
    
    # figure out number of rows
    genotype <- fread(matrix_file,header=TRUE, select = 1)
    total_genos <- dim(genotype)[1]
    if (!is.null(gene_end) && gene_end > total_genos) {
        gene_end <- total_genos
    }
    if ((is.null(gene_start) && !is.null(gene_end)) || (!is.null(gene_start) && is.null(gene_end))) {
        stop("geno block not defined correctly")
    }
    
    #header
    header <- fread(matrix_file,header=TRUE, nrows = 0)
    if (is.null(gene_start)) {
        geno_block <- fread(matrix_file,header=FALSE,  skip = 1) 
    } else {
        nrows <- gene_end - gene_start + 1
        geno_block <- fread(matrix_file,header=FALSE, nrows = nrows, skip = gene_start) 
    }
    
    genes <- unlist(geno_block[,1])
    genes <- gsub("\"","", genes)
    geno_block <- as.data.frame(t(geno_block[,-1]))
    geno_block[geno_block >0] <- 1
    geno_block <- cbind(names(header)[-1],geno_block)
    names(geno_block) <- c("userID", genes)
    row.names(geno_block) <- names(header)[-1]
    geno_block$userID <- as.numeric(as.character(geno_block$userID))
    return(geno_block)
}

read.pheno.matrix <- function(matrix_file, pheno_start = NULL, pheno_end = NULL ) {
    if (!file.exists(matrix_file)) {
        stop("pheno matrix not found.")
    }
    
    header <- fread(matrix_file,header=TRUE, nrows = 0)
    total_phenos <- length(header) - 1
    
    if ((is.null(pheno_start) && !is.null(pheno_end)) || (!is.null(pheno_start) && is.null(pheno_end))) {
        stop("pheno block not defined correctly")
    }
    
    if (is.null(pheno_start)) {
        pheno_block <- fread(matrix_file,header=TRUE) 
    } else {
        pheno_block <- fread(matrix_file,header=TRUE, select = c(1,(pheno_start+1):(pheno_end+1))) 
    }
    
    row.names(pheno_block) <- pheno_block$userID
    return(pheno_block)
    
}