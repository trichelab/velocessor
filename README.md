# velocessor

[![Build Status](https://travis-ci.org/trichelab/velocessor.png?branch=master)](https://travis-ci.org/trichelab/velocessor)  [![codecov](https://codecov.io/gh/trichelab/velocessor/branch/master/graph/badge.svg)](https://codecov.io/gh/trichelab/velocessor)

The pre-release version of the package can be pulled from GitHub using the [devtools](https://github.com/hadley/devtools) package:

    install.packages("BiocManager") # from CRAN
    BiocManager::install("trichelab/velocessor")

## Super quickstart

### 0. Ram everything through Salmon/Alevin and tximeta

See the [Alevin velocity tutorial](https://combine-lab.github.io/alevin-tutorial/2020/alevin-velocity/) for more details on this.

### 1. Load annotations

    library(AnnotationHub)
    ah <- AnnotationHub()
    # query(ah, c("ensdb", "Ensembl 99", "Homo sapiens") )
    ens99 <- ah[["AH78783"]]
    gencode33 <- genes(ens99, columns=c("symbol", "gene_biotype"))

### 2. Load quantifications 

    library(velocessor)
    alevin_suffix <- "_alevin"
    runs <- list.files(pattern=paste0(alevin_suffix, "$"))
    names(runs) <- sub(alevin_suffix, "", runs)
    qm <- "alevin/quants_mat.gz"
    txstub <- "gencode.v33.annotation.withintrons.expanded"
    txis <- process_velo_txis(runs, txstub, anno=gencode33, 
                              HARMONY=TRUE, SCVELO=TRUE, 
                              BPPARAM=MulticoreParam(3))

### 3. Plot the results (optionally on rescaled axes)

    txis <- rescale_dimred(txis, "UMAP") 
    plot_velo(txis, embed="scaled_UMAP", sizeref=3)

### 4. Explore the results

plot_velo generates an [interactive plot colored by velocity pseudotime](https://trichelab.github.io/CHLA9_CHLA10_MSCs_pseudotime/) by default (above) .

Try zooming in and out with the middle mouse button, hovering, and rotating.



## For developers

The repository includes a Makefile to facilitate some common tasks.

### Running tests

`$ make test`. Requires the [testthat](https://github.com/hadley/testthat) package. You can also specify a specific test file or files to run by adding a "file=" argument, like `$ make test file=logging`. `test_package` will do a regular-expression pattern match within the file names. See its documentation in the `testthat` package.

### Updating documentation

`$ make doc`. Requires the [roxygen2](https://github.com/klutometis/roxygen) package.
