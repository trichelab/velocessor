# velocessor

[![Build Status](https://travis-ci.org/trichelab/velocessor.png?branch=master)](https://travis-ci.org/trichelab/velocessor)  [![codecov](https://codecov.io/gh/trichelab/velocessor/branch/master/graph/badge.svg)](https://codecov.io/gh/trichelab/velocessor)

## How to finish setting up your new package

Now that you've got a working package skeleton, there are a few steps to finish setting up all the integrations:

### 1. Git(Hub)

Go to https://github.com/trichelab and create a new repository. Then, in the directory where this package is, create your git repository from the command line, add the files, and push it to GitHub:

    git init
    git add --all
    git commit -m "Initial commit of package skeleton"
    git remote add origin git@github.com:trichelab/velocessor.git
    git push -u origin master

### 2. Travis

Now you can go to [Travis](https://travis-ci.org/profile/trichelab) and turn on continuous integration for your new package. You may need to click the "Sync account" button to get your new package to show up in the list.

If you have a codecov.io account, running your tests on Travis will trigger the code coverage job. No additional configuration is necessary

## Installing

<!-- If you're putting `velocessor` on CRAN, it can be installed with

    install.packages("velocessor") -->

The pre-release version of the package can be pulled from GitHub using the [devtools](https://github.com/hadley/devtools) package:

    # install.packages("devtools")
    devtools::install_github("trichelab/velocessor", build_vignettes=TRUE)

## For developers

The repository includes a Makefile to facilitate some common tasks.

### Running tests

`$ make test`. Requires the [testthat](https://github.com/hadley/testthat) package. You can also specify a specific test file or files to run by adding a "file=" argument, like `$ make test file=logging`. `test_package` will do a regular-expression pattern match within the file names. See its documentation in the `testthat` package.

### Updating documentation

`$ make doc`. Requires the [roxygen2](https://github.com/klutometis/roxygen) package.
