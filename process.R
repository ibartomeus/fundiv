#what I do to create a package, based in :https://github.com/jtleek/rpackages and
  #http://adv-r.had.co.nz/Package-basics.html

## Setup
#install.packages(c("devtools", "roxygen2", "knitr"))

## Load the libraries
library("devtools")
library("roxygen2")
library("knitr")

##To read roxigen comments and create the man folder
document("fundiv")

##to do some unit testing
#http://ivory.idyll.org/blog/automated-testing-and-research-software.html
#http://www.johnmyleswhite.com/notebook/2010/08/17/unit-testing-in-r-the-bare-minimum/

library(testthat)
library(fundiv)
test_package("fundiv")
# R CMD check

test_file("test-defunction.R")

#test2
source('sample.R')
test_dir('tests', reporter = 'Summary')
