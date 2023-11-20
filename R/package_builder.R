library(Rcpp)
library(RcppArmadillo)
library(microbenchmark)
library(inline)
library(BH)
library(doParallel)
library(parallel)


sourceCpp("D:/OneDrive - rwth-aachen.de/Studium/Master/Projekt ISW/rcpp/package/misc.cpp")


# Install package with rcppArmadillo input


RcppArmadillo.package.skeleton(name = "EMcensoring")

# Then copy cpp file in Package/source directory

setwd("C:/Users/Cquec/Documents/EMcensoring")
Rcpp::compileAttributes()

# Edit DESCRIPTION file, add BH to Imports and Linking to

# I think this is optional with C++
# tools::package_native_routine_registration_skeleton(dir = "C:/Users/Cquec/Documents", character_only = TRUE)

# Set wd back to default
setwd("C:/Users/Cquec/Documents")

# Execute in Terminal
# R CMD build PackageDirectory/PackageSkeletonName 

install.packages("EMcensoring", repos = NULL, type="source")