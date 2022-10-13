#### Implememting penalized B-spline method with aghq:
library(fda)

### Functions to create the given set of B spline basis:
## region: the interval of interest
## k: the number of basis functions
## p: the order of the IWP
create_Osullivan_basis <- function(region, k, p){
  B <- suppressWarnings(create.bspline.basis(rangeval = c(min(region),max(region)),
                                              nbasis = k,
                                              norder = (2*p)))
  B
}


### Functions to create design matrix:
## basis: output from create_Osullivan_basis
## x: the vector of observed locations
create_Osullivan_design <- function(basis, x){
  B <- basis
  Bmatrix <- eval.basis(x, B, Lfdobj=0, returnMatrix=T)
  colnames(Bmatrix) <- NULL
  as(Bmatrix, "dgTMatrix")
}


### Functions to create penalty matrix:
## basis: output from create_Osullivan_basis
## p: the order of the IWP
create_Osullivan_Q <- function(basis, p){
  B <- basis
  region <- basis$rangeval
  eval.penalty(basisobj = B, Lfdobj=int2Lfd(p), rng=region)
}


### Function to obtain the sample path of derivatives
## basis: output from create_Osullivan_basis
## weights: samples of regression weights, each column is one sample.
## x: the vector of observed locations
sample_Osullivan <- function(basis, weights, q = 0, x){
  B <- basis
  Bmatrix <- eval.basis(x, B, Lfdobj=q, returnMatrix=T)
  dim(Bmatrix)
  Bmatrix %*% weights
}


