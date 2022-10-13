###### Suppress cat output:
quiet <- function(x) { 
  sink(tempfile()) 
  on.exit(sink()) 
  invisible(force(x)) 
} 


############################################
###########################################
### For precision/design construction:
compute_Dp <- function(p,dif){
  n <- length(dif)
  D <- matrix(0, nrow = n, ncol = n)
  for (j in 1:n) {
    for (i in 1:j) {
      if(i == j) {D[j,i] <- (dif[i]^p)/(factorial(p))}
      else {
        k <- 1:p
        D[j,i]<- sum((dif[i]^k)*(sum(dif[(i+1):j])^(p-k))/(factorial(k)*factorial(p-k)))
      }
    }
  }
  as(D,"dtCMatrix")
}
compute_cov_piece_cons_RWp <- function(svec, glob_prec, sd_Wt = 1){
  p <- length(glob_prec)
  d <- diff(svec)
  D <- compute_Dp(p,d)
  T_matrix <- NULL
  for (i in 0:(p-1)) {
    T_matrix <- cbind(T_matrix,svec[-1]^i)
  }
  T_matrix <- cbind(T_matrix,D)
  Sigweights <- diag(c(1/glob_prec,sd_Wt/d))
  T_matrix %*% Sigweights %*% t(T_matrix)
}

### For Precision matrix of the exact method: actual space
compute_Wp_cov <- function(svec, p = 2){
  n <- length(svec)
  result <- matrix(0, ncol = n, nrow = n)
  for (i in 1:n) {
    for (j in i:n) {
      ele <- 0:(p - 1)
      elements <- ((-1)^(p - 1 - ele)) * choose((2*p - 1),ele) * (svec[i]^(2*p - 1 - ele)) * (svec[j]^ele)
      term <- sum(elements)/factorial(2*p - 1)
      result[i,j] <- term 
    }
  }
  result <- Matrix::forceSymmetric(result)
  result
}

### For Design/Precision of Proposed O-Spline
compute_weights_precision <- function(svec){
  d <- diff(svec)
  Precweights <- diag(d)
  Precweights
}
compute_weights_design <- function(svec,p){
  d <- diff(svec)
  D <- compute_Dp(p,d)
  T_matrix <- NULL
  for (i in 0:(p-1)) {
    T_matrix <- cbind(T_matrix,svec[-1]^i)
  }
  T_matrix <- cbind(T_matrix,D)
  T_matrix
}


##### When there are different initial conditions:
##### Implementing the middle-condition inference using O-Spline


### Rotating a matrix 90 degrees:
rotate <- function(x) t(apply(x, 2, rev))


### Obtain local design matrix:
condition_local_design <- function(xvec, s = NULL, p){
  if(is.null(s)){
    s <- ceiling(length(xvec)/2)
  }
  vec1 <- xvec[s] - xvec[(s):1]
  vec2 <- xvec[(s):length(xvec)] - xvec[s]
  B1 <- as.matrix(compute_Dp(dif = c(0,diff(sort(vec1))), p = p))
  B1 <- rotate(rotate(B1))
  B2 <- compute_Dp(dif = diff(sort(vec2)), p = p)
  B <- bdiag(B1,B2)
  B[,-s]
}

## Example:
# x <- c(1,2,3,4,5)
# condition_local_design(x,s = 3,2)



### Obtain global basis matrix:
condition_global_design <- function(xvec, s = NULL, p){
  if(is.null(s)){
    s <- ceiling(length(xvec)/2)
  }
  vec1 <- xvec[1:(s-1)] - xvec[s]
  vec2 <- xvec[(s+1):length(xvec)] - xvec[s]
  T_matrix <- NULL
  for (i in 0:(p-1)) {
    T_matrix <- cbind(T_matrix,c(vec1,0,vec2)^i)
  }
  T_matrix
}
# ## Example:
# x <- c(1,2,3,4,5)
# condition_global_design(x,s = 3,2)




### Evaluate the local basis at a refined grid with particular condition:

condition_local_poly <- function(knots, refined_x, p, s = NULL){
  if(is.null(s)){
    s <- ceiling(length(knots)/2)
  }
  knots1 <- knots[s] - knots[s:1]
  knots2 <- knots[(s):length(knots)] - knots[s]
  refined_x_part1 <- sort((knots[s] - refined_x)[(knots[s]-refined_x) >= 0])
  refined_x_part2 <- sort((refined_x - knots[s])[(knots[s]-refined_x) < 0])
  
  Part1 <- local_poly(knots1, refined_x = refined_x_part1, p = p)
  Part1 <- rotate(rotate(Part1))
  
  Part2 <- local_poly(knots2, refined_x = refined_x_part2, p = p)
  
  all_result <- bdiag(Part1,Part2)
  all_result
}

condition_local_poly_for_deriv <- function(knots, refined_x, p, s = NULL, degree){
  if(is.null(s)){
    s <- ceiling(length(knots)/2)
  }
  knots1 <- knots[s] - knots[s:1]
  knots2 <- knots[(s):length(knots)] - knots[s]
  refined_x_part1 <- sort((knots[s] - refined_x)[(knots[s]-refined_x) >= 0])
  refined_x_part2 <- sort((refined_x - knots[s])[(knots[s]-refined_x) < 0])
  
  Part1 <- local_poly(knots1, refined_x = refined_x_part1, p = p)
  Part1 <- ((-1)^degree)* rotate(rotate(Part1))
  
  Part2 <- local_poly(knots2, refined_x = refined_x_part2, p = p)
  
  all_result <- bdiag(Part1,Part2)
  all_result
}

condition_global_poly <- function(knots, refined_x, p, s = NULL){
  if(is.null(s)){
    s <- ceiling(length(knots)/2)
  }
  condition_value <- knots[s]
  s <- which.min(abs(refined_x - knots[s]))
  condition_global_design(refined_x, s = s, p = p)
}

# condition_local_poly(knots = c(1,2,3,4,5), refined_x = seq(1,5,by = 0.1), s = 3, p = 2)




### Obtain the precision matrix:
condition_precision <- function(xvec, s = NULL, p){
  if(is.null(s)){
    s <- ceiling(length(xvec)/2)
  }
  vec1 <- xvec[s] - xvec[s:1] 
  vec2 <- xvec[s:length(xvec)] - xvec[s]
  Q1 <- compute_weights_precision(vec1)
  Q1 <- rotate(rotate(Q1))
  Q2 <- compute_weights_precision(vec2)
  bdiag(Q1,Q2)
}


### For Design/Precision of Proposed B-Spline
get_knots_Bspline <- function(q, min, max, m){
  boundary <- c(min,max)
  inner <- seq(boundary[1],boundary[2], length.out = (q - m - 1 + 2))
  knots <- sort(c(boundary, inner[-c(1,(q - m - 1 + 2))]))
  trival <- rep(boundary, each = m)
  all_knots <- sort(c(trival, knots))
}
get_design_Bspline <- function(x,q,min,max,p){
  m <- 2*p - 1
  all_knots <- get_knots_Bspline(q, min, max, m = m)
  bs3 <- splineDesign(x = x,knots = all_knots, ord = (m+1), outer.ok = F) 
  ### Removing those basis functions that won't be used
  bs3[,-c(1:p)]
}
get_precision_Bspline <- function(q,min,max,p){
  m <- 2*p - 1
  ord <- m + 1
  all_knots <- get_knots_Bspline(q, min, max, m = m)
  inner_knots <- all_knots[ord:(length(all_knots) - (ord-1))]
  object <- s(x,k = q, bs = "bs",  m = c(m,p))
  data <- data.frame(x = inner_knots)
  SS <- smoothCon(object = object, data = data, knots = data.frame(x = all_knots))
  Q <- SS[[1]]$S[[1]][-c(1:p),-c(1:p)]
  Q
}

### For Precision matrix using RW2 (2008)
compute_H_rue <- function(d,n){
  H <- matrix(data = 0, nrow = n, ncol = n)
  for (i in 2:(nrow(H)-1)) {
    H[i,i] <- -(1/d[i-1]+1/d[i])
    H[i,i-1] <- 1/d[i-1]
    H[i,i+1] <- 1/d[i]
  }
  H
}
compute_B <- function(d,n){
  B <-matrix(0, nrow = n, ncol = n)
  B[1,1] <- d[1]/3
  B[1,2] <- d[1]/6
  B[n,n-1] <- d[n-1]/6
  B[n,n] <- d[n-1]/3
  for (i in 2:(nrow(B)-1)) {
    B[i,i-1] <- d[i-1]/6
    B[i,i] <- (d[i-1]+d[i])/3
    B[i,i+1] <- d[i]/6
  }
  B
}
compute_A <- function(d,n){
  A <-matrix(0, nrow = n, ncol = n)
  A[1,1] <- d[1]/2
  A[n,n] <- d[n-1]/2
  for (i in 2:(nrow(A)-1)) {
    A[i,i] <- (d[i-1]+d[i])/2
  }
  A
}

### For Precision matrix of the exact method: augmented space
Compute_Ti <- function(svec,p = 2,i){
  Ti <- matrix(0,nrow = p, ncol = p)
  delta <- diff(c(0,svec))
  denom <- factorial(c(0:(p-1)))
  numr <- delta[i+1]^(0:(p-1))
  Ti[1,] <- numr/denom
  for (i in 2:p) {
    Ti[i,] <- c(rep(0,(i-1)),Ti[(i-1),((i-1):(p-1))])
  }
  Ti
}
Compute_Ci <- function(svec, p = 2, i, is.cov = FALSE){
  delta <- diff(c(0,svec))
  Result <- matrix(0,nrow = p, ncol = p)
  index <- i+1
  for (i in 1:p) {
    for (j in i:p) {
      Result[i,j] <- (delta[index]^(2*p + 1 - i - j))/((2*p + 1 - i - j)*factorial(p-i)*factorial(p-j))
    }
  }
  Result <- Matrix::forceSymmetric(Result)
  if(is.cov == T){
    return(Result)
  }
  else{
    round(solve(Result),digits = 5)
  }
}
Compute_Ai <- function(svec, p = 2, i){
  Ci <- Compute_Ci(svec,p,i)
  Ti <- Compute_Ti(svec,p,i)
  Ai <- t(Ti) %*% Ci
  Ai <- Ai %*% Ti
  Ai
}
Compute_Bi <- function(svec, p = 2, i){
  Ci <- Compute_Ci(svec,p,i)
  Ti <- Compute_Ti(svec,p,i)
  Bi <- -t(Ti) %*% Ci
  Bi
}
Compute_Aug_Wp_Prec <- function(svec, p = 2) {
  n <- length(svec)
  Blist <- list()
  AClist <- list()
  for (i in 1:(n - 1)) {
    AClist[[i]] <- Compute_Ai(svec = svec, i = i , p = p) + Compute_Ci(svec = svec, i = (i-1), p = p)
  }
  AClist[[n]] <- Compute_Ci(svec = svec, i = n - 1, p = p)
  for (i in 1:(n - 1)) {
    Blist[[i]] <- Compute_Bi(svec = svec, i = i, p = p)
  }
  Mlist <- list()
  M <- matrix(0, nrow = 0, ncol = n*p)
  for (i in 1:(n-1)) {
    Mlist[[i]] <- cbind(matrix(0,nrow = p, ncol = p * (i-1)), AClist[[i]], Blist[[i]], matrix(0,nrow = p, ncol = (p * (n-i-1))) )
    M <- rbind(M,Mlist[[i]])
  }
  M <- rbind(M,cbind(matrix(0,nrow = p, ncol = p * (n-1)), AClist[[n]]))
  M <- Matrix::forceSymmetric(M)
  as(as.matrix(M), "dgTMatrix")
}
Compute_design_Aug <- function(svec, p = 2){
  Design <- Diagonal((p * length(svec)), x = 0)
  diag(Design)[seq(1,nrow(Design), by = p)] <- 1
  as(as.matrix(Design[seq(1,nrow(Design), by = p),]), "dgTMatrix")
}
Compute_design_Aug_deriv <- function(svec, p = 2, degree){
  Design <- Diagonal((p * length(svec)), x = 0)
  diag(Design)[(seq(1,nrow(Design), by = p) + degree)] <- 1
  as(as.matrix(Design[(seq(1,nrow(Design), by = p) + degree),]), "dgTMatrix")
}








### For posterior summary:
simulate_refined_resolution <- function(x, refined_x, gx_samps){
  global_poly <- function(x){
    cbind(rep(1, length(x)), x)
  }
  local_poly <- function(x, refined_x){
    dif <- diff(x)
    nn <- length(refined_x)
    n <- length(x)
    D <- matrix(0, nrow = nn, ncol = (n-1))
    for (j in 1:nn) {
      for (i in 1:(n-1)) {
        if(refined_x[j] <= x[i]){
          D[j, i] <- 0
        }
        else if (refined_x[j] <= x[i+1] & refined_x[j] >= x[i]){
          D[j, i] <- 0.5 * (refined_x[j] - x[i])^2
        }
        else{
          D[j, i] <- 0.5 * dif[i]^2 + dif[i] * (refined_x[j] - x[i+1])
        }
      }
    }
    D
  }
  global_sample_refined <- global_poly(refined_x) %*% gx_samps[length(x):nrow(gx_samps),]
  local_sample_refined <- local_poly(x, refined_x) %*% gx_samps[1:(length(x)-1),]
  overall_refined <- local_sample_refined + global_sample_refined
  overall_refined
}
compute_numeric_deriv <- function(gx, h, degree = 1){
  diff(gx, differences = degree)/(h^degree)
}


### For high-resolution representation:
global_poly <- function(x, p = 2){
  result <- NULL
  for (i in 1:p) {
    result <- cbind(result, x^(i-1))
  }
  result
}
local_poly <- function(x, refined_x, p = 2){
  dif <- diff(x)
  nn <- length(refined_x)
  n <- length(x)
  D <- matrix(0, nrow = nn, ncol = (n-1))
  for (j in 1:nn) {
    for (i in 1:(n-1)) {
      if(refined_x[j] <= x[i]){
        D[j, i] <- 0
      }
      else if (refined_x[j] <= x[i+1] & refined_x[j] >= x[i]){
        D[j, i] <- (1/factorial(p)) * (refined_x[j] - x[i])^p
      }
      else{
        k <- 1:p
        D[j,i] <- sum((dif[i]^k)*((refined_x[j] - x[i+1])^(p-k))/(factorial(k)*factorial(p-k)))
      }
    }
  }
  D
}
knots_RW2 <- function(x, refined_x){
  dif <- diff(x)
  nn <- length(refined_x)
  n <- length(x)
  D <- matrix(0, nrow = nn, ncol = n)
  for (j in 1:nn) {
    for (i in 2:(n-1)) {
      if (refined_x[j] >= x[i-1] & refined_x[j] <= x[i+1] ){
        if(refined_x[j] >= x[i-1] & refined_x[j] <= x[i]){
          D[j, i] <- (refined_x[j] - x[i-1])/dif[i-1]
        }
        else{
          D[j, i] <- 1 - (refined_x[j] - x[i])/dif[i]
        }
      }
      else{
        D[j,i] <- 0
      }
    }
    if(refined_x[j] <= x[2]){
      D[j, 1] <- 1 - (refined_x[j] - x[1])/dif[1]
    }
    else if(refined_x[j] >= x[n-1]){
      D[j, n] <- (refined_x[j] - x[n-1])/dif[n-1]
    }
  }
  D
}

# try <- knots_RW2(x = c(1,2,3,4,5),refined_x = seq(0.5,5.5,by = 0.1)) 
# plot(try[,1]~ seq(0.5,5.5,by = 0.1), type = 'l')
# lines(try[,2]~ seq(0.5,5.5,by = 0.1), type = 'l')
# lines(try[,3]~ seq(0.5,5.5,by = 0.1), type = 'l')
# lines(try[,4]~ seq(0.5,5.5,by = 0.1), type = 'l')
# lines(try[,5]~ seq(0.5,5.5,by = 0.1), type = 'l')








##########################################################################################
##########################################################################################
### Implement the smoothing for Gaussian likelihood: Report the AGHQ object
## The proposed method with Overlapping Spline
Imple_BayesRegression <- function(y,x,knots = NULL, p = 2, prior = NULL, aghq_k = 7, likelihood = "Gaussian", prior.type = "PC", overdispersion = F, Xf = NULL){
  if(prior.type == "PC"){
    if(likelihood == "Gaussian"){
      if(is.null(prior)){
        u1 = 1
        alpha1 = 0.5
        u2 = 1
        alpha2 = 0.5
        betaprec = 10^(-6)
      }
      else{
        u1 = prior$u1
        alpha1 = prior$alpha1
        u2 = prior$u2
        alpha2 = prior$alpha2
        betaprec = prior$betaprec
      }
      
      if(is.null(knots)){
        P_proposed <- compute_weights_precision(x)
        D_proposed <- compute_weights_design(x,p)
        ## Design matrix for the global polynomial
        X = rbind(c(1,rep(0,p-1)), D_proposed[,1:p])
        X = cbind(X, Xf)
        ## Design matrix for the spline basis weights
        B = rbind(rep(0, (length(x)-1)), D_proposed[,(p+1):ncol(D_proposed)])
        tmbdat <- list(
          # Design matrix
          X = X,
          B = B,
          P = as(P_proposed,'dgTMatrix'),
          logPdet = as.numeric(determinant(P_proposed,logarithm = T)$modulus),
          # Response
          y = y,
          # PC Prior params
          u1 = u1,
          alpha1 = alpha1,
          u2 = u2,
          alpha2 = alpha2,
          betaprec = betaprec
        )
        tmbparams <- list(
          W = c(rep(0, (ncol(X) + ncol(B)))), # W = c(U,beta); U = B-Spline coefficients
          theta1 = 0, # -2log(sigma)
          theta2 = 0
        )
        dyn.load(dynlib("00_Gaussian_Smoothing_PC"))
        ff <- TMB::MakeADFun(
          data = tmbdat,
          parameters = tmbparams,
          random = "W",
          DLL = "00_Gaussian_Smoothing_PC",
          silent = TRUE
        )
        # Hessian not implemented for RE models
        ff$he <- function(w) numDeriv::jacobian(ff$gr,w)
      }
      else{
        P_proposed <- compute_weights_precision(knots)
        D_proposed <- compute_weights_design(x,p)
        ## Design matrix for the global polynomial
        X = rbind(global_poly(x[1], p), D_proposed[,1:p])
        X = cbind(X, Xf)
        ## Design matrix for the spline basis weights
        B = as(local_poly(x = knots, refined_x = x, p = p),"dgTMatrix")
        tmbdat <- list(
          # Design matrix
          X = X,
          B = B,
          P = as(P_proposed,'dgTMatrix'),
          logPdet = as.numeric(determinant(P_proposed,logarithm = T)$modulus),
          # Response
          y = y,
          # PC Prior params
          u1 = u1,
          alpha1 = alpha1,
          u2 = u2,
          alpha2 = alpha2,
          betaprec = betaprec
        )
        tmbparams <- list(
          W = c(rep(0, (ncol(X) + ncol(B)))), # W = c(U,beta); U = B-Spline coefficients
          theta1 = 0, # -2log(sigma)
          theta2 = 0
        )
        dyn.load(dynlib("00_Gaussian_Smoothing_PC"))
        ff <- TMB::MakeADFun(
          data = tmbdat,
          parameters = tmbparams,
          random = "W",
          DLL = "00_Gaussian_Smoothing_PC",
          silent = TRUE
        )
        # Hessian not implemented for RE models
        ff$he <- function(w) numDeriv::jacobian(ff$gr,w)
      }
      
      aghq::marginal_laplace_tmb(ff,aghq_k,c(0,0))
    }
    else if(likelihood == "Poisson"){
      if(is.null(prior)){
        u1 = 1
        alpha1 = 0.5
        betaprec = 10^(-6)
        if(overdispersion){
          u2 = 1
          alpha2 = 0.5
        }
      }
      else{
        u1 = prior$u1
        alpha1 = prior$alpha1
        betaprec = prior$betaprec
        if(overdispersion){
          u2 = prior$u2
          alpha2 = prior$alpha2
        }
      }
      if(overdispersion){
        if(is.null(knots)){
          P_proposed <- compute_weights_precision(x)
          D_proposed <- compute_weights_design(x,p)
          ## Design matrix for the global polynomial
          X = rbind(c(1,rep(0,p-1)), D_proposed[,1:p])
          X = cbind(X, Xf)
          ## Design matrix for the spline basis weights
          B = rbind(rep(0, (length(x)-1)), D_proposed[,(p+1):ncol(D_proposed)])
          I = diag(nrow = length(y), ncol = length(y))
          tmbdat <- list(
            # Design matrix
            X = X,
            B = B,
            I = as(I, "dgTMatrix"),
            P = as(P_proposed,'dgTMatrix'),
            logPdet = as.numeric(determinant(P_proposed,logarithm = T)$modulus),
            # Response
            y = y,
            # PC Prior params
            u1 = u1,
            alpha1 = alpha1,
            u2 = u2,
            alpha2 = alpha2,
            betaprec = betaprec
          )
          tmbparams <- list(
            W = c(rep(0, (ncol(X) + ncol(B) + ncol(I)))), # W = c(U,beta); U = B-Spline coefficients
            theta1 = 0,
            theta2 = 0
          )
          dyn.load(dynlib("00_Poisson_Smoothing_PC_overdisp"))
          ff <- TMB::MakeADFun(
            data = tmbdat,
            parameters = tmbparams,
            random = "W",
            DLL = "00_Poisson_Smoothing_PC_overdisp",
            silent = TRUE
          )
          # Hessian not implemented for RE models
          ff$he <- function(w) numDeriv::jacobian(ff$gr,w)
        }
        else{
          P_proposed <- compute_weights_precision(knots)
          D_proposed <- compute_weights_design(x,p)
          ## Design matrix for the global polynomial
          X = rbind(global_poly(x[1], p), D_proposed[,1:p])
          X = cbind(X, Xf)
          ## Design matrix for the spline basis weights
          B = as(local_poly(x = knots, refined_x = x, p = p),"dgTMatrix")
          I = diag(nrow = length(y), ncol = length(y))
          tmbdat <- list(
            # Design matrix
            X = X,
            B = B,
            I = as(I, "dgTMatrix"),
            P = as(P_proposed,'dgTMatrix'),
            logPdet = as.numeric(determinant(P_proposed,logarithm = T)$modulus),
            # Response
            y = y,
            # PC Prior params
            u1 = u1,
            alpha1 = alpha1,
            u2 = u2,
            alpha2 = alpha2,
            betaprec = betaprec
          )
          tmbparams <- list(
            W = c(rep(0, (ncol(X) + ncol(B) + ncol(I)))), # W = c(U,beta); U = B-Spline coefficients
            theta1 = 0,
            theta2 = 0
          )
          dyn.load(dynlib("00_Poisson_Smoothing_PC_overdisp"))
          ff <- TMB::MakeADFun(
            data = tmbdat,
            parameters = tmbparams,
            random = "W",
            DLL = "00_Poisson_Smoothing_PC_overdisp",
            silent = TRUE
          )
          # Hessian not implemented for RE models
          ff$he <- function(w) numDeriv::jacobian(ff$gr,w)
        }
        aghq::marginal_laplace_tmb(ff,aghq_k,c(0,0))
      }
      else{
        if(is.null(knots)){
          P_proposed <- compute_weights_precision(x)
          D_proposed <- compute_weights_design(x,p)
          ## Design matrix for the global polynomial
          X = rbind(c(1,rep(0,p-1)), D_proposed[,1:p])
          X = cbind(X, Xf)
          ## Design matrix for the spline basis weights
          B = rbind(rep(0, (length(x)-1)), D_proposed[,(p+1):ncol(D_proposed)])
          tmbdat <- list(
            # Design matrix
            X = X,
            B = B,
            P = as(P_proposed,'dgTMatrix'),
            logPdet = as.numeric(determinant(P_proposed,logarithm = T)$modulus),
            # Response
            y = y,
            # PC Prior params
            u1 = u1,
            alpha1 = alpha1,
            betaprec = betaprec
          )
          tmbparams <- list(
            W = c(rep(0, ncol(D_proposed))), # W = c(U,beta); U = B-Spline coefficients
            theta1 = 0
          )
          dyn.load(dynlib("00_Poisson_Smoothing_PC"))
          ff <- TMB::MakeADFun(
            data = tmbdat,
            parameters = tmbparams,
            random = "W",
            DLL = "00_Poisson_Smoothing_PC",
            silent = TRUE
          )
          # Hessian not implemented for RE models
          ff$he <- function(w) numDeriv::jacobian(ff$gr,w)
        }
        else{
          P_proposed <- compute_weights_precision(knots)
          D_proposed <- compute_weights_design(x,p)
          ## Design matrix for the global polynomial
          X = rbind(global_poly(x[1], p), D_proposed[,1:p])
          X = cbind(X, Xf)
          ## Design matrix for the spline basis weights
          B = as(local_poly(x = knots, refined_x = x, p = p),"dgTMatrix")
          tmbdat <- list(
            # Design matrix
            X = X,
            B = B,
            P = as(P_proposed,'dgTMatrix'),
            logPdet = as.numeric(determinant(P_proposed,logarithm = T)$modulus),
            # Response
            y = y,
            # PC Prior params
            u1 = u1,
            alpha1 = alpha1,
            betaprec = betaprec
          )
          tmbparams <- list(
            W = c(rep(0, (ncol(X) + ncol(B)))), # W = c(U,beta); U = B-Spline coefficients
            theta1 = 0
          )
          dyn.load(dynlib("00_Poisson_Smoothing_PC"))
          ff <- TMB::MakeADFun(
            data = tmbdat,
            parameters = tmbparams,
            random = "W",
            DLL = "00_Poisson_Smoothing_PC",
            silent = TRUE
          )
          # Hessian not implemented for RE models
          ff$he <- function(w) numDeriv::jacobian(ff$gr,w)
        }
        aghq::marginal_laplace_tmb(ff,aghq_k,c(0))
      }
    }
  }
  else if(prior.type == "log-gamma"){
    if(likelihood == "Gaussian"){
      if(is.null(prior)){
        # a is shape, b is inverse-scale, mean is a/b, variance is a/(b^2)
        a = 1
        b = 5e-5
        u2 = 1
        alpha2 = 0.5
        betaprec = 10^(-6)
      }
      else{
        a = prior$a
        b = prior$b
        u2 = prior$u2
        alpha2 = prior$alpha2
        betaprec = prior$betaprec
      }
      
      if(is.null(knots)){
        P_proposed <- compute_weights_precision(x)
        D_proposed <- compute_weights_design(x,p)
        ## Design matrix for the global polynomial
        X = rbind(c(1,rep(0,p-1)), D_proposed[,1:p])
        X = cbind(X, Xf)
        ## Design matrix for the spline basis weights
        B = rbind(rep(0, (length(x)-1)), D_proposed[,(p+1):ncol(D_proposed)])
        tmbdat <- list(
          # Design matrix
          X = X,
          B = B,
          P = as(P_proposed,'dgTMatrix'),
          logPdet = as.numeric(determinant(P_proposed,logarithm = T)$modulus),
          # Response
          y = y,
          # PC Prior params
          u1 = u1,
          alpha1 = alpha1,
          # PC Prior params
          u2 = u2,
          alpha2 = alpha2,
          betaprec = betaprec
        )
        tmbparams <- list(
          W = c(rep(0, ncol(D_proposed))), # W = c(U,beta); U = B-Spline coefficients
          theta1 = 0, # -2log(sigma)
          theta2 = 0
        )
        dyn.load(dynlib("00_Gaussian_Smoothing_logGamma"))
        ff <- TMB::MakeADFun(
          data = tmbdat,
          parameters = tmbparams,
          random = "W",
          DLL = "00_Gaussian_Smoothing_logGamma",
          silent = TRUE
        )
        # Hessian not implemented for RE models
        ff$he <- function(w) numDeriv::jacobian(ff$gr,w)
      }
      else{
        P_proposed <- compute_weights_precision(knots)
        D_proposed <- compute_weights_design(x,p)
        ## Design matrix for the global polynomial
        X = rbind(global_poly(x[1], p), D_proposed[,1:p])
        X = cbind(X, Xf)
        ## Design matrix for the spline basis weights
        B = as(local_poly(x = knots, refined_x = x, p = p),"dgTMatrix")
        tmbdat <- list(
          # Design matrix
          X = X,
          B = B,
          P = as(P_proposed,'dgTMatrix'),
          logPdet = as.numeric(determinant(P_proposed,logarithm = T)$modulus),
          # Response
          y = y,
          # Log Gamma:
          a = a,
          b = b,
          # PC prior
          u2 = u2,
          alpha2 = alpha2,
          betaprec = betaprec
        )
        tmbparams <- list(
          W = c(rep(0, (ncol(X) + ncol(B)))), # W = c(U,beta); U = B-Spline coefficients
          theta1 = 0, # -2log(sigma)
          theta2 = 0
        )
        dyn.load(dynlib("00_Gaussian_Smoothing_logGamma"))
        ff <- TMB::MakeADFun(
          data = tmbdat,
          parameters = tmbparams,
          random = "W",
          DLL = "00_Gaussian_Smoothing_logGamma",
          silent = TRUE
        )
        # Hessian not implemented for RE models
        ff$he <- function(w) numDeriv::jacobian(ff$gr,w)
      }
      
      aghq::marginal_laplace_tmb(ff,aghq_k,c(0,0))
    }
    else if(likelihood == "Poisson"){
      if(is.null(prior)){
        a <- 1
        b <- 5e-5
        betaprec = 10^(-6)
      }
      else{
        a = prior$a
        b = prior$b
        betaprec = prior$betaprec
      }
      
      if(is.null(knots)){
        P_proposed <- compute_weights_precision(x)
        D_proposed <- compute_weights_design(x,p)
        ## Design matrix for the global polynomial
        X = rbind(c(1,rep(0,p-1)), D_proposed[,1:p])
        X = cbind(X, Xf)
        ## Design matrix for the spline basis weights
        B = rbind(rep(0, (length(x)-1)), D_proposed[,(p+1):ncol(D_proposed)])
        tmbdat <- list(
          # Design matrix
          X = X,
          B = B,
          P = as(P_proposed,'dgTMatrix'),
          logPdet = as.numeric(determinant(P_proposed,logarithm = T)$modulus),
          # Response
          y = y,
          # PC Prior params
          u1 = u1,
          alpha1 = alpha1,
          betaprec = betaprec
        )
        tmbparams <- list(
          W = c(rep(0, ncol(D_proposed))), # W = c(U,beta); U = B-Spline coefficients
          theta1 = 0
        )
        dyn.load(dynlib("00_Poisson_Smoothing_logGamma"))
        ff <- TMB::MakeADFun(
          data = tmbdat,
          parameters = tmbparams,
          random = "W",
          DLL = "00_Poisson_Smoothing_logGamma",
          silent = TRUE
        )
        # Hessian not implemented for RE models
        ff$he <- function(w) numDeriv::jacobian(ff$gr,w)
      }
      else{
        P_proposed <- compute_weights_precision(knots)
        D_proposed <- compute_weights_design(x,p)
        ## Design matrix for the global polynomial
        X = rbind(global_poly(x[1], p), D_proposed[,1:p])
        X = cbind(X, Xf)
        ## Design matrix for the spline basis weights
        B = as(local_poly(x = knots, refined_x = x, p = p),"dgTMatrix")
        tmbdat <- list(
          # Design matrix
          X = X,
          B = B,
          P = as(P_proposed,'dgTMatrix'),
          logPdet = as.numeric(determinant(P_proposed,logarithm = T)$modulus),
          # Response
          y = y,
          # log gamma prior
          a = a,
          b = b,
          betaprec = betaprec
        )
        tmbparams <- list(
          W = c(rep(0, (ncol(X) + ncol(B)))), # W = c(U,beta); U = B-Spline coefficients
          theta1 = 0
        )
        dyn.load(dynlib("00_Poisson_Smoothing_logGamma"))
        ff <- TMB::MakeADFun(
          data = tmbdat,
          parameters = tmbparams,
          random = "W",
          DLL = "00_Poisson_Smoothing_logGamma",
          silent = TRUE
        )
        # Hessian not implemented for RE models
        ff$he <- function(w) numDeriv::jacobian(ff$gr,w)
      }
      
      aghq::marginal_laplace_tmb(ff,aghq_k,c(0))
    }
  }
  else {
    return(stop("Not supported prior type input, please check."))
  }
  
}

## The proposed method with Overlapping Spline (any starting point)
Imple_BayesRegression_constr <- function(y, x, knots = NULL,p = 2,prior = NULL, aghq_k = 7, likelihood = "Gaussian", s = NULL, prior.type = "PC"){
  if(prior.type == "PC"){
    if(is.null(s)){
      if(is.null(knots)){
        s <- ceiling(length(x)/2)
      }
      else{
        s <- ceiling(length(knots)/2)
      }
    }
    if(likelihood == "Gaussian"){
      if(is.null(prior)){
        u1 = 1
        alpha1 = 0.5
        u2 = 1
        alpha2 = 0.5
        betaprec = 10^(-6)
      }
      else{
        u1 = prior$u1
        alpha1 = prior$alpha1
        u2 = prior$u2
        alpha2 = prior$alpha2
        betaprec = prior$betaprec
      }
      cat("start setting up the proposed designs \n")
      if(is.null(knots)){
        P_proposed <- condition_precision(x, s = s, p = p)
        D_proposed <- condition_local_design(x, s = s, p = p)
        ## Design matrix for the global polynomial
        X = condition_global_design(x, s = s, p = p)
        X = as(X,'dgTMatrix')
        
        ## Design matrix for the spline basis weights
        B = as(D_proposed,'dgTMatrix')
        tmbdat <- list(
          # Design matrix
          X = X,
          B = B,
          P = as(P_proposed,'dgTMatrix'),
          logPdet = as.numeric(determinant(P_proposed,logarithm = T)$modulus),
          # Response
          y = y,
          # PC Prior params
          u1 = u1,
          alpha1 = alpha1,
          u2 = u2,
          alpha2 = alpha2,
          betaprec = betaprec
        )
        tmbparams <- list(
          W = c(rep(0, (ncol(B) + ncol(X)))), # W = c(U,beta); U = B-Spline coefficients
          theta1 = 0, # -2log(sigma)
          theta2 = 0
        )
        dyn.load(dynlib("00_Gaussian_Smoothing_PC"))
        ff <- TMB::MakeADFun(
          data = tmbdat,
          parameters = tmbparams,
          random = "W",
          DLL = "00_Gaussian_Smoothing_PC",
          silent = TRUE
        )
        # Hessian not implemented for RE models
        ff$he <- function(w) numDeriv::jacobian(ff$gr,w)
      }
      else{
        P_proposed <- condition_precision(knots, s = s, p = p)
        D_proposed <- condition_local_poly(knots, refined_x = x, s = s, p = p)
        ## Design matrix for the global polynomial
        X = condition_global_poly(knots = knots, refined_x = x, s = s, p = p)
        X = as(X,'dgTMatrix')
        ## Design matrix for the spline basis weights
        B = as(D_proposed,'dgTMatrix')
        
        tmbdat <- list(
          # Design matrix
          X = X,
          B = B,
          P = as(P_proposed,'dgTMatrix'),
          logPdet = as.numeric(determinant(P_proposed,logarithm = T)$modulus),
          # Response
          y = y,
          # PC Prior params
          u1 = u1,
          alpha1 = alpha1,
          u2 = u2,
          alpha2 = alpha2,
          betaprec = betaprec
        )
        tmbparams <- list(
          W = c(rep(0, (ncol(X) + ncol(B)))), # W = c(U,beta); U = B-Spline coefficients
          theta1 = 0, # -2log(sigma)
          theta2 = 0
        )
        dyn.load(dynlib("00_Gaussian_Smoothing_PC"))
        ff <- TMB::MakeADFun(
          data = tmbdat,
          parameters = tmbparams,
          random = "W",
          DLL = "00_Gaussian_Smoothing_PC",
          silent = TRUE
        )
        # Hessian not implemented for RE models
        ff$he <- function(w) numDeriv::jacobian(ff$gr,w)
      }
      cat("start the inferences using aghq \n")
      aghq::marginal_laplace_tmb(ff,aghq_k,c(0,0))
    }
    else if(likelihood == "Poisson"){
      if(is.null(prior)){
        u1 = 1
        alpha1 = 0.5
        betaprec = 10^(-6)
      }
      else{
        u1 = prior$u1
        alpha1 = prior$alpha1
        betaprec = prior$betaprec
      }
      cat("start setting up the proposed designs \n")
      if(is.null(knots)){
        P_proposed <- condition_precision(x, s = s, p = p)
        D_proposed <- condition_local_design(x, s = s, p = p)
        ## Design matrix for the global polynomial
        X = condition_global_design(x, s = s, p = p)
        X = as(X,'dgTMatrix')
        ## Design matrix for the spline basis weights
        B = as(D_proposed,'dgTMatrix')
        tmbdat <- list(
          # Design matrix
          X = X,
          B = B,
          P = as(P_proposed,'dgTMatrix'),
          logPdet = as.numeric(determinant(P_proposed,logarithm = T)$modulus),
          # Response
          y = y,
          # PC Prior params
          u1 = u1,
          alpha1 = alpha1,
          betaprec = betaprec
        )
        tmbparams <- list(
          W = c(rep(0, (ncol(B) + ncol(X)))), # W = c(U,beta); U = B-Spline coefficients
          theta1 = 0 # -2log(sigma)
        )
        dyn.load(dynlib("00_Poisson_Smoothing_PC"))
        ff <- TMB::MakeADFun(
          data = tmbdat,
          parameters = tmbparams,
          random = "W",
          DLL = "00_Poisson_Smoothing_PC",
          silent = TRUE
        )
        # Hessian not implemented for RE models
        ff$he <- function(w) numDeriv::jacobian(ff$gr,w)
      }
      else{
        P_proposed <- condition_precision(knots, s = s, p = p)
        D_proposed <- condition_local_poly(knots, refined_x = x, s = s, p = p)
        ## Design matrix for the global polynomial
        X = condition_global_poly(knots = knots, refined_x = x, s = s, p = p)
        
        X = as(X,'dgTMatrix')
        ## Design matrix for the spline basis weights
        B = as(D_proposed,'dgTMatrix')
        
        tmbdat <- list(
          # Design matrix
          X = X,
          B = B,
          P = as(P_proposed,'dgTMatrix'),
          logPdet = as.numeric(determinant(P_proposed,logarithm = T)$modulus),
          # Response
          y = y,
          # PC Prior params
          u1 = u1,
          alpha1 = alpha1,
          betaprec = betaprec
        )
        tmbparams <- list(
          W = c(rep(0, (ncol(X) + ncol(B)))), # W = c(U,beta); U = B-Spline coefficients
          theta1 = 0
        )
        dyn.load(dynlib("00_Poisson_Smoothing_PC"))
        ff <- TMB::MakeADFun(
          data = tmbdat,
          parameters = tmbparams,
          random = "W",
          DLL = "00_Poisson_Smoothing_PC",
          silent = TRUE
        )
        # Hessian not implemented for RE models
        ff$he <- function(w) numDeriv::jacobian(ff$gr,w)
      }
      cat("start the inferences using aghq \n")
      aghq::marginal_laplace_tmb(ff,aghq_k,c(0))
    }
  }
  else if(prior.type == "log-Gamma"){
    if(is.null(s)){
      if(is.null(knots)){
        s <- ceiling(length(x)/2)
      }
      else{
        s <- ceiling(length(knots)/2)
      }
    }
    if(likelihood == "Gaussian"){
      if(is.null(prior)){
        a = 1
        b = 5e-5
        u2 = 1
        alpha2 = 0.5
        betaprec = 10^(-6)
      }
      else{
        a = prior$a
        b = prior$b
        u2 = prior$u2
        alpha2 = prior$alpha2
        betaprec = prior$betaprec
      }
      cat("start setting up the proposed designs \n")
      if(is.null(knots)){
        P_proposed <- condition_precision(x, s = s, p = p)
        D_proposed <- condition_local_design(x, s = s, p = p)
        ## Design matrix for the global polynomial
        X = condition_global_design(x, s = s, p = p)
        X = as(X,'dgTMatrix')
        
        ## Design matrix for the spline basis weights
        B = as(D_proposed,'dgTMatrix')
        tmbdat <- list(
          # Design matrix
          X = X,
          B = B,
          P = as(P_proposed,'dgTMatrix'),
          logPdet = as.numeric(determinant(P_proposed,logarithm = T)$modulus),
          # Response
          y = y,
          # log Gamma Prior params
          a = a,
          b = b,
          # PC Prior params
          u2 = u2,
          alpha2 = alpha2,
          betaprec = betaprec
        )
        tmbparams <- list(
          W = c(rep(0, (ncol(B) + ncol(X)))), # W = c(U,beta); U = B-Spline coefficients
          theta1 = 0, # -2log(sigma)
          theta2 = 0
        )
        dyn.load(dynlib("00_Gaussian_Smoothing_logGamma"))
        ff <- TMB::MakeADFun(
          data = tmbdat,
          parameters = tmbparams,
          random = "W",
          DLL = "00_Gaussian_Smoothing_logGamma",
          silent = TRUE
        )
        # Hessian not implemented for RE models
        ff$he <- function(w) numDeriv::jacobian(ff$gr,w)
      }
      else{
        P_proposed <- condition_precision(knots, s = s, p = p)
        D_proposed <- condition_local_poly(knots, refined_x = x, s = s, p = p)
        ## Design matrix for the global polynomial
        X = condition_global_poly(knots = knots, refined_x = x, s = s, p = p)
        X = as(X,'dgTMatrix')
        ## Design matrix for the spline basis weights
        B = as(D_proposed,'dgTMatrix')
        
        tmbdat <- list(
          # Design matrix
          X = X,
          B = B,
          P = as(P_proposed,'dgTMatrix'),
          logPdet = as.numeric(determinant(P_proposed,logarithm = T)$modulus),
          # Response
          y = y,
          # log Gamma Prior params
          a = a,
          b = b,
          # PC Prior params
          u2 = u2,
          alpha2 = alpha2,
          betaprec = betaprec
        )
        tmbparams <- list(
          W = c(rep(0, (ncol(X) + ncol(B)))), # W = c(U,beta); U = B-Spline coefficients
          theta1 = 0, # -2log(sigma)
          theta2 = 0
        )
        dyn.load(dynlib("00_Gaussian_Smoothing_logGamma"))
        ff <- TMB::MakeADFun(
          data = tmbdat,
          parameters = tmbparams,
          random = "W",
          DLL = "00_Gaussian_Smoothing_logGamma",
          silent = TRUE
        )
        # Hessian not implemented for RE models
        ff$he <- function(w) numDeriv::jacobian(ff$gr,w)
      }
      cat("start the inferences using aghq \n")
      aghq::marginal_laplace_tmb(ff,aghq_k,c(0,0))
    }
    else if(likelihood == "Poisson"){
      if(is.null(prior)){
        a = 1
        b = 5e-5
        betaprec = 10^(-6)
      }
      else{
        a = prior$a
        b = prior$b
        betaprec = prior$betaprec
      }
      cat("start setting up the proposed designs \n")
      if(is.null(knots)){
        P_proposed <- condition_precision(x, s = s, p = p)
        D_proposed <- condition_local_design(x, s = s, p = p)
        ## Design matrix for the global polynomial
        X = condition_global_design(x, s = s, p = p)
        X = as(X,'dgTMatrix')
        ## Design matrix for the spline basis weights
        B = as(D_proposed,'dgTMatrix')
        tmbdat <- list(
          # Design matrix
          X = X,
          B = B,
          P = as(P_proposed,'dgTMatrix'),
          logPdet = as.numeric(determinant(P_proposed,logarithm = T)$modulus),
          # Response
          y = y,
          # log Gamma Prior params
          a = a,
          b = b,
          betaprec = betaprec
        )
        tmbparams <- list(
          W = c(rep(0, (ncol(B) + ncol(X)))), # W = c(U,beta); U = B-Spline coefficients
          theta1 = 0 # -2log(sigma)
        )
        dyn.load(dynlib("00_Poisson_Smoothing_logGamma"))
        ff <- TMB::MakeADFun(
          data = tmbdat,
          parameters = tmbparams,
          random = "W",
          DLL = "00_Poisson_Smoothing_logGamma",
          silent = TRUE
        )
        # Hessian not implemented for RE models
        ff$he <- function(w) numDeriv::jacobian(ff$gr,w)
      }
      else{
        P_proposed <- condition_precision(knots, s = s, p = p)
        D_proposed <- condition_local_poly(knots, refined_x = x, s = s, p = p)
        ## Design matrix for the global polynomial
        X = condition_global_poly(knots = knots, refined_x = x, s = s, p = p)
        
        X = as(X,'dgTMatrix')
        ## Design matrix for the spline basis weights
        B = as(D_proposed,'dgTMatrix')
        
        tmbdat <- list(
          # Design matrix
          X = X,
          B = B,
          P = as(P_proposed,'dgTMatrix'),
          logPdet = as.numeric(determinant(P_proposed,logarithm = T)$modulus),
          # Response
          y = y,
          # log Gamma Prior params
          a = a,
          b = b,
          betaprec = betaprec
        )
        tmbparams <- list(
          W = c(rep(0, (ncol(X) + ncol(B)))), # W = c(U,beta); U = B-Spline coefficients
          theta1 = 0
        )
        dyn.load(dynlib("00_Poisson_Smoothing_logGamma"))
        ff <- TMB::MakeADFun(
          data = tmbdat,
          parameters = tmbparams,
          random = "W",
          DLL = "00_Poisson_Smoothing_logGamma",
          silent = TRUE
        )
        # Hessian not implemented for RE models
        ff$he <- function(w) numDeriv::jacobian(ff$gr,w)
      }
      cat("start the inferences using aghq \n")
      aghq::marginal_laplace_tmb(ff,aghq_k,c(0))
    }
  }
  else {
    return(stop("Not supported prior type input, please check."))
  }
}

## The RW2 method (only PC prior allowed currently)
Imple_RW2 <- function(y,x,knots = NULL,prior = NULL, aghq_k = 7, diagonal = 1e-9, likelihood = "Gaussian"){
  if(likelihood == "Gaussian"){
    if(is.null(prior)){
      u1 = 1
      alpha1 = 0.5
      u2 = 1
      alpha2 = 0.5
      betaprec = 10^(-6)
    }
    else{
      u1 = prior$u1
      alpha1 = prior$alpha1
      u2 = prior$u2
      alpha2 = prior$alpha2
      betaprec = prior$betaprec
    }
    if(is.null(knots)){
      d <- diff(x)
      H <- compute_H_rue(d,n = length(x))
      A <- compute_A(d,n = length(x))
      n <- length(x)
      QRW2 <- as(t(H) %*% solve(A) %*% H, 'dgTMatrix')
      QRW2 <- as(QRW2 + Diagonal(n, x = diagonal), "dgTMatrix")
      designB <- as(Diagonal(n,1), "matrix")
      ### Implement exact SDE method:
      tmbdat <- list(
        X = as(designB,"dgTMatrix"),
        # Penalty(Precision) matrix
        P = as(as.matrix(QRW2),"dgTMatrix"),
        # Log determinant of penalty matrix (without the sigma part)
        # logPdet = as.numeric(sum(log(sort(all_values)[-(1:2)]))),
        logPdet = as.numeric(determinant(QRW2,logarithm = T)$modulus),
        # Response
        y = y,  # PC Prior params
        u1 = u1,
        alpha1 = alpha1,
        u2 = u2,
        alpha2 = alpha2,
        betaprec = betaprec
      )
      tmbparams <- list(
        W = c(rep(0, n)), # W = c(U,beta); U = B-Spline coefficients
        theta1 = 0, # -2log(sigma)
        theta2 = 0
      )
      dyn.load(dynlib("02_RW2Comparison"))
      ff <- TMB::MakeADFun(
        data = tmbdat,
        parameters = tmbparams,
        random = "W",
        DLL = "02_RW2Comparison",
        silent = TRUE
      )
      # Hessian not implemented for RE models
      ff$he <- function(w) numDeriv::jacobian(ff$gr,w)
    }
    else{
      designB_rw2_reduced <- knots_RW2(x = knots, refined_x =  x)
      d <- diff(knots)
      H <- compute_H_rue(d,n = length(knots))
      A <- compute_A(d,n = length(knots))
      n <- length(knots)
      Q_rw2_reduced <- as(t(H) %*% solve(A) %*% H, 'dgTMatrix')
      Q_rw2_reduced <- as(Q_rw2_reduced + Diagonal(n, x = diagonal), "dgTMatrix")
      tmbdat <- list(
        X = as(designB_rw2_reduced,"dgTMatrix"),
        # Penalty(Precision) matrix
        P = as(as.matrix(Q_rw2_reduced),"dgTMatrix"),
        # Log determinant of penalty matrix (without the sigma part)
        # logPdet = as.numeric(sum(log(sort(all_values)[-(1:2)]))),
        logPdet = as.numeric(determinant(Q_rw2_reduced,logarithm = T)$modulus),
        # Response
        y = y,  # PC Prior params
        u1 = u1,
        alpha1 = alpha1,
        u2 = u2,
        alpha2 = alpha2,
        betaprec = betaprec
      )
      tmbparams <- list(
        W = c(rep(0, n)), # W = c(U,beta); U = B-Spline coefficients
        theta1 = 0, # -2log(sigma)
        theta2 = 0
      )
      dyn.load(dynlib("02_RW2Comparison"))
      ff <- TMB::MakeADFun(
        data = tmbdat,
        parameters = tmbparams,
        random = "W",
        DLL = "02_RW2Comparison",
        silent = TRUE
      )
      # Hessian not implemented for RE models
      ff$he <- function(w) numDeriv::jacobian(ff$gr,w)
      
    }
    aghq::marginal_laplace_tmb(ff,aghq_k,c(0,0))
  }
  else if(likelihood == "Poisson"){
    if(is.null(prior)){
      u = 1
      alpha = 0.5
      betaprec = 10^(-6)
    }
    else{
      u = prior$u
      alpha = prior$alpha
      betaprec = prior$betaprec
    }
    if(is.null(knots)){
      d <- diff(x)
      H <- compute_H_rue(d,n = length(x))
      A <- compute_A(d,n = length(x))
      n <- length(x)
      QRW2 <- as(t(H) %*% solve(A) %*% H, 'dgTMatrix')
      QRW2 <- as(QRW2 + Diagonal(n, x = diagonal), "dgTMatrix")
      designB <- as(Diagonal(n,1), "matrix")
      ### Implement exact SDE method:
      tmbdat <- list(
        X = as(designB,"dgTMatrix"),
        # Penalty(Precision) matrix
        P = as(as.matrix(QRW2),"dgTMatrix"),
        # Log determinant of penalty matrix (without the sigma part)
        # logPdet = as.numeric(sum(log(sort(all_values)[-(1:2)]))),
        logPdet = as.numeric(determinant(QRW2,logarithm = T)$modulus),
        # Response
        y = y,  # PC Prior params
        u = u,
        alpha = alpha,
        betaprec = betaprec
      )
      tmbparams <- list(
        W = c(rep(0, n)), # W = c(U,beta); U = B-Spline coefficients
        theta = 0)
      
      dyn.load(dynlib("00_Poisson_Smoothing_PC"))
      ff <- TMB::MakeADFun(
        data = tmbdat,
        parameters = tmbparams,
        random = "W",
        DLL = "00_Poisson_Smoothing_PC",
        silent = TRUE
      )
      # Hessian not implemented for RE models
      ff$he <- function(w) numDeriv::jacobian(ff$gr,w)
    }
    else{
      designB_rw2_reduced <- knots_RW2(x = knots, refined_x =  x)
      d <- diff(knots)
      H <- compute_H_rue(d,n = length(knots))
      A <- compute_A(d,n = length(knots))
      n <- length(knots)
      Q_rw2_reduced <- as(t(H) %*% solve(A) %*% H, 'dgTMatrix')
      Q_rw2_reduced <- as(Q_rw2_reduced + Diagonal(n, x = diagonal), "dgTMatrix")
      tmbdat <- list(
        X = as(designB_rw2_reduced,"dgTMatrix"),
        # Penalty(Precision) matrix
        P = as(as.matrix(Q_rw2_reduced),"dgTMatrix"),
        # Log determinant of penalty matrix (without the sigma part)
        # logPdet = as.numeric(sum(log(sort(all_values)[-(1:2)]))),
        logPdet = as.numeric(determinant(Q_rw2_reduced,logarithm = T)$modulus),
        # Response
        y = y,  # PC Prior params
        u = u,
        alpha = alpha,
        betaprec = betaprec
      )
      tmbparams <- list(
        W = c(rep(0, n)), # W = c(U,beta); U = B-Spline coefficients
        theta = 0
      )
      dyn.load(dynlib("00_Poisson_Smoothing_PC"))
      ff <- TMB::MakeADFun(
        data = tmbdat,
        parameters = tmbparams,
        random = "W",
        DLL = "00_Poisson_Smoothing_PC",
        silent = TRUE
      )
      # Hessian not implemented for RE models
      ff$he <- function(w) numDeriv::jacobian(ff$gr,w)
    }
    aghq::marginal_laplace_tmb(ff,aghq_k,c(0))
  }
}

## The proposed method with B spline (only PC prior allowed currently)
Imple_BsplineSmooth <- function(y,x, min = NULL, max = NULL,q, p = 2, prior = NULL, aghq_k = 7, likelihood = "Gaussian"){
  m <- 2*p - 1
  if(likelihood == "Gaussian"){
    if(is.null(prior)){
      u1 = 1
      alpha1 = 0.5
      u2 = 1
      alpha2 = 0.5
      betaprec = 10^(-6)
    }
    else{
      u1 = prior$u1
      alpha1 = prior$alpha1
      u2 = prior$u2
      alpha2 = prior$alpha2
      betaprec = prior$betaprec
    }
    if(is.null(min)|is.null(max)){
      min = range(x)[1]
      max = range(x)[2]
    }
    
    Q <- suppressWarnings(get_precision_Bspline(q,min,max,p))
    B <- get_design_Bspline(x,q,min,max,p)
    ## Design matrix for the global polynomial
    X = global_poly(x, p)
    tmbdat <- list(
      # Design matrix
      X = as(X,'dgTMatrix'),
      B = as(B,'dgTMatrix'),
      P = as(Q,'dgTMatrix'),
      logPdet = as.numeric(determinant(Q,logarithm = T)$modulus),
      # Response
      y = y,
      # PC Prior params
      u1 = u1,
      alpha1 = alpha1,
      u2 = u2,
      alpha2 = alpha2,
      betaprec = betaprec
    )
    tmbparams <- list(
      W = c(rep(0, (ncol(X) + ncol(B)))), # W = c(U,beta); U = B-Spline coefficients
      theta1 = 0, # -2log(sigma)
      theta2 = 0
    )
    dyn.load(dynlib("00_Gaussian_Smoothing_PC"))
    ff <- TMB::MakeADFun(
      data = tmbdat,
      parameters = tmbparams,
      random = "W",
      DLL = "00_Gaussian_Smoothing_PC",
      silent = TRUE
    )
    # Hessian not implemented for RE models
    ff$he <- function(w) numDeriv::jacobian(ff$gr,w)
    
    aghq::marginal_laplace_tmb(ff,aghq_k,c(0,0))
  }
  else if(likelihood == "Poisson"){
    if(is.null(prior)){
      u1 = 1
      alpha1 = 0.5
      betaprec = 10^(-6)
    }
    else{
      u1 = prior$u1
      alpha1 = prior$alpha1
      betaprec = prior$betaprec
    }
    
    if(is.null(min)|is.null(max)){
      min = range(x)[1]
      max = range(x)[2]
    }
    
    Q <- suppressWarnings(get_precision_Bspline(q,min,max,p))
    B <- get_design_Bspline(x,q,min,max,p)
    ## Design matrix for the global polynomial
    X = global_poly(x, p)
    tmbdat <- list(
      # Design matrix
      X = X,
      B = B,
      P = as(Q,'dgTMatrix'),
      logPdet = as.numeric(determinant(Q,logarithm = T)$modulus),
      # Response
      y = y,
      # PC Prior params
      u1 = u1,
      alpha1 = alpha1,
      u2 = u2,
      alpha2 = alpha2,
      betaprec = betaprec
    )
    tmbparams <- list(
      W = c(rep(0, (ncol(X) + ncol(B)))), # W = c(U,beta); U = B-Spline coefficients
      theta1 = 0, # -2log(sigma)
      theta2 = 0
    )
    dyn.load(dynlib("00_Poisson_Smoothing_PC"))
    ff <- TMB::MakeADFun(
      data = tmbdat,
      parameters = tmbparams,
      random = "W",
      DLL = "00_Poisson_Smoothing_PC",
      silent = TRUE
    )
    # Hessian not implemented for RE models
    ff$he <- function(w) numDeriv::jacobian(ff$gr,w)
    
    aghq::marginal_laplace_tmb(ff,aghq_k,c(0))
  }
}
#### Polynomial regression:
Imple_PolyReg <- function(y,x, p, prior = NULL, aghq_k = 7, likelihood = "Gaussian"){
  if(likelihood == "Gaussian"){
    if(is.null(prior)){
      u1 = 1
      alpha1 = 0.5
      betaprec = 10^(-6)
    }
    else{
      u1 = prior$u1
      alpha1 = prior$alpha1
      betaprec = prior$betaprec
    }
    
    
    ## Design matrix for the global polynomial
    X = global_poly(x, (p+1))
    tmbdat <- list(
      # Design matrix
      X = as(X,'dgTMatrix'),
      # Response
      y = y,
      # PC Prior params
      u1 = u1,
      alpha1 = alpha1,
      betaprec = betaprec
    )
    tmbparams <- list(
      W = c(rep(0,(p+1))), # W = c(U,beta); U = B-Spline coefficients
      theta1 = 0, # -2log(sigma)
      theta2 = 0
    )
    dyn.load(dynlib("00_Poly_Reg"))
    ff <- TMB::MakeADFun(
      data = tmbdat,
      parameters = tmbparams,
      random = "W",
      DLL = "00_Poly_Reg",
      silent = TRUE
    )
    # Hessian not implemented for RE models
    ff$he <- function(w) numDeriv::jacobian(ff$gr,w)
    
    aghq::marginal_laplace_tmb(ff,aghq_k,c(0))
  }
  else{
    print("Model not implemented yet")
  }
}


#### Exact WP regression with derivative:
Imple_WP <- function(y, x, p, prior = NULL, aghq_k = 7, likelihood = "Gaussian", prior.type = "PC"){
  if(prior.type == "PC"){
    x <- x - x[1] ### normalize with x[1]
    if(likelihood == "Gaussian"){
      if(is.null(prior)){
        u1 = 1
        alpha1 = 0.5
        u2 = 1
        alpha2 = 0.5
        betaprec = 10^(-6)
      }
      else{
        u1 = prior$u1
        alpha1 = prior$alpha1
        u2 = prior$u2
        alpha2 = prior$alpha2
        betaprec = prior$betaprec
      }
      
      svec <- x[-1] ## assume x[1] = 0
      B <- Compute_design_Aug(svec = svec, p = p)
      Q <- Compute_Aug_Wp_Prec(svec = svec,p = p)
      D_proposed <- compute_weights_design(x,p)
      ## Design matrix for the global polynomial
      X = rbind(c(1,rep(0,p-1)), D_proposed[,1:p])
      B <- rbind(rep(0, ncol(B)), B)
      tmbdat <- list(
        # Design matrix
        X = X,
        B = B,
        P = as(Q,'dgTMatrix'),
        logPdet = as.numeric(determinant(Q,logarithm = T)$modulus),
        # Response
        y = y,
        # PC Prior params
        u1 = u1,
        alpha1 = alpha1,
        u2 = u2,
        alpha2 = alpha2,
        betaprec = betaprec
      )
      tmbparams <- list(
        W = c(rep(0, (ncol(Q) + ncol(X)))), 
        theta1 = 0, # -2log(sigma)
        theta2 = 0
      )
      dyn.load(dynlib("00_Gaussian_Smoothing_PC"))
      ff <- TMB::MakeADFun(
        data = tmbdat,
        parameters = tmbparams,
        random = "W",
        DLL = "00_Gaussian_Smoothing_PC",
        silent = TRUE
      )
      # Hessian not implemented for RE models
      ff$he <- function(w) numDeriv::jacobian(ff$gr,w)
      
      aghq::marginal_laplace_tmb(ff,aghq_k,c(0,0))
    }
    else if(likelihood == "Poisson"){
      if(is.null(prior)){
        u1 = 1
        alpha1 = 0.5
        betaprec = 10^(-6)
      }
      else{
        u1 = prior$u1
        alpha1 = prior$alpha1
        betaprec = prior$betaprec
      }
      
      svec <- x[-1] ## assume x[1] = 0
      B <- Compute_design_Aug(svec = svec, p = p)
      Q <- Compute_Aug_Wp_Prec(svec = svec,p = p)
      D_proposed <- compute_weights_design(x,p)
      ## Design matrix for the global polynomial
      X = rbind(c(1,rep(0,p-1)), D_proposed[,1:p])
      tmbdat <- list(
        # Design matrix
        X = X,
        B = B,
        P = as(Q,'dgTMatrix'),
        logPdet = as.numeric(determinant(Q,logarithm = T)$modulus),
        # Response
        y = y,
        # PC Prior params
        u1 = u1,
        alpha1 = alpha1,
        betaprec = betaprec
      )
      tmbparams <- list(
        W = c(rep(0, (ncol(Q) + ncol(X)))), 
        theta1 = 0
      )
      dyn.load(dynlib("01_Spline_Smoothing"))
      ff <- TMB::MakeADFun(
        data = tmbdat,
        parameters = tmbparams,
        random = "W",
        DLL = "01_Spline_Smoothing",
        silent = TRUE
      )
      # Hessian not implemented for RE models
      ff$he <- function(w) numDeriv::jacobian(ff$gr,w)
      
      aghq::marginal_laplace_tmb(ff,aghq_k,c(0))
    }
  }
  else if(prior.type == "log-Gamma"){
    x <- x - x[1] ### normalize with x[1]
    if(likelihood == "Gaussian"){
      if(is.null(prior)){
        a = 1
        b = 5e-5
        u2 = 1
        alpha2 = 0.5
        betaprec = 10^(-6)
      }
      else{
        a = prior$a
        b = prior$b
        u2 = prior$u2
        alpha2 = prior$alpha2
        betaprec = prior$betaprec
      }
      svec <- x[-1] ## assume x[1] = 0
      B <- Compute_design_Aug(svec = svec, p = p)
      Q <- Compute_Aug_Wp_Prec(svec = svec,p = p)
      D_proposed <- compute_weights_design(x,p)
      ## Design matrix for the global polynomial
      X = rbind(c(1,rep(0,p-1)), D_proposed[,1:p])
      B <- rbind(rep(0, ncol(B)), B)
      tmbdat <- list(
        # Design matrix
        X = X,
        B = B,
        P = as(Q,'dgTMatrix'),
        logPdet = as.numeric(determinant(Q,logarithm = T)$modulus),
        # Response
        y = y,
        # log gamma Prior params
        a = a,
        b = b,
        u2 = u2,
        alpha2 = alpha2,
        betaprec = betaprec
      )
      tmbparams <- list(
        W = c(rep(0, (ncol(Q) + ncol(X)))), 
        theta1 = 0, # -2log(sigma)
        theta2 = 0
      )
      dyn.load(dynlib("00_Gaussian_Smoothing_logGamma"))
      ff <- TMB::MakeADFun(
        data = tmbdat,
        parameters = tmbparams,
        random = "W",
        DLL = "00_Gaussian_Smoothing_logGamma",
        silent = TRUE
      )
      # Hessian not implemented for RE models
      ff$he <- function(w) numDeriv::jacobian(ff$gr,w)
      
      aghq::marginal_laplace_tmb(ff,aghq_k,c(0,0))
    }
    else if(likelihood == "Poisson"){
      if(is.null(prior)){
        a = 1
        b = 5e-5
        betaprec = 10^(-6)
      }
      else{
        a = prior$a
        b = prior$b
        betaprec = prior$betaprec
      }
      
      svec <- x[-1] ## assume x[1] = 0
      B <- Compute_design_Aug(svec = svec, p = p)
      Q <- Compute_Aug_Wp_Prec(svec = svec,p = p)
      D_proposed <- compute_weights_design(x,p)
      ## Design matrix for the global polynomial
      X = rbind(c(1,rep(0,p-1)), D_proposed[,1:p])
      tmbdat <- list(
        # Design matrix
        X = X,
        B = B,
        P = as(Q,'dgTMatrix'),
        logPdet = as.numeric(determinant(Q,logarithm = T)$modulus),
        # Response
        y = y,
        # log gamma Prior params
        a = a,
        b = b,
        betaprec = betaprec
      )
      tmbparams <- list(
        W = c(rep(0, (ncol(Q) + ncol(X)))), 
        theta1 = 0
      )
      dyn.load(dynlib("00_Poisson_Smoothing_logGamma"))
      ff <- TMB::MakeADFun(
        data = tmbdat,
        parameters = tmbparams,
        random = "W",
        DLL = "00_Poisson_Smoothing_logGamma",
        silent = TRUE
      )
      # Hessian not implemented for RE models
      ff$he <- function(w) numDeriv::jacobian(ff$gr,w)
      
      aghq::marginal_laplace_tmb(ff,aghq_k,c(0))
    } 
  }
  else {
    return(stop("Not supported prior type input, please check."))
  }
  
}


#### Fit exact method, with known parameters, for Gaussian Data
### beta: known regression coefficients
### sigmaE: known SD parameter
### theta: known smoothing parameter (log prec)
### x: normed or unnormed x values (size should be n)
### y: observed data with size n
Imple_WP_known <- function(beta, sigmaE, theta, x, y, p = 3, optimizer = F){
  x_norm <- x - x[1]
  svec <- x_norm[-1] ## assume x[1] = 0
  B <- Compute_design_Aug(svec = svec, p = p)
  Q <- Compute_Aug_Wp_Prec(svec = svec,p = p)
  D_proposed <- compute_weights_design(x_norm, p)
  ## Design matrix for the global polynomial
  X = rbind(c(1,rep(0,p-1)), D_proposed[,1:p])
  B <- rbind(rep(0, ncol(B)), B)
  tmbdat <- list(
    # Design matrix
    X = X,
    B = B,
    P = as(Q,'dgTMatrix'),
    logPdet = as.numeric(determinant(Q,logarithm = T)$modulus),
    # Response
    y = y,
    beta = beta,
    sigmaE = sigmaE,
    theta1 = theta
  )
  tmbparams <- list(
    W = c(rep(0, ncol(Q)))
  )
  
  ff <- TMB::MakeADFun(
    data = tmbdat,
    parameters = tmbparams,
    # random = "W",
    DLL = "00_Gaussian_Smoothing_known",
    silent = TRUE
  )
  # Hessian not implemented for RE models
  ff$he <- function(w) numDeriv::jacobian(ff$gr,w)
  
  ## posterior mode and posterior precision matrix
  opt <- nlminb(start = ff$par, objective = ff$fn, gradient = ff$gr, hessian = ff$he, 
                control = list(eval.max = 200, iter.max = 200))
  prec_matrix <- forceSymmetric(ff$he(opt$par))
  
  if (optimizer == T){
    return(list(mean = opt$par, prec = as.matrix(prec_matrix), opt = opt))
  }
  
  ## Return posterior mean and posterior precision matrix
  return(list(mean = opt$par, prec = as.matrix(prec_matrix)))
}


#### Fit OS method, with known parameters, for Gaussian Data
Imple_OS_known <- function(beta, sigmaE, theta, x, k, y, p = 3, optimizer = F, nsamps = 10000){
  if(k == length(x)){
    knots <- x
  }
  else{
    knots <- seq(0, 10, length.out = k)
  }
  P_proposed <- compute_weights_precision(knots)
  D_proposed <- compute_weights_design(x,p)
  ## Design matrix for the global polynomial
  X = rbind(global_poly(x[1], p), D_proposed[,1:p])
  ## Design matrix for the spline basis weights
  B = as(local_poly(x = knots, refined_x = x, p = p),"dgTMatrix")
  
  tmbdat <- list(
    # Design matrix
    X = X,
    B = B,
    P = as(P_proposed,'dgTMatrix'),
    logPdet = as.numeric(determinant(P_proposed,logarithm = T)$modulus),
    # Response
    y = y,
    beta = beta,
    sigmaE = sigmaE,
    theta1 = theta
  )
  tmbparams <- list(
    W = c(rep(0, ncol(P_proposed)))
  )
  
  ff <- TMB::MakeADFun(
    data = tmbdat,
    parameters = tmbparams,
    # random = "W",
    DLL = "00_Gaussian_Smoothing_known",
    silent = TRUE
  )
  # Hessian not implemented for RE models
  ff$he <- function(w) numDeriv::jacobian(ff$gr,w)
  
  ## posterior mode and posterior precision matrix
  opt <- nlminb(start = ff$par, objective = ff$fn, gradient = ff$gr, hessian = ff$he, 
                control = list(eval.max = 200, iter.max = 200))
  mean = opt$par
  prec_matrix <- forceSymmetric(ff$he(opt$par))
  
  sample_weights <- rmvnp(n = nsamps, mu = mean, Omega = as.matrix(prec_matrix))
  samples_fun <- data.frame(x = x)
  samples_fun <- cbind(samples_fun, as.data.frame(as.matrix(B %*% t(sample_weights))))
  Blower <- as(local_poly(x = knots, refined_x = x, p = (p-1)),"dgTMatrix")
  
  sample_deriv <- data.frame(x = x)
  sample_deriv <- cbind(sample_deriv, as.data.frame(as.matrix(Blower %*% t(sample_weights))))
  if (optimizer == T){
    return(list(samples_fun = samples_fun, sample_deriv = sample_deriv, opt = opt))
  }
  
  ## Return posterior mean and posterior precision matrix
  return(list(samples_fun = samples_fun, sample_deriv = sample_deriv, mean = mean, prec = as.matrix(prec_matrix)))
}


#### Fit RW2 method, with known parameters, for Gaussian Data
Imple_RW2_known <- function(beta, sigmaE, theta, x, k, y, p = 3, optimizer = F, nsamps = 10000, diagonal = 1e-9){
  if(k < length(x)){
    knots <- seq(0, 10, length.out = k)
  }
  else {
    knots <- x
  }
  designB_rw2_reduced <- knots_RW2(x = knots, refined_x = x)
  d <- diff(knots)
  H <- compute_H_rue(d,n = length(knots))
  A <- compute_A(d,n = length(knots))
  Q_rw2_reduced <- as(t(H) %*% solve(A) %*% H, 'dgTMatrix')
  Q_rw2_reduced <- as(Q_rw2_reduced + Diagonal(k, x = diagonal), "dgTMatrix")
  D_proposed <- compute_weights_design(x,p)
  X = rbind(global_poly(x[1], p), D_proposed[,1:p])
  
  tmbdat <- list(
    X = as(X,"dgTMatrix"),
    B = as(designB_rw2_reduced,"dgTMatrix"),
    # Penalty(Precision) matrix
    P = as(as.matrix(Q_rw2_reduced),"dgTMatrix"),
    # Log determinant of penalty matrix (without the sigma part)
    # logPdet = as.numeric(sum(log(sort(all_values)[-(1:2)]))),
    logPdet = as.numeric(determinant(Q_rw2_reduced,logarithm = T)$modulus),
    # Response
    y = y,
    beta = beta,
    sigmaE = sigmaE,
    theta1 = theta
  )
  tmbparams <- list(
    W = c(rep(0, k))
  )
  ff <- TMB::MakeADFun(
    data = tmbdat,
    parameters = tmbparams,
    # random = "W",
    DLL = "00_Gaussian_Smoothing_known",
    silent = TRUE
  )
  # Hessian not implemented for RE models
  ff$he <- function(w) numDeriv::jacobian(ff$gr,w)
  ## posterior mode and posterior precision matrix
  opt <- nlminb(start = ff$par, objective = ff$fn, gradient = ff$gr, hessian = ff$he, 
                control = list(eval.max = 200, iter.max = 200))
  mean = opt$par
  prec_matrix <- forceSymmetric(ff$he(opt$par))
  
  sample_weights <- rmvnp(n = nsamps, mu = mean, Omega = as.matrix(prec_matrix))
  samples_fun <- data.frame(x = x)
  samples_fun <- cbind(samples_fun, as.data.frame(as.matrix(designB_rw2_reduced %*% t(sample_weights))))

  sample_deriv_raw <- samples_fun[,-1] %>% apply(MARGIN = 2, compute_numeric_deriv, h = mean(diff(x)), degree = 1)
  sample_deriv_raw <- rbind(rep(0,nsamps), sample_deriv_raw)
  sample_deriv <- data.frame(x = x)
  sample_deriv <- cbind(sample_deriv, sample_deriv_raw)
  if (optimizer == T){
    return(list(samples_fun = samples_fun, sample_deriv = sample_deriv, opt = opt))
  }
  
  ## Return posterior mean and posterior precision matrix
  return(list(samples_fun = samples_fun, sample_deriv = sample_deriv, mean = mean, prec = as.matrix(prec_matrix)))
}




##########################################################################################
##########################################################################################
### Posterior Inferences for the function:

### Construct posterior inference given an aghq-object using O-Spline: (for the observed grid)
extract_mean_interval <- function(quad, n_samp = 3000, level = 0.95, p = 2, x, type = "overall"){
  pos_samps <- sample_marginal(quad, n_samp)
  D_proposed <- compute_weights_design(x,p)
  X = rbind(c(1,rep(0,p-1)), D_proposed[,1:p])
  ## Design matrix for the spline basis weights
  B = rbind(rep(0, (length(x)-1)), D_proposed[,(p+1):ncol(D_proposed)])
  a <- (1 - level)/2
  if(type == "global"){
    x_samps <- pos_samps$samps[length(x):(nrow(pos_samps$samps)),]
    para_samps <- X %*% x_samps
    mean_x  <- apply(para_samps,1, mean)
    upper_x  <- apply(para_samps,1, quantile, p = (level + a))
    lower_x  <- apply(para_samps,1, quantile, p = a)
    return(data.frame(x = x, mean = mean_x, lower = lower_x, upper = upper_x))
  }
  else if(type == "local"){
    w_samps <- pos_samps$samps[1:(length(x) - 1),]
    nonpar_samps <- B %*% w_samps
    mean_w  <- apply(nonpar_samps,1, mean) 
    upper_w  <- apply(nonpar_samps,1, quantile, p = (level + a))
    lower_w  <- apply(nonpar_samps,1, quantile, p = a)
    return(data.frame(x = x, mean = mean_w, lower = lower_w, upper = upper_w))
    
  }
  else{
    fitted_samps <- cbind(B,X) %*% pos_samps$samps 
    mean_f  <- apply(fitted_samps,1, mean) 
    upper_f  <- apply(fitted_samps,1, quantile, p = (level + a))
    lower_f  <- apply(fitted_samps,1, quantile, p = a)
    return(data.frame(x = x, mean = mean_f, lower = lower_f, upper = upper_f))
  }
}

### Construct posterior inference given samples
extract_mean_interval_given_samps <- function(samps, n_samp = 3000, level = 0.95){
    x <- samps[,1]
    samples <- samps[,-1]
    result <- data.frame(x = x)
    alpha <- 1-level
    result$plower <- as.numeric(apply(samples, MARGIN = 1, quantile, p = (alpha/2)))
    result$pupper <- as.numeric(apply(samples, MARGIN = 1, quantile, p = (level + (alpha/2))))
    result$mean <- as.numeric(apply(samples, MARGIN = 1, mean))
    result
}



### Construct posterior inference given an aghq-object using O-Spline: (for the refined grid), given overdispersion
extract_mean_interval_refined_overd <- function(quad, n_samp = 3000, level = 0.95, p = 2, x, refined_x, type = 'overall', Xf = NULL){
  pos_samps <- sample_marginal(quad, n_samp)
  D_proposed <- compute_weights_design(refined_x,p)
  X = rbind(global_poly(refined_x[1], p = p), D_proposed[,1:p])
  X = cbind(X, Xf)
  ## Design matrix for the spline basis weights
  B = as(local_poly(x, refined_x = refined_x, p = p),"dgTMatrix")
  I = diag(1, nrow = length(refined_x), ncol = length(refined_x))
  
  ## Overall matrix:
  M <- cbind(B,X,I)
  
  a <- (1 - level)/2
  if(type == "overall"){
    fitted_samps_PTR3 <- M %*% pos_samps$samps
    mean_f  <- apply(fitted_samps_PTR3,1, mean) 
    upper_f  <- apply(fitted_samps_PTR3,1, quantile, p = (level + a))
    lower_f  <- apply(fitted_samps_PTR3,1, quantile, p = a)
    return(data.frame(x = refined_x, mean = mean_f, lower = lower_f, upper = upper_f))
  }
  else if(type == "function"){
    fitted_samps_PTR3 <- cbind(B,X[,1:p]) %*% pos_samps$samps[-((p + ncol(B) + 1): nrow(pos_samps$samps)), ]
    mean_f  <- apply(fitted_samps_PTR3,1, mean) 
    upper_f  <- apply(fitted_samps_PTR3,1, quantile, p = (level + a))
    lower_f  <- apply(fitted_samps_PTR3,1, quantile, p = a)
    return(data.frame(x = refined_x, mean = mean_f, lower = lower_f, upper = upper_f))
  }
  else if(type == "local"){
    fitted_samps_PTR3 <- B %*% pos_samps$samps[1:(length(x) - 1),]
    mean_f  <- apply(fitted_samps_PTR3,1, mean) 
    upper_f  <- apply(fitted_samps_PTR3,1, quantile, p = (level + a))
    lower_f  <- apply(fitted_samps_PTR3,1, quantile, p = a)
    return(data.frame(x = refined_x, mean = mean_f, lower = lower_f, upper = upper_f))
  }
  else{
    fitted_samps_PTR3 <- X[,1:p] %*% pos_samps$samps[(ncol(B)+1):(ncol(B) + p),]
    mean_f  <- apply(fitted_samps_PTR3,1, mean) 
    upper_f  <- apply(fitted_samps_PTR3,1, quantile, p = (level + a))
    lower_f  <- apply(fitted_samps_PTR3,1, quantile, p = a)
    return(data.frame(x = refined_x, mean = mean_f, lower = lower_f, upper = upper_f))
  }
}
extract_samples_Ospline_refined_overd <- function(quad, n_samp = 3000, p = 2, x, refined_x, type = 'overall', Xf = NULL){
  pos_samps <- sample_marginal(quad, n_samp)
  D_proposed <- compute_weights_design(refined_x,p)
  X = cbind(rbind(global_poly(refined_x[1], p = p), D_proposed[,1:p]), Xf)
  ## Design matrix for the spline basis weights
  B = as(local_poly(x, refined_x = refined_x, p = p),"dgTMatrix")
  I = diag(1, nrow = length(refined_x), ncol = length(refined_x))
  result <- list()
  if(type == "overall"){
    fitted_samps_PTR3 <- cbind(B,X,I) %*% pos_samps$samps
    result <- cbind(x = refined_x, data.frame(as.matrix(fitted_samps_PTR3)))
    return(result)
  }
  else if(type == "function"){
    fitted_samps_PTR3 <- cbind(B,X) %*% pos_samps$samps[-((ncol(X) + ncol(B) + 1): nrow(pos_samps$samps)), ]
    result <- cbind(x = refined_x, data.frame(as.matrix(fitted_samps_PTR3)))
    return(result)
  }
  else{
    fitted_samps_PTR3 <- X %*% pos_samps$samps[length(x):(length(x) + (p-1)),]
    result <- cbind(x = refined_x, data.frame(as.matrix(fitted_samps_PTR3)))
    return(result)
  }
}
extract_deriv_samples_OSpline_overd <- function(quad, n_samp = 3000, p = 2, x, refined_x, degree){
  pos_samps <- sample_marginal(quad, n_samp)
  D_proposed <- compute_weights_design(refined_x,p)
  X = rbind(global_poly(refined_x[1], p = p), D_proposed[,1:p])
  X <- as.matrix(X[,1:(p-degree)])
  for (i in 1:ncol(X)) {
    X[,i] <- (factorial(i + degree - 1)/factorial(i-1))* X[,i]
  }
  ## Design matrix for the spline basis weights
  B = as(local_poly(x, refined_x = refined_x, p = (p-degree)),"dgTMatrix")
  fitted_samps_deriv <- X %*% pos_samps$samps[(length(x) + degree):(length(x) + p - 1), ] + B %*% pos_samps$samps[1:(length(x)-1), ]
  result <- cbind(x = refined_x, data.frame(as.matrix(fitted_samps_deriv)))
  result
}


### Construct posterior inference of RW2 given an aghq-object: (for the observed grid)
extract_mean_interval_rw2 <- function(quad, n_samp = 3000, level = 0.95, x){
  pos_samps_RW2F <- sample_marginal(quad, n_samp)
  a <- (1-level)/2
  n <- length(x)
  designB <- as(Diagonal(n,1), "matrix")
  fitted_samps_RW2F <- designB %*% pos_samps_RW2F$samps
  mean_fitted_RW2F  <- apply(fitted_samps_RW2F,1, mean)
  upper_fitted_RW2F  <- apply(fitted_samps_RW2F,1, quantile, p = level + a)
  lower_fitted_RW2F  <- apply(fitted_samps_RW2F,1, quantile, p = a)
  return(data.frame(x = x, mean = mean_fitted_RW2F, lower = lower_fitted_RW2F, upper = upper_fitted_RW2F))
}

### Construct posterior inference of RW2 given an aghq-object: (for the refined grid)
extract_mean_interval_rw2_refined <- function(quad, n_samp = 3000, level = 0.95, x, refined_x){
  designB <- knots_RW2(x, refined_x)
  pos_samps_RW2F <- sample_marginal(quad, n_samp)
  a <- (1-level)/2
  fitted_samps_RW2F <- designB %*% pos_samps_RW2F$samps
  mean_fitted_RW2F  <- apply(fitted_samps_RW2F,1, mean)
  upper_fitted_RW2F  <- apply(fitted_samps_RW2F,1, quantile, p = level + a)
  lower_fitted_RW2F  <- apply(fitted_samps_RW2F,1, quantile, p = a)
  return(data.frame(x = refined_x, mean = mean_fitted_RW2F, lower = lower_fitted_RW2F, upper = upper_fitted_RW2F))
}

#### Construct posterior inference of Proposed B-spline method given an aghq-object:
extract_mean_interval_Bspline <- function(mod, n_samp = 3000, level = 0.95, x, p = 2, min, max){
  m <- 2*p - 1
  pos_samps <- sample_marginal(mod, n_samp)
  q <- nrow(pos_samps$samps)
  if(is.null(min)|is.null(max)){
    min = range(x)[1]
    max = range(x)[2]
  }
  B <- get_design_Bspline(x,q,min = min, max = max, p = p)
  X <- global_poly(x, p = p)
  fitted_samps <- cbind(B,X) %*% pos_samps$samps
  mean_fitted  <- apply(fitted_samps,1, mean)
  a <- (1-level)/2
  upper_fitted  <- apply(fitted_samps,1, quantile, p = level + a)
  lower_fitted  <- apply(fitted_samps,1, quantile, p = a)
  return(data.frame(x = x, mean = mean_fitted, lower = lower_fitted, upper = upper_fitted))
}

#### Construct posterior inference of polynomial method given an aghq-object:
extract_mean_interval_Poly <- function(mod, n_samp = 3000, level = 0.95, x, p = 2){
  pos_samps <- sample_marginal(mod, n_samp)
  q <- nrow(pos_samps$samps)
  X <- global_poly(x, p = (p+1))
  fitted_samps <- X %*% pos_samps$samps
  mean_fitted  <- apply(fitted_samps,1, mean)
  a <- (1-level)/2
  upper_fitted  <- apply(fitted_samps,1, quantile, p = level + a)
  lower_fitted  <- apply(fitted_samps,1, quantile, p = a)
  return(data.frame(x = x, mean = mean_fitted, lower = lower_fitted, upper = upper_fitted))
}

#### Construct posterior inference of exact method given an aghq-object:
extract_mean_interval_WP <- function(quad, n_samp = 3000, level = 0.95, x, p = 2){
  pos_samps_WP <- sample_marginal(quad, n_samp)
  a <- (1-level)/2
  n <- length(x)
  designB <- Compute_design_Aug(svec = x[-1],p)
  designB <- rbind(rep(0, ncol(designB)), designB)
  x_standardized <- x - x[1]
  D_proposed <- compute_weights_design(x_standardized,p)
  ## Design matrix for the global polynomial
  designX <- rbind(c(1,rep(0,p-1)), D_proposed[,1:p])
  fitted_samps <- cbind(designB,designX) %*% pos_samps_WP$samps 
  mean_fitted  <- apply(fitted_samps,1, mean)
  upper_fitted  <- apply(fitted_samps,1, quantile, p = level + a)
  lower_fitted  <- apply(fitted_samps,1, quantile, p = a)
  return(data.frame(x = x, mean = mean_fitted, lower = lower_fitted, upper = upper_fitted))
}

#### Extract samples from posterior fitted from Imple_WP_known
extract_samples_WP_known <- function(result, beta, x, nsam = 3000, p = 3, target = "function"){
  samps <- rmvnp(n = nsam, mu = as.numeric(result$mean), Omega = as.matrix(result$prec))
  samps <- cbind(matrix(0,nrow = nsam, ncol = p), samps)
  if(target == "function"){
    samps_func <- data.frame(x = x)
    indx <- seq(1,(p*length(x)), by = p)
    samps_func <- cbind(samps_func, t(samps[,indx]))
    samps_func
  }
  else if(target == "derivative"){
    samps_deriv <- data.frame(x = x)
    indx <- seq(1,(p*length(x)), by = p) + 1
    samps_deriv <- cbind(samps_deriv, t(samps[,indx]))
    samps_deriv
  }
  else{
    stop(errorCondition("No selected target identified."))
  }
}



##############################################################################################
##############################################################################################
#### Extract Samples and global region:


### Extract samples from the fitted O-Spline object: (with refined resolution)
extract_samples_Ospline_refined <- function(quad, n_samp = 3000, p = 2, x, refined_x, type = 'overall'){
  pos_samps <- sample_marginal(quad, n_samp)
  D_proposed <- compute_weights_design(refined_x,p)
  X = rbind(global_poly(refined_x[1], p = p), D_proposed[,1:p])
  ## Design matrix for the spline basis weights
  B = as(local_poly(x, refined_x = refined_x, p = p),"dgTMatrix")
  result <- list()
  if(type == "overall"){
    fitted_samps_PTR3 <- cbind(B,X) %*% pos_samps$samps
    result <- cbind(x = refined_x, data.frame(as.matrix(fitted_samps_PTR3)))
    return(result)
  }
  else if(type == "local"){
    fitted_samps_PTR3 <- B %*% pos_samps$samps[1:(length(x) - 1),]
    result <- cbind(x = refined_x, data.frame(as.matrix(fitted_samps_PTR3)))
    return(result)
  }
  else{
    fitted_samps_PTR3 <- X %*% pos_samps$samps[length(x):(nrow(pos_samps$samps)),]
    result <- cbind(x = refined_x, data.frame(as.matrix(fitted_samps_PTR3)))
    return(result)
  }
}
### Extract samples of derivatives from the fitted O-Spline object: (with refined resolution)
extract_deriv_samples_OSpline <- function(quad, n_samp = 3000, p = 2, x, refined_x, degree){
  pos_samps <- sample_marginal(quad, n_samp)
  D_proposed <- compute_weights_design(refined_x,p)
  X = rbind(global_poly(refined_x[1], p = p), D_proposed[,1:p])
  X <- as.matrix(X[,1:(p-degree)])
  for (i in 1:ncol(X)) {
    X[,i] <- (factorial(i + degree - 1)/factorial(i-1))* X[,i]
  }
  ## Design matrix for the spline basis weights
  B = as(local_poly(x, refined_x = refined_x, p = (p-degree)),"dgTMatrix")
  fitted_samps_deriv <- X %*% pos_samps$samps[(length(x) + degree):nrow(pos_samps$samps), ] + B %*% pos_samps$samps[1:(length(x)-1), ]
  result <- cbind(x = refined_x, data.frame(as.matrix(fitted_samps_deriv)))
  result
}

### Extract samples of RW
extract_RW_samples <- function(quad, n_samp = 3000, x, refined_x){
  n <- length(x)
  pos_samps_RW2F <- sample_marginal(quad, n_samp)
  designB <- knots_RW2(x, refined_x)
  pos_samps_RW2F <- sample_marginal(quad, n_samp)
  fitted_samps_RW2F <- designB %*% pos_samps_RW2F$samps
  return(cbind(refined_x, fitted_samps_RW2F))
}

####### When there are different initial conditions

### Construct posterior inference given an aghq-object using O-Spline with condition: (for the observed grid)
extract_mean_interval_cond <- function(quad, n_samp = 3000, level = 0.95, p = 2, x, type = "overall", s = NULL){
  if(is.null(s)){
    s <- ceiling(length(x)/2)
  }
  pos_samps <- sample_marginal(quad, n_samp)
  
  ## Design matrix for the spline basis weights
  D_proposed <- condition_local_design(x, s = s, p = p)
  ## Design matrix for the global polynomial
  X = condition_global_design(x, s = s, p = p)
  B = as(D_proposed,'dgTMatrix')
  
  
  a <- (1 - level)/2
  if(type == "global"){
    x_samps <- pos_samps$samps[length(x):(nrow(pos_samps$samps)),]
    para_samps <- X %*% x_samps
    mean_x  <- apply(para_samps,1, mean)
    upper_x  <- apply(para_samps,1, quantile, p = (level + a))
    lower_x  <- apply(para_samps,1, quantile, p = a)
    return(data.frame(x = x, mean = mean_x, lower = lower_x, upper = upper_x))
  }
  else if(type == "local"){
    w_samps <- pos_samps$samps[1:(length(x) - 1),]
    nonpar_samps <- B %*% w_samps
    mean_w  <- apply(nonpar_samps,1, mean) 
    upper_w  <- apply(nonpar_samps,1, quantile, p = (level + a))
    lower_w  <- apply(nonpar_samps,1, quantile, p = a)
    return(data.frame(x = x, mean = mean_w, lower = lower_w, upper = upper_w))
    
  }
  else{
    fitted_samps <- cbind(B,X) %*% pos_samps$samps 
    mean_f  <- apply(fitted_samps,1, mean) 
    upper_f  <- apply(fitted_samps,1, quantile, p = (level + a))
    lower_f  <- apply(fitted_samps,1, quantile, p = a)
    return(data.frame(x = x, mean = mean_f, lower = lower_f, upper = upper_f))
  }
}

### Construct posterior inference given an aghq-object using O-Spline: (for the refined grid)
extract_mean_interval_refined_cond <- function(quad, n_samp = 3000, level = 0.95, p = 2, x, refined_x, type = 'overall', s = NULL){
  if(is.null(s)){
    s <- ceiling(length(x)/2)
  }
  pos_samps <- sample_marginal(quad, n_samp)
  D_proposed <- condition_local_poly(knots = x, refined_x = refined_x, s = s, p = p)
  ## Design matrix for the global polynomial
  X = condition_global_poly(knots = x, refined_x = refined_x, s = s, p = p)
  X = as(X,'dgTMatrix')
  ## Design matrix for the spline basis weights
  B = as(D_proposed,'dgTMatrix')
  
  a <- (1 - level)/2
  if(type == "overall"){
    fitted_samps_PTR3 <- cbind(B,X) %*% pos_samps$samps
    mean_f  <- apply(fitted_samps_PTR3,1, mean) 
    upper_f  <- apply(fitted_samps_PTR3,1, quantile, p = (level + a))
    lower_f  <- apply(fitted_samps_PTR3,1, quantile, p = a)
    return(data.frame(x = refined_x, mean = mean_f, lower = lower_f, upper = upper_f))
  }
  else if(type == "local"){
    fitted_samps_PTR3 <- B %*% pos_samps$samps[1:(length(x) - 1),]
    mean_f  <- apply(fitted_samps_PTR3,1, mean) 
    upper_f  <- apply(fitted_samps_PTR3,1, quantile, p = (level + a))
    lower_f  <- apply(fitted_samps_PTR3,1, quantile, p = a)
    return(data.frame(x = refined_x, mean = mean_f, lower = lower_f, upper = upper_f))
  }
  else{
    fitted_samps_PTR3 <- X %*% pos_samps$samps[length(x):(nrow(pos_samps$samps)),]
    mean_f  <- apply(fitted_samps_PTR3,1, mean) 
    upper_f  <- apply(fitted_samps_PTR3,1, quantile, p = (level + a))
    lower_f  <- apply(fitted_samps_PTR3,1, quantile, p = a)
    return(data.frame(x = refined_x, mean = mean_f, lower = lower_f, upper = upper_f))
  }
}

### Construct exact estimate for the derivative, with additional condition:
extract_deriv_OSpline_cond <- function(quad, n_samp = 3000, level = 0.95, p = 2, x, refined_x, degree, s = NULL){
  if(is.null(s)){
    s <- ceiling(length(x)/2)
  }
  pos_samps <- sample_marginal(quad, n_samp)
  D_proposed <- condition_local_poly_for_deriv(knots = x, refined_x = refined_x, s = s, p = (p-degree), degree = degree)
  ## Design matrix for the global polynomial
  X = condition_global_poly(knots = x, refined_x = refined_x, s = s, p = p)
  X = as(X,'dgTMatrix')
  ## Design matrix for the spline basis weights
  B = as(D_proposed,'dgTMatrix')
  X <- as.matrix(X[,1:(p-degree)])
  for (i in 1:ncol(X)) {
    X[,i] <- (factorial(i + degree - 1)/factorial(i-1))* X[,i]
  }
  ## Design matrix for the spline basis weights
  a <- (1 - level)/2
  fitted_samps_deriv <- X %*% pos_samps$samps[(length(x) + degree):nrow(pos_samps$samps), ] + B %*% pos_samps$samps[1:(length(x)-1), ]
  mean_fitted  <- apply(fitted_samps_deriv,1, mean)
  upper_fitted  <- apply(fitted_samps_deriv,1, quantile, p = level + a)
  lower_fitted  <- apply(fitted_samps_deriv,1, quantile, p = a)
  return(data.frame(x = refined_x, mean = mean_fitted, lower = lower_fitted, upper = upper_fitted))
}


## Obtain samples for the function
extract_samples_Ospline_refined_cond <- function(quad, n_samp = 3000, p = 2, x, refined_x, type = 'overall', s = NULL){
  if(is.null(s)){
    s <- ceiling(length(x)/2)
  }
  pos_samps <- sample_marginal(quad, n_samp)
  D_proposed <- condition_local_poly(knots = x, refined_x = refined_x, s = s, p = p)
  ## Design matrix for the global polynomial
  X = condition_global_poly(knots = x, refined_x = refined_x, s = s, p = p)
  X = as(X,'dgTMatrix')
  ## Design matrix for the spline basis weights
  B = as(D_proposed,'dgTMatrix')
  if(type == "overall"){
    fitted_samps_PTR3 <- cbind(B,X) %*% pos_samps$samps
    result <- cbind(x = refined_x, data.frame(as.matrix(fitted_samps_PTR3)))
    return(result)
  }
  else if(type == "local"){
    fitted_samps_PTR3 <- B %*% pos_samps$samps[1:(length(x) - 1),]
    result <- cbind(x = refined_x, data.frame(as.matrix(fitted_samps_PTR3)))
    return(result)
  }
  else{
    fitted_samps_PTR3 <- X %*% pos_samps$samps[length(x):(nrow(pos_samps$samps)),]
    result <- cbind(x = refined_x, data.frame(as.matrix(fitted_samps_PTR3)))
    return(result)
  }
}


## Obtain samples for the derivatives
extract_deriv_samples_OSpline_cond <- function(quad, n_samp = 3000, p = 2, x, refined_x, degree, s = NULL){
  if(is.null(s)){
    s <- ceiling(length(x)/2)
  }
  pos_samps <- sample_marginal(quad, n_samp)
  D_proposed <- condition_local_poly_for_deriv(knots = x, refined_x = refined_x, s = s, p = (p-degree), degree = degree)
  ## Design matrix for the global polynomial
  X = condition_global_poly(knots = x, refined_x = refined_x, s = s, p = p)
  X = as(X,'dgTMatrix')
  ## Design matrix for the spline basis weights
  B = as(D_proposed,'dgTMatrix')
  X <- as.matrix(X[,1:(p-degree)])
  for (i in 1:ncol(X)) {
    X[,i] <- (factorial(i + degree - 1)/factorial(i-1))* X[,i]
  }
  ## Design matrix for the spline basis weights
  fitted_samps_deriv <- X %*% pos_samps$samps[(length(x) + degree):nrow(pos_samps$samps), ] + B %*% pos_samps$samps[1:(length(x)-1), ]
  result <- cbind(x = refined_x, data.frame(as.matrix(fitted_samps_deriv)))
  result
}




### Extract samples from the fitted B-Spline object: (with refined resolution)
extract_samples_Bspline <- function(mod, n_samp = 3000, x, p = 2, min, max){
  m <- 2*p - 1
  pos_samps <- sample_marginal(mod, n_samp)
  q <- nrow(pos_samps$samps)
  if(is.null(min)|is.null(max)){
    min = range(x)[1]
    max = range(x)[2]
  }
  B <- get_design_Bspline(x,q,min = min, max = max, p = p)
  X <- global_poly(x, p = p)
  fitted_samps <- cbind(B,X) %*% pos_samps$samps
  result <- cbind(x = x, data.frame(as.matrix(fitted_samps)))
  result
}

### Extract samples of derivatives from the fitted B-Spline object: (with refined resolution)
extract_deriv_samples_Bspline <- function(mod, n_samp = 3000, x, p = 2, min, max, degree){
  m <- 2*p - 1
  pos_samps <- sample_marginal(mod, n_samp)
  q <- nrow(pos_samps$samps)
  if(is.null(min)|is.null(max)){
    min = range(x)[1]
    max = range(x)[2]
  }
  B <- get_design_Bspline(x,q,min = min, max = max, p = p)
  X <- global_poly(x, p = p)
  fitted_samps <- cbind(B,X) %*% pos_samps$samps
  h <- mean(diff(x))
  fitted_samps_deriv <- apply(fitted_samps,2, compute_numeric_deriv, h = h, degree = degree)
  result <- cbind(x = x[-c(1:degree)], data.frame(as.matrix(fitted_samps_deriv)))
  result
}

### Extract samples from the fitted WP object:
extract_samples_WP <- function(quad, n_samp = 3000, x, p = 2){
  pos_samps_WP <- sample_marginal(quad, n_samp)
  n <- length(x)
  designB <- Compute_design_Aug(svec = x[-1],p)
  designB <- rbind(rep(0, ncol(designB)), designB)
  x_standardized <- x - x[1]
  D_proposed <- compute_weights_design(x_standardized,p)
  ## Design matrix for the global polynomial
  designX <- rbind(c(1,rep(0,p-1)), D_proposed[,1:p])
  fitted_samps <- cbind(designB,designX) %*% pos_samps_WP$samps 
  result <- cbind(x = x, data.frame(as.matrix(fitted_samps)))
  result
}

### Extract samples of derivatives from the fitted WP object:
extract_deriv_samples_WP <- function(quad, n_samp = 3000, x, degree, p = 2){
  pos_samps <- sample_marginal(quad, n_samp)
  B <- Compute_design_Aug_deriv(svec = x[-1],p, degree = degree)
  B <- rbind(rep(0, ncol(B)), B)
  x_standardized <- x - x[1]
  D_proposed <- compute_weights_design(x_standardized,p)
  X = rbind(global_poly(0, p = p), D_proposed[,1:p])
  X <- as.matrix(X[,1:(p-degree)])
  for (i in 1:ncol(X)) {
    X[,i] <- (factorial(i + degree - 1)/factorial(i-1))* X[,i]
  }
  fitted_samps_deriv <- X %*% pos_samps$samps[(ncol(B) + 1 + degree):nrow(pos_samps$samps), ] + B %*% pos_samps$samps[1:(ncol(B)), ]
  result <- cbind(x = x, data.frame(as.matrix(fitted_samps_deriv)))
  result
}

### Construct posterior global credible region (GCR) given an aghq-object of any class:
extract_GCR <- function(fitted_samples, pointwise = F, level = 0.95, method = "excursion", GET = 'erl'){
  if(method == "excursion"){
    x <- fitted_samples$x
    samples <- fitted_samples[,-1]
    CR <- quiet(simconf.mc(samples = samples, alpha = (1-level)))
    result <- data.frame(x = x, lower = CR$a, upper = CR$b)
    if(pointwise == T){
      result$plower <- CR$a.marginal
      result$pupper <- CR$b.marginal
    }
    result 
  }
  else if(method == "GET"){
    x <- fitted_samples$x
    samples <- fitted_samples[,-1]
    cset_g <- create_curve_set(list(r = x, obs = as.matrix(samples)))
    CR <- central_region(curve_sets = cset_g, type = GET, coverage = level)
 
    result <- data.frame(x = x, lower = CR$lo, upper = CR$hi)
    if(pointwise == T){
      alpha <- 1-level
      result$plower <- as.numeric(apply(samples, MARGIN = 1, quantile, p = (alpha/2)))
      result$pupper <- as.numeric(apply(samples, MARGIN = 1, quantile, p = (level + (alpha/2))))
    }
    result 
  }
  else{
    stop(errorCondition("Unknown method"))
  }
}


### Construct posterior global credible region (GCR) given an aghq-object of any class:
extract_GCR_known <- function(mu, Q, x, pointwise = F, level = 0.95){
  CR <- quiet(simconf(mu = mu, Q = Q, alpha = (1-level)))
  result <- data.frame(x = x, lower = CR$a, upper = CR$b)
  if(pointwise == T){
    result$plower <- CR$a.marginal
    result$pupper <- CR$b.marginal
  }
  result
}




##########################################################################################
##########################################################################################
### Posterior Inferences for the function (numeric) derivative:
extract_deriv_OSpline_numeric <- function(quad, n_samp = 3000, level = 0.95, p = 2, x, refined_x, degree){
  pos_samps <- sample_marginal(quad, n_samp)
  D_proposed <- compute_weights_design(refined_x,p)
  X = rbind(global_poly(refined_x[1], p = p), D_proposed[,1:p])
  ## Design matrix for the spline basis weights
  B = as(local_poly(x, refined_x = refined_x, p = p),"dgTMatrix")
  a <- (1 - level)/2
  fitted_samps <- cbind(B,X) %*% pos_samps$samps
  h <- mean(diff(refined_x))
  fitted_samps_deriv <- apply(fitted_samps,2, compute_numeric_deriv, h = h, degree = degree)
  mean_fitted  <- apply(fitted_samps_deriv,1, mean)
  upper_fitted  <- apply(fitted_samps_deriv,1, quantile, p = level + a)
  lower_fitted  <- apply(fitted_samps_deriv,1, quantile, p = a)
  return(data.frame(x = refined_x[-c(1:degree)], mean = mean_fitted, lower = lower_fitted, upper = upper_fitted))
}
extract_deriv_OSpline <- function(quad, n_samp = 3000, level = 0.95, p = 2, x, refined_x, degree){
  pos_samps <- sample_marginal(quad, n_samp)
  D_proposed <- compute_weights_design(refined_x,p)
  X = rbind(global_poly(refined_x[1], p = p), D_proposed[,1:p])
  X <- as.matrix(X[,1:(p-degree)])
  for (i in 1:ncol(X)) {
    X[,i] <- (factorial(i + degree - 1)/factorial(i-1))* X[,i]
  }
  ## Design matrix for the spline basis weights
  B = as(local_poly(x, refined_x = refined_x, p = (p-degree)),"dgTMatrix")
  a <- (1 - level)/2
  fitted_samps_deriv <- X %*% pos_samps$samps[(length(x) + degree):nrow(pos_samps$samps), ] + B %*% pos_samps$samps[1:(length(x)-1), ]
  mean_fitted  <- apply(fitted_samps_deriv,1, mean)
  upper_fitted  <- apply(fitted_samps_deriv,1, quantile, p = level + a)
  lower_fitted  <- apply(fitted_samps_deriv,1, quantile, p = a)
  return(data.frame(x = refined_x, mean = mean_fitted, lower = lower_fitted, upper = upper_fitted))
}
extract_deriv_BSpline_numeric <- function(mod, n_samp = 3000, level = 0.95, x, p = 2, min, max, degree){
  m <- 2*p - 1
  a <- (1 - level)/2
  pos_samps <- sample_marginal(mod, n_samp)
  q <- nrow(pos_samps$samps)
  if(is.null(min)|is.null(max)){
    min = range(x)[1]
    max = range(x)[2]
  }
  B <- get_design_Bspline(x,q,min = min, max = max, p = p)
  X <- global_poly(x, p = p)
  fitted_samps <- cbind(B,X) %*% pos_samps$samps
  h <- mean(diff(x))
  fitted_samps_deriv <- apply(fitted_samps,2, compute_numeric_deriv, h = h, degree = degree)
  mean_fitted  <- apply(fitted_samps_deriv,1, mean)
  upper_fitted  <- apply(fitted_samps_deriv,1, quantile, p = level + a)
  lower_fitted  <- apply(fitted_samps_deriv,1, quantile, p = a)
  return(data.frame(x = x[-c(1:degree)], mean = mean_fitted, lower = lower_fitted, upper = upper_fitted))
}
extract_deriv_RW2_numeric <- function(quad, n_samp = 3000, level = 0.95, x, refined_x, degree){
  designB <- knots_RW2(x, refined_x)
  pos_samps_RW2F <- sample_marginal(quad, n_samp)
  a <- (1-level)/2
  fitted_samps_RW2F <- designB %*% pos_samps_RW2F$samps
  h <- mean(diff(refined_x))
  fitted_samps_RW2F <- apply(fitted_samps_RW2F,2, compute_numeric_deriv, h = h, degree = degree)
  mean_fitted_RW2F  <- apply(fitted_samps_RW2F,1, mean)
  upper_fitted_RW2F  <- apply(fitted_samps_RW2F,1, quantile, p = level + a)
  lower_fitted_RW2F  <- apply(fitted_samps_RW2F,1, quantile, p = a)
  return(data.frame(x = refined_x[-c(1:degree)], mean = mean_fitted_RW2F, lower = lower_fitted_RW2F, upper = upper_fitted_RW2F))
}
extract_deriv_WP <- function(quad, n_samp = 3000, level = 0.95, x, degree, p = 2){
  pos_samps <- sample_marginal(quad, n_samp)
  B <- Compute_design_Aug_deriv(svec = x[-1],p, degree = degree)
  B <- rbind(rep(0, ncol(B)), B)
  x_standardized <- x - x[1]
  D_proposed <- compute_weights_design(x_standardized,p)
  X = rbind(global_poly(0, p = p), D_proposed[,1:p])
  X <- as.matrix(X[,1:(p-degree)])
  for (i in 1:ncol(X)) {
    X[,i] <- (factorial(i + degree - 1)/factorial(i-1))* X[,i]
  }
  a <- (1 - level)/2
  fitted_samps_deriv <- X %*% pos_samps$samps[(ncol(B) + 1 + degree):nrow(pos_samps$samps), ] + B %*% pos_samps$samps[1:(ncol(B)), ]
  mean_fitted  <- apply(fitted_samps_deriv,1, mean)
  upper_fitted  <- apply(fitted_samps_deriv,1, quantile, p = level + a)
  lower_fitted  <- apply(fitted_samps_deriv,1, quantile, p = a)
  return(data.frame(x = x, mean = mean_fitted, lower = lower_fitted, upper = upper_fitted))
}





##################################################################################################
#### Model Selection:
Compute_log_Bayes_Factor <- function(M1,M2){
  get_log_normconst(M1) - get_log_normconst(M2)
}





###################################################################################################
### Functions to simulate from GP:
sim_BM <- function(mesh_size, max_t){
  t <- seq(0, max_t, by = mesh_size)
  n <- length(t)
  noises <- rnorm((n-1), mean = 0, sd = 1)
  ft <- sqrt(mesh_size)*cumsum(c(0,noises))
  data.frame(t = t, ft = ft)
}
sim_Wp <- function(mesh_size, max_t, p){
  Wp <- numeric()
  t <- seq(0, max_t, by = mesh_size)
  n <- length(t)
  noises <- rnorm((n-1), mean = 0, sd = 1)
  if(p == 1){
    return(sim_BM(mesh_size, max_t))
  }
  else{
    ft <- sqrt(mesh_size)*cumsum(c(0,noises))
    for (i in 2:p) {
      n <- length(length(ft))
      ft <- mesh_size*cumsum(c(0,ft[-n]))
    } 
    data.frame(t = t, ft = ft)
  }
}

sim_Wp_Var <- function(t = NULL, mesh_size = 0.01, max_t = 10, p, sigmaS = 1){
  if(is.null(t)){
    t <- seq(0, max_t, by = mesh_size)
  }
  n <- length(t) - 1
  matT <- Compute_Ti(svec = t[-1], p = p, i = 1)
  Sig <- Compute_Ci(svec = t[-1], p = p, i = 1, is.cov = T)
  sample_path <- sigmaS * tsDyn::VAR.sim(B = matT, lag = 1, n = n, starting = NULL, varcov = Sig, include = "none", returnStarting = T)
  result <- data.frame(t = t, sample_path)
  for (i in 2:ncol(result)) {
    names(result)[i] <- paste0("IW",(p-i+2))
  }
  result
}


sigdigs <- function(x){
  x <- as.character(x)
  x <- sub("\\.", "", x)
  x <- gsub("(^0+|0+$)", "", x)
  nchar(x)
}
evaluate_GP <- function(x, GP){
  signif1 <- max(sigdigs(x))
  signif2 <- max(sigdigs(GP$t))
  signif3 <- max(signif1, signif2)
  if(all(round(x,signif3) %in% (round(GP$t,signif3))) == T){
    return(GP$ft[sort(which((round(GP$t,signif3) %in% round(x,signif3))))])
  }
  else{
    print("Error, X is not included in the discretization")
  }
}
evaluate_GP_Var <- function(x, GP){
  signif1 <- max(sigdigs(x))
  signif2 <- max(sigdigs(GP$t))
  signif3 <- max(signif1, signif2)
  if(all(round(x,signif3) %in% (round(GP$t,signif3))) == T){
    return(GP[sort(which((round(GP$t,signif3) %in% round(x,signif3)))), ])
  }
  else{
    print("Error, X is not included in the discretization")
  }
}

true_ploy <- function(x, folds = 2){
  func <- 0
  for (i in 1:folds) {
    func <- func + rnorm(mean = 0, sd = 1, n = 1)*x^(i-1)
  }
  func
}
true_deviation <- function(x, max_t, folds = 2, weight = 0.01, mesh_size = 0.0001){
  GP <- sim_Wp(mesh_size = mesh_size, p = folds, max_t = max_t)
  GP$ft <- weight * GP$ft
  evaluate_GP(x, GP)
}


##### Functions for simulation to compare knot effect, reporting a sequence of posterior mean estimate for different knots number
sim_PT_knots <- function(y,x,k_vec,min = 0, max = 100,p, placement = "equal", likelihood = "Gaussian", aghq_k = 7){
  PTRmod <- Imple_BayesRegression(y = y,x = x, p = p, likelihood = likelihood)
  PTRpos <- extract_mean_interval(PTRmod, x = x, p = p)
  result <- data.frame(k = n, x = x ,mean = PTRpos$mean)
  if(placement == "equal"){
    for (i in 1:length(k_vec)) {
      k <- k_vec[i]
      selected_x <- seq(min,max, length.out = k)
      PTRmod <- Imple_BayesRegression(y = y,knots = selected_x, x = x, p = p, likelihood = likelihood)
      PTRpos <- extract_mean_interval_refined(PTRmod, x = selected_x, refined_x = x, p = p)
      result <- rbind(result,data.frame(k = k, x = x ,mean = PTRpos$mean))
    } 
  }
  else if(placement == "quantile"){
    for (i in 1:length(k_vec)) {
      k <- k_vec[i]
      selected_x <- quantile(x, p = seq(0,1,length.out = k))
      PTRmod <- Imple_BayesRegression(y = y,knots = selected_x, x = x, p = p, likelihood = likelihood)
      PTRpos <- extract_mean_interval_refined(PTRmod, x = selected_x, refined_x = x, p = p)
      result <- rbind(result,data.frame(k = k, x = x ,mean = PTRpos$mean))
    } 
  }
  else{ 
    ## assuming consistent knots from the observed points, this is just to make sure 
    ## the knot placement across different knot value will be consistent
    not_selected_x <- x
    selected_x <- NULL
    for (i in 1:length(k_vec)) {
      k <- k_vec[i]
      selected_x <- sort(c(selected_x, sample(not_selected_x, size = (k-length(selected_x)))))
      not_selected_x <- not_selected_x[!not_selected_x %in% selected_x]
      PTRmod <- Imple_BayesRegression(y = y,knots = selected_x, x = x, p = p, likelihood = likelihood)
      PTRpos <- extract_mean_interval_refined(PTRmod, x = selected_x, refined_x = x, p = p)
      result <- rbind(result,data.frame(k = k, x = x ,mean = PTRpos$mean))
    } 
  }
  result %>% arrange(k) %>% mutate(k = as.factor(k))
}
sim_RW2_knots <- function(y,x,k_vec,min = 0, max = 100, placement = "equal", diagonal = 1e-9, likelihood = "Gaussian", aghq_k = 7){
  mod <- Imple_RW2(y = y,x = x, diagonal = diagonal, aghq_k = aghq_k, likelihood = likelihood)
  pos <- extract_mean_interval_rw2(mod, x = x)
  result <- data.frame(k = n, x = x ,mean = pos$mean)
  if(placement == "equal"){
    for (i in 1:length(k_vec)) {
      k <- k_vec[i]
      selected_x <- seq(min,max, length.out = k)
      mod <- Imple_RW2(y = y, knots = selected_x, x = x, likelihood = likelihood)
      pos <- extract_mean_interval_rw2_refined(mod, x = selected_x, refined_x = x)
      result <- rbind(result,data.frame(k = k, x = x ,mean = pos$mean))
    } 
  }
  else if(placement == "quantile"){
    for (i in 1:length(k_vec)) {
      k <- k_vec[i]
      selected_x <- quantile(x, p = seq(0,1,length.out = k))
      mod <- Imple_RW2(y = y, knots = selected_x, x = x, likelihood = likelihood)
      pos <- extract_mean_interval_rw2_refined(mod, x = selected_x, refined_x = x)
      result <- rbind(result,data.frame(k = k, x = x ,mean = pos$mean))
    }
  }
  else{
    not_selected_x <- x
    selected_x <- NULL
    for (i in 1:length(k_vec)) {
      k <- k_vec[i]
      selected_x <- sort(c(selected_x, sample(not_selected_x, size = (k-length(selected_x)))))
      not_selected_x <- not_selected_x[!not_selected_x %in% selected_x]
      mod <- Imple_RW2(y = y, knots = selected_x, x = x, likelihood = likelihood)
      pos <- extract_mean_interval_rw2_refined(mod, x = selected_x, refined_x = x)
      result <- rbind(result,data.frame(k = k, x = x ,mean = pos$mean))
    }
  }
  result %>% arrange(k) %>% mutate(k = as.factor(k))
}


##### Functions for simulation for MSE, MCW with different functions:
sim_Accurary_compare <- function(x, num_replic = 300, sigmaS, prior.used = NULL, k = NULL){
  result_gx <- data.frame()
  result_deriv <- data.frame()
  if(is.null(prior.used) == T){
    prior.used <- list(
      u1 = 1,
      alpha1 = 0.5,
      u2 = 1,
      alpha2 = 0.5,
      betaprec = 10^(-6)
    )
  }
  n <- length(x)
  for (i in 1:num_replic) {
    true_para <- true_ploy(x)
    true_nonpara <- true_deviation(x, folds = 2, weight = sigmaS, max_t = max(x), mesh_size = 0.0001)
    gx <- true_para + true_nonpara
    y <- gx + rnorm(n, sd = 0.2)
    if(is.null(k) == T){
      k <- n
    }
    knots <- seq(min(x),max(x), length.out = k)
    ### Fit Proposed O-Spline:
    O2_mod <- Imple_BayesRegression(y,x,knots = knots, p = 2, prior = prior.used)
    O2_result <- extract_mean_interval_refined(O2_mod, refined_x = x, x = knots, p = 2)
    ### Fit Proposed B-Spline:
    B2_mod <- Imple_BsplineSmooth(y,x, q = (k + 2), p = 2, prior = prior.used)
    B2_result <- extract_mean_interval_Bspline(mod = B2_mod, x = x, p = 2, min = min(x), max = max(x))
    ### Fit RW2 method:
    RW2_mod <- Imple_RW2(y,x, knots = knots, prior = prior.used)
    RW2_result <- extract_mean_interval_rw2_refined(quad = RW2_mod, x = knots, refined_x = x)
    ### Fit WP result:
    WP_mod <- Imple_WP(y = y, x = x, p = 2, prior = prior.used)
    WP_result <- extract_mean_interval_WP(quad = WP_mod, x = x)
    MCW <- c()
    MCW[1] <- mean(O2_result$upper - O2_result$lower)
    MCW[2] <- mean(B2_result$upper - B2_result$lower)
    MCW[3] <- mean(RW2_result$upper - RW2_result$lower)
    MCW[4] <- mean(WP_result$upper - WP_result$lower)
    ### MSE with the true function:
    MSE <- c()
    MSE[1] <- mean((O2_result$mean - gx)^2)
    MSE[2] <- mean((B2_result$mean - gx)^2)
    MSE[3] <- mean((RW2_result$mean - gx)^2)
    MSE[4] <- mean((WP_result$mean - gx)^2)
    ### Display result:
    result_gx <- rbind(result_gx, data.frame(method = c("O2","B2","RW2", "WP"), MSE = MSE, MCW = MCW))
    O2_1st <- extract_deriv_OSpline(O2_mod, refined_x = x, x = knots, degree = 1, p = 2)
    B2_1st <- extract_deriv_BSpline_numeric(B2_mod, x = x, degree = 1, p = 2, min = min(x), max = max(x))
    RW2_1st <- extract_deriv_RW2_numeric(quad = RW2_mod, refined_x = x, x = knots, degree = 1)
    WP_1st <- extract_deriv_WP(WP_mod,x = x, degree = 1, p = 2)
    MCW <- c()
    MCW[1] <- mean(O2_1st$upper - O2_1st$lower)
    MCW[2] <- mean(B2_1st$upper - B2_1st$lower)
    MCW[3] <- mean(RW2_1st$upper - RW2_1st$lower)
    MCW[4] <- mean(WP_1st$upper - WP_1st$lower)
    ### MSE with the true function:
    MSE <- c()
    MSE[1] <- mean((O2_1st$mean[-1] - diff(gx)/mean(diff(x)))^2)
    MSE[2] <- mean((B2_1st$mean - diff(gx)/mean(diff(x)))^2)
    MSE[3] <- mean((RW2_1st$mean - diff(gx)/mean(diff(x)))^2)
    MSE[4] <- mean((WP_1st$mean[-1] - diff(gx)/mean(diff(x)))^2)
    ### Display result:
    result_deriv <- rbind(result_deriv, data.frame(method = c("O2","B2","RW2", "WP"), MSE = MSE, MCW = MCW)) 
  }
  list(result_gx = result_gx, result_deriv = result_deriv)
}

#### Simulation using exact method with known parameters, for replication:
### with sample based GCR from Excursion
sim_WP_known <- function(B = 100, sigmaE, theta, x, p){
  ptw_cov_g <- c()
  cov_g <- c()
  ptw_cov_g1st <- c()
  cov_g1st <- c()
  pb <- progress::progress_bar$new(total = B)
  for (i in 1:B) {
    RandomFunc <- sim_Wp(mesh_size = 0.001, max_t = 10, p = p)
    RandomDeriv <- as.data.frame(apply(RandomFunc, 2, compute_numeric_deriv, h = 0.001))
    RandomDeriv$t <- RandomFunc$t[-1]
    RandomDeriv <- rbind(c(0,0), RandomDeriv)
    
    sigmaS <- sqrt(1/exp(theta)) ## convert theta to sigma
    gx <- sigmaS * evaluate_GP(x = x, GP = RandomFunc)
    gx1st <- sigmaS * evaluate_GP(x = x, GP = RandomDeriv)
    
    sigmaE <- 1
    y <- rnorm(n = length(x), mean = gx, sd = sigmaE)
    exact_result <- Imple_WP_known(beta = rep(0,p), sigmaE = 1, theta = 0, p = p, x = x, y = y, optimizer = T)
    cat(paste("\n", exact_result$opt$message, "\n"))
    
    if (exact_result$opt$convergence != 0) {
      cat("\n Convergence fails at this replication: this dataset is skipped. \n")
    }
    
    else{
      
      ### Convert to samples of functions
      samps_func <- extract_samples_WP_known(result = exact_result, beta = rep(0,p), x = x, target = "function")
      samps_deriv <- extract_samples_WP_known(result = exact_result, beta = rep(0,p), x = x, target = "derivative")
      
      ### Compute interval
      exact_GCR <- extract_GCR(samps_func, level = 0.8, pointwise = T)
      exact_GCR_deriv <- extract_GCR(samps_deriv, level = 0.8, pointwise = T)
      
      ### Compute coverage:
      ptw_cov_g <- c(ptw_cov_g, mean(exact_GCR$pupper >= (gx) & exact_GCR$plower <= (gx)))
      cov_g <- c(cov_g, all(exact_GCR$upper >= (gx) & exact_GCR$lower <= (gx)))
      ptw_cov_g1st <- c(ptw_cov_g1st, mean(exact_GCR_deriv$pupper >= (gx1st) & exact_GCR_deriv$plower <= (gx1st)))
      cov_g1st <- c(cov_g1st, all(exact_GCR_deriv$upper >= (gx1st) & exact_GCR_deriv$lower <= (gx1st)))
      
    }
    pb$tick()
  }
  
  list(ptw_cov_g = ptw_cov_g, cov_g = cov_g, ptw_cov_g1st = ptw_cov_g1st, cov_g1st = cov_g1st)
  
}
#### Simulation using exact method with known parameters, for replication:
### with precision based GCR from Excursion
sim_WP_known_v2 <- function(B = 100, sigmaE, theta, x, p, level = 0.8){
  ptw_cov_g <- c()
  cov_g <- c()
  ptw_cov_g1st <- c()
  cov_g1st <- c()
  pb <- progress::progress_bar$new(total = B)
  for (i in 1:B) {
    RandomFunc <- sim_Wp(mesh_size = 0.001, max_t = 10, p = p)
    RandomDeriv <- as.data.frame(apply(RandomFunc, 2, compute_numeric_deriv, h = 0.001))
    RandomDeriv$t <- RandomFunc$t[-1]
    RandomDeriv <- rbind(c(0,0), RandomDeriv)
    
    sigmaS <- sqrt(1/exp(theta)) ## convert theta to sigma
    gx <- sigmaS * evaluate_GP(x = x, GP = RandomFunc)
    gx1st <- sigmaS * evaluate_GP(x = x, GP = RandomDeriv)
    
    sigmaE <- 1
    y <- rnorm(n = length(x), mean = gx, sd = sigmaE)
    exact_result <- Imple_WP_known(beta = rep(0,p), sigmaE = sigmaE, theta = theta, p = p, x = x, y = y, optimizer = T)
    cat(paste("\n", exact_result$opt$message, "\n"))
    
    if (exact_result$opt$convergence != 0) {
      cat("\n Convergence fails at this replication: this dataset is skipped. \n")
    }
    
    else{
      ### Partition the mean vector and Precision matrix
      mu <- exact_result$opt$par
      Q <- exact_result$prec
      q <- length(x) - 1
      indx_g <- seq(1, (q*p), by = p)
      indx_g1st <- seq(1, (q*p), by = p) + 1
      Qg <- forceSymmetric(solve(solve(Q)[indx_g, indx_g]))
      Qg1st <- forceSymmetric(solve(solve(Q)[indx_g1st, indx_g1st]))
      mu_g <- mu[indx_g]
      mu_g1st <- mu[indx_g1st]
      ### Compute interval
      exact_GCR <- simconf(alpha = (1-level), mu = mu_g, Q = Qg)
      exact_GCR_deriv <- simconf(alpha = (1-level), mu = mu_g1st, Q = Qg1st)
      
      ### Compute coverage:
      ### There is a bug in simconf: a and b mean opposite things for marginal and global interval
      
      if(exact_GCR$b[1] > exact_GCR$a[1]){
        ptw_cov_g <- c(ptw_cov_g, mean(c(0,exact_GCR$a.marginal) >= (gx) & c(0, exact_GCR$b.marginal) <= (gx)))
        cov_g <- c(cov_g, all(c(0, exact_GCR$b) >= (gx) & c(0, exact_GCR$a) <= (gx)))
        ptw_cov_g1st <- c(ptw_cov_g1st, mean(c(0, exact_GCR_deriv$a.marginal) >= (gx1st) & c(0, exact_GCR_deriv$b.marginal) <= (gx1st)))
        cov_g1st <- c(cov_g1st, all(c(0, exact_GCR_deriv$b) >= (gx1st) & c(0, exact_GCR_deriv$a) <= (gx1st)))
      }
      else{
        ptw_cov_g <- c(ptw_cov_g, mean(c(0,exact_GCR$a.marginal) >= (gx) & c(0, exact_GCR$b.marginal) <= (gx)))
        cov_g <- c(cov_g, all(c(0, exact_GCR$a) >= (gx) & c(0, exact_GCR$b) <= (gx)))
        ptw_cov_g1st <- c(ptw_cov_g1st, mean(c(0, exact_GCR_deriv$a.marginal) >= (gx1st) & c(0, exact_GCR_deriv$b.marginal) <= (gx1st)))
        cov_g1st <- c(cov_g1st, all(c(0, exact_GCR_deriv$a) >= (gx1st) & c(0, exact_GCR_deriv$b) <= (gx1st)))
      }

    }
    pb$tick()
  }
  
  list(ptw_cov_g = ptw_cov_g, cov_g = cov_g, ptw_cov_g1st = ptw_cov_g1st, cov_g1st = cov_g1st)
  
}


#### Simulation using exact method with known parameters, for replication:
### with precision based GCR from Excursion
### Generate from IWP using VAR method
sim_WP_known_v3 <- function(B = 100, sigmaE, theta, x, p, level = 0.8){
  ptw_cov_g <- c()
  cov_g <- c()
  ptw_cov_g1st <- c()
  cov_g1st <- c()
  pb <- progress::progress_bar$new(total = B)
  for (i in 1:B) {
    sigmaS <- sqrt(1/exp(theta)) ## convert theta to sigma
    
    RandomFunc <- sim_Wp_Var(mesh_size = 0.001, max_t = 10, p = p)
    RandomFunc[,-1] <- sigmaS * RandomFunc[,-1]
    
    gx_data <- evaluate_GP_Var(x, RandomFunc)[,1:2]
    gx1st_data <- evaluate_GP_Var(x, RandomFunc)[,c(1,3)]
    
    gx <- gx_data[,2]
    gx1st <- gx1st_data[,2]
    
    
    sigmaE <- 1
    y <- rnorm(n = length(x), mean = gx, sd = sigmaE)
    exact_result <- Imple_WP_known(beta = rep(0,p), sigmaE = sigmaE, theta = theta, p = p, x = x, y = y, optimizer = T)
    cat(paste("\n", exact_result$opt$message, "\n"))
    
    if (exact_result$opt$convergence != 0) {
      cat("\n Convergence fails at this replication: this dataset is skipped. \n")
    }
    
    else{
      ### Partition the mean vector and Precision matrix
      mu <- exact_result$opt$par
      Q <- exact_result$prec
      q <- length(x) - 1
      indx_g <- seq(1, (q*p), by = p)
      indx_g1st <- seq(1, (q*p), by = p) + 1
      Qg <- forceSymmetric(solve(solve(Q)[indx_g, indx_g]))
      Qg1st <- forceSymmetric(solve(solve(Q)[indx_g1st, indx_g1st]))
      mu_g <- mu[indx_g]
      mu_g1st <- mu[indx_g1st]
      ### Compute interval
      exact_GCR <- simconf(alpha = (1-level), mu = mu_g, Q = Qg)
      exact_GCR_deriv <- simconf(alpha = (1-level), mu = mu_g1st, Q = Qg1st)
      
      ### Compute coverage:
      ### There is a bug in simconf: a and b mean opposite things for marginal and global interval
      
      if(exact_GCR$b[1] > exact_GCR$a[1]){
        ptw_cov_g <- c(ptw_cov_g, mean(c(0,exact_GCR$a.marginal) >= (gx) & c(0, exact_GCR$b.marginal) <= (gx)))
        cov_g <- c(cov_g, all(c(0, exact_GCR$b) >= (gx) & c(0, exact_GCR$a) <= (gx)))
        ptw_cov_g1st <- c(ptw_cov_g1st, mean(c(0, exact_GCR_deriv$a.marginal) >= (gx1st) & c(0, exact_GCR_deriv$b.marginal) <= (gx1st)))
        cov_g1st <- c(cov_g1st, all(c(0, exact_GCR_deriv$b) >= (gx1st) & c(0, exact_GCR_deriv$a) <= (gx1st)))
      }
      else{
        ptw_cov_g <- c(ptw_cov_g, mean(c(0,exact_GCR$a.marginal) >= (gx) & c(0, exact_GCR$b.marginal) <= (gx)))
        cov_g <- c(cov_g, all(c(0, exact_GCR$a) >= (gx) & c(0, exact_GCR$b) <= (gx)))
        ptw_cov_g1st <- c(ptw_cov_g1st, mean(c(0, exact_GCR_deriv$a.marginal) >= (gx1st) & c(0, exact_GCR_deriv$b.marginal) <= (gx1st)))
        cov_g1st <- c(cov_g1st, all(c(0, exact_GCR_deriv$a) >= (gx1st) & c(0, exact_GCR_deriv$b) <= (gx1st)))
      }
      
    }
    pb$tick()
  }
  
  list(ptw_cov_g = ptw_cov_g, cov_g = cov_g, ptw_cov_g1st = ptw_cov_g1st, cov_g1st = cov_g1st)
  
}





#### Simulation using OS method with known parameters, for replication:
### with sample based GCR from Excursion or GET
sim_OS_known <- function(B = 100, sigmaE, theta, x, p, k, method = "excursion", GET = "erl", level = 0.8){
  ptw_cov_g <- c()
  cov_g <- c()
  ptw_cov_g1st <- c()
  cov_g1st <- c()
  pb <- progress::progress_bar$new(total = B)
  for (i in 1:B) {
    RandomFunc <- sim_Wp(mesh_size = 0.001, max_t = 10, p = p)
    RandomDeriv <- as.data.frame(apply(RandomFunc, 2, compute_numeric_deriv, h = 0.001))
    RandomDeriv$t <- RandomFunc$t[-1]
    RandomDeriv <- rbind(c(0,0), RandomDeriv)
    
    sigmaS <- sqrt(1/exp(theta)) ## convert theta to sigma
    gx <- sigmaS * evaluate_GP(x = x, GP = RandomFunc)
    gx1st <- sigmaS * evaluate_GP(x = x, GP = RandomDeriv)
    
    sigmaE <- 1
    y <- rnorm(n = length(x), mean = gx, sd = sigmaE)
    OSmod <- Imple_OS_known(beta = rep(0,p), sigmaE = sigmaE, theta = theta, x = x, k = k, y = y, p = p, optimizer = T)
    
    cat(paste("\n", OSmod$opt$message, "\n"))
    
    if (OSmod$opt$convergence != 0) {
      cat("\n Convergence fails at this replication: this dataset is skipped. \n")
    }
    
    else{
      GCR_g <- extract_GCR(OSmod$samples_fun, level = level, pointwise = T, method = method, GET = GET)
      GCR_g1st <- extract_GCR(OSmod$sample_deriv, level = level, pointwise = T, method = method, GET = GET)
      
      ### Compute coverage:
      ptw_cov_g <- c(ptw_cov_g, mean(GCR_g$pupper >= (gx) & GCR_g$plower <= (gx)))
      cov_g <- c(cov_g, all(GCR_g$upper >= (gx) & GCR_g$lower <= (gx)))
      ptw_cov_g1st <- c(ptw_cov_g1st, mean(GCR_g1st$pupper >= (gx1st) & GCR_g1st$plower <= (gx1st)))
      cov_g1st <- c(cov_g1st, all(GCR_g1st$upper >= (gx1st) & GCR_g1st$lower <= (gx1st)))

    }
    pb$tick()
  }
  
  list(ptw_cov_g = ptw_cov_g, cov_g = cov_g, ptw_cov_g1st = ptw_cov_g1st, cov_g1st = cov_g1st)
  
}


#### Simulation using RW2 method with known parameters, for replication:
### with sample based GCR from Excursion or GET
sim_RW2_known <- function(B = 100, sigmaE, theta, x, p, k, method = "excursion", GET = "erl", level = 0.8, diagonal = 1e-9){
  ptw_cov_g <- c()
  cov_g <- c()
  ptw_cov_g1st <- c()
  cov_g1st <- c()
  pb <- progress::progress_bar$new(total = B)
  for (i in 1:B) {
    RandomFunc <- sim_Wp(mesh_size = 0.001, max_t = 10, p = p)
    RandomDeriv <- as.data.frame(apply(RandomFunc, 2, compute_numeric_deriv, h = 0.001))
    RandomDeriv$t <- RandomFunc$t[-1]
    RandomDeriv <- rbind(c(0,0), RandomDeriv)
    
    sigmaS <- sqrt(1/exp(theta)) ## convert theta to sigma
    gx <- sigmaS * evaluate_GP(x = x, GP = RandomFunc)
    gx1st <- sigmaS * evaluate_GP(x = x, GP = RandomDeriv)
    
    sigmaE <- 1
    y <- rnorm(n = length(x), mean = gx, sd = sigmaE)
    RW2mod <- Imple_RW2_known(beta = rep(0,p), sigmaE = sigmaE, theta = theta, x = x, k = k, y = y, p = p, optimizer = T, diagonal = diagonal)
    
    cat(paste("\n", RW2mod$opt$message, "\n"))
    
    if (RW2mod$opt$convergence != 0) {
      cat("\n Convergence fails at this replication: this dataset is skipped. \n")
    }
    
    else{
      GCR_g <- extract_GCR(RW2mod$samples_fun, level = level, pointwise = T, method = method, GET = GET)
      GCR_g1st <- extract_GCR(RW2mod$sample_deriv, level = level, pointwise = T, method = method, GET = GET)
      
      ### Compute coverage:
      ptw_cov_g <- c(ptw_cov_g, mean(GCR_g$pupper >= (gx) & GCR_g$plower <= (gx)))
      cov_g <- c(cov_g, all(GCR_g$upper >= (gx) & GCR_g$lower <= (gx)))
      ptw_cov_g1st <- c(ptw_cov_g1st, mean(GCR_g1st$pupper >= (gx1st) & GCR_g1st$plower <= (gx1st)))
      cov_g1st <- c(cov_g1st, all(GCR_g1st$upper >= (gx1st) & GCR_g1st$lower <= (gx1st)))
      
    }
    pb$tick()
  }
  
  list(ptw_cov_g = ptw_cov_g, cov_g = cov_g, ptw_cov_g1st = ptw_cov_g1st, cov_g1st = cov_g1st)
  
}



## samps: output from extract_samps
## intervals: output from extract_GCR
compute_COV_samps <- function(samps, intervals, indx){
  cov <- apply(samps[,-1], 2, function(x)  all(x[indx]<= intervals$upper[indx]) &  all(x [indx] >= intervals$lower[indx]))
  mean(cov)
}

## samps: output from extract_samps
## intervals: output from extract_GCR with pointwise = T
compute_pCOV_samps <- function(samps, intervals, indx){
  cov <- apply(samps[,-1], 2, function(x)  mean(x[indx] <= intervals$pupper[indx] & x[indx] >= intervals$plower[indx]))
  mean(cov)
}


#### In this version, implement the new idea of coverage of posterior sample
#### pathes from the exact method:

### Input:
# x: The location vector on which the data will be simulated
# x_compare: The location vector on which to evaluate the comparison,
#            should be choose to be far away from 0 to avoid boundary effect.
#            This have to be a subset of x vector.
# p: What is the smoothing degree to compare
# kvec: How many knots will be used in splines-based approximation (can be a vector)
# sigmaE: The residual SD in the simulation
# theta: The log precision of IWP to consider
# B: The number of replications
# nsamps: The number of samples to be drawn from each method

### Output: A list of two elements: overall_result_g and overall_result_g1st
### each element is a dataframe, with variable: Method, k, COV

compare_all_three_methods_B_times_sample_cov <- function(x, x_compare, p, kvec, sigmaE, theta, seed = 12345,
                                                         control = list(method = "GET", GET = "rank", level = 0.8), 
                                                         nsamps = 8000, parallel = FALSE){
  method = control$method
  if(is.null(kvec)){
    kvec <- c(length(x))
  }
  GET = control$GET
  level = control$level
  indx <- which(x %in% x_compare)
  sigmaS <- sqrt(1/exp(theta))
  g_result <- data.frame()
  g1st_result <- data.frame()
  
  if(parallel == T){
    set.seed(seed)
    RandomFunc <- sim_Wp_Var(mesh_size = 0.001, max_t = 10, p = p)
    RandomFunc[,-1] <- sigmaS * RandomFunc[,-1]
    
    # The function
    RandomFunc_IWPp <- RandomFunc[,1:2]
    RandomFunc_IWPp$IWp <- RandomFunc_IWPp$IWp
    # The first derivative, so q = p - 1
    RandomFunc_IWPq <- RandomFunc[,c(1,3)]
    RandomFunc_IWPq$IWq <- RandomFunc_IWPq$IWq
    
    n <- length(x)
    gx <- evaluate_GP_Var(x, RandomFunc)[,1:2]
    names(gx) <- c("t", "IWp")
    gx1st <- evaluate_GP_Var(x, RandomFunc)[,c(1,3)]
    names(gx1st) <- c("t", "IWq")
    
    y <- gx$IWp + rnorm(n = n, mean = 0, sd = sigmaE)
    
    
    ## Fit the exact method
    exact_result <- Imple_WP_known(beta = rep(0, p), sigmaE = sigmaE, theta = theta, p = p, x = x, y = y, optimizer = T)
    
    ### Convert to samples of functions
    samps_func <- extract_samples_WP_known(result = exact_result, beta = rep(0, p), x = x, target = "function", p = p, nsam = nsamps)
    samps_deriv <- extract_samples_WP_known(result = exact_result, beta = rep(0, p), x = x, target = "derivative", p = p, nsam = nsamps)
    exact_g_mean <- samps_func[,c(-1)] %>% apply(1, mean)
    exact_g1st_mean <- samps_deriv[,c(-1)] %>% apply(1, mean)
    
    
    do_once <- function(seed,k){
      set.seed(seed)
      OSmod <- Imple_OS_known(beta = rep(0,p), sigmaE = sigmaE, theta = theta , x = x, k = (k+1), y = y, p = p, optimizer = T, nsamps = nsamps)
      OS_GCR_g <- extract_GCR(OSmod$samples_fun, level = level, pointwise = T, method = method, GET = GET)
      OS_GCR_g1st <- extract_GCR(OSmod$sample_deriv, level = level, pointwise = T, method = method, GET = GET)
      
      diagonal <- 1e-10
      RW2mod <- Imple_RW2_known(beta = rep(0,2), sigmaE = sigmaE, theta = theta, x = x, k = k, y = y, p = 2, optimizer = T, diagonal = diagonal, nsamps = nsamps)
      names(RW2mod$samples_fun)[1] <- "x"
      RW2_GCR_g <- extract_GCR(RW2mod$samples_fun, level = level, pointwise = T, method = method, GET = GET)
      RW2_GCR_g1st <- extract_GCR(RW2mod$sample_deriv, level = level, pointwise = T, method = method, GET = GET)
      
      
      data <- data.frame(x = x, y = y)
      newdata <- data.frame(x = x)
      
      ##### Fit P-Spline:
      ### cubic basis and second order penalty
      lambda <- (sigmaE^2)/((sigmaS^2))
      control_list <- list(scalePenalty = FALSE)
      BSmod <- gam(formula = y~s(x,bs="bs",m=c(3,p), k = k, pc = c(0), sp = lambda), data = data, family = gaussian(), scale = (sigmaE^2), drop.intercept = T, control = control_list)
      
      ### Samples from the fitted GAM:
      gam_coef_samples <- t(rmvn(nsamps,coef(BSmod),vcov(BSmod)))
      gam_design <- predict.gam(BSmod, newdata = newdata, type = 'lpmatrix')
      gam_samples <- gam_design %*% gam_coef_samples
      gam_samples <- cbind(x, gam_samples)
      gam_samples <- as.data.frame(gam_samples)
      gam_conf_simu <- extract_GCR(gam_samples, level = level, pointwise = T, method = method, GET = GET)
      
      gam_conf_simu$upper[1] <- 0
      gam_conf_simu$lower[1] <- 0
      gam_conf_simu$pupper[1] <- 0
      gam_conf_simu$plower[1] <- 0
      
      gam_newdata <- data.frame(x = x)
      
      gam_conf <- confint(BSmod, parm = "s(x)", type = "confidence", shift = F, level = level, newdata = newdata)
      
      gam_deriv <- derivatives(
        BSmod,
        term = "s(x)",
        newdata = gam_newdata,
        order = 1L,
        type = c("backward"),
        eps = 1e-10,
        interval = c("simultaneous"),
        n_sim = nsamps,
        level = 0.8,
        unconditional = F,
        ncores = 1
      )
      gam_deriv_pt <- derivatives(
        BSmod,
        term = "s(x)",
        newdata = gam_newdata,
        order = 1L,
        type = c("backward"),
        eps = 1e-10,
        interval = c("confidence"),
        n_sim = nsamps,
        level = 0.8,
        unconditional = F,
        ncores = 1
      )
      
      BS_GCR_g1st <- data.frame(x = x, upper = gam_deriv$upper, lower = gam_deriv$lower, pupper = gam_deriv_pt$upper, plower = gam_deriv_pt$lower)
      
      
      OS_g_mean <- apply(OSmod$samples_fun[,-1], 1, mean)
      RW2_g_mean <- apply(RW2mod$samples_fun[,-1], 1, mean)
      OS_g1st_mean <- apply(OSmod$sample_deriv[,-1], 1, mean)
      RW2_g1st_mean <- apply(RW2mod$sample_deriv[,-1], 1, mean)
      
      g_RMSE_OS_compared <- sqrt(mean(((OS_g_mean-exact_g_mean)^2)[indx])) # by definition
      g1st_RMSE_OS_compared <- sqrt(mean(((OS_g1st_mean-exact_g1st_mean)^2)[indx])) # by definition
      
      g_RMSE_RW_compared <- sqrt(mean(((RW2_g_mean-exact_g_mean)^2)[indx])) # by definition
      g1st_RMSE_RW_compared <- sqrt(mean(((RW2_g1st_mean-exact_g1st_mean)^2)[indx])) # by definition
      
      g_RMSE_BS_compared <- sqrt(mean(((gam_conf$est-exact_g_mean)^2)[indx])) # by definition
      g1st_RMSE_BS_compared <- sqrt(mean(((gam_deriv$derivative-exact_g1st_mean)^2)[indx])) # by definition
      
      
      g_OS_cov <- compute_COV_samps(samps = samps_func, intervals = OS_GCR_g, indx = indx)
      g1st_OS_cov <- compute_COV_samps(samps = samps_deriv, intervals = OS_GCR_g1st, indx = indx)
      g_OS_pcov <- compute_pCOV_samps(samps = samps_func, intervals = OS_GCR_g, indx = indx)
      g1st_OS_pcov <- compute_pCOV_samps(samps = samps_deriv, intervals = OS_GCR_g1st, indx = indx)
      
      g_result_OS <- data.frame(Method = "OS", k = k, COV = g_OS_cov, PCOV = g_OS_pcov, Type = "g", rMSE = g_RMSE_OS_compared)
      g1st_result_OS <- data.frame(Method = "OS", k = k, COV = g1st_OS_cov, PCOV = g1st_OS_pcov, Type = "g1st", rMSE = g1st_RMSE_OS_compared)
      result_OS <- rbind(g_result_OS, g1st_result_OS)
      
      g_RW2_cov <- compute_COV_samps(samps = samps_func, intervals = RW2_GCR_g, indx = indx)
      g1st_RW2_cov <- compute_COV_samps(samps = samps_deriv, intervals = RW2_GCR_g1st, indx = indx)
      g_RW2_pcov <- compute_pCOV_samps(samps = samps_func, intervals = RW2_GCR_g, indx = indx)
      g1st_RW2_pcov <- compute_pCOV_samps(samps = samps_deriv, intervals = RW2_GCR_g1st, indx = indx)
      
      g_result_RW <- data.frame(Method = "RW", k = k, COV = g_RW2_cov, PCOV = g_RW2_pcov, Type = "g", rMSE = g_RMSE_RW_compared)
      g1st_result_RW <- data.frame(Method = "RW", k = k, COV = g1st_RW2_cov, PCOV = g1st_RW2_pcov, Type = "g1st", rMSE = g1st_RMSE_RW_compared)
      result_RW <- rbind(g_result_RW, g1st_result_RW)
      
      g_BS_cov <- compute_COV_samps(samps = samps_func, intervals = gam_conf_simu, indx = indx)
      g1st_BS_cov <- compute_COV_samps(samps = samps_deriv, intervals = BS_GCR_g1st, indx = indx)
      g_BS_pcov <- compute_pCOV_samps(samps = samps_func, intervals = gam_conf_simu, indx = indx)
      g1st_BS_pcov <- compute_pCOV_samps(samps = samps_deriv, intervals = BS_GCR_g1st, indx = indx)
      
      g_result_BS <- data.frame(Method = "BS", k = k, COV = g_BS_cov, PCOV = g_BS_pcov, Type = "g", rMSE = g_RMSE_BS_compared)
      g1st_result_BS <- data.frame(Method = "BS", k = k, COV = g1st_BS_cov, PCOV = g1st_BS_pcov, Type = "g1st", rMSE = g1st_RMSE_BS_compared)
      result_BS <- rbind(g_result_BS, g1st_result_BS)
      
      overall_result <- rbind(result_OS, result_RW, result_BS)
      return(overall_result)
    }
    result_all <- foreach(i = 1:length(kvec), .combine = "rbind", .packages = c("mgcv", "gratia", "TMB", "aghq", "Matrix", "tidyverse", "GET", "tsDyn", "LaplacesDemon")) %dopar% {
      source("14_functions_defined.R") ## Only on Windows
      dyn.load(dynlib("00_Gaussian_Smoothing_known")) ## Only on Windows
      do_once(seed = i, k = kvec[i])
    }
    g_result <- result_all %>% filter(Type == "g") %>% select(-Type)
    g1st_result <- result_all %>% filter(Type == "g1st") %>% select(-Type)
  }
  else{
    set.seed(seed)
    pb <- progress_bar$new(total = length(kvec))
    overall_result <- data.frame()
    RandomFunc <- sim_Wp_Var(mesh_size = 0.001, max_t = 10, p = p)
    RandomFunc[,-1] <- sigmaS * RandomFunc[,-1]
    RandomFunc_IWPp <- RandomFunc[,1:2]
    RandomFunc_IWPp$IWp <- RandomFunc_IWPp$IWp
    RandomFunc_IWPq <- RandomFunc[,c(1,3)]
    RandomFunc_IWPq$IWq <- RandomFunc_IWPq$IWq
    n <- length(x)
    gx <- evaluate_GP_Var(x, RandomFunc)[,1:2]
    names(gx) <- c("t", "IWp")
    gx1st <- evaluate_GP_Var(x, RandomFunc)[,c(1,3)]
    names(gx1st) <- c("t", "IWq")
    y <- gx$IWp + rnorm(n = n, mean = 0, sd = sigmaE)
    
    ## Fit the exact method
    exact_result <- Imple_WP_known(beta = rep(0, p), sigmaE = sigmaE, theta = theta, p = p, x = x, y = y, optimizer = T)
    
    ### Convert to samples of functions
    samps_func <- extract_samples_WP_known(result = exact_result, beta = rep(0, p), x = x, target = "function", p = p, nsam = nsamps)
    samps_deriv <- extract_samples_WP_known(result = exact_result, beta = rep(0, p), x = x, target = "derivative", p = p, nsam = nsamps)
    exact_g_mean <- samps_func[,c(-1)] %>% apply(1, mean)
    exact_g1st_mean <- samps_deriv[,c(-1)] %>% apply(1, mean)
    
    overall_result <- data.frame()
    for (i in 1:length(kvec)) {
      k <- kvec[i]
      pb$tick()
      set.seed(i)
      
      OSmod <- Imple_OS_known(beta = rep(0,p), sigmaE = sigmaE, theta = theta , x = x, k = (k+1), y = y, p = p, optimizer = T)
      OS_GCR_g <- extract_GCR(OSmod$samples_fun, level = level, pointwise = T, method = method, GET = GET)
      OS_GCR_g1st <- extract_GCR(OSmod$sample_deriv, level = level, pointwise = T, method = method, GET = GET)
      
      diagonal <- 1e-10
      RW2mod <- Imple_RW2_known(beta = rep(0,2), sigmaE = sigmaE, theta = theta, x = x, k = k, y = y, p = 2, optimizer = T, diagonal = diagonal)
      names(RW2mod$samples_fun)[1] <- "x"
      RW2_GCR_g <- extract_GCR(RW2mod$samples_fun, level = level, pointwise = T, method = method, GET = GET)
      RW2_GCR_g1st <- extract_GCR(RW2mod$sample_deriv, level = level, pointwise = T, method = method, GET = GET)
      
      
      data <- data.frame(x = x, y = y)
      newdata <- data.frame(x = x)
      ##### Fit P-Spline:
      ### cubic basis and second order penalty
      lambda <- (sigmaE^2)/((sigmaS^2))
      control_list <- list(scalePenalty = FALSE)
      BSmod <- gam(formula = y~s(x,bs="bs",m=c(3,p), k = (k), pc = c(0), sp = lambda), data = data, family = gaussian(), scale = (sigmaE^2), drop.intercept = T, control = control_list)
      
      ### Samples from the fitted GAM:
      gam_coef_samples <- t(rmvn(nsamps,coef(BSmod),vcov(BSmod)))
      gam_design <- predict.gam(BSmod, newdata = newdata, type = 'lpmatrix')
      gam_samples <- gam_design %*% gam_coef_samples
      gam_samples <- cbind(x, gam_samples)
      gam_samples <- as.data.frame(gam_samples)
      gam_conf_simu <- extract_GCR(gam_samples, level = level, pointwise = T, method = method, GET = GET)
      
      gam_conf_simu$upper[1] <- 0
      gam_conf_simu$lower[1] <- 0
      gam_conf_simu$pupper[1] <- 0
      gam_conf_simu$plower[1] <- 0
      
      gam_newdata <- data.frame(x = x)
      
      gam_conf <- confint(BSmod, parm = "s(x)", type = "confidence", shift = F, level = level, newdata = newdata)
      
      gam_deriv <- derivatives(
        BSmod,
        term = "s(x)",
        newdata = gam_newdata,
        order = 1L,
        type = c("backward"),
        eps = 1e-10,
        interval = c("simultaneous"),
        n_sim = nsamps,
        level = 0.8,
        unconditional = F,
        ncores = 1
      )
      gam_deriv_pt <- derivatives(
        BSmod,
        term = "s(x)",
        newdata = gam_newdata,
        order = 1L,
        type = c("backward"),
        eps = 1e-10,
        interval = c("confidence"),
        n_sim = nsamps,
        level = 0.8,
        unconditional = F,
        ncores = 1
      )
      
      BS_GCR_g1st <- data.frame(x = x, upper = gam_deriv$upper, lower = gam_deriv$lower, pupper = gam_deriv_pt$upper, plower = gam_deriv_pt$lower)
      
      
      OS_g_mean <- apply(OSmod$samples_fun[,-1], 1, mean)
      RW2_g_mean <- apply(RW2mod$samples_fun[,-1], 1, mean)
      OS_g1st_mean <- apply(OSmod$sample_deriv[,-1], 1, mean)
      RW2_g1st_mean <- apply(RW2mod$sample_deriv[,-1], 1, mean)
      
      g_RMSE_OS_compared <- sqrt(mean(((OS_g_mean-exact_g_mean)^2)[indx])) # by definition
      g1st_RMSE_OS_compared <- sqrt(mean(((OS_g1st_mean-exact_g1st_mean)^2)[indx])) # by definition
      
      g_RMSE_RW_compared <- sqrt(mean(((RW2_g_mean-exact_g_mean)^2)[indx])) # by definition
      g1st_RMSE_RW_compared <- sqrt(mean(((RW2_g1st_mean-exact_g1st_mean)^2)[indx])) # by definition
      
      g_RMSE_BS_compared <- sqrt(mean(((gam_conf$est-exact_g_mean)^2)[indx])) # by definition
      g1st_RMSE_BS_compared <- sqrt(mean(((gam_deriv$derivative-exact_g1st_mean)^2)[indx])) # by definition
      
      
      g_OS_cov <- compute_COV_samps(samps = samps_func, intervals = OS_GCR_g, indx = indx)
      g1st_OS_cov <- compute_COV_samps(samps = samps_deriv, intervals = OS_GCR_g1st, indx = indx)
      g_OS_pcov <- compute_pCOV_samps(samps = samps_func, intervals = OS_GCR_g, indx = indx)
      g1st_OS_pcov <- compute_pCOV_samps(samps = samps_deriv, intervals = OS_GCR_g1st, indx = indx)
      
      g_result_OS <- data.frame(Method = "OS", k = k, COV = g_OS_cov, PCOV = g_OS_pcov, Type = "g", rMSE = g_RMSE_OS_compared)
      g1st_result_OS <- data.frame(Method = "OS", k = k, COV = g1st_OS_cov, PCOV = g1st_OS_pcov, Type = "g1st", rMSE = g1st_RMSE_OS_compared)
      result_OS <- rbind(g_result_OS, g1st_result_OS)
      
      g_RW2_cov <- compute_COV_samps(samps = samps_func, intervals = RW2_GCR_g, indx = indx)
      g1st_RW2_cov <- compute_COV_samps(samps = samps_deriv, intervals = RW2_GCR_g1st, indx = indx)
      g_RW2_pcov <- compute_pCOV_samps(samps = samps_func, intervals = RW2_GCR_g, indx = indx)
      g1st_RW2_pcov <- compute_pCOV_samps(samps = samps_deriv, intervals = RW2_GCR_g1st, indx = indx)
      
      g_result_RW <- data.frame(Method = "RW", k = k, COV = g_RW2_cov, PCOV = g_RW2_pcov, Type = "g", rMSE = g_RMSE_RW_compared)
      g1st_result_RW <- data.frame(Method = "RW", k = k, COV = g1st_RW2_cov, PCOV = g1st_RW2_pcov, Type = "g1st", rMSE = g1st_RMSE_RW_compared)
      result_RW <- rbind(g_result_RW, g1st_result_RW)
      
      g_BS_cov <- compute_COV_samps(samps = samps_func, intervals = gam_conf_simu, indx = indx)
      g1st_BS_cov <- compute_COV_samps(samps = samps_deriv, intervals = BS_GCR_g1st, indx = indx)
      g_BS_pcov <- compute_pCOV_samps(samps = samps_func, intervals = gam_conf_simu, indx = indx)
      g1st_BS_pcov <- compute_pCOV_samps(samps = samps_deriv, intervals = BS_GCR_g1st, indx = indx)
      
      g_result_BS <- data.frame(Method = "BS", k = k, COV = g_BS_cov, PCOV = g_BS_pcov, Type = "g", rMSE = g_RMSE_BS_compared)
      g1st_result_BS <- data.frame(Method = "BS", k = k, COV = g1st_BS_cov, PCOV = g1st_BS_pcov, Type = "g1st", rMSE = g1st_RMSE_BS_compared)
      result_BS <- rbind(g_result_BS, g1st_result_BS)
      
      overall_result <- rbind(overall_result, result_OS, result_RW, result_BS)
      
    }
    g_result <- overall_result %>% filter(Type == "g") %>% select(-Type)
    g1st_result <- overall_result %>% filter(Type == "g1st") %>% select(-Type)
  }
  return(list(g = g_result, g1st = g1st_result))
}



#### Summarize the result from the simulation


#### In this version, implement the new idea of coverage of posterior sample
#### pathes from the exact method (with diffuse boundary): no longer need to compute
#### the global interval for computational efficiency reason.

#### Also compute the point-wise coverage rate at the maximum of g and g'

### Input:
# x: The location vector on which the data will be simulated
# x_compare: The location vector on which to evaluate the comparison,
#            should be choose to be far away from 0 to avoid boundary effect.
#            This have to be a subset of x vector.
# p: What is the smoothing degree to compare
# kvec: How many knots will be used in splines-based approximation (can be a vector)
# sigmaE: The residual SD in the simulation
# theta: The log precision of IWP to consider
# B: The number of replications
# nsamps: The number of samples to be drawn from each method

### Output: A list of two elements: overall_result_g and overall_result_g1st
### each element is a dataframe, with variable: Method, k, COV
compare_all_three_methods_B_times_sample_cov_diffuse_v3 <- function(x, x_compare = NULL, p, kvec, sigmaE, theta, seed = 12345,
                                                                    control = list(method = "GET", GET = "rank", level = 0.8), 
                                                                    nsamps = 8000, parallel = FALSE){
  if(is.null(x_compare)){
    x_compare <- x
  }
  method = control$method
  if(is.null(kvec)){
    kvec <- c(length(x))
  }
  GET = control$GET
  level = control$level
  indx <- which(x %in% x_compare)
  sigmaS <- sqrt(1/exp(theta))
  g_result <- data.frame()
  g1st_result <- data.frame()
  set.seed(seed)
  if(parallel == T){
    n <- length(x)
    PD_condition <- FALSE
    while(!PD_condition){
      RandomFunc <- sim_Wp_Var(mesh_size = 0.001, max_t = 10, p = p)
      RandomFunc[,-1] <- sigmaS * RandomFunc[,-1]
      
      # The function
      RandomFunc_IWPp <- RandomFunc[,1:2]
      # The first derivative, so q = p - 1
      RandomFunc_IWPq <- RandomFunc[,c(1,3)]
      
      # Plus the polynomial parts:
      beta <- rnorm(p, mean = 0, sd = 1/(factorial((1:p)-1)))
      poly_parts <- matrix(x, nrow = length(x), ncol = p)
      for (i in 1:p) {
        poly_parts[,i] <- poly_parts[,i]^(i-1)
      }
      poly_parts_gx <- poly_parts %*% beta
      poly_parts_1st <- matrix(poly_parts[,-1],ncol = (p-1))
      for (i in 2:p) {
        poly_parts_1st[,(i-1)] <- (i-1)*poly_parts[,(i-1)]
      }
      poly_parts_gx1st <- poly_parts_1st %*% beta[-1]
      gx <- evaluate_GP_Var(x, RandomFunc)[,1:2]
      names(gx) <- c("t", "IWp")
      gx$IWp <- gx$IWp + poly_parts_gx
      gx1st <- evaluate_GP_Var(x, RandomFunc)[,c(1,3)]
      names(gx1st) <- c("t", "IWq")
      gx1st$IWq <- gx1st$IWq + poly_parts_gx1st
      
      y <- gx$IWp + rnorm(n = n, mean = 0, sd = sigmaE)
      ## Fit the exact method
      exact_result <- Imple_WP_diff_bound(sigmaE = sigmaE, theta = theta, p = p, x = x, y = y, optimizer = T)
      PD_condition <- is.positive.definite(exact_result$prec)
    }
    
    ### Convert to samples of functions
    samps_func <- extract_samples_WP_rand_bound(result = exact_result, x = x, target = "function", p = p, nsam = nsamps)
    samps_deriv <- extract_samples_WP_rand_bound(result = exact_result, x = x, target = "derivative", p = p, nsam = nsamps)
    exact_g_mean <- samps_func[,c(-1)] %>% apply(1, mean)
    exact_g1st_mean <- samps_deriv[,c(-1)] %>% apply(1, mean)
    
    
    do_once <- function(k){
      OSmod <- Imple_OS_diffuse_bound(sigmaE = sigmaE, theta = theta , x = x, k = (k+1), y = y, p = p, optimizer = T, nsamps = nsamps)
      OS_GCR_g <- extract_mean_interval_given_samps(OSmod$samples_fun, level = level)
      OS_GCR_g1st <- extract_mean_interval_given_samps(OSmod$sample_deriv, level = level)
      
      diagonal <- 0
      RW2mod <- Imple_RW_diff_bound(sigmaE = sigmaE, theta = theta, x = x, k = k, y = y, p = 2, optimizer = T, diagonal = diagonal, nsamps = nsamps)
      names(RW2mod$samples_fun)[1] <- "x"
      RW2_GCR_g <- extract_mean_interval_given_samps(RW2mod$samples_fun, level = level)
      RW2_GCR_g1st <- extract_mean_interval_given_samps(RW2mod$sample_deriv, level = level)
      
      ##### Fit P-Spline:
      ### cubic basis and second order penalty
      BSmod <- Imple_MGCV_diffuse_bound(sigmaE = sigmaE, theta = theta, x = x, k = k, y = y, p = p)
      
      ### Samples from the fitted GAM:
      BS_g_samps <- extract_mgcv_samps(BSmod = BSmod, x = x, nsamps = nsamps)
      BS_g1st_samps <- extract_mgcv_deriv_samps(BSmod = BSmod, x = x, nsamps = nsamps)
      BS_g1st_samps <- cbind(x,BS_g1st_samps)
      
      
      BS_GCR_g <- extract_mean_interval_given_samps(BS_g_samps, level = level)
      BS_GCR_g1st <- extract_mean_interval_given_samps(BS_g1st_samps, level = level)
      
      
      OS_g_mean <- apply(OSmod$samples_fun[,-1], 1, mean)
      RW2_g_mean <- apply(RW2mod$samples_fun[,-1], 1, mean)
      OS_g1st_mean <- apply(OSmod$sample_deriv[,-1], 1, mean)
      RW2_g1st_mean <- apply(RW2mod$sample_deriv[,-1], 1, mean)
      
      g_RMSE_OS_compared <- sqrt(mean(((OS_g_mean-exact_g_mean)^2)[indx])) # by definition
      g1st_RMSE_OS_compared <- sqrt(mean(((OS_g1st_mean-exact_g1st_mean)^2)[indx])) # by definition
      
      g_RMSE_RW_compared <- sqrt(mean(((RW2_g_mean-exact_g_mean)^2)[indx])) # by definition
      g1st_RMSE_RW_compared <- sqrt(mean(((RW2_g1st_mean-exact_g1st_mean)^2)[indx])) # by definition
      
      g_RMSE_BS_compared <- sqrt(mean(((BS_GCR_g$mean-exact_g_mean)^2)[indx])) # by definition
      g1st_RMSE_BS_compared <- sqrt(mean(((BS_GCR_g1st$mean-exact_g1st_mean)^2)[indx])) # by definition
      
      g_OS_pcov <- compute_pCOV_samps(samps = samps_func, intervals = OS_GCR_g, indx = indx)
      g1st_OS_pcov <- compute_pCOV_samps(samps = samps_deriv, intervals = OS_GCR_g1st, indx = indx)
      
      g_result_OS <- data.frame(Method = "OS", k = k, COV = g_OS_cov, PCOV = g_OS_pcov, Type = "g", rMSE = g_RMSE_OS_compared)
      g1st_result_OS <- data.frame(Method = "OS", k = k, COV = g1st_OS_cov, PCOV = g1st_OS_pcov, Type = "g1st", rMSE = g1st_RMSE_OS_compared)
      result_OS <- rbind(g_result_OS, g1st_result_OS)
      
      g_RW2_pcov <- compute_pCOV_samps(samps = samps_func, intervals = RW2_GCR_g, indx = indx)
      g1st_RW2_pcov <- compute_pCOV_samps(samps = samps_deriv, intervals = RW2_GCR_g1st, indx = indx)
      
      g_result_RW <- data.frame(Method = "RW", k = k, PCOV = g_RW2_pcov, Type = "g", rMSE = g_RMSE_RW_compared)
      g1st_result_RW <- data.frame(Method = "RW", k = k, PCOV = g1st_RW2_pcov, Type = "g1st", rMSE = g1st_RMSE_RW_compared)
      result_RW <- rbind(g_result_RW, g1st_result_RW)
      
      g_BS_pcov <- compute_pCOV_samps(samps = samps_func, intervals = BS_GCR_g, indx = indx)
      g1st_BS_pcov <- compute_pCOV_samps(samps = samps_deriv, intervals = BS_GCR_g1st, indx = indx)
      
      g_result_BS <- data.frame(Method = "BS", k = k, PCOV = g_BS_pcov, Type = "g", rMSE = g_RMSE_BS_compared)
      g1st_result_BS <- data.frame(Method = "BS", k = k, PCOV = g1st_BS_pcov, Type = "g1st", rMSE = g1st_RMSE_BS_compared)
      result_BS <- rbind(g_result_BS, g1st_result_BS)
      
      overall_result <- rbind(result_OS, result_RW, result_BS)
      return(overall_result)
    }
    result_all <- foreach(i = 1:length(kvec), .combine = "rbind", .packages = c("mgcv", "gratia", "TMB", "aghq", "Matrix", "tidyverse", "GET", "tsDyn", "LaplacesDemon")) %dopar% {
      source("16_functions_defined.R") ## Only on Windows
      dyn.load(dynlib("00_Gaussian_Smoothing_known")) ## Only on Windows
      do_once(k = kvec[i])
    }
    g_result <- result_all %>% filter(Type == "g") %>% select(-Type)
    g1st_result <- result_all %>% filter(Type == "g1st") %>% select(-Type)
  }
  else{
    pb <- progress_bar$new(total = length(kvec))
    overall_result <- data.frame()
    n <- length(x)
    PD_condition <- FALSE
    while(!PD_condition){
      RandomFunc <- sim_Wp_Var(mesh_size = 0.001, max_t = 10, p = p)
      RandomFunc[,-1] <- sigmaS * RandomFunc[,-1]
      
      # The function
      RandomFunc_IWPp <- RandomFunc[,1:2]
      RandomFunc_IWPp$IWp <- RandomFunc_IWPp$IWp
      # The first derivative, so q = p - 1
      RandomFunc_IWPq <- RandomFunc[,c(1,3)]
      RandomFunc_IWPq$IWq <- RandomFunc_IWPq$IWq
      
      gx <- evaluate_GP_Var(x, RandomFunc)[,1:2]
      names(gx) <- c("t", "IWp")
      gx1st <- evaluate_GP_Var(x, RandomFunc)[,c(1,3)]
      names(gx1st) <- c("t", "IWq")
      y <- gx$IWp + rnorm(n = n, mean = 0, sd = sigmaE)
      ## Fit the exact method
      exact_result <- Imple_WP_diff_bound(sigmaE = sigmaE, theta = theta, p = p, x = x, y = y, optimizer = T)
      PD_condition <- is.positive.definite(exact_result$prec)
    }
    
    ### Convert to samples of functions
    samps_func <- extract_samples_WP_rand_bound(result = exact_result, x = x, target = "function", p = p, nsam = nsamps)
    samps_deriv <- extract_samples_WP_rand_bound(result = exact_result, x = x, target = "derivative", p = p, nsam = nsamps)
    exact_g_mean <- samps_func[,c(-1)] %>% apply(1, mean)
    exact_g1st_mean <- samps_deriv[,c(-1)] %>% apply(1, mean)
    
    overall_result <- data.frame()
    for (i in 1:length(kvec)) {
      k <- kvec[i]
      pb$tick()
      OSmod <- Imple_OS_diffuse_bound(sigmaE = sigmaE, theta = theta , x = x, k = (k+1), y = y, p = p, optimizer = T)
      OS_GCR_g <- extract_mean_interval_given_samps(OSmod$samples_fun, level = level)
      OS_GCR_g1st <- extract_mean_interval_given_samps(OSmod$sample_deriv, level = level)
      
      diagonal <- 0
      RW2mod <- Imple_RW_diff_bound(sigmaE = sigmaE, theta = theta, x = x, k = k, y = y, p = 2, optimizer = T, diagonal = diagonal)
      names(RW2mod$samples_fun)[1] <- "x"
      RW2_GCR_g <- extract_mean_interval_given_samps(RW2mod$samples_fun, level = level)
      RW2_GCR_g1st <- extract_mean_interval_given_samps(RW2mod$sample_deriv, level = level)
      
      
      
      ##### Fit P-Spline:
      ### cubic basis and second order penalty
      BSmod <- Imple_MGCV_diffuse_bound(sigmaE = sigmaE, theta = theta, x = x, k = k, y = y, p = p)
      
      ### Samples from the fitted GAM:
      BS_g_samps <- extract_mgcv_samps(BSmod = BSmod, x = x, nsamps = nsamps)
      BS_g1st_samps <- extract_mgcv_deriv_samps(BSmod = BSmod, x = x, nsamps = nsamps)
      BS_g1st_samps <- cbind(x,BS_g1st_samps)
      
      
      BS_GCR_g <- extract_mean_interval_given_samps(BS_g_samps, level = level)
      BS_GCR_g1st <- extract_mean_interval_given_samps(BS_g1st_samps, level = level)
      
      
      OS_g_mean <- apply(OSmod$samples_fun[,-1], 1, mean)
      RW2_g_mean <- apply(RW2mod$samples_fun[,-1], 1, mean)
      OS_g1st_mean <- apply(OSmod$sample_deriv[,-1], 1, mean)
      RW2_g1st_mean <- apply(RW2mod$sample_deriv[,-1], 1, mean)
      
      g_RMSE_OS_compared <- sqrt(mean(((OS_g_mean-exact_g_mean)^2)[indx])) # by definition
      g1st_RMSE_OS_compared <- sqrt(mean(((OS_g1st_mean-exact_g1st_mean)^2)[indx])) # by definition
      
      g_RMSE_RW_compared <- sqrt(mean(((RW2_g_mean-exact_g_mean)^2)[indx])) # by definition
      g1st_RMSE_RW_compared <- sqrt(mean(((RW2_g1st_mean-exact_g1st_mean)^2)[indx])) # by definition
      
      g_RMSE_BS_compared <- sqrt(mean(((BS_GCR_g$mean-exact_g_mean)^2)[indx])) # by definition
      g1st_RMSE_BS_compared <- sqrt(mean(((BS_GCR_g1st$mean-exact_g1st_mean)^2)[indx])) # by definition
      
      
      g_OS_pcov <- compute_pCOV_samps(samps = samps_func, intervals = OS_GCR_g, indx = indx)
      g1st_OS_pcov <- compute_pCOV_samps(samps = samps_deriv, intervals = OS_GCR_g1st, indx = indx)
      
      g_result_OS <- data.frame(Method = "OS", k = k, PCOV = g_OS_pcov, Type = "g", rMSE = g_RMSE_OS_compared)
      g1st_result_OS <- data.frame(Method = "OS", k = k, PCOV = g1st_OS_pcov, Type = "g1st", rMSE = g1st_RMSE_OS_compared)
      result_OS <- rbind(g_result_OS, g1st_result_OS)
      
      g_RW2_pcov <- compute_pCOV_samps(samps = samps_func, intervals = RW2_GCR_g, indx = indx)
      g1st_RW2_pcov <- compute_pCOV_samps(samps = samps_deriv, intervals = RW2_GCR_g1st, indx = indx)
      
      g_result_RW <- data.frame(Method = "RW", k = k, PCOV = g_RW2_pcov, Type = "g", rMSE = g_RMSE_RW_compared)
      g1st_result_RW <- data.frame(Method = "RW", k = k, PCOV = g1st_RW2_pcov, Type = "g1st", rMSE = g1st_RMSE_RW_compared)
      result_RW <- rbind(g_result_RW, g1st_result_RW)
      
      g_BS_pcov <- compute_pCOV_samps(samps = samps_func, intervals = BS_GCR_g, indx = indx)
      g1st_BS_pcov <- compute_pCOV_samps(samps = samps_deriv, intervals = BS_GCR_g1st, indx = indx)
      
      g_result_BS <- data.frame(Method = "BS", k = k, PCOV = g_BS_pcov, Type = "g", rMSE = g_RMSE_BS_compared)
      g1st_result_BS <- data.frame(Method = "BS", k = k, PCOV = g1st_BS_pcov, Type = "g1st", rMSE = g1st_RMSE_BS_compared)
      result_BS <- rbind(g_result_BS, g1st_result_BS)
      
      overall_result <- rbind(overall_result, result_OS, result_RW, result_BS)
      
    }
    g_result <- overall_result %>% filter(Type == "g") %>% select(-Type)
    g1st_result <- overall_result %>% filter(Type == "g1st") %>% select(-Type)
  }
  return(list(g = g_result, g1st = g1st_result))
}



#### Aggregate the above simulation through B replications
agg_compare_diff_v3 <- function(B = 300, x, x_compare, p, kvec, sigmaE, theta,
                                control = list(method = "GET", GET = "rank", level = 0.8), 
                                nsamps = 8000, parallel = FALSE){
  ### Aggregate the result B times:
  
  if(parallel == FALSE){
    result_cov_boundary_g_agg <- data.frame()
    result_cov_boundary_g1st_agg <- data.frame()
    
    for (i in 1:B) {
      result_cov_boundary_new <- compare_all_three_methods_B_times_sample_cov_diffuse_v3(x = x, x_compare = x_compare, p = p, kvec = kvec,
                                                                                         sigmaE = sigmaE, theta = theta, seed = i,
                                                                                         control = control, 
                                                                                         nsamps = nsamps, parallel = FALSE)
      result_cov_boundary_new$g$Rep = (i+1)
      result_cov_boundary_new$g1st$Rep = (i+1)
      result_cov_boundary_g_agg <- rbind(result_cov_boundary_g_agg, result_cov_boundary_new$g)
      result_cov_boundary_g1st_agg <- rbind(result_cov_boundary_g1st_agg, result_cov_boundary_new$g1st)
    }
  }
  else{
    do_once <- function(seed){
      result_cov_boundary_new <- compare_all_three_methods_B_times_sample_cov_diffuse_v3(x = x, x_compare = x_compare, p = p, kvec = kvec,
                                                                                         sigmaE = sigmaE, theta = theta, seed = seed,
                                                                                         control = control, 
                                                                                         nsamps = nsamps, parallel = FALSE)
      result_cov_boundary_new$g$Type <- "g"
      result_cov_boundary_new$g1st$Type <- "g1st"
      rbind(result_cov_boundary_new$g, result_cov_boundary_new$g1st)
    }
    result_all <- foreach(i = 1:B, .combine = "rbind", .packages = c("progress","mgcv", "gratia", "TMB", "aghq", "Matrix", "tidyverse", "GET", "tsDyn", "LaplacesDemon")) %dopar% {
      ## Only on Windows
      source("16_functions_defined.R")
      source("02_functions_random_boundary.R")
      dyn.load(dynlib("00_Gaussian_Smoothing_unknown_bound"))
      dyn.load(dynlib("00_Gaussian_Smoothing_diffuse_bound"))
      dyn.load(dynlib("00_RW2_Smoothing_diff_bound"))
      do_once(seed = i)
    }
    result_cov_boundary_g_agg <- result_all %>% filter(Type == "g") %>% select(-Type)
    result_cov_boundary_g1st_agg <- result_all %>% filter(Type == "g1st") %>% select(-Type)
  }
  
  return(list(g_result = result_cov_boundary_g_agg, g1st_result = result_cov_boundary_g1st_agg))
  
}




