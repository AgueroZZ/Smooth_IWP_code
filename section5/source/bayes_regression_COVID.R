### Implement Poission regression for the COVID example:
# compile("00_Poisson_Smoothing_PC_overdisp_covid.cpp")
# dyn.load(dynlib("00_Poisson_Smoothing_PC_overdisp_covid"))
library(mgcv)

### The proposed method:
Imple_BayesRegression_COVID <- function(y, x, knots = NULL, p = 2, prior = NULL, aghq_k = 7, Xf = NULL) {
  if (is.null(prior)) {
    u1 = 1
    alpha1 = 0.5
    betaprec = 10^(-6)
    Xfprec = 10^(-6)
    u2 = 1
    alpha2 = 0.5
  }
  else{
    u1 = prior$u1
    alpha1 = prior$alpha1
    betaprec = prior$betaprec
    Xfprec = prior$Xfprec
    u2 = prior$u2
    alpha2 = prior$alpha2
  }
    if (is.null(knots)) {
      P_proposed <- compute_weights_precision(x)
      D_proposed <- compute_weights_design(x, p)
      ## Design matrix for the global polynomial
      X_global = rbind(c(1, rep(0, p - 1)), D_proposed[, 1:p])
      ## Design matrix for the spline basis weights
      B = rbind(rep(0, (length(x) - 1)), D_proposed[, (p + 1):ncol(D_proposed)])
      I = diag(nrow = length(y), ncol = length(y))
      tmbdat <- list(
        # Design matrix
        X_global = as(X_global, "dgTMatrix"),
        Xf = as(Xf, "dgTMatrix"),
        B = B,
        I = as(I, "dgTMatrix"),
        P = as(P_proposed, 'dgTMatrix'),
        logPdet = as.numeric(determinant(P_proposed, logarithm = T)$modulus),
        # Response
        y = y,
        # PC Prior params
        u1 = u1,
        alpha1 = alpha1,
        u2 = u2,
        alpha2 = alpha2,
        betaprec = betaprec,
        Xfprec = Xfprec
      )
      tmbparams <- list(W = c(rep(0, (
        ncol(X_global) + ncol(Xf) + ncol(B) + ncol(I)
      ))),
      # W = c(U,beta); U = B-Spline coefficients
      theta1 = 0,
      theta2 = 0)
      # dyn.load(dynlib("00_Poisson_Smoothing_PC_overdisp_covid"))
      ff <- TMB::MakeADFun(
        data = tmbdat,
        parameters = tmbparams,
        random = "W",
        DLL = "00_Poisson_Smoothing_PC_overdisp_covid",
        silent = TRUE
      )
      # Hessian not implemented for RE models
      ff$he <- function(w)
        numDeriv::jacobian(ff$gr, w)
    }
    else{
      P_proposed <- compute_weights_precision(knots)
      D_proposed <- compute_weights_design(x, p)
      ## Design matrix for the global polynomial
      X_global = rbind(global_poly(x[1], p), D_proposed[, 1:p])
      ## Design matrix for the spline basis weights
      B = as(local_poly(
        x = knots,
        refined_x = x,
        p = p
      ), "dgTMatrix")
      I = diag(nrow = length(y), ncol = length(y))
      tmbdat <- list(
        # Design matrix
        X_global = as(X_global, "dgTMatrix"),
        Xf = as(Xf, "dgTMatrix"),
        B = B,
        I = as(I, "dgTMatrix"),
        P = as(P_proposed, 'dgTMatrix'),
        logPdet = as.numeric(determinant(P_proposed, logarithm = T)$modulus),
        # Response
        y = y,
        # PC Prior params
        u1 = u1,
        alpha1 = alpha1,
        u2 = u2,
        alpha2 = alpha2,
        betaprec = betaprec,
        Xfprec = Xfprec
      )
      tmbparams <- list(W = c(rep(0, (
        ncol(X_global) + ncol(Xf) + ncol(B) + ncol(I)
      ))),
      # W = c(U,beta); U = B-Spline coefficients
      theta1 = 0,
      theta2 = 0)
      # dyn.load(dynlib("00_Poisson_Smoothing_PC_overdisp_covid"))
      ff <- TMB::MakeADFun(
        data = tmbdat,
        parameters = tmbparams,
        random = "W",
        DLL = "00_Poisson_Smoothing_PC_overdisp_covid",
        silent = TRUE
      )
      # Hessian not implemented for RE models
      ff$he <- function(w)
        numDeriv::jacobian(ff$gr, w)
    }
    aghq::marginal_laplace_tmb(ff, aghq_k, c(0, 0))
}



### The O-Sullivan method:
Imple_OSullivan_COVID <- function(region = NULL, y, x, k, p = 2, prior = NULL, aghq_k = 7, Xf = NULL){
  if (is.null(prior)) {
    u1 = 1
    alpha1 = 0.5
    betaprec = 10^(-6)
    Xfprec = 10^(-6)
    u2 = 1
    alpha2 = 0.5
  }
  else{
    u1 = prior$u1
    alpha1 = prior$alpha1
    betaprec = prior$betaprec
    Xfprec = prior$Xfprec
    u2 = prior$u2
    alpha2 = prior$alpha2
  }
  if (is.null(region)) {
    region <- c(min(x), max(x))
    basis <- create_Osullivan_basis(region = region, k = k, p = p)
    P <- as(create_Osullivan_Q(basis = basis, p = p), "dgTMatrix")
    B <- as(create_Osullivan_design(basis = basis, x = x), "dgTMatrix")

    ## Design matrix for the spline basis weights
    I = diag(nrow = length(y), ncol = length(y))
    tmbdat <- list(
      # Design matrix
      Xf = as(Xf, "dgTMatrix"),
      B = B,
      I = as(I, "dgTMatrix"),
      P = as(P, 'dgTMatrix'),
      logPdet = sum(log(eigen(P, only.values = T)$values[1:(k-p)])),
      # Response
      y = y,
      # PC Prior params
      u1 = u1,
      alpha1 = alpha1,
      u2 = u2,
      alpha2 = alpha2,
      Xfprec = Xfprec
    )
    tmbparams <- list(W = c(rep(0, (
      ncol(Xf) + ncol(B) + ncol(I)
    ))),
    # W = c(U,beta); U = B-Spline coefficients
    theta1 = 0,
    theta2 = 0)
    # dyn.load(dynlib("00_Poisson_Smoothing_PC_overdisp_covid_osullivan"))
    ff <- TMB::MakeADFun(
      data = tmbdat,
      parameters = tmbparams,
      random = "W",
      DLL = "00_Poisson_Smoothing_PC_overdisp_covid_osullivan",
      silent = TRUE
    )
    # Hessian not implemented for RE models
    ff$he <- function(w)
      numDeriv::jacobian(ff$gr, w)
  }
  else{
    basis <- create_Osullivan_basis(region = region, k = k, p = p)
    P <- as(create_Osullivan_Q(basis = basis, p = p), "dgTMatrix")
    B <- as(create_Osullivan_design(basis = basis, x = x), "dgTMatrix")
    
    ## Design matrix for the spline basis weights
    I = diag(nrow = length(y), ncol = length(y))
    tmbdat <- list(
      # Design matrix
      Xf = as(Xf, "dgTMatrix"),
      B = B,
      I = as(I, "dgTMatrix"),
      P = as(P, 'dgTMatrix'),
      logPdet = sum(log(eigen(P, only.values = T)$values[1:(k-p)])),
      # Response
      y = y,
      # PC Prior params
      u1 = u1,
      alpha1 = alpha1,
      u2 = u2,
      alpha2 = alpha2,
      Xfprec = Xfprec
    )
    tmbparams <- list(W = c(rep(0, (
      ncol(Xf) + ncol(B) + ncol(I)
    ))),
    # W = c(U,beta); U = B-Spline coefficients
    theta1 = 0,
    theta2 = 0)
    # dyn.load(dynlib("00_Poisson_Smoothing_PC_overdisp_covid_osullivan"))
    ff <- TMB::MakeADFun(
      data = tmbdat,
      parameters = tmbparams,
      random = "W",
      DLL = "00_Poisson_Smoothing_PC_overdisp_covid_osullivan",
      silent = TRUE
    )
    # Hessian not implemented for RE models
    ff$he <- function(w)
      numDeriv::jacobian(ff$gr, w)
  }
  aghq::marginal_laplace_tmb(ff, aghq_k, c(0, 0))
}




#### Method based on mgcv:
Imple_mgcv_COVID <- function(y, x, k = 100, p = 3, Xf = NULL){
  data <- data.frame(y = y, x = x, ID = factor(c(1:nrow(Xf))))
  data <- cbind(data, Xf)
  mod <- mgcv::gam(formula = y ~ -1 + weekdays1 + weekdays2 + weekdays3 + weekdays4 + weekdays5 + weekdays6 + s(x, k = k, bs="bs",m=c(5,3)),
                   family = nb(), data = data, control = list(scalePenalty=F), method = "ML")
  mod
}


Sample_mgcv_COVID <- function(mod, nsamps){
  data <- mod$model
  all <- predict.gam(mod, newdata = data, type = 'lpmatrix')
  bs.design <- as.matrix(all[,-c(1:6)])
  bs.weights.samples <- t(rmvn(n = nsamps, mu = mod$coefficients[-c(1:6)], Sigma = vcov(mod)[-c(1:6), -c(1:6)]))
  bs.fit <- bs.design %*% bs.weights.samples
  bs.fit
}
