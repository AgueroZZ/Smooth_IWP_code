#include <TMB.hpp>                                // Links in the TMB libraries
//#include <fenv.h>

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(y); //response variable
  DATA_SPARSE_MATRIX(X_global); // Design matrix (for fixed global effects)
  DATA_SPARSE_MATRIX(Xf); // Design matrix (for fixed covariate effects)
  DATA_SPARSE_MATRIX(B); // Design matrix (for random effects)
  DATA_SPARSE_MATRIX(P); // Penalty matrix
  DATA_SPARSE_MATRIX(I); // Design matrix (for OLRE to model over-dispersion)
  
  int d = P.cols(); // Number of B-Spline coefficients
  DATA_SCALAR(logPdet); // Determinant of (fixed) penalty matrix
  DATA_SCALAR(u1); // pc prior, u param
  DATA_SCALAR(alpha1); // pc prior, alpha param
  DATA_SCALAR(u2); // pc prior, u param (for eps)
  DATA_SCALAR(alpha2); // pc prior, alpha param (for eps)
  DATA_SCALAR(betaprec); // beta ~iid N(0,1/betaprec)
  DATA_SCALAR(Xfprec); // beta ~iid N(0,1/betaprec)


  // Parameter
  PARAMETER_VECTOR(W); // W = c(U,beta,eps), eta = B * U + X * beta
  int Wdim = W.size();
  int n = y.size();
  int beta1dim = X_global.cols();
  int beta2dim = Xf.cols();

  vector<Type> U(d);
  vector<Type> beta1(beta1dim);
  vector<Type> beta2(beta2dim);
  vector<Type> eps(n);
  for (int i=0;i<d;i++) U(i) = W(i);
  for (int i=0;i<beta1dim;i++) beta1(i) = W(i+d);
  for (int i=0;i<beta2dim;i++) beta2(i) = W(i+beta1dim+d);
  for (int i=0;i<n;i++) eps(i) = W(i+d+beta1dim+beta2dim);
  PARAMETER(theta1); // theta1 = -2log(sigma1)
  PARAMETER(theta2); // theta2 = -2log(sigma2)


  // Transformations
  vector<Type> eta = X_global * beta1 + Xf * beta2 + B * U + I * eps;
  Type sigma1 = exp(-0.5*theta1);
  Type sigma2 = exp(-0.5*theta2);
  REPORT(sigma1);
  REPORT(sigma2);


  // Log likelihood
  Type ll = 0;
  ll = sum(dpois(y, exp(eta), TRUE));
  REPORT(ll);
  
  
  // Log prior on W
  Type lpW = 0;

  // Cross product
  vector<Type> PU = P*U;
  Type UPU = (U * PU).sum();
  lpW += -0.5 * exp(theta1) * UPU; // U part
  Type bb = (beta1 * beta1).sum();
  lpW += -0.5 * betaprec * bb; // Beta part
  Type bb2 = (beta2 * beta2).sum();
  lpW += -0.5 * Xfprec * bb2; // Beta part
  Type ee = (eps * eps).sum();
  lpW += -0.5 * exp(theta2) * ee; // eps part

  // Log determinant
  Type logdet1 = d * theta1 + logPdet;
  lpW += 0.5 * logdet1; // P part
  Type logdet2 = n * theta2;
  lpW += 0.5 * logdet2; // over-dispersion part
  REPORT(logdet1);
  REPORT(lpW);

  
  // Log prior for theta
  Type lpT = 0;
  Type phi1 = -log(alpha1) / u1;
  lpT += log(0.5 * phi1) - phi1*exp(-0.5*theta1) - 0.5*theta1;
  Type phi2 = -log(alpha2) / u2;
  lpT += log(0.5 * phi2) - phi2*exp(-0.5*theta2) - 0.5*theta2;
  REPORT(lpT);
  
  // Final result!
  Type logpost = -1 * (ll + lpW + lpT);
  
  return logpost;
}