#include <Rcpp.h>
using namespace Rcpp;

//' @title Model error ME in AMA example
 //' @description
 //' Compute the model error
 //' \eqn{ME(\hat\beta) = (\hat\beta - \beta_0)^\top \Sigma (\hat\beta - \beta_0)}
 //' used in the AMA Example 1 simulation, implemented in C++ for efficiency.
 //'
 //' @param beta_hat Numeric vector of estimated coefficients \eqn{\hat\beta}.
 //' @param beta0 Numeric vector of true coefficients \eqn{\beta_0}.
 //' @param Sigma Numeric covariance matrix \eqn{\Sigma} (p by p).
 //'
 //' @return A scalar numeric value giving the model error.
 //'
 //' @examples
 //' \dontrun{
 //'   p <- 5
 //'   beta0    <- rep(1, p)
 //'   beta_hat <- beta0 + rnorm(p, 0, 0.1)
 //'   Sigma    <- diag(p)
 //'   compute_ME_cpp(beta_hat, beta0, Sigma)
 //' }
 //'
 //' @export
 // [[Rcpp::export]]
 double compute_ME_cpp(const NumericVector& beta_hat,
                       const NumericVector& beta0,
                       const NumericMatrix& Sigma) {
   int p = beta_hat.size();
   if (beta0.size() != p)
     stop("beta_hat and beta0 must have the same length");
   if (Sigma.nrow() != p || Sigma.ncol() != p)
     stop("Sigma must be a p x p matrix.");
   
   // diff = beta_hat - beta0
   NumericVector diff(p);
   for (int j = 0; j < p; ++j) {
     diff[j] = beta_hat[j] - beta0[j];
   }
   
   // tmp = Sigma * diff
   NumericVector tmp(p);
   for (int i = 0; i < p; ++i) {
     double s = 0.0;
     for (int j = 0; j < p; ++j) {
       s += Sigma(i, j) * diff[j];
     }
     tmp[i] = s;
   }
   
   // ME = diff' * tmp
   double me = 0.0;
   for (int j = 0; j < p; ++j) {
     me += diff[j] * tmp[j];
   }
   return me;
 }

//' @title Sparsity measures C and IC in AMA example
 //' @description
 //' Compute the sparsity measures used in Zhang (2022) AMA Example 1:
 //' \itemize{
 //'   \item \code{C}: number of truly zero coefficients estimated as zero;
 //'   \item \code{IC}: number of truly non-zero coefficients incorrectly estimated as zero.
 //' }
 //' Implemented in C++ and used inside the simulation grid.
 //'
 //' @param beta_hat Numeric vector of estimated coefficients.
 //' @param beta0 Numeric vector of true coefficients.
 //' @param tol Numeric tolerance to decide whether a coefficient is treated as zero
 //'   (default \code{1e-6}).
 //'
 //' @return A named numeric vector with elements \code{C} and \code{IC}.
 //'
 //' @examples
 //' \dontrun{
 //'   beta0    <- c(1, 0, 0.5, 0, 0)
 //'   beta_hat <- c(0.9, 0, 0, 0, 0.1)
 //'   compute_C_IC_cpp(beta_hat, beta0)
 //' }
 //'
 //' @export
 // [[Rcpp::export]]
 NumericVector compute_C_IC_cpp(const NumericVector& beta_hat,
                                const NumericVector& beta0,
                                double tol = 1e-6) {
   int p = beta_hat.size();
   if (beta0.size() != p)
     stop("beta_hat and beta0 must have the same length");
   
   int C = 0;
   int IC = 0;
   
   for (int j = 0; j < p; ++j) {
     bool true_zero = std::fabs(beta0[j]) < tol;
     bool est_zero  = std::fabs(beta_hat[j]) < tol;
     
     if (true_zero && est_zero)  C++;
     if (!true_zero && est_zero) IC++;
   }
   
   NumericVector out(2);
   out[0] = C;
   out[1] = IC;
   out.names() = CharacterVector::create("C", "IC");
   return out;
 }
