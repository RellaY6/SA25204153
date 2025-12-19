#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double compute_ME_cpp(const NumericVector& beta_hat,
                      const NumericVector& beta0,
                      const NumericMatrix& Sigma) {
  int p = beta_hat.size();
  if (beta0.size() != p) stop("beta_hat and beta0 must have the same length.");
  if (Sigma.nrow() != p || Sigma.ncol() != p) stop("Sigma must be a p x p matrix.");
  
  double me = 0.0;
  for (int i = 0; i < p; ++i) {
    double di = beta_hat[i] - beta0[i];
    double s  = 0.0;
    for (int j = 0; j < p; ++j) {
      s += Sigma(i, j) * (beta_hat[j] - beta0[j]);
    }
    me += di * s;
  }
  return me;
}

// [[Rcpp::export]]
NumericVector compute_C_IC_cpp(const NumericVector& beta_hat,
                               const NumericVector& beta0,
                               double tol = 1e-6) {
  int p = beta_hat.size();
  if (beta0.size() != p) stop("beta_hat and beta0 must have the same length.");
  
  int C = 0, IC = 0;
  for (int j = 0; j < p; ++j) {
    bool true_zero = std::fabs(beta0[j]) < tol;
    bool est_zero  = std::fabs(beta_hat[j]) < tol;
    if (true_zero && est_zero)  ++C;
    if (!true_zero && est_zero) ++IC;
  }
  
  NumericVector out = NumericVector::create(_["C"] = C, _["IC"] = IC);
  return out;
}
