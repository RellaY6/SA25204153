#' Model error (ME) in AMA Example 1 (C++ implementation)
#'
#' Compute the model error
#' \deqn{ME(\hat\beta) = (\hat\beta - \beta_0)^\top \Sigma (\hat\beta - \beta_0)}
#' used in the AMA Example 1 simulation.
#'
#' @param beta_hat Numeric vector of estimated coefficients \eqn{\hat\beta}.
#' @param beta0 Numeric vector of true coefficients \eqn{\beta_0}.
#' @param Sigma Numeric covariance matrix \eqn{\Sigma} (p by p).
#'
#' @return A scalar numeric value giving the model error.
#'
#' @examples
#' \dontrun{
#' p <- 5
#' beta0    <- rep(1, p)
#' beta_hat <- beta0 + rnorm(p, 0, 0.1)
#' Sigma    <- diag(p)
#' compute_ME_cpp(beta_hat, beta0, Sigma)
#' }
#'
#' @name compute_ME_cpp
#' @rdname compute_ME_cpp
#' @export
NULL


#' Sparsity measures C and IC in AMA Example 1 (C++ implementation)
#'
#' Compute the sparsity measures used in AMA Example 1:
#' \itemize{
#'   \item \code{C}: number of truly zero coefficients estimated as zero;
#'   \item \code{IC}: number of truly non-zero coefficients incorrectly estimated as zero.
#' }
#'
#' @param beta_hat Numeric vector of estimated coefficients.
#' @param beta0 Numeric vector of true coefficients.
#' @param tol Numeric tolerance for deciding whether a coefficient is treated as zero
#'   (default \code{1e-6}).
#'
#' @return A named numeric vector with elements \code{C} and \code{IC}.
#'
#' @examples
#' \dontrun{
#' beta0    <- c(1, 0, 0.5, 0, 0)
#' beta_hat <- c(0.9, 0, 0, 0, 0.1)
#' compute_C_IC_cpp(beta_hat, beta0)
#' }
#'
#' @name compute_C_IC_cpp
#' @rdname compute_C_IC_cpp
#' @export
NULL
