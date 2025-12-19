
gen_example1_data <- function(n, rho = 0.5, R2 = 0.4) {
  ## 1. 维度：p_n = [4 n^{1/2}] - 5
  p <- 4 * floor(sqrt(n)) - 5
  
  ## 2. q = [p_n / 9], 非零个数 |A| = 3q
  q  <- floor(p / 9)
  nonzero <- 3 * q
  stopifnot(nonzero <= p)
  
  ## 3. 协方差矩阵
  idx   <- 1:p
  Sigma <- rho^abs(outer(idx, idx, "-"))
  
  ## 4. 真系数
  beta0 <- c(rep(3, nonzero), rep(0, p - nonzero))
  
  ## 5. 生成 X ~ N(0, Sigma)
  X <- MASS::mvrnorm(n = n, mu = rep(0, p), Sigma = Sigma)
  
  ## 6. 根据 R2 决定噪声方差
  var_signal <- as.numeric(t(beta0) %*% Sigma %*% beta0)
  sigma2     <- var_signal * (1 - R2) / R2
  eps        <- rnorm(n, mean = 0, sd = sqrt(sigma2))
  
  ## 7. 生成 y
  y <- as.numeric(X %*% beta0 + eps)
  
  list(
    X       = X,
    y       = y,
    beta0   = beta0,
    Sigma   = Sigma,
    sigma2  = sigma2,
    p       = p,
    q       = q,
    nonzero = nonzero
  )
}