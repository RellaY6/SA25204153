#' @importFrom MASS mvrnorm
#' @importFrom lars lars
#' @importFrom quadprog solve.QP
NULL
## ==========================================================
## Adaptive Model Averaging(AMA) Numerical Example 1: 内部辅助函数
## ==========================================================

ama_step1_alasso_lars <- function(X, y, gamma_alasso = 1, tol_active = 1e-8) {
  # X: n x p 设计矩阵
  # y: n 维响应向量
  # gamma_alasso: ALASSO 权重中的指数，默认 gamma_alasso = 1
  # tol_active: 判断“非零系数”的阈值
  
  stopifnot(is.matrix(X), length(y) == nrow(X))
  n <- nrow(X); p <- ncol(X)
  
  if (p >= n) {
    warning("当前 p >= n，不是 AMA 设定的 p < n 场景，OLS 初始估计可能不可用或不稳定。")
  }
  
  ## 0. 预处理：标准化 X，中心化 y
  X_scaled   <- scale(X)
  y_centered <- as.numeric(scale(y, center = TRUE, scale = FALSE))
  
  ## 1. 初始估计 beta^(init)：OLS
  lm_init   <- lm(y_centered ~ X_scaled)
  beta_init <- coef(lm_init)[-1]  # 去掉截距
  
  ## 2. 构造 ALASSO 权重
  eps     <- 1e-6
  weights <- 1 / (abs(beta_init) + eps)^gamma_alasso
  
  ## 3. 重新缩放设计矩阵: X*_j = X_scaled_j / w_j
  X_star <- sweep(X_scaled, 2, weights, FUN = "/")
  
  ## 4. 在 (X*, y_centered) 上跑 LARS-LASSO 路径
  if (!requireNamespace("lars", quietly = TRUE)) {
    stop("Package 'lars' is required. Please install it first.")
  }
  
  lars_fit <- lars::lars(
    x          = X_star,
    y          = y_centered,
    type       = "lasso",
    normalize  = FALSE,
    intercept  = TRUE
  )
  
  beta_lars_path   <- lars_fit$beta  # K x p
  beta_alasso_path <- sweep(beta_lars_path, 2, weights, FUN = "/")
  
  ## 5. 提取每一步 active set
  K <- nrow(beta_alasso_path)
  active_sets <- vector("list", K)
  for (k in seq_len(K)) {
    bk <- beta_alasso_path[k, ]
    active_sets[[k]] <- which(abs(bk) > tol_active)
  }
  
  list(
    X_scaled         = X_scaled,
    y_centered       = y_centered,
    beta_init        = beta_init,
    weights          = weights,
    lars_fit         = lars_fit,
    beta_alasso_path = beta_alasso_path,
    active_sets      = active_sets
  )
}

ama_step2_select_H_mBIC <- function(X, y, step1_res,
                                    use_scaled = TRUE) {
  n <- nrow(X); p <- ncol(X)
  
  if (use_scaled) {
    X_use <- step1_res$X_scaled
    y_use <- step1_res$y_centered
  } else {
    X_use <- scale(X)
    y_use <- as.numeric(scale(y, center = TRUE, scale = FALSE))
  }
  
  active_sets <- step1_res$active_sets
  
  ## 去重
  keys <- sapply(active_sets, function(idx) paste(sort(idx), collapse = ","))
  keep <- !duplicated(keys)
  active_sets_uniq <- active_sets[keep]
  K <- length(active_sets_uniq)
  
  mBIC_values <- numeric(K)
  sigma2_hat  <- numeric(K)
  df_vec      <- integer(K)
  
  for (k in seq_len(K)) {
    S   <- active_sets_uniq[[k]]
    d_s <- length(S)
    df_vec[k] <- d_s
    
    if (d_s == 0) {
      resid <- y_use
    } else {
      Xk    <- X_use[, S, drop = FALSE]
      fit_k <- lm(y_use ~ Xk - 1)
      resid <- residuals(fit_k)
    }
    
    sigma2_hat[k] <- mean(resid^2)
    mBIC_values[k] <- n * log(sigma2_hat[k]) +
      d_s * log(n) * log(log(p))
  }
  
  best_k <- which.min(mBIC_values)
  H      <- active_sets_uniq[[best_k]]
  
  ## 如果整体最优是空模型，再在非空模型里选一次
  if (length(H) == 0) {
    nonempty <- which(df_vec > 0)
    if (length(nonempty) == 0) {
      stop("LARS 路径中没有非空 active set。")
    }
    best_k <- nonempty[which.min(mBIC_values[nonempty])]
    H      <- active_sets_uniq[[best_k]]
    message("mBIC 最优为空模型，改用最佳非空模型 k = ", best_k)
  }
  
  list(
    H               = H,
    best_index      = best_k,
    mBIC_values     = mBIC_values,
    df              = df_vec,
    sigma2_hat      = sigma2_hat,
    active_sets_uniq = active_sets_uniq
  )
}

ama_step34_build_candidates <- function(X, y,
                                        step1_res,
                                        step2_res) {
  n <- nrow(X); p <- ncol(X)
  
  active_sets <- step1_res$active_sets
  K <- length(active_sets)
  
  H <- sort(unique(step2_res$H))
  q_star <- length(H)
  if (q_star == 0) {
    stop("mBIC 选出的 H 为空，无法继续。")
  }
  
  ## A. 每个变量的“进入时间”
  entry_step <- rep(Inf, p)
  for (k in seq_len(K)) {
    idx <- active_sets[[k]]
    if (length(idx) == 0) next
    for (j in idx) {
      if (is.infinite(entry_step[j])) {
        entry_step[j] <- k
      }
    }
  }
  
  steps_H   <- entry_step[H]
  ord_H     <- order(steps_H)
  H_ordered <- H[ord_H]   # He1,...,He_{q*}
  
  ## B. Step 2 的 q* 个嵌套模型
  nested_sets <- vector("list", q_star)
  for (j in seq_len(q_star)) {
    nested_sets[[j]] <- sort(H_ordered[1:j])
  }
  
  ## C. Step 3: “包含 H 的更大模型”
  k_all <- max(steps_H)
  extra_sets <- list()
  if (is.finite(k_all)) {
    for (k in seq(from = k_all, to = K)) {
      S <- active_sets[[k]]
      if (length(S) == 0) next
      if (!all(H %in% S)) next
      extra_sets[[length(extra_sets) + 1]] <- sort(S)
    }
  }
  
  ## D. 合并去重
  all_sets <- c(nested_sets, extra_sets)
  keys     <- sapply(all_sets, function(idx) paste(idx, collapse = ","))
  keep     <- !duplicated(keys)
  candidate_sets <- all_sets[keep]
  q_n           <- length(candidate_sets)
  
  ## E. 对每个候选模型做一次 OLS
  beta_hat_mat <- matrix(0, nrow = q_n, ncol = p)
  X_list       <- vector("list", q_n)
  df_vec       <- integer(q_n)
  
  for (s in seq_len(q_n)) {
    S <- candidate_sets[[s]]
    df_vec[s] <- length(S)
    if (length(S) == 0) {
      X_list[[s]] <- NULL
      next
    }
    Xs   <- X[, S, drop = FALSE]
    XtX  <- crossprod(Xs)
    Xty  <- crossprod(Xs, y)
    beta_sub  <- solve(XtX, Xty)
    beta_full <- rep(0, p)
    beta_full[S] <- beta_sub
    beta_hat_mat[s, ] <- beta_full
    X_list[[s]]       <- Xs
  }
  
  list(
    H_ordered      = H_ordered,
    nested_sets    = nested_sets,
    candidate_sets = candidate_sets,
    q_n            = q_n,
    beta_hat_mat   = beta_hat_mat,
    X_list         = X_list,
    df             = df_vec
  )
}

ama_weights_mu_for_y <- function(X, y,
                                 step34_res,
                                 phi_n,
                                 sigma2_hat,
                                 ridge_eps = 1e-4) {
  if (!requireNamespace("quadprog", quietly = TRUE)) {
    stop("Package 'quadprog' is required. Please install it first.")
  }
  
  n  <- nrow(X)
  qn <- step34_res$q_n
  
  candidate_sets <- step34_res$candidate_sets
  df_vec         <- step34_res$df
  
  ## 1. 构造 M(y): n x q_n
  M <- matrix(0, nrow = n, ncol = qn)
  for (s in seq_len(qn)) {
    S <- candidate_sets[[s]]
    if (length(S) == 0) {
      M[, s] <- 0
    } else {
      Xs   <- X[, S, drop = FALSE]
      XtX  <- crossprod(Xs)
      Xty  <- crossprod(Xs, y)
      beta_sub <- solve(XtX, Xty)
      M[, s]   <- as.numeric(Xs %*% beta_sub)
    }
  }
  
  ## 2. S_n(w) 的二次型部分
  A <- crossprod(M)
  b <- crossprod(M, y)
  v <- df_vec
  
  Dmat <- 2 * A
  Dmat <- (Dmat + t(Dmat)) / 2
  eig      <- eigen(Dmat, symmetric = TRUE, only.values = TRUE)$values
  min_eig  <- min(eig)
  if (min_eig <= 0) {
    Dmat <- Dmat + (abs(min_eig) + ridge_eps) * diag(qn)
  } else {
    Dmat <- Dmat + ridge_eps * diag(qn)
  }
  
  dvec <- as.numeric(2 * b - phi_n * sigma2_hat * v)
  
  ## 3. 约束: sum w = 1, w >= 0
  Amat <- cbind(rep(1, qn), diag(qn))
  bvec <- c(1, rep(0, qn))
  meq  <- 1
  
  sol   <- quadprog::solve.QP(Dmat, dvec, Amat, bvec, meq = meq)
  w_hat <- sol$solution
  w_hat[w_hat < 0] <- 0
  w_hat <- w_hat / sum(w_hat)
  
  mu_hat <- as.numeric(M %*% w_hat)
  Sn_val <- sum((y - mu_hat)^2) + phi_n * sigma2_hat * sum(v * w_hat)
  
  list(
    w_hat    = w_hat,
    mu_hat   = mu_hat,
    Sn_value = Sn_val
  )
}

ama_step5_6_section3 <- function(X, y,
                                 step2_res,
                                 step34_res,
                                 n_phi_grid = 20,
                                 B_mc       = 50,
                                 seed       = NULL) {
  n <- nrow(X)
  p <- ncol(X)
  
  ## 1. 用 H 模型估计 sigma^2
  H <- step2_res$H
  if (length(H) == 0) stop("mBIC 选出的 H 为空，无法估计 sigma^2。")
  X_H   <- X[, H, drop = FALSE]
  fit_H <- lm(y ~ X_H - 1)
  resid_H <- residuals(fit_H)
  df_H    <- length(H)
  sigma2_hat <- sum(resid_H^2) / (n - df_H)
  
  ## 2. φ_n 搜索区间
  low  <- log(log(p))
  high <- log(log(p)) * log(n)
  if (!is.finite(low) || low <= 0) low <- 1
  if (high <= low) high <- low + 1
  phi_grid <- seq(from = low, to = high, length.out = n_phi_grid)
  
  ## 3. MC 设置
  tau <- 0.5 * sqrt(sigma2_hat)
  if (!is.null(seed)) set.seed(seed)
  
  qn     <- step34_res$q_n
  df_vec <- step34_res$df
  
  C_vals  <- numeric(n_phi_grid)
  Sn_vals <- numeric(n_phi_grid)
  g0_vals <- numeric(n_phi_grid)
  w_list  <- vector("list", n_phi_grid)
  mu_list <- vector("list", n_phi_grid)
  
  for (i in seq_along(phi_grid)) {
    phi <- phi_grid[i]
    
    base_res <- ama_weights_mu_for_y(
      X, y, step34_res,
      phi_n      = phi,
      sigma2_hat = sigma2_hat
    )
    w_hat  <- base_res$w_hat
    mu_hat <- base_res$mu_hat
    
    acc <- 0
    for (b_idx in seq_len(B_mc)) {
      delta   <- rnorm(n, mean = 0, sd = tau)
      y_tilde <- y + delta
      res_tilde <- ama_weights_mu_for_y(
        X, y_tilde, step34_res,
        phi_n      = phi,
        sigma2_hat = sigma2_hat
      )
      mu_tilde <- res_tilde$mu_hat
      acc <- acc + sum(delta * mu_tilde)
    }
    g0_hat <- 2 * acc / (tau^2 * B_mc)
    
    C_vals[i]  <- sum((mu_hat - y)^2) + sigma2_hat * g0_hat
    Sn_vals[i] <- base_res$Sn_value
    
    w_list[[i]]  <- w_hat
    mu_list[[i]] <- mu_hat
    g0_vals[i]   <- g0_hat
  }
  
  best_i  <- which.min(C_vals)
  phi_hat <- phi_grid[best_i]
  w_hat   <- w_list[[best_i]]
  mu_hat  <- mu_list[[best_i]]
  
  beta_ma <- as.numeric(w_hat %*% step34_res$beta_hat_mat)
  
  list(
    phi_grid   = phi_grid,
    C_vals     = C_vals,
    Sn_vals    = Sn_vals,
    g0_hat     = g0_vals,
    best_index = best_i,
    phi_hat    = phi_hat,
    w_hat      = w_hat,
    mu_hat     = mu_hat,
    beta_ma    = beta_ma,
    sigma2_hat = sigma2_hat,
    tau        = tau
  )
}

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

select_H_criterion <- function(X, y, step1_res,
                               use_scaled = TRUE,
                               criterion = c("mBIC", "eBIC")) {
  criterion <- match.arg(criterion)
  
  n <- nrow(X); p <- ncol(X)
  
  if (use_scaled) {
    X_use <- step1_res$X_scaled
    y_use <- step1_res$y_centered
  } else {
    X_use <- scale(X)
    y_use <- as.numeric(scale(y, center = TRUE, scale = FALSE))
  }
  
  active_sets <- step1_res$active_sets
  
  keys <- sapply(active_sets, function(idx) paste(sort(idx), collapse = ","))
  keep <- !duplicated(keys)
  active_sets_uniq <- active_sets[keep]
  K <- length(active_sets_uniq)
  
  crit_values <- numeric(K)
  sigma2_hat  <- numeric(K)
  df_vec      <- integer(K)
  
  gamma_eBIC <- 1 - 2 * log(n) / log(p)
  
  for (k in seq_len(K)) {
    S   <- active_sets_uniq[[k]]
    d_s <- length(S)
    df_vec[k] <- d_s
    
    if (d_s == 0) {
      resid <- y_use
    } else {
      Xk    <- X_use[, S, drop = FALSE]
      fit_k <- lm(y_use ~ Xk - 1)
      resid <- residuals(fit_k)
    }
    
    sigma2_hat[k] <- mean(resid^2)
    
    if (criterion == "mBIC") {
      crit_values[k] <- n * log(sigma2_hat[k]) +
        d_s * log(n) * log(log(p))
    } else {
      crit_values[k] <- n * log(sigma2_hat[k]) +
        d_s * log(n) +
        2 * gamma_eBIC * lchoose(p, d_s)
    }
  }
  
  best_k <- which.min(crit_values)
  H      <- active_sets_uniq[[best_k]]
  
  list(
    H                = sort(unique(H)),
    best_index       = best_k,
    crit_values      = crit_values,
    df               = df_vec,
    sigma2_hat       = sigma2_hat,
    active_sets_uniq = active_sets_uniq,
    gamma_eBIC       = gamma_eBIC
  )
}

alasso_select_final <- function(X, y, step1_res, crit_res) {
  p <- ncol(X)
  S <- crit_res$H
  
  beta_final <- rep(0, p)
  if (length(S) > 0) {
    Xs <- X[, S, drop = FALSE]
    beta_sub <- solve(crossprod(Xs), crossprod(Xs, y))
    beta_final[S] <- beta_sub
  }
  beta_final
}

sim_example1_one_setting <- function(
    n,
    rho,
    R2,
    B_rep      = 100,
    n_phi_grid = 20,
    B_mc       = 50,
    seed       = 2025
) {
  set.seed(seed)
  
  RME_ama   <- numeric(B_rep)
  RME_mBIC  <- numeric(B_rep)
  RME_eBIC  <- numeric(B_rep)
  
  C_ama  <- IC_ama  <- numeric(B_rep)
  C_mBIC <- IC_mBIC <- numeric(B_rep)
  C_eBIC <- IC_eBIC <- numeric(B_rep)
  
  for (b in seq_len(B_rep)) {
    dat   <- gen_example1_data(n = n, rho = rho, R2 = R2)
    X     <- dat$X
    y     <- dat$y
    beta0 <- dat$beta0
    Sigma <- dat$Sigma
    
    beta_ls <- solve(crossprod(X), crossprod(X, y))
    ME_ls   <- compute_ME_cpp(beta_ls, beta0, Sigma)
    
    step1 <- ama_step1_alasso_lars(X, y, gamma_alasso = 1)
    
    step2_mBIC <- ama_step2_select_H_mBIC(X, y, step1)
    step34     <- ama_step34_build_candidates(X, y, step1, step2_mBIC)
    step56     <- ama_step5_6_section3(
      X, y,
      step2_res  = step2_mBIC,
      step34_res = step34,
      n_phi_grid = n_phi_grid,
      B_mc       = B_mc,
      seed       = b
    )
    beta_ama <- step56$beta_ma
    
    crit_mBIC <- select_H_criterion(X, y, step1, criterion = "mBIC")
    beta_mBIC <- alasso_select_final(X, y, step1, crit_mBIC)
    
    crit_eBIC <- select_H_criterion(X, y, step1, criterion = "eBIC")
    beta_eBIC <- alasso_select_final(X, y, step1, crit_eBIC)
    
    ME_ama  <- compute_ME_cpp(beta_ama,  beta0, Sigma)
    ME_mBIC <- compute_ME_cpp(beta_mBIC, beta0, Sigma)
    ME_eBIC <- compute_ME_cpp(beta_eBIC, beta0, Sigma)
    
    RME_ama[b]  <- ME_ama  / ME_ls
    RME_mBIC[b] <- ME_mBIC / ME_ls
    RME_eBIC[b] <- ME_eBIC / ME_ls
    
    CIC_ama  <- compute_C_IC_cpp(beta_ama,  beta0)
    CIC_mBIC <- compute_C_IC_cpp(beta_mBIC, beta0)
    CIC_eBIC <- compute_C_IC_cpp(beta_eBIC, beta0)
    
    C_ama[b]  <- CIC_ama["C"];   IC_ama[b]  <- CIC_ama["IC"]
    C_mBIC[b] <- CIC_mBIC["C"];  IC_mBIC[b] <- CIC_mBIC["IC"]
    C_eBIC[b] <- CIC_eBIC["C"];  IC_eBIC[b] <- CIC_eBIC["IC"]
    
    if (b %% 10 == 0) cat("rep", b, "done\n")
  }
  
  list(
    MRME = c(
      AMA         = median(RME_ama),
      ALASSO_mBIC = median(RME_mBIC),
      ALASSO_eBIC = median(RME_eBIC)
    ),
    C = c(
      AMA         = mean(C_ama),
      ALASSO_mBIC = mean(C_mBIC),
      ALASSO_eBIC = mean(C_eBIC)
    ),
    IC = c(
      AMA         = mean(IC_ama),
      ALASSO_mBIC = mean(IC_mBIC),
      ALASSO_eBIC = mean(IC_eBIC)
    )
  )
}

## ==========================================================
## 下面两个是“用户级函数”，加 roxygen comments 并 @export
## ==========================================================

#' @title Adaptive Model Averaging 拟合函数（Example 1 设计）
#'
#' @description
#' 在 AMA 论文
#' Zhang (2022), "Adaptive Model Averaging Estimation with a Diverging Number of Parameters"
#' 的 Example 1 设定下, 完成自适应模型平均 (AMA) 的整个拟合流程。
#' 给定设计矩阵 \code{X} 和响应向量 \code{y}，
#' 按照自适应模型平均 (Adaptive Model Averaging, AMA)
#' (ALASSO 初始模型选择、mBIC 选取集合 \eqn{H}, 构造候选模型集合,
#' 以及 Cp-type 准则选择正则化参数 \eqn{\phi_n}, 最终给出 AMA 模型平均估计。)
#' 文献 Example 1 的设定，完成：
#' \enumerate{
#'   \item Step 1：ALASSO via LARS，得到一条自适应 LASSO 解路径；
#'   \item Step 2：基于 mBIC 从路径中选出候选变量集合 \eqn{H}；
#'   \item Step 3–4：构造候选模型集合并对每个模型做 OLS 估计；
#'   \item Step 5–6：在 \eqn{\phi_n} 网格上用 Cp-type 准则选取
#'         \eqn{\hat\phi_n}，得到最终 AMA 模型平均估计。
#' }
#' 函数返回最终 AMA 估计以及各步的中间结果，便于进一步分析和作图。
#'
#' @references
#' Zhang, X. (2022).
#' \emph{Adaptive Model Averaging Estimation With A Diverging Number of Parameters}.
#'
#' @param X 数值矩阵，维度 \eqn{n \times p}，每列对应一个协变量。
#' @param y 数值向量，长度 \eqn{n}，响应变量。
#' @param n_phi_grid 搜索 \eqn{\phi_n} 时网格上的点数，默认 20。
#' @param B_mc Monte Carlo 重复次数，用于估计 \eqn{g_0(\phi_n)}，默认 50。
#' @param seed 随机种子，方便重复结果。若为 \code{NULL} 则不设定。
#'
#' @return 一个列表，包含：
#' \itemize{
#'   \item \code{beta_ma} 最终 AMA 模型平均估计 \eqn{\hat\beta(\hat\phi_n)}。
#'   \item \code{phi_hat} 选出的 \eqn{\hat\phi_n}。
#'   \item \code{w_hat} 最终的模型平均权重向量。
#'   \item \code{sigma2_hat} 使用集合 \eqn{H} 拟合得到的 \eqn{\hat\sigma^2}。
#'   \item \code{H} mBIC 选出的候选变量集合 \eqn{H}。
#'   \item \code{candidate_sets} Step 2+3 构造出的所有候选模型集合。
#'   \item \code{step1, step2, step34, step56} 四个中间步骤的完整输出。
#' }
#'
#' @examples
#' \dontrun{
#'   set.seed(1)
#'   dat <- gen_example1_data(n = 100, rho = 0.5, R2 = 0.4)
#'   fit <- ama_fit_example1(dat$X, dat$y,
#'                           n_phi_grid = 10, B_mc = 20)
#'   str(fit$beta_ma)
#' }
#'
#' @export
ama_fit_example1 <- function(X, y,
                             n_phi_grid = 20,
                             B_mc       = 50,
                             seed       = 2025) {
  step1 <- ama_step1_alasso_lars(X, y, gamma_alasso = 1)
  step2 <- ama_step2_select_H_mBIC(X, y, step1)
  step34 <- ama_step34_build_candidates(X, y, step1, step2)
  step56 <- ama_step5_6_section3(
    X, y,
    step2_res  = step2,
    step34_res = step34,
    n_phi_grid = n_phi_grid,
    B_mc       = B_mc,
    seed       = seed
  )
  
  list(
    beta_ma        = step56$beta_ma,
    phi_hat        = step56$phi_hat,
    w_hat          = step56$w_hat,
    sigma2_hat     = step56$sigma2_hat,
    H              = step2$H,
    candidate_sets = step34$candidate_sets,
    step1          = step1,
    step2          = step2,
    step34         = step34,
    step56         = step56
  )
}

#' @title AMA Numerical Example 1 模拟研究
#'
#' @description
#' 在 AMA 文献 Numerical Example 1 的设计下，对一系列
#' \code{(n, rho, R2)} 组合重复模拟 \code{B_rep} 次，
#' 比较三种方法的性能：
#' \itemize{
#'   \item AMA（自适应模型平均）；
#'   \item 基于 ALASSO 的单模型选择 + mBIC；
#'   \item 基于 ALASSO 的单模型选择 + eBIC。
#' }
#' 输出各组合下模型误差相对 LS 的中位数（MRME），以及
#' 稀疏性指标 \code{C} / \code{IC} 的平均值。
#'
#' @references
#' Zhang, X. (2022).
#' \emph{Adaptive Model Averaging Estimation With A Diverging Number of Parameters}.
#'
#' @param n_vec 数值向量，备选样本量，例如 \code{200} 或 \code{c(200, 400)}。
#' @param rho_vec 数值向量，相关系数 \eqn{\rho} 的备选值。
#' @param R2_vec 数值向量，信号强度 \eqn{R^2} 的备选值。
#' @param B_rep 每个设定下的重复次数，默认 100。
#' @param n_phi_grid \eqn{\phi_n} 网格点个数，默认 20。
#' @param B_mc 估计 \eqn{g_0(\phi_n)} 的 Monte Carlo 次数，默认 50。
#' @param seed 随机种子，默认 2025。
#'
#' @return 一个 \code{data.frame}，每一行对应一个
#' \code{(n, rho, R2, method)} 组合，包含列：
#' \itemize{
#'   \item \code{MRME}：模型误差相对 LS 的中位数；
#'   \item \code{C}：真正为零且被正确估计为零的系数个数；
#'   \item \code{IC}：真正非零却被估计为零的系数个数。
#' }
#'
#' @examples
#' \dontrun{
#'   res <- sim_example1_grid(
#'     n_vec   = 200,
#'     rho_vec = c(0.5, 0.75),
#'     R2_vec  = c(0.2, 0.4),
#'     B_rep   = 10,
#'     n_phi_grid = 10,
#'     B_mc    = 20
#'   )
#'   head(res)
#' }
#'
#' @export
sim_example1_grid <- function(
    n_vec      = 200,
    rho_vec    = c(0.5, 0.75),
    R2_vec     = c(0.2, 0.4, 0.6, 0.8),
    B_rep      = 100,
    n_phi_grid = 20,
    B_mc       = 50,
    seed       = 2025
) {
  res_all <- data.frame(
    n      = integer(0),
    rho    = numeric(0),
    R2     = numeric(0),
    method = character(0),
    MRME   = numeric(0),
    C      = numeric(0),
    IC     = numeric(0),
    stringsAsFactors = FALSE
  )
  
  for (n in n_vec) {
    for (rho in rho_vec) {
      for (R2 in R2_vec) {
        cat("Running n =", n, ", rho =", rho, ", R2 =", R2, "...\n")
        
        res_one <- sim_example1_one_setting(
          n          = n,
          rho        = rho,
          R2         = R2,
          B_rep      = B_rep,
          n_phi_grid = n_phi_grid,
          B_mc       = B_mc,
          seed       = seed
        )
        
        methods <- names(res_one$MRME)
        
        df_tmp <- data.frame(
          n      = rep(n,    length(methods)),
          rho    = rep(rho,  length(methods)),
          R2     = rep(R2,   length(methods)),
          method = methods,
          MRME   = as.numeric(res_one$MRME),
          C      = as.numeric(res_one$C),
          IC     = as.numeric(res_one$IC),
          stringsAsFactors = FALSE
        )
        
        res_all <- rbind(res_all, df_tmp)
      }
    }
  }
  
  res_all
}
