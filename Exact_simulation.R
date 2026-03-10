# Implementation of Algorithm 1 from the paper "Exact simulation of max-stable
# processes" (2016) by Dombry, Engelke and Oesting

Ornstein_Uhlenbeck_cov <- function(s, t, a) {
  stopifnot(is.numeric(a), length(a) == 1L, is.finite(a), a > 0)
  s <- as.numeric(s); t <- as.numeric(t)
  stopifnot(length(s) == length(t), length(s) > 0L)
  stopifnot(all(is.finite(s)), all(is.finite(t)))
  d <- sqrt(drop(crossprod(s-t)))
  return(exp(-d/a))
}

Levy_cov <- function(s, t, sigma2) {
  stopifnot(is.numeric(sigma2), length(sigma2) == 1L, is.finite(sigma2), sigma2 > 0)
  s <- as.numeric(s); t <- as.numeric(t)
  stopifnot(length(s) == length(t), length(s) > 0L)
  stopifnot(all(is.finite(s)), all(is.finite(t)))
  0.5 * sigma2 * (sqrt(drop(crossprod(s))) + sqrt(drop(crossprod(t))) - sqrt(drop(crossprod(s-t))))
}

# Build covariance matrix
cov_matrix <- function(coords, cov_type = c("Levy", "OU"), sigma2 = NULL, a = NULL) {
  stopifnot(is.matrix(coords), is.numeric(coords), nrow(coords) > 0L, ncol(coords) > 0L)
  stopifnot(all(is.finite(coords)))
  N <- nrow(coords)
  Sigma <- matrix(0, N, N)
  cov_type <- match.arg(cov_type)
  
  if (cov_type == "Levy" && is.null(sigma2)) stop("Parameter 'sigma2' must be provided for Lévy covariance.")
  if (cov_type == "OU" && is.null(a)) stop("Parameter 'a' must be provided for OU covariance.")
  
  kfun <- if (cov_type == "Levy") {
    function(i, j) Levy_cov(s = coords[i, ], t = coords[j, ], sigma2 = sigma2)
  } else {
    function(i, j) Ornstein_Uhlenbeck_cov(s = coords[i, ], t = coords[j, ], a = a)
  }
  
  for (i in seq_len(N)) {
    Sigma[i, i] <- kfun(i, i)
    if (i < N) {
      for (j in (i+1L):N) {
        cij <- kfun(i, j)
        Sigma[i, j] <- cij
        Sigma[j, i] <- cij
      }
    }
  }
  
  return(Sigma)
}


safe_chol <- function(Sigma) {
  R <- try(chol(Sigma), silent = TRUE)
  if (inherits(R, "try-error")) stop("Error: Sigma has no Cholesky decomposition")
  return(R)
}

vdiff_mat <- function(Sigma) {
  d <- diag(Sigma)
  outer(d, d, "+") - 2*Sigma
}


sample_W <- function(R) {
  stopifnot(is.matrix(R), is.numeric(R), nrow(R) > 0L, nrow(R) == ncol(R))
  stopifnot(all(is.finite(R)))
  N <- nrow(R)
  W <- t(R) %*% rnorm(N) # Gaussian vector with covariance Sigma = R^T * R
  return(as.vector(W))
}

anchored_Y_from_W <- function(W, j0, vdiff_mat_) {
  stopifnot(is.numeric(W), all(is.finite(W)))
  stopifnot(is.matrix(vdiff_mat_), is.numeric(vdiff_mat_), all(is.finite(vdiff_mat_)))
  
  N <- length(W)
  stopifnot(nrow(vdiff_mat_) == N, ncol(vdiff_mat_) == N)
  stopifnot(is.integer(j0), length(j0) == 1L, j0 >= 1L, j0 <= N)
  
  vdiff_ <- vdiff_mat_[, j0]
  exp(W - W[j0] - 0.5*vdiff_)
}


# Sampler for Smith process (Proposition 3)
r_Px_smith_1d <- function(j0, coords, sigma2_for_smith) {
  stopifnot(is.matrix(coords), is.numeric(coords), all(is.finite(coords)))
  stopifnot(nrow(coords) > 0L, ncol(coords) == 1L)
  stopifnot(is.numeric(sigma2_for_smith), is.finite(sigma2_for_smith), 
            length(sigma2_for_smith) == 1L, sigma2_for_smith > 0)
  
  N <- nrow(coords)
  stopifnot(is.integer(j0), length(j0) == 1L, j0 >= 1L, j0 <= N)
  
  x0 <- drop(coords[j0, ])
  sd <- sqrt(sigma2_for_smith)
  chi_var <- rnorm(1L, sd = sd)
  num <- dnorm(coords + chi_var - x0, sd = sd)
  den <- dnorm(chi_var, sd = sd)
  return(drop(num/den))
}

r_Px_smith_2d <- function(j0, coords, Sigma_for_smith, Sigma_inv = NULL, R_chol = NULL) {
  stopifnot(is.matrix(coords), is.numeric(coords), all(is.finite(coords)))
  stopifnot(nrow(coords) > 0L, ncol(coords) > 1L)
  stopifnot(is.matrix(Sigma_for_smith), is.numeric(Sigma_for_smith), 
            all(is.finite(Sigma_for_smith)), 
            nrow(Sigma_for_smith) == ncol(Sigma_for_smith))
  N <- nrow(coords)
  stopifnot(is.integer(j0), length(j0) == 1L, j0 >= 1L, j0 <= N)
  
  d <- ncol(coords)
  stopifnot(nrow(Sigma_for_smith) == d)
  
  x0 <- coords[j0, ]
  
  if (is.null(R_chol)) R_chol <- safe_chol(Sigma = Sigma_for_smith)
  if (is.null(Sigma_inv)) Sigma_inv <- chol2inv(R_chol)
  
  chi <- drop(t(R_chol) %*% rnorm(d))
  
  Delta <- sweep(coords, 2L, x0, "-")
  AD <- Delta %*% Sigma_inv
  q <- rowSums(AD * Delta)
  b <- as.vector(AD %*% chi)
  
  exp(-0.5 * q - b)
}


rt_process_from_mu_scale <- function(K_eg, R_eg, j0) {
  stopifnot(is.matrix(K_eg), is.numeric(K_eg), all(is.finite(K_eg)))
  stopifnot(is.matrix(R_eg), is.numeric(R_eg), all(is.finite(R_eg)))
  stopifnot(nrow(K_eg) == ncol(K_eg), nrow(R_eg) == ncol(R_eg), nrow(K_eg) == nrow(R_eg))
  
  N <- nrow(K_eg)
  stopifnot(is.integer(j0), length(j0) == 1L, j0 >= 1L, j0 <= N)
  
  r_j0 <- R_eg[, j0]
  mu <- K_eg[, j0]
  z <- rnorm(N)
  y <- z - r_j0 * as.numeric(crossprod(r_j0, z))
  X <- drop(crossprod(R_eg, y))
  U <- rchisq(1L, df = 2L)
  return(mu + X / sqrt(U))
}

# Algorithm 1 requires to draw a process Y ~ P_x0. The next function 
# simulates Y using Propositions 3-5
getY <- function(coords, type = c("brownresnick", "smith", "extremalgaussian"), 
                 j0, R_br = NULL, vdiff_mat_ = NULL, K_eg = NULL, R_eg = NULL, 
                 sigma2_for_smith = NULL, Sigma_for_smith = NULL) {
  type <- match.arg(type)
  
  if (type == "brownresnick") {
    if (is.null(R_br) || is.null(vdiff_mat_)) stop("'R_br' and 'vdiff_mat_' must be provided for Brown-Resnick.")
    W <- sample_W(R = R_br)
    Y <- anchored_Y_from_W(W = W, j0 = j0, vdiff_mat_ = vdiff_mat_)
    return(Y)
  } else if (type == "smith") {
    if (ncol(coords) == 1L) {
      if (is.null(sigma2_for_smith)) stop("'sigma2_for_smith' must be provided for Smith 1D.")
      Y <- r_Px_smith_1d(j0 = j0, coords = coords, sigma2_for_smith = sigma2_for_smith)
      return(Y)
    } else {
      if (is.null(Sigma_for_smith)) stop("'Sigma_for_smith' must be provided for Smith 2D.")
      R_smith <- safe_chol(Sigma_for_smith)
      SigInv <- chol2inv(R_smith)
      Y <- r_Px_smith_2d(j0 = j0, coords = coords, Sigma_for_smith = Sigma_for_smith, 
                         Sigma_inv = SigInv, R_chol = R_smith)
      return(Y)
    }
  } else {
    if (is.null(K_eg) || is.null(R_eg)) stop("'K_eg' and 'R_eg' must be provided for Extremal Gaussian.")
    T_ <- rt_process_from_mu_scale(K_eg = K_eg, R_eg = R_eg, j0 = j0)
    Y <- pmax(T_, 0)
    return(Y)
  }
}

# Algorithm 1
algo1_exact <- function(coords, type = c("brownresnick", "smith", "extremalgaussian"), ...) {
  stopifnot(is.matrix(coords), is.numeric(coords), all(is.finite(coords)))
  stopifnot(ncol(coords) > 0L)
  stopifnot(nrow(coords) > 1L)
  type <- match.arg(type)
  N <- nrow(coords)
  xi_inv <- rexp(1L); xi <- 1/xi_inv
  Y <- getY(coords = coords, type = type, j0 = 1L, ...)
  
  Z <- xi * Y
  
  for (n in 2L:N) {
    prev <- 1L:(n-1L)
    xi_inv <- rexp(1L); xi <- 1/xi_inv
    
    while (xi > Z[n]) {
      Y <- getY(coords = coords, type = type, j0 = n, ...)
      cand <- xi * Y
      if (all(cand[prev] < Z[prev])) {
        idx <- cand > Z
        if (any(idx)) Z[idx] <- cand[idx]
      }
      xi_inv <- xi_inv + rexp(1L); xi <- 1/xi_inv
    }
  }
  
  return(as.vector(Z))
}

# draw multiple processes and store them row wise in a matrix M
draw_multiple_1d_processes <- function(coords, type = c("brownresnick", "smith", "extremalgaussian"), numb, ...) {
  stopifnot(is.matrix(coords), is.numeric(coords), all(is.finite(coords)))
  stopifnot(ncol(coords) == 1L)
  stopifnot(nrow(coords) > 1L)
  stopifnot(is.integer(numb), length(numb) == 1L, numb >= 1L)
  
  type <- match.arg(type)
  N <- nrow(coords)
  
  M <- matrix(0, numb, N)
  for (i in seq_len(numb)) {
    M[i, ] <- algo1_exact(coords = coords, type = type, ...)
    if (i %% 100L == 0L) message(sprintf("process %d done", i))
  }
  
  return(M)
}

################################################################################
################################################################################
################################################################################


# A) 1D Simulation

get_br_parameter <- function(theta) {
  2*qnorm(theta/2)
}

get_sm_parameter <- function(theta) {
  1 / (2*qnorm(theta/2))
}

get_eg_parameter <- function(theta) {
  - 1 / log(1-2*(theta-1)^2)
}


theta_vec <- seq(1.1, 1.9, by = 0.1)

br_pars <- get_br_parameter(theta = theta_vec)
names(br_pars) <- as.character(theta_vec)

sm_pars <- get_sm_parameter(theta = theta_vec)
names(sm_pars) <- as.character(theta_vec)

theta_vec_eg <- seq(1.1, 1.7, by = 0.1)

eg_pars <- get_eg_parameter(theta = theta_vec_eg)
names(eg_pars) <- as.character(theta_vec_eg)

################################################################################


simulate_1d_data <- function(N, numb, dir_name, seed_number, 
                             br_pars = NULL, sm_pars = NULL, eg_pars = NULL) {
  stopifnot(is.integer(N), length(N) == 1L, N > 1L)
  stopifnot(is.integer(numb), length(numb) == 1L, numb >= 1L)
  
  coords <- matrix(1L:N, ncol = 1L)
  dir.create(dir_name, recursive = TRUE, showWarnings = FALSE)
  
  set.seed(seed_number)
  # A1) Brown-Resnick
  if (!is.null(br_pars)) {
    for (k in seq_along(br_pars)) {
      Sigma <- cov_matrix(coords = coords, cov_type = "Levy", sigma2 = br_pars[k]^2)
      R_br <- safe_chol(Sigma = Sigma)
      vdiff_mat_ <- vdiff_mat(Sigma = Sigma)
      
      M <- draw_multiple_1d_processes(coords = coords, type = "brownresnick", 
                                      numb = numb, R_br = R_br, vdiff_mat_ = vdiff_mat_)
      saveRDS(object = M, file = paste0(dir_name, "/M_br_", names(br_pars)[k], ".rds"))
    }
  }
  
  # A2) Smith
  if (!is.null(sm_pars)) {
    for (k in seq_along(sm_pars)) {
      M <- draw_multiple_1d_processes(coords = coords, type = "smith", numb = numb, 
                                      sigma2_for_smith = sm_pars[k]^2)
      saveRDS(object = M, file = paste0(dir_name, "/M_sm_", names(sm_pars)[k], ".rds"))
    }
  }
  
  # A3) Extremal Gaussian
  if (!is.null(eg_pars)) {
    for (k in seq_along(eg_pars)) {
      K_eg <- cov_matrix(coords = coords, cov_type = "OU", a = eg_pars[k])
      R_eg <- safe_chol(Sigma = K_eg)
      
      M <- draw_multiple_1d_processes(coords = coords, type = "extremalgaussian", 
                                      numb = numb, K_eg = K_eg, R_eg = R_eg)
      saveRDS(object = M, file = paste0(dir_name, "/M_eg_", names(eg_pars)[k], ".rds"))
    }
  }
}


# simulate 1d data for the "pre-experiments"
simulate_1d_data(N = 203L, 
                 numb = 1000L, 
                 dir_name = "Data/pre_experiments/1", 
                 seed_number = 100, 
                 br_pars = br_pars, 
                 sm_pars = sm_pars, 
                 eg_pars = eg_pars)

simulate_1d_data(N = 203L, 
                 numb = 1000L, 
                 dir_name = "Data/pre_experiments/2", 
                 seed_number = 200, 
                 br_pars = br_pars, 
                 sm_pars = sm_pars, 
                 eg_pars = eg_pars)

simulate_1d_data(N = 203L, 
                 numb = 1000L, 
                 dir_name = "Data/pre_experiments/3", 
                 seed_number = 300, 
                 br_pars = br_pars, 
                 sm_pars = sm_pars, 
                 eg_pars = eg_pars)

simulate_1d_data(N = 203L, 
                 numb = 10000L, 
                 dir_name = "Data/pre_experiments", 
                 seed_number = 1000, 
                 br_pars = br_pars, 
                 sm_pars = sm_pars, 
                 eg_pars = eg_pars)


# simulate data for the main experiment (20 step forecast)
simulate_1d_data(N = 2141, numb = 1000, dir_name = "Data/20step_prediction", seed_number = 100, 
                 br_pars = br_pars, sm_pars = sm_pars, eg_pars = eg_pars)





