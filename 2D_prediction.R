get_2d_forecast_and_learning_samples <- function(Z, target, k) {
  stopifnot(is.matrix(Z), nrow(Z) == ncol(Z))
  n <- nrow(Z)
  
  stopifnot(length(target) == 2L)
  i0 <- as.integer(target[1])
  j0 <- as.integer(target[2])
  
  k <- as.integer(k)
  stopifnot(k >= 1L)
  
  cell_id <- function(i, j) {
    i + (j - 1L) * n
  }
  
  # forecast sample = k closest points inside {1..n}^2 to the target
  gridd <- as.matrix(expand.grid(i = 1:n, j = 1:n))
  
  dx <- gridd[, "i"] - i0
  dy <- gridd[, "j"] - j0
  d2 <- dx*dx + dy*dy
  ang <- atan2(dy, dx)
  
  ord <- order(d2, ang, gridd[, "i"], gridd[, "j"])
  
  d2_k <- d2[ord][k]
  cand <- ord[d2[ord] <= d2_k]
  if (length(cand) > k) {
    cand <- cand[order(d2[cand], ang[cand], gridd[cand, "i"], gridd[cand, "j"])][seq_len(k)]
  }
  
  nn_idx <- cand
  forecast_points <- gridd[nn_idx, ]
  forecast_sample <- Z[cbind(forecast_points[, 1], forecast_points[, 2])]
  
  offsets <- cbind(di = forecast_points[, 1] - i0, dj = forecast_points[, 2] - j0)
  
  # forecast sample cells (these must not be used by any learning sample)
  forecast_ids <- cell_id(forecast_points[, 1], forecast_points[, 2])
  
  # feasible shifted targets inside the grid
  centers <- as.matrix(expand.grid(it = 1:n, jt = 1:n))
  N0 <- nrow(centers)
  
  feasible <- logical(N0)
  for (m in seq_len(N0)) {
    it <- centers[m, "it"]
    jt <- centers[m, "jt"]
    ii <- it + offsets[, "di"]
    jj <- jt + offsets[, "dj"]
    feasible[m] <- all(ii >= 1 & ii <= n & jj >= 1 & jj <= n)
  }
  
  centers <- centers[feasible, ]
  if (!nrow(centers)) stop("No feasible learning samples with this target/k on this grid.")
  
  footprints <- vector("list", nrow(centers))
  for (m in seq_len(nrow(centers))) {
    it <- centers[m, "it"]
    jt <- centers[m, "jt"]
    ii <- it + offsets[, "di"]
    jj <- jt + offsets[, "dj"]
    footprints[[m]] <- sort(unique(cell_id(ii, jj)))
  }
  
  # greedy maximal set of non-overlapping learning samples
  used <- rep(FALSE, n*n)
  used[forecast_ids] <- TRUE
  
  chosen <- logical(length(footprints))
  for (m in seq_along(footprints)) {
    ids <- footprints[[m]]
    if (!any(used[ids])) {
      chosen[m] <- TRUE
      used[ids] <- TRUE
    }
  }
  
  centers_chosen <- centers[chosen, ]
  if (!nrow(centers_chosen)) stop("No non-overlapping learning samples found under these constraints.")
  
  # build learning_samples + observed_values
  N <- nrow(centers_chosen)
  learning_samples <- matrix(0, N, k)
  observed_values <- numeric(N)
  
  for (m in seq_len(N)) {
    it <- centers_chosen[m, "it"]
    jt <- centers_chosen[m, "jt"]
    ii <- it + offsets[, "di"]
    jj <- jt + offsets[, "dj"]
    
    learning_samples[m, ] <- Z[cbind(ii, jj)]
    observed_values[m] <- Z[it, jt]
  }
  
  return(list(forecast_sample = forecast_sample, 
              learning_samples = learning_samples, 
              observed_values = observed_values))
}


make_ring_targets <- function(n, r) {
  N <- n + r
  rbind(cbind(i = 1:(N-1L), j = N), cbind(i = N, j = 1:N))
}

make_targets_rings <- function(n, R) {
  do.call(rbind, lapply(1:R, function(r) make_ring_targets(n, r)))
}


pred_one_target <- function(Zobs, i0, j0, k_nn, gam, max_iter = 10000L, 
                            stag_pat = 1000L, bs = 1L, eta = 0.1, 
                            optimizer = "adam", reparam = TRUE, eval_every = 1L, 
                            plot_conv = FALSE, full_grad = FALSE) {
  res <- get_2d_forecast_and_learning_samples(Z = Zobs, target = c(i0, j0), k = k_nn)
  
  fit <- gradient_descent(lambda_start = rep(1, k_nn), 
                          learning_samples = res$learning_samples, 
                          observed_values = res$observed_values, 
                          gam = gam, 
                          max_iter = max_iter, stag_pat = stag_pat, 
                          bootstrap = FALSE, X_tilde = NULL, 
                          plot_conv = plot_conv, bs = bs, eta = eta, 
                          optimizer = optimizer, reparam = reparam, 
                          eval_every = eval_every, full_grad = full_grad)
  
  max(res$forecast_sample * fit$best_lambda)
}

extrapolate_rings_2d <- function(Zobs, R, k_nn, gam, max_iter = 10000L, 
                                 stag_pat = 1000L, bs = 1L, eta = 0.1, 
                                 optimizer = "adam", reparam = TRUE, 
                                 eval_every = 1L, plot_conv = FALSE, 
                                 full_grad = FALSE) {
  stopifnot(is.matrix(Zobs), nrow(Zobs) == ncol(Zobs))
  n <- nrow(Zobs)
  N <- n + R
  
  targets <- make_targets_rings(n, R)
  L <- nrow(targets)
  
  Z_hat <- matrix(0, N, N)
  Z_hat[1:n, 1:n] <- Zobs
  
  preds <- numeric(L)
  for (m in seq_len(L)) {
    i0 <- targets[m, "i"]
    j0 <- targets[m, "j"]
    preds[m] <- pred_one_target(Zobs = Zobs, i0 = i0, j0 = j0, k_nn = k_nn, 
                                gam = gam, max_iter = max_iter, stag_pat = stag_pat, 
                                bs = bs, eta = eta, optimizer = optimizer, 
                                reparam = reparam, eval_every = eval_every, 
                                plot_conv = plot_conv, full_grad = full_grad)
    Z_hat[i0, j0] <- preds[m]
  }
  
  list(Z_hat = Z_hat, targets = targets, preds = preds)
}


plot_extrapolated_surface <- function(Z_hat, n_obs, R, 
                                      show_pred_points = FALSE) {
  N <- n_obs + R
  stopifnot(nrow(Z_hat) == N, ncol(Z_hat) == N)
  
  cols_facet <- matrix("grey80", N-1, N-1)
  
  cols_facet[n_obs:(N-1), ] <- "mistyrose"
  cols_facet[, n_obs:(N-1)] <- "mistyrose"
  
  pmat <- persp(x = 1:N, y = 1:N, z = Z_hat, theta = 135, phi = 25, expand = 0.6, 
                xlab = "", ylab = "", zlab = "", ticktype = "detailed", 
                axes = TRUE, box = TRUE, main = "", col = cols_facet, 
                border = NA, shade = 0.4)
  
  if (isTRUE(show_pred_points)) {
    targets <- make_targets_rings(n_obs, R)
    xs <- targets[, "i"]
    ys <- targets[, "j"]
    zs <- Z_hat[cbind(xs, ys)]
    pp <- trans3d(xs, ys, zs, pmat)
    points(pp, pch = 16, col = "red", cex = 0.6)
  }
}

################################################################################
################################################################################

n_obs <- 50
R <- 10
k_nn <- 11

coords <- as.matrix(expand.grid(1:n_obs, 1:n_obs))

################################################################################

S <- cov_matrix(coords = coords, cov_type = "Levy", sigma2 = br_pars["1.6"])
R_br <- safe_chol(Sigma = S)
vdiff_mat_ <- vdiff_mat(Sigma = S)
Z <- algo1_exact(coords = coords, type = "brownresnick", R_br = R_br, vdiff_mat_ = vdiff_mat_)
Zobs <- matrix(as.numeric(Z), nrow = n_obs, ncol = n_obs)


Z <- algo1_exact(coords = coords, type = "smith", Sigma_for_smith = sm_pars["1.6"]*as.matrix(diag(2)))
Zobs <- matrix(as.numeric(Z), nrow = n_obs, ncol = n_obs)


K_eg <- cov_matrix(coords = coords, cov_type = "OU", a = eg_pars["1.6"])
R_eg <- safe_chol(Sigma = K_eg)
Z <- algo1_exact(coords = coords, type = "extremalgaussian", K_eg = K_eg, R_eg = R_eg)
Zobs <- matrix(as.numeric(Z), nrow = n_obs, ncol = n_obs)

################################################################################

out_rings <- extrapolate_rings_2d(Zobs = Zobs, 
                                  R = R, 
                                  k_nn = k_nn, 
                                  gam = 1.0, 
                                  max_iter = 10000L, 
                                  stag_pat = 1000L, 
                                  bs = 1L, 
                                  eta = 0.1, 
                                  optimizer = "adam", 
                                  reparam = TRUE, 
                                  eval_every = 1L, 
                                  plot_conv = TRUE, 
                                  full_grad = FALSE)

Z_hat_plot <- out_rings$Z_hat

dir.create("Results/2d_experiments/eg/1_6", recursive = TRUE, showWarnings = FALSE)
write.csv(Z_hat_plot, file = "Results/2d_experiments/eg/1_6/Z_hat_plot.csv", row.names = FALSE)


Z_hat_plot <- as.matrix(read.csv("Results/2d_experiments/eg/1_6/Z_hat_plot.csv"))
Z_hat_plot <- log(Z_hat_plot)

plot_extrapolated_surface(Z_hat = Z_hat_plot, n_obs = n_obs, R = R)




