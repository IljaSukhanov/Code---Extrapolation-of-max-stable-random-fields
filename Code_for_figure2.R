# compute mse like in remark 5 of the paper
mse <- function(X_hat, M = 1000L, alpha = 1, loc = 0, sca = 1) {
  U_hat <- H_cdf(X_hat, alpha = alpha, loc = loc, sca = sca)
  u_grid <- (1:M) / (M + 1)
  F_emp <- ecdf(U_hat)(u_grid)
  mean((F_emp - u_grid)^2)
}

# evaluate objective function on a grid
make_grid <- function(Z, t_end, FH, len, overlap = 0, gam, bootstrap = FALSE, 
                      X_tilde = NULL, alpha = 1, loc = 0, sca = 1, 
                      lambda_from = 0.01, lambda_to = 2, lambda_by = 0.01) {
  
  if (isTRUE(bootstrap) && is.null(X_tilde)) {
    stop("Error: for the bootstrap version 'X_tilde' must be provided.")
  }
  
  ls_pre <- get_learning_samples(Z = Z[1:t_end], FH = FH, len = len, 
                                 overlap = overlap)
  
  learning_samples <- ls_pre$learning_samples
  observed_values <- ls_pre$observed_values
  
  lambda_1 <- seq(lambda_from, lambda_to, by = lambda_by)
  lambda_2 <- seq(lambda_from, lambda_to, by = lambda_by)
  
  Z_grid <- matrix(0, nrow = length(lambda_1), ncol = length(lambda_2))
  
  for (i in seq_along(lambda_1)) {
    for (j in seq_along(lambda_2)) {
      Z_grid[i, j] <- if (isTRUE(bootstrap)) {
        phi(lambda = c(lambda_1[i], lambda_2[j]), 
            learning_samples = learning_samples, 
            observed_values = observed_values, X_tilde = X_tilde, gam = gam, 
            alpha = alpha, loc = loc, sca = sca)$total
      } else {
        phi_variant(lambda = c(lambda_1[i], lambda_2[j]), 
                    learning_samples = learning_samples, 
                    observed_values = observed_values, gam = gam, 
                    alpha = alpha, loc = loc, sca = sca)$total
      }
    }
  }
  
  return(list(
    lambda_1 = lambda_1, lambda_2 = lambda_2, Z_grid = Z_grid
  ))
}

# create 2d surface plot of the objective function for figure 2
plot_phi_with_gd_plotly <- function(lambda_1, lambda_2, Z_grid, 
                                    show_true_min = TRUE, show_contours = TRUE, 
                                    contour_line_freq = 60L) {
  nx <- length(lambda_1)
  ny <- length(lambda_2)
  dz <- dim(Z_grid)
  
  if (!setequal(dz, c(nx, ny))) {
    stop(sprintf(
      "Dimension mismatch: dim(Z_grid) = (%d,%d) but length(lambda_1)=%d, length(lambda_2)=%d", 
      dz[1], dz[2], nx, ny
    ))
  }
  
  contour_spec <- if (isTRUE(show_contours)) {
    rng <- range(Z_grid)
    stepp <- diff(rng) / contour_line_freq
    list(
      z = list(
        show = TRUE, 
        coloring = "lines", 
        start = rng[1], 
        end = rng[2], 
        size = stepp, 
        width = 1.2
      )
    )
  } else {
    list()
  }
  
  p <- plotly::plot_ly() |>
    plotly::add_surface(
      x = lambda_1, y = lambda_2, z = Z_grid, colorscale = "Viridis", 
      showscale = FALSE, opacity = 0.8, contours = contour_spec
    )
  
  if (isTRUE(show_true_min)) {
    idx <- arrayInd(which.min(Z_grid), dim(Z_grid))
    r_ <- idx[1]
    c_ <- idx[2]
    
    iy <- if (dz[1] == ny && dz[2] == nx) {
      r_
    } else {
      c_
    }
    
    ix <- if (dz[1] == ny && dz[2] == nx) {
      c_
    } else {
      r_
    }
    
    lam_true <- c(lambda_1[ix], lambda_2[iy])
    z_true <- if (dz[1] == ny && dz[2] == nx) {
      Z_grid[iy, ix]
    } else {
      Z_grid[ix, iy]
    }
    
    p <- p |>
      plotly::add_markers(
        x = lam_true[1], y = lam_true[2], z = z_true + 1e-6*diff(range(Z_grid)), 
        marker = list(size = 5, color = "red", line = list(width = 2, color = "white")), 
        name = "", showlegend = FALSE
      )
  }
  
  p |> 
    plotly::layout(
      scene = list(
        xaxis = list(title = list(text = "λ₂")), 
        yaxis = list(title = list(text = "λ₁")), 
        zaxis = list(title = list(text = "Φ(λ)")), 
        aspectmode = "cube"
      )
    )
}

################################################################################
################################################################################

# pre-experiments

# A1) Brown-Resnick Process

M_br_01 <- readRDS("Data/pre_experiments/1/M_br_1.3.rds")
M_br_02 <- readRDS("Data/pre_experiments/1/M_br_1.6.rds")
M_br_03 <- readRDS("Data/pre_experiments/1/M_br_1.7.rds")


# A2) Smith Process

M_sm_01 <- readRDS("Data/pre_experiments/1/M_sm_1.3.rds")
M_sm_02 <- readRDS("Data/pre_experiments/1/M_sm_1.6.rds")
M_sm_03 <- readRDS("Data/pre_experiments/1/M_sm_1.7.rds")


# A3) Extremal Gaussian Process

M_eg_01 <- readRDS("Data/pre_experiments/1/M_eg_1.3.rds")
M_eg_02 <- readRDS("Data/pre_experiments/1/M_eg_1.6.rds")
M_eg_03 <- readRDS("Data/pre_experiments/1/M_eg_1.7.rds")



# plot functions
Z_grid_br <- make_grid(Z = M_br_02[1, ], t_end = 202, FH = 1, len = 2, 
                       gam = 100, lambda_from = 0.01, lambda_to = 2.0, 
                       lambda_by = 0.01, overlap = 0, bootstrap = FALSE, 
                       X_tilde = NULL)

plot_phi_with_gd_plotly(lambda_1 = Z_grid_br$lambda_1, 
                        lambda_2 = Z_grid_br$lambda_2, 
                        Z_grid = Z_grid_br$Z_grid, 
                        show_true_min = TRUE, 
                        show_contours = TRUE, 
                        contour_line_freq = 60L)

################################################################################
################################################################################

n_runs <- 10
best_phis <- numeric(n_runs)
best_lambdas <- matrix(0, n_runs, 2)

for (i in seq_len(n_runs)) {
  set.seed(1711 + i)
  result_br <- run_extrapolation(Z = M_br_02[1, ], t_end = 202, FH_end = 1, 
                                 len = 2, overlap = 0, gam = 100, plot_ = FALSE, 
                                 max_iter = 10000L, stag_pat = 1000L, bs = 1L, 
                                 eta = 0.1, plot_conv = TRUE, bootstrap = FALSE, 
                                 optimizer = "adam", reparam = TRUE, 
                                 eval_every = 1L, full_grad = FALSE)
  best_phis[i] <- result_br$best_phis
  best_lambdas[i, ] <- result_br$lambdas
}

mean(best_phis)
sd(best_phis)
min(best_phis)
which.min(best_phis)

set.seed(1711 + which.min(best_phis))
result_br <- run_extrapolation(Z = M_br_02[1, ], t_end = 202, FH_end = 1, 
                               len = 2, overlap = 0, gam = 100, plot_ = FALSE, 
                               max_iter = 10000L, stag_pat = 1000L, bs = 1L, 
                               eta = 0.1, plot_conv = TRUE, bootstrap = FALSE, 
                               optimizer = "adam", reparam = TRUE, 
                               eval_every = 1L, full_grad = FALSE)

result_br$best_phis
result_br$lambdas

################################################################################

n_runs <- 10
best_phis <- numeric(n_runs)
best_lambdas <- matrix(0, n_runs, 2)

for (i in seq_len(n_runs)) {
  set.seed(1171 + i)
  result_br <- run_extrapolation(Z = M_br_02[1, ], t_end = 202, FH_end = 1, 
                                 len = 2, overlap = 0, gam = 100, plot_ = FALSE, 
                                 max_iter = 10000L, stag_pat = 1000L, bs = 1L, 
                                 eta = 0.01, plot_conv = TRUE, bootstrap = FALSE, 
                                 optimizer = "sgd", reparam = TRUE, 
                                 eval_every = 1L, full_grad = FALSE)
  best_phis[i] <- result_br$best_phis
  best_lambdas[i, ] <- result_br$lambdas
}

mean(best_phis)
sd(best_phis)
min(best_phis)
which.min(best_phis)

set.seed(1171 + which.min(best_phis))
result_br <- run_extrapolation(Z = M_br_02[1, ], t_end = 202, FH_end = 1, 
                               len = 2, overlap = 0, gam = 100, plot_ = FALSE, 
                               max_iter = 10000L, stag_pat = 1000L, bs = 1L, 
                               eta = 0.01, plot_conv = TRUE, bootstrap = FALSE, 
                               optimizer = "sgd", reparam = TRUE, 
                               eval_every = 1L, full_grad = FALSE)

result_br$best_phis
result_br$lambdas

################################################################################

result_br <- run_extrapolation(Z = M_br_02[1, ], t_end = 202, FH_end = 1, 
                               len = 2, overlap = 0, gam = 100, plot_ = FALSE, 
                               max_iter = 10000L, stag_pat = 1000L, bs = 1L, 
                               eta = 0.1, plot_conv = TRUE, bootstrap = FALSE, 
                               optimizer = "nlminb", reparam = FALSE, 
                               eval_every = 1L, full_grad = FALSE)

result_br$best_phis
result_br$lambdas

################################################################################
################################################################################

# function to carry out "pre-experiment", i.e. make 1000 1-step predictions and 
# compute empirical versions of the excursion metric and the mean squared error
# to find "best" penalty weight gamma
pre_experiment <- function(read_in_path, out_dir, gamma_vals = seq(2, 20, by = 2)) {
  M_test <- readRDS(read_in_path)
  numb <- nrow(M_test)
  len <- 2L
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  X_hats_all <- matrix(NA_real_, length(gamma_vals), numb)
  em_curves <- matrix(NA_real_, length(gamma_vals), numb)
  mse_curves <- matrix(NA_real_, length(gamma_vals), numb)
  
  lambdas_all <- array(NA_real_, dim = c(length(gamma_vals), numb, len))
  
  truth <- M_test[, 203]
  
  for (k in seq_along(gamma_vals)) {
    for (i in seq_len(numb)) {
      test_result <- extrapolation(Z = M_test[i, ], t_end = 202, len = len, 
                                   overlap = 0, gam = gamma_vals[k], FH = 1, 
                                   max_iter = 10000L, stag_pat = 1000L, 
                                   bootstrap = FALSE, eta = 0.1, bs = 1L, 
                                   plot_conv = FALSE, reparam = FALSE, 
                                   eval_every = 1L, optimizer = "nlminb", 
                                   full_grad = FALSE)
      X_hats_all[k, i] <- test_result$X_hat
      em_curves[k, i] <- excursion_metric(X = truth[1:i], X_hat = X_hats_all[k, 1:i])
      mse_curves[k, i] <- mse(X_hat = X_hats_all[k, 1:i])
      
      lambdas_all[k, i, ] <- test_result$lambda_end
    }
    message(sprintf("[DONE] gamma=%d", gamma_vals[k]))
  }
  write.csv(X_hats_all, file = paste0(out_dir, "/X_hats_all.csv"), row.names = FALSE)
  write.csv(em_curves, file = paste0(out_dir, "/em_curves.csv"), row.names = FALSE)
  write.csv(mse_curves, file = paste0(out_dir, "/mse_curves.csv"), row.names = FALSE)
  
  lam_df <- data.frame(
    gamma_idx = rep(seq_along(gamma_vals), each = numb), 
    gamma = rep(gamma_vals, each = numb), 
    i = rep(seq_len(numb), times = length(gamma_vals)), 
    lambda1 = as.vector(lambdas_all[, , 1]), 
    lambda2 = as.vector(lambdas_all[, , 2])
  )
  
  write.csv(lam_df, file = paste0(out_dir, "/lambdas_all_long.csv"), row.names = FALSE)
}

################################################################################
################################################################################

# pre-experiments #1 (1000 processes)

pre_experiment(read_in_path = "Data/pre_experiments/1/M_br_1.1.rds", 
               out_dir = "Results/pre_experiments/1/br_01", 
               gamma_vals = seq(0, 200, by = 1))

pre_experiment(read_in_path = "Data/pre_experiments/1/M_br_1.2.rds", 
               out_dir = "Results/pre_experiments/1/br_02", 
               gamma_vals = seq(0, 200, by = 1))

pre_experiment(read_in_path = "Data/pre_experiments/1/M_br_1.3.rds", 
               out_dir = "Results/pre_experiments/1/br_03", 
               gamma_vals = seq(0, 200, by = 1))

pre_experiment(read_in_path = "Data/pre_experiments/1/M_br_1.4.rds", 
               out_dir = "Results/pre_experiments/1/br_04", 
               gamma_vals = seq(0, 200, by = 1))

pre_experiment(read_in_path = "Data/pre_experiments/1/M_br_1.5.rds", 
               out_dir = "Results/pre_experiments/1/br_05", 
               gamma_vals = seq(0, 200, by = 1))

pre_experiment(read_in_path = "Data/pre_experiments/1/M_br_1.6.rds", 
               out_dir = "Results/pre_experiments/1/br_06", 
               gamma_vals = seq(0, 200, by = 1))

pre_experiment(read_in_path = "Data/pre_experiments/1/M_br_1.7.rds", 
               out_dir = "Results/pre_experiments/1/br_07", 
               gamma_vals = seq(0, 200, by = 1))

pre_experiment(read_in_path = "Data/pre_experiments/1/M_br_1.8.rds", 
               out_dir = "Results/pre_experiments/1/br_08", 
               gamma_vals = seq(0, 200, by = 1))

pre_experiment(read_in_path = "Data/pre_experiments/1/M_br_1.9.rds", 
               out_dir = "Results/pre_experiments/1/br_09", 
               gamma_vals = seq(0, 200, by = 1))

################################################################################

pre_experiment(read_in_path = "Data/pre_experiments/1/M_sm_1.1.rds", 
               out_dir = "Results/pre_experiments/1/sm_01", 
               gamma_vals = seq(0, 200, by = 1))

pre_experiment(read_in_path = "Data/pre_experiments/1/M_sm_1.2.rds", 
               out_dir = "Results/pre_experiments/1/sm_02", 
               gamma_vals = seq(0, 200, by = 1))

pre_experiment(read_in_path = "Data/pre_experiments/1/M_sm_1.3.rds", 
               out_dir = "Results/pre_experiments/1/sm_03", 
               gamma_vals = seq(0, 200, by = 1))

pre_experiment(read_in_path = "Data/pre_experiments/1/M_sm_1.4.rds", 
               out_dir = "Results/pre_experiments/1/sm_04", 
               gamma_vals = seq(0, 200, by = 1))

pre_experiment(read_in_path = "Data/pre_experiments/1/M_sm_1.5.rds", 
               out_dir = "Results/pre_experiments/1/sm_05", 
               gamma_vals = seq(0, 200, by = 1))

pre_experiment(read_in_path = "Data/pre_experiments/1/M_sm_1.6.rds", 
               out_dir = "Results/pre_experiments/1/sm_06", 
               gamma_vals = seq(0, 200, by = 1))

pre_experiment(read_in_path = "Data/pre_experiments/1/M_sm_1.7.rds", 
               out_dir = "Results/pre_experiments/1/sm_07", 
               gamma_vals = seq(0, 200, by = 1))

pre_experiment(read_in_path = "Data/pre_experiments/1/M_sm_1.8.rds", 
               out_dir = "Results/pre_experiments/1/sm_08", 
               gamma_vals = seq(0, 200, by = 1))

pre_experiment(read_in_path = "Data/pre_experiments/1/M_sm_1.9.rds", 
               out_dir = "Results/pre_experiments/1/sm_09", 
               gamma_vals = seq(0, 200, by = 1))

################################################################################

pre_experiment(read_in_path = "Data/pre_experiments/1/M_eg_1.1.rds", 
               out_dir = "Results/pre_experiments/1/eg_01", 
               gamma_vals = seq(0, 200, by = 1))

pre_experiment(read_in_path = "Data/pre_experiments/1/M_eg_1.2.rds", 
               out_dir = "Results/pre_experiments/1/eg_02", 
               gamma_vals = seq(0, 200, by = 1))

pre_experiment(read_in_path = "Data/pre_experiments/1/M_eg_1.3.rds", 
               out_dir = "Results/pre_experiments/1/eg_03", 
               gamma_vals = seq(0, 200, by = 1))

pre_experiment(read_in_path = "Data/pre_experiments/1/M_eg_1.4.rds", 
               out_dir = "Results/pre_experiments/1/eg_04", 
               gamma_vals = seq(0, 200, by = 1))

pre_experiment(read_in_path = "Data/pre_experiments/1/M_eg_1.5.rds", 
               out_dir = "Results/pre_experiments/1/eg_05", 
               gamma_vals = seq(0, 200, by = 1))

pre_experiment(read_in_path = "Data/pre_experiments/1/M_eg_1.6.rds", 
               out_dir = "Results/pre_experiments/1/eg_06", 
               gamma_vals = seq(0, 200, by = 1))

pre_experiment(read_in_path = "Data/pre_experiments/1/M_eg_1.7.rds", 
               out_dir = "Results/pre_experiments/1/eg_07", 
               gamma_vals = seq(0, 200, by = 1))

################################################################################
################################################################################

# pre-experiments #2 (1000 processes, different seed)

pre_experiment(read_in_path = "Data/pre_experiments/2/M_br_1.1.rds", 
               out_dir = "Results/pre_experiments/2/br_01", 
               gamma_vals = seq(0, 200, by = 1))

pre_experiment(read_in_path = "Data/pre_experiments/2/M_br_1.2.rds", 
               out_dir = "Results/pre_experiments/2/br_02", 
               gamma_vals = seq(0, 200, by = 1))

pre_experiment(read_in_path = "Data/pre_experiments/2/M_br_1.3.rds", 
               out_dir = "Results/pre_experiments/2/br_03", 
               gamma_vals = seq(0, 200, by = 1))

pre_experiment(read_in_path = "Data/pre_experiments/2/M_br_1.4.rds", 
               out_dir = "Results/pre_experiments/2/br_04", 
               gamma_vals = seq(0, 200, by = 1))

pre_experiment(read_in_path = "Data/pre_experiments/2/M_br_1.5.rds", 
               out_dir = "Results/pre_experiments/2/br_05", 
               gamma_vals = seq(0, 200, by = 1))

pre_experiment(read_in_path = "Data/pre_experiments/2/M_br_1.6.rds", 
               out_dir = "Results/pre_experiments/2/br_06", 
               gamma_vals = seq(0, 200, by = 1))

pre_experiment(read_in_path = "Data/pre_experiments/2/M_br_1.7.rds", 
               out_dir = "Results/pre_experiments/2/br_07", 
               gamma_vals = seq(0, 200, by = 1))

pre_experiment(read_in_path = "Data/pre_experiments/2/M_br_1.8.rds", 
               out_dir = "Results/pre_experiments/2/br_08", 
               gamma_vals = seq(0, 200, by = 1))

pre_experiment(read_in_path = "Data/pre_experiments/2/M_br_1.9.rds", 
               out_dir = "Results/pre_experiments/2/br_09", 
               gamma_vals = seq(0, 200, by = 1))

################################################################################

pre_experiment(read_in_path = "Data/pre_experiments/2/M_sm_1.1.rds", 
               out_dir = "Results/pre_experiments/2/sm_01", 
               gamma_vals = seq(0, 200, by = 1))

pre_experiment(read_in_path = "Data/pre_experiments/2/M_sm_1.2.rds", 
               out_dir = "Results/pre_experiments/2/sm_02", 
               gamma_vals = seq(0, 200, by = 1))

pre_experiment(read_in_path = "Data/pre_experiments/2/M_sm_1.3.rds", 
               out_dir = "Results/pre_experiments/2/sm_03", 
               gamma_vals = seq(0, 200, by = 1))

pre_experiment(read_in_path = "Data/pre_experiments/2/M_sm_1.4.rds", 
               out_dir = "Results/pre_experiments/2/sm_04", 
               gamma_vals = seq(0, 200, by = 1))

pre_experiment(read_in_path = "Data/pre_experiments/2/M_sm_1.5.rds", 
               out_dir = "Results/pre_experiments/2/sm_05", 
               gamma_vals = seq(0, 200, by = 1))

pre_experiment(read_in_path = "Data/pre_experiments/2/M_sm_1.6.rds", 
               out_dir = "Results/pre_experiments/2/sm_06", 
               gamma_vals = seq(0, 200, by = 1))

pre_experiment(read_in_path = "Data/pre_experiments/2/M_sm_1.7.rds", 
               out_dir = "Results/pre_experiments/2/sm_07", 
               gamma_vals = seq(0, 200, by = 1))

pre_experiment(read_in_path = "Data/pre_experiments/2/M_sm_1.8.rds", 
               out_dir = "Results/pre_experiments/2/sm_08", 
               gamma_vals = seq(0, 200, by = 1))

pre_experiment(read_in_path = "Data/pre_experiments/2/M_sm_1.9.rds", 
               out_dir = "Results/pre_experiments/2/sm_09", 
               gamma_vals = seq(0, 200, by = 1))

################################################################################

pre_experiment(read_in_path = "Data/pre_experiments/2/M_eg_1.1.rds", 
               out_dir = "Results/pre_experiments/2/eg_01", 
               gamma_vals = seq(0, 200, by = 1))

pre_experiment(read_in_path = "Data/pre_experiments/2/M_eg_1.2.rds", 
               out_dir = "Results/pre_experiments/2/eg_02", 
               gamma_vals = seq(0, 200, by = 1))

pre_experiment(read_in_path = "Data/pre_experiments/2/M_eg_1.3.rds", 
               out_dir = "Results/pre_experiments/2/eg_03", 
               gamma_vals = seq(0, 200, by = 1))

pre_experiment(read_in_path = "Data/pre_experiments/2/M_eg_1.4.rds", 
               out_dir = "Results/pre_experiments/2/eg_04", 
               gamma_vals = seq(0, 200, by = 1))

pre_experiment(read_in_path = "Data/pre_experiments/2/M_eg_1.5.rds", 
               out_dir = "Results/pre_experiments/2/eg_05", 
               gamma_vals = seq(0, 200, by = 1))

pre_experiment(read_in_path = "Data/pre_experiments/2/M_eg_1.6.rds", 
               out_dir = "Results/pre_experiments/2/eg_06", 
               gamma_vals = seq(0, 200, by = 1))

pre_experiment(read_in_path = "Data/pre_experiments/2/M_eg_1.7.rds", 
               out_dir = "Results/pre_experiments/2/eg_07", 
               gamma_vals = seq(0, 200, by = 1))

################################################################################
################################################################################

# pre-experiments #3 (10000 processes)

pre_experiment(read_in_path = "Data/pre_experiments/M_br_1.1.rds", 
               out_dir = "Results/pre_experiments/br_01", 
               gamma_vals = seq(0, 200, by = 10))

pre_experiment(read_in_path = "Data/pre_experiments/M_br_1.2.rds", 
               out_dir = "Results/pre_experiments/br_02", 
               gamma_vals = seq(0, 200, by = 10))

pre_experiment(read_in_path = "Data/pre_experiments/M_br_1.3.rds", 
               out_dir = "Results/pre_experiments/br_03", 
               gamma_vals = seq(0, 200, by = 10))

pre_experiment(read_in_path = "Data/pre_experiments/M_br_1.4.rds", 
               out_dir = "Results/pre_experiments/br_04", 
               gamma_vals = seq(0, 200, by = 10))

pre_experiment(read_in_path = "Data/pre_experiments/M_br_1.5.rds", 
               out_dir = "Results/pre_experiments/br_05", 
               gamma_vals = seq(0, 200, by = 10))

pre_experiment(read_in_path = "Data/pre_experiments/M_br_1.6.rds", 
               out_dir = "Results/pre_experiments/br_06", 
               gamma_vals = seq(0, 200, by = 10))

pre_experiment(read_in_path = "Data/pre_experiments/M_br_1.7.rds", 
               out_dir = "Results/pre_experiments/br_07", 
               gamma_vals = seq(0, 200, by = 10))

pre_experiment(read_in_path = "Data/pre_experiments/M_br_1.8.rds", 
               out_dir = "Results/pre_experiments/br_08", 
               gamma_vals = seq(0, 200, by = 10))

pre_experiment(read_in_path = "Data/pre_experiments/M_br_1.9.rds", 
               out_dir = "Results/pre_experiments/br_09", 
               gamma_vals = seq(0, 200, by = 10))

################################################################################

pre_experiment(read_in_path = "Data/pre_experiments/M_sm_1.1.rds", 
               out_dir = "Results/pre_experiments/sm_01", 
               gamma_vals = seq(0, 200, by = 10))

pre_experiment(read_in_path = "Data/pre_experiments/M_sm_1.2.rds", 
               out_dir = "Results/pre_experiments/sm_02", 
               gamma_vals = seq(0, 200, by = 10))

pre_experiment(read_in_path = "Data/pre_experiments/M_sm_1.3.rds", 
               out_dir = "Results/pre_experiments/sm_03", 
               gamma_vals = seq(0, 200, by = 10))

pre_experiment(read_in_path = "Data/pre_experiments/M_sm_1.4.rds", 
               out_dir = "Results/pre_experiments/sm_04", 
               gamma_vals = seq(0, 200, by = 10))

pre_experiment(read_in_path = "Data/pre_experiments/M_sm_1.5.rds", 
               out_dir = "Results/pre_experiments/sm_05", 
               gamma_vals = seq(0, 200, by = 10))

pre_experiment(read_in_path = "Data/pre_experiments/M_sm_1.6.rds", 
               out_dir = "Results/pre_experiments/sm_06", 
               gamma_vals = seq(0, 200, by = 10))

pre_experiment(read_in_path = "Data/pre_experiments/M_sm_1.7.rds", 
               out_dir = "Results/pre_experiments/sm_07", 
               gamma_vals = seq(0, 200, by = 10))

pre_experiment(read_in_path = "Data/pre_experiments/M_sm_1.8.rds", 
               out_dir = "Results/pre_experiments/sm_08", 
               gamma_vals = seq(0, 200, by = 10))

pre_experiment(read_in_path = "Data/pre_experiments/M_sm_1.9.rds", 
               out_dir = "Results/pre_experiments/sm_09", 
               gamma_vals = seq(0, 200, by = 10))

################################################################################

pre_experiment(read_in_path = "Data/pre_experiments/M_eg_1.1.rds", 
               out_dir = "Results/pre_experiments/eg_01", 
               gamma_vals = seq(0, 200, by = 10))

pre_experiment(read_in_path = "Data/pre_experiments/M_eg_1.2.rds", 
               out_dir = "Results/pre_experiments/eg_02", 
               gamma_vals = seq(0, 200, by = 10))

pre_experiment(read_in_path = "Data/pre_experiments/M_eg_1.3.rds", 
               out_dir = "Results/pre_experiments/eg_03", 
               gamma_vals = seq(0, 200, by = 10))

pre_experiment(read_in_path = "Data/pre_experiments/M_eg_1.4.rds", 
               out_dir = "Results/pre_experiments/eg_04", 
               gamma_vals = seq(0, 200, by = 10))

pre_experiment(read_in_path = "Data/pre_experiments/M_eg_1.5.rds", 
               out_dir = "Results/pre_experiments/eg_05", 
               gamma_vals = seq(0, 200, by = 10))

pre_experiment(read_in_path = "Data/pre_experiments/M_eg_1.6.rds", 
               out_dir = "Results/pre_experiments/eg_06", 
               gamma_vals = seq(0, 200, by = 10))

pre_experiment(read_in_path = "Data/pre_experiments/M_eg_1.7.rds", 
               out_dir = "Results/pre_experiments/eg_07", 
               gamma_vals = seq(0, 200, by = 10))

################################################################################
################################################################################

# pre-experiments #3 again but with finer grid search

pre_experiment(read_in_path = "Data/pre_experiments/M_br_1.1.rds", 
               out_dir = "Results/pre_experiments/br_01/new", 
               gamma_vals = seq(0, 20, 1))

pre_experiment(read_in_path = "Data/pre_experiments/M_br_1.2.rds", 
               out_dir = "Results/pre_experiments/br_02/new", 
               gamma_vals = seq(0, 20, 1))

pre_experiment(read_in_path = "Data/pre_experiments/M_br_1.3.rds", 
               out_dir = "Results/pre_experiments/br_03/new", 
               gamma_vals = seq(0, 20, 1))

pre_experiment(read_in_path = "Data/pre_experiments/M_br_1.4.rds", 
               out_dir = "Results/pre_experiments/br_04/new", 
               gamma_vals = seq(0, 20, 1))

pre_experiment(read_in_path = "Data/pre_experiments/M_br_1.5.rds", 
               out_dir = "Results/pre_experiments/br_05/new", 
               gamma_vals = seq(0, 20, 1))

pre_experiment(read_in_path = "Data/pre_experiments/M_br_1.6.rds", 
               out_dir = "Results/pre_experiments/br_06/new", 
               gamma_vals = seq(0, 20, 1))

pre_experiment(read_in_path = "Data/pre_experiments/M_br_1.7.rds", 
               out_dir = "Results/pre_experiments/br_07/new", 
               gamma_vals = seq(0, 20, 1))

pre_experiment(read_in_path = "Data/pre_experiments/M_br_1.8.rds", 
               out_dir = "Results/pre_experiments/br_08/new", 
               gamma_vals = seq(0, 20, 1))

pre_experiment(read_in_path = "Data/pre_experiments/M_br_1.9.rds", 
               out_dir = "Results/pre_experiments/br_09/new", 
               gamma_vals = seq(0, 20, 1))

################################################################################

pre_experiment(read_in_path = "Data/pre_experiments/M_sm_1.1.rds", 
               out_dir = "Results/pre_experiments/sm_01/new", 
               gamma_vals = seq(0, 20, 1))

pre_experiment(read_in_path = "Data/pre_experiments/M_sm_1.2.rds", 
               out_dir = "Results/pre_experiments/sm_02/new", 
               gamma_vals = seq(0, 20, 1))

pre_experiment(read_in_path = "Data/pre_experiments/M_sm_1.3.rds", 
               out_dir = "Results/pre_experiments/sm_03/new", 
               gamma_vals = seq(0, 20, 1))

pre_experiment(read_in_path = "Data/pre_experiments/M_sm_1.4.rds", 
               out_dir = "Results/pre_experiments/sm_04/new", 
               gamma_vals = seq(0, 20, 1))

pre_experiment(read_in_path = "Data/pre_experiments/M_sm_1.5.rds", 
               out_dir = "Results/pre_experiments/sm_05/new", 
               gamma_vals = seq(0, 20, 1))

pre_experiment(read_in_path = "Data/pre_experiments/M_sm_1.6.rds", 
               out_dir = "Results/pre_experiments/sm_06/new", 
               gamma_vals = seq(0, 20, 1))

pre_experiment(read_in_path = "Data/pre_experiments/M_sm_1.7.rds", 
               out_dir = "Results/pre_experiments/sm_07/new", 
               gamma_vals = seq(0, 20, 1))

pre_experiment(read_in_path = "Data/pre_experiments/M_sm_1.8.rds", 
               out_dir = "Results/pre_experiments/sm_08/new", 
               gamma_vals = seq(0, 20, 1))

pre_experiment(read_in_path = "Data/pre_experiments/M_sm_1.9.rds", 
               out_dir = "Results/pre_experiments/sm_09/new", 
               gamma_vals = seq(0, 20, 1))

################################################################################

pre_experiment(read_in_path = "Data/pre_experiments/M_eg_1.1.rds", 
               out_dir = "Results/pre_experiments/eg_01/new", 
               gamma_vals = seq(0, 20, 1))

pre_experiment(read_in_path = "Data/pre_experiments/M_eg_1.2.rds", 
               out_dir = "Results/pre_experiments/eg_02/new", 
               gamma_vals = seq(0, 20, 1))

pre_experiment(read_in_path = "Data/pre_experiments/M_eg_1.3.rds", 
               out_dir = "Results/pre_experiments/eg_03/new", 
               gamma_vals = seq(0, 20, 1))

pre_experiment(read_in_path = "Data/pre_experiments/M_eg_1.4.rds", 
               out_dir = "Results/pre_experiments/eg_04/new", 
               gamma_vals = seq(0, 20, 1))

pre_experiment(read_in_path = "Data/pre_experiments/M_eg_1.5.rds", 
               out_dir = "Results/pre_experiments/eg_05/new", 
               gamma_vals = seq(0, 20, 1))

pre_experiment(read_in_path = "Data/pre_experiments/M_eg_1.6.rds", 
               out_dir = "Results/pre_experiments/eg_06/new", 
               gamma_vals = seq(0, 20, 1))

pre_experiment(read_in_path = "Data/pre_experiments/M_eg_1.7.rds", 
               out_dir = "Results/pre_experiments/eg_07/new", 
               gamma_vals = seq(0, 20, 1))

################################################################################
################################################################################

# evaluate results of the pre-experiments

gamma_vals <- seq(0, 200, by = 10)
gamma_vals <- seq(0, 20, by = 1)


read_in_results <- function(sub_dir, gamma_vals) {
  em_curves <- read.csv(file = paste0("Results/pre_experiments/", sub_dir, "/em_curves.csv"))
  mse_curves <- read.csv(file = paste0("Results/pre_experiments/", sub_dir, "/mse_curves.csv"))
  X_hat_curves <- read.csv(file = paste0("Results/pre_experiments/", sub_dir, "/X_hats_all.csv"))
  
  em_vec <- numeric(length(gamma_vals))
  mse_vec <- numeric(length(gamma_vals))
  
  for (k in seq_along(gamma_vals)) {
    em_vec[k] <- mean(tail(as.numeric(em_curves[k, ]), 1))
    mse_vec[k] <- mean(tail(as.numeric(mse_curves[k, ]), 1))
  }
  
  em_vec_norm <- (em_vec - min(em_vec)) / (max(em_vec) - min(em_vec))
  mse_vec_norm <- (mse_vec - min(mse_vec)) / (max(mse_vec) - min(mse_vec))
  
  return(list(em_vec = em_vec, em_vec_norm = em_vec_norm, 
              mse_vec = mse_vec, mse_vec_norm = mse_vec_norm))
}

results <- read_in_results(sub_dir = "sm_07", gamma_vals = gamma_vals)

em_vec <- results$em_vec
em_vec_norm <- results$em_vec_norm

mse_vec <- results$mse_vec
mse_vec_norm <- results$mse_vec_norm


plot(gamma_vals, em_vec, type = "b", xlab = expression(gamma), ylab = "em", 
     main = "", lwd = 1.5)

plot(gamma_vals, mse_vec, type = "b", xlab = expression(gamma), ylab = "mse", 
     main = "", lwd = 1.5)

plot(gamma_vals, em_vec_norm, type = "b", col = "red", xlab = expression(gamma), ylab = "", lwd = 2.0)
lines(gamma_vals, mse_vec_norm, type = "b", col = "blue", lwd = 2.0)

################################################################################

# this code chunk is used to create the figure 3

yl <- range(c(em_vec_norm, mse_vec_norm))
yr_em <- range(em_vec)
yr_mse <- range(mse_vec)
to_em <- function(z) yr_em[1] + (z - yl[1]) * diff(yr_em) / diff(yl)
to_mse <- function(z) yr_mse[1] + (z - yl[1]) * diff(yr_mse) / diff(yl)
plot(gamma_vals, em_vec_norm, type = "l", col = "red", 
     xlab = expression(gamma), ylab = "", lwd = 2, ylim = yl, yaxt = "n", 
     main = expression(atop("Extremal Gaussian", theta == 1.7)))
lines(gamma_vals, mse_vec_norm, col = "blue", lwd = 2)
at_norm <- pretty(yl)
lab_left <- format(round(to_em(at_norm), 3), nsmall = 2, trim = TRUE)
lab_right <- format(round(to_mse(at_norm), 3), nsmall = 3, trim = TRUE)
axis(side = 2, at = at_norm, labels = lab_left, col = "red", col.axis = "red", 
     col.ticks = "red")
mtext("", side = 2, line = 3, col = "red")
axis(side = 4, at = at_norm, labels = lab_right, col = "blue", col.axis = "blue", 
     col.ticks = "blue")
mtext("", side = 4, line = 3, col = "blue")
################################################################################

goal <- function(em_vec_norm, mse_vec_norm, sum_ = TRUE) {
  if (isTRUE(sum_)) em_vec_norm + mse_vec_norm else pmax(em_vec_norm, mse_vec_norm)
}


goal_vec <- goal(em_vec = em_vec_norm, mse_vec = mse_vec_norm, sum_ = FALSE)
plot(gamma_vals, goal_vec, type = "b", xlab = expression(gamma), ylab = "", 
     main = "", lwd = 1.5)
gamma_min <- which.min(goal_vec)
gamma_vals[which.min(goal_vec)]
points(gamma_vals[gamma_min], goal_vec[gamma_min], col = "red", pch = 19, cex = 1.5)
legend("bottom", legend = bquote(gamma[min] == .(gamma_vals[gamma_min])), 
       bty = "n", text.col = "red")

################################################################################
################################################################################

gamma_vals <- seq(0, 200, by = 10)
gamma_vals <- seq(0, 20, by = 1)
gamma_vals <- seq(0, 200, by = 1)

# read in and create plots

out_dir <- file.path("Images", "gamma_plots", "br", "br_10000", "new")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

save_gamma_plot <- function(x, y, ylab, filename, highlight_min = FALSE, 
                            main_title = "") {
  png(filename = filename, width = 1200, height = 900, res = 160)
  on.exit(dev.off(), add = TRUE)
  
  plot(x, y, type = "b", xlab = expression(gamma), ylab = ylab, main = main_title, lwd = 1.5)
  
  if (isTRUE(highlight_min)) {
    i_min <- which.min(y)
    points(x[i_min], y[i_min], col = "red", pch = 19, cex = 1.5)
    legend("bottom", legend = bquote(gamma[min] == .(x[i_min])), 
           bty = "n", text.col = "red")
  }
}

sub_dirs <- sprintf("br_%02d/new", 1:9)

for (l in seq_along(sub_dirs)) {
  message("Processing: ", sub_dirs[l])
  
  results <- read_in_results(sub_dir = sub_dirs[l], gamma_vals = gamma_vals)
  
  em_vec <- results$em_vec
  em_vec_norm <- results$em_vec_norm
  
  mse_vec <- results$mse_vec
  mse_vec_norm <- results$mse_vec_norm
  
  suffix <- sprintf("%02d", l)
  
  # EM plot
  save_gamma_plot(gamma_vals, em_vec, ylab = "excursion metric", highlight_min = FALSE, 
                  filename = file.path(out_dir, sprintf("gamma_vs_em_%s.png", suffix)))
  
  # MSE plot
  save_gamma_plot(gamma_vals, mse_vec, ylab = "MSE", highlight_min = FALSE, 
                  filename = file.path(out_dir, sprintf("gamma_vs_mse_%s.png", suffix)))
  
  # Goal function
  goal_vec <- goal(em_vec_norm = em_vec_norm, mse_vec_norm = mse_vec_norm, sum_ = FALSE)
  
  save_gamma_plot(gamma_vals, goal_vec, ylab = "", highlight_min = TRUE, 
                  filename = file.path(out_dir, sprintf("gamma_vs_goal_%s.png", suffix)))
}



