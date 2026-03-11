# cumulative distribution function of the Fréchet distribution with location 'loc', 
# scale 'sca' and shape 'alpha'
H_cdf <- function(x, loc = 0, sca = 1, alpha = 1) {
  z <- (x - loc) / sca
  out <- numeric(length(x))
  ok <- is.finite(z) & z > 0
  out[ok] <- exp(-z[ok]^(-alpha))
  return(out)
}

# probability density function of the Fréchet distribution with location 'loc', 
# scale 'sca' and shape 'alpha'
h_pdf <- function(x, loc = 0, sca = 1, alpha = 1) {
  z <- (x - loc) / sca
  out <- numeric(length(x))
  ok <- is.finite(z) & z > 0
  out[ok] <- (alpha/sca) * z[ok]^(-alpha-1) * exp(-z[ok]^(-alpha))
  return(out)
}

# function to construct 1d learning samples
get_learning_samples <- function(Z, FH, len, overlap = 0L) {
  if (FH >= len) stop("'len' must be strictly greater than 'FH'")
  T_ <- length(Z)
  stepp <- len - overlap
  if (stepp < 1L) stop("Invalid overlap")
  last_start <- T_ - 2*len + 1L
  if (last_start < 1L) stop("Time series too short for the requested len")
  starts <- seq.int(from = last_start, to = 1L, by = -stepp)
  ends <- starts + len - 1L
  targets <- ends + FH
  
  learning_samples <- t(vapply(starts, function(s) Z[s:(s+len-1L)], numeric(len)))
  observed_values <- Z[targets]
  
  return(list(learning_samples = learning_samples, observed_values = observed_values))
}


# this function implements \Phi (see equation (32)) using the non-bootstrap version (formulation (30))
phi_variant <- function(lambda, learning_samples, observed_values, gam, 
                        alpha = 1, loc = 0, sca = 1) {
  N_ls <- nrow(learning_samples)
  
  sprod <- t(t(learning_samples)*lambda)
  M_all <- matrixStats::rowMaxs(sprod)
  
  H_all <- H_cdf(M_all, alpha = alpha, loc = loc, sca = sca)
  
  H_max <- H_cdf(pmax(observed_values, M_all), alpha = alpha, loc = loc, sca = sca)
  
  r <- rank(H_all, ties.method = "first")
  tri_sum <- (r - 1) * H_all
  
  em <- mean(2*H_max - H_all) - 0.5
  ws <- mean(H_all^2 - (1/N_ls)*(H_all + 2*tri_sum)) + (1/3)
  y <- em + gam * ws
  
  return(list(total = y, em = em, ws = ws))
}


# this function implements \Phi (see equation (32)) using the bootstrap version (formulation (29))
phi <- function(lambda, learning_samples, observed_values, X_tilde, gam, 
                alpha = 1, loc = 0, sca = 1) {
  sprod <- t(t(learning_samples)*lambda)
  sprod_tilde <- t(t(X_tilde)*lambda)
  
  M_all <- matrixStats::rowMaxs(sprod)
  M_all_tilde <- matrixStats::rowMaxs(sprod_tilde)
  
  H_all <- H_cdf(M_all, alpha = alpha, loc = loc, sca = sca)
  H_all_tilde <- H_cdf(M_all_tilde, alpha = alpha, loc = loc, sca = sca)
  
  H_max <- H_cdf(pmax(observed_values, M_all), alpha = alpha, loc = loc, sca = sca)
  H_max_boot <- pmax(H_all, H_all_tilde)
  
  em <- mean(2*H_max - H_all) - 0.5
  ws <- mean(H_all^2 - H_max_boot) + (1/3)
  y <- em + gam * ws
  
  return(list(total = y, em = em, ws = ws))
}


# this function implements the gradient of Q_j (see section 8) for the non-bootstrap version (formulation (30))
get_gradient_variant <- function(j, learning_samples, observed_values, lambda, 
                                 gam, alpha = 1, loc = 0, sca = 1) {
  N_ls <- nrow(learning_samples)
  
  X_j <- observed_values[j]
  sub_sample <- learning_samples[j, ]
  sprod_j <- lambda * sub_sample
  X_hat_j <- max(sprod_j)
  
  if (abs(X_j - X_hat_j) < 1e-12) return(NULL)
  
  is_max <- abs(sprod_j - X_hat_j) < 1e-12
  
  if (sum(is_max) > 1L) return(NULL)
  
  i_star <- which(is_max)[1L]
  
  X_hat_sum_p <- numeric(length(lambda))
  X_hat_sum_q <- 0
  
  if (j > 1L) {
    S <- t(t(learning_samples[1L:(j-1L), ])*lambda)
    X_hat_m <- matrixStats::rowMaxs(S)
    if (any(abs(X_hat_m - X_hat_j) < 1e-12)) return(NULL)
    if (any(rowSums(abs(S - X_hat_m) < 1e-12) > 1L)) return(NULL)
    X_hat_sum_q <- sum(X_hat_m < X_hat_j)
    idx <- which(X_hat_m > X_hat_j)
    if (length(idx) > 0L) {
      w <- max.col(m = S[idx, ], ties.method = "first")
      xm <- learning_samples[cbind(idx, w)]
      contrib <- h_pdf(X_hat_m[idx], alpha = alpha, loc = loc, sca = sca) * xm
      rs <- rowsum(x = contrib, group = w, reorder = FALSE)
      X_hat_sum_p[as.integer(rownames(rs))] <- rs[, 1L]
    }
  }
  
  dQ <- -(2*gam/N_ls) * X_hat_sum_p
  red <- h_pdf(X_hat_j, alpha = alpha, loc = loc, sca = sca) * sub_sample[i_star]
  
  brac <- (2*(X_j < X_hat_j) - 1 + 2*gam*H_cdf(X_hat_j, alpha = alpha, loc = loc, sca = sca) - (gam/N_ls) - (2*gam/N_ls)*X_hat_sum_q)
  
  dQ[i_star] <- dQ[i_star] + red * brac
  
  return(dQ)
}


full_grad_cached <- function(learning_samples, observed_values, gam, 
                             M_all, w_all, tie_all, red_all, H_all) {
  N_ls <- length(M_all)
  g <- numeric(ncol(learning_samples))
  
  non_diff <- 0
  
  pref_tie <- logical(N_ls) # FALSE ... FALSE
  if (N_ls >= 2L) {
    pref_tie[2L:N_ls] <- cummax(tie_all[1L:(N_ls-1L)])
  }
  
  for (j in seq_len(N_ls)) {
    X_j <- observed_values[j]
    X_hat_j <- M_all[j]
    
    if (abs(X_j - X_hat_j) < 1e-12 || tie_all[j]) {
      non_diff <- non_diff + 1L
      next
    }
    if (j > 1L) {
      if (pref_tie[j]) {
        non_diff <- non_diff + 1L
        next
      }
      if (any(abs(M_all[1L:(j-1L)] - X_hat_j) < 1e-12)) {
        non_diff <- non_diff + 1L
        next
      }
    }
    
    X_hat_sum_q <- if (j == 1L) 0 else sum(M_all[1L:(j-1L)] < X_hat_j)
    
    X_hat_sum_p <- numeric(length(g))
    if (j > 1L) {
      idx <- which(M_all[1L:(j-1L)] > X_hat_j)
      if (length(idx) > 0L) {
        rs <- rowsum(x = red_all[idx], group = w_all[idx], reorder = FALSE)
        X_hat_sum_p[as.integer(rownames(rs))] <- rs[, 1L]
      }
    }
    
    dQ <- -(2*gam/N_ls) * X_hat_sum_p
    i_star <- w_all[j]
    brac <- (2*(X_j < X_hat_j) - 1 + 2*gam*H_all[j] - (gam/N_ls) - (2*gam/N_ls)*X_hat_sum_q)
    dQ[i_star] <- dQ[i_star] + red_all[j] * brac
    g <- g + dQ
  }
  
  return(g/(N_ls - non_diff))
}


# this function implements the gradient of Q_j (see section 8) for the bootstrap version (formulation (29))
get_gradient <- function(j, learning_samples, observed_values, lambda, gam, X_tilde, 
                         alpha = 1, loc = 0, sca = 1) {
  dQ <- numeric(length(lambda))
  
  X_j <- observed_values[j]
  
  sub_sample <- learning_samples[j, ]
  sub_sample_tilde <- X_tilde[j, ]
  
  sprod_j <- lambda * sub_sample
  sprod_j_tilde <- lambda * sub_sample_tilde
  
  X_hat_j <- max(sprod_j)
  X_hat_j_tilde <- max(sprod_j_tilde)
  
  is_max <- abs(sprod_j - X_hat_j) < 1e-12
  is_max_tilde <- abs(sprod_j_tilde - X_hat_j_tilde) < 1e-12
  
  if (sum(is_max) > 1L) return(NULL)
  if (sum(is_max_tilde) > 1L) return(NULL)
  
  if (abs(X_j - X_hat_j) < 1e-12) return(NULL)
  if (abs(X_hat_j - X_hat_j_tilde) < 1e-12) return(NULL)
  
  i_star <- which(is_max)[1L]
  i_star_tilde <- which(is_max_tilde)[1L]
  
  dQ[i_star] <- (2*(X_j < X_hat_j) - 1 - gam*(X_hat_j_tilde < X_hat_j) + 2*gam*H_cdf(X_hat_j, alpha = alpha, loc = loc, sca = sca)) * h_pdf(X_hat_j, alpha = alpha, loc = loc, sca = sca) * sub_sample[i_star]
  
  dQ[i_star_tilde] <- dQ[i_star_tilde] - gam*(X_hat_j_tilde > X_hat_j) * h_pdf(X_hat_j_tilde, alpha = alpha, loc = loc, sca = sca) * sub_sample_tilde[i_star_tilde]
  
  return(dQ)
}


# functions for minimization and extrapolation

# stochastic gradient descent
gradient_descent <- function(lambda_start, learning_samples, observed_values, 
                             gam, max_iter = 10000L, stag_pat = 1000L, 
                             bootstrap = FALSE, X_tilde = NULL, plot_conv = FALSE, 
                             alpha = 1, loc = 0, sca = 1, bs = 1L, eta = 1e-1, 
                             optimizer = c("sgd", "adam"), reparam = FALSE, 
                             eval_every = 1L, full_grad = FALSE) {
  N_ls <- nrow(learning_samples)
  p <- length(lambda_start)
  
  if (isTRUE(bootstrap) && is.null(X_tilde)) {
    stop("Error: for the bootstrap version 'X_tilde' must be provided.")
  }
  
  k <- 0L
  
  theta_k <- if (isTRUE(reparam)) log(pmax(lambda_start, 1e-12)) else NULL
  
  lambda_k <- if (isTRUE(reparam)) exp(theta_k) else lambda_start
  
  # Function evaluation
  phi0 <- if (isTRUE(bootstrap)) {
    phi(lambda = lambda_k, learning_samples = learning_samples, 
        observed_values = observed_values, X_tilde = X_tilde, gam = gam, 
        alpha = alpha, loc = loc, sca = sca)
  } else {
    phi_variant(lambda = lambda_k, learning_samples = learning_samples, 
                observed_values = observed_values, gam = gam, alpha = alpha, 
                loc = loc, sca = sca)
  }
  
  phi_last_eval <- phi0
  
  # initialize history for plotting
  iter_hist <- integer(max_iter + 1L)
  phi_hist <- numeric(max_iter + 1L)
  em_hist <- numeric(max_iter + 1L)
  ws_hist <- numeric(max_iter + 1L)
  used <- 0L
  
  # initialize best lambda and best phi value
  best_lambda <- lambda_k
  best_phi <- phi0$total
  best_em <- phi0$em
  best_ws <- phi0$ws
  
  # track the last iteration index where best_phi improved
  last_improve_iter <- 0L
  
  used <- used + 1L
  iter_hist[used] <- 0L
  phi_hist[used] <- phi0$total
  em_hist[used] <- phi0$em
  ws_hist[used] <- phi0$ws
  
  optimizer <- match.arg(optimizer)
  
  beta1 <- if (optimizer == "adam") 0.9 else NULL
  beta2 <- if (optimizer == "adam") 0.999 else NULL
  eps <- if (optimizer == "adam") 1e-8 else NULL
  m <- if (optimizer == "adam") rep(0, p) else NULL
  v <- if (optimizer == "adam") rep(0, p) else NULL
  t_ <- if (optimizer == "adam") 0L else NULL
  
  while (k < max_iter) {
    k <- k + 1L
    
    g <- numeric(p)
    if (!isTRUE(bootstrap) && isTRUE(full_grad)) {
      S_all <- t(t(learning_samples)*lambda_k)
      M_all <- matrixStats::rowMaxs(S_all)
      tie_all <- rowSums(abs(S_all - M_all) < 1e-12) > 1L
      w_all <- max.col(S_all, ties.method = "first")
      xm_win <- learning_samples[cbind(seq_len(N_ls), w_all)]
      red_all <- h_pdf(M_all, alpha = alpha, loc = loc, sca = sca) * xm_win
      H_all <- H_cdf(M_all, alpha = alpha, loc = loc, sca = sca)
      
      g <- full_grad_cached(learning_samples = learning_samples, 
                            observed_values = observed_values, gam = gam, 
                            M_all = M_all, w_all = w_all, tie_all = tie_all, 
                            red_all = red_all, H_all = H_all)
    } else {
      good <- 0L
      avail <- seq_len(N_ls)
      
      while(good < bs) {
        if (length(avail) == 0L) {
          if (good == 0L) stop("Error: Non-differentiable")
          warning(sprintf("Iter %d: only %d/%d differentiable grads; using smaller batch.", 
                          k, good, bs), call. = FALSE)
          break
        }
        jj <- sample(avail, size = 1L)
        avail <- avail[avail != jj]
        
        gl_ <- if (isTRUE(bootstrap)) {
          get_gradient(j = jj, learning_samples = learning_samples, 
                       observed_values = observed_values, lambda = lambda_k, 
                       gam = gam, X_tilde = X_tilde, alpha = alpha, loc = loc, 
                       sca = sca)
        } else {
          get_gradient_variant(j = jj, learning_samples = learning_samples, 
                               observed_values = observed_values, 
                               lambda = lambda_k, gam = gam, alpha = alpha, 
                               loc = loc, sca = sca)
        }
        
        if (is.null(gl_)) next
        
        g <- g + gl_
        good <- good + 1L
      }
      g <- g / good
    }
    
    if (isTRUE(reparam)) g <- g * lambda_k
    
    gg <- if (optimizer == "adam") {
      t_ <- t_ + 1L
      m <- beta1 * m + (1 - beta1) * g
      v <- beta2 * v + (1 - beta2) * (g * g)
      mhat <- m / (1 - beta1^t_)
      vhat <- v / (1 - beta2^t_)
      mhat / (sqrt(vhat) + eps)
    } else g
    
    theta_new <- if (isTRUE(reparam)) theta_k - eta * gg else NULL
    
    lambda_new <- if (isTRUE(reparam)) exp(theta_new) else lambda_k - eta * gg
    
    # Projection
    if (!isTRUE(reparam)) lambda_new <- pmax(lambda_new, 1e-12)
    
    phi_now <- if (k %% eval_every == 0L) {
      if (isTRUE(bootstrap)) {
        phi(lambda = lambda_new, learning_samples = learning_samples, 
            observed_values = observed_values, X_tilde = X_tilde, gam = gam, 
            alpha = alpha, loc = loc, sca = sca)
      } else {
        phi_variant(lambda = lambda_new, learning_samples = learning_samples, 
                    observed_values = observed_values, gam = gam, alpha = alpha, 
                    loc = loc, sca = sca)
      }
    } else phi_last_eval
    
    phi_last_eval <- phi_now
    
    # check if phi_now is smaller than best_phi
    if (phi_now$total < best_phi) {
      best_phi <- phi_now$total
      best_em <- phi_now$em
      best_ws <- phi_now$ws
      best_lambda <- lambda_new
      last_improve_iter <- k
    }
    
    # update history
    used <- used + 1L
    iter_hist[used] <- k
    phi_hist[used] <- phi_now$total
    em_hist[used] <- phi_now$em
    ws_hist[used] <- phi_now$ws
    
    # update lambda
    lambda_k <- lambda_new
    theta_k <- theta_new
    
    # stopping criterion
    if (k - last_improve_iter >= stag_pat) break
  }
  
  if (used < length(iter_hist)) {
    iter_hist <- iter_hist[1L:used]
    phi_hist <- phi_hist[1L:used]
    em_hist <- em_hist[1L:used]
    ws_hist <- ws_hist[1L:used]
  }
  
  if (isTRUE(plot_conv) && length(iter_hist) == length(phi_hist)) {
    plot_convergence(iter_hist = iter_hist, phi_hist = phi_hist, 
                     last_improve_iter = last_improve_iter, best_phi = best_phi)
  }
  
  return(list(best_lambda = best_lambda, 
              best_phi = best_phi, phi_hist = phi_hist, 
              best_em = best_em, em_hist = em_hist, 
              best_ws = best_ws, ws_hist = ws_hist, 
              last_improve_iter = last_improve_iter))
}


extrapolation <- function(Z, t_end, len, overlap, gam, FH, max_iter = 10000L, 
                          bootstrap = FALSE, eta = 0.1, bs = 1L, 
                          stag_pat = 1000L, alpha = 1, loc = 0, sca = 1, 
                          plot_conv = FALSE, reparam = FALSE, eval_every = 1L, 
                          optimizer = c("sgd", "adam", "nlminb", "optim_bfgs"), 
                          full_grad = FALSE) {
  optimizer <- match.arg(optimizer)
  
  # Initialize forecast sample
  forecast_sample <- Z[seq(t_end - len + 1, t_end)]
  
  # Get learning samples and the values corresponding to forecast horizon
  list_learning_samples <- get_learning_samples(Z = Z[seq(t_end)], FH = FH, 
                                                len = len, overlap = overlap)
  learning_samples <- list_learning_samples$learning_samples
  observed_values <- list_learning_samples$observed_values
  
  X_tilde <- NULL
  
  if (isTRUE(bootstrap)) {
    N_ls <- nrow(learning_samples)
    idx <- sample.int(N_ls, size = N_ls, replace = TRUE)
    X_tilde <- learning_samples[idx, ]
  }
  
  # Initial guess; must be same size as forecast sample
  lambda00 <- rep(1, len)
  
  theta00 <- if (isTRUE(reparam)) rep(0, len) else NULL
  
  # Define objective function for standard solvers (optim, nlminb)
  obj_box <- if (optimizer == "nlminb" || optimizer == "optim_bfgs") {
    if (isTRUE(bootstrap)) {
      if (isTRUE(reparam)) {
        function(theta) {
          phi(lambda = exp(theta), learning_samples = learning_samples, 
              observed_values = observed_values, X_tilde = X_tilde, gam = gam, 
              alpha = alpha, loc = loc, sca = sca)$total
        }
      } else {
        function(lambda) {
          phi(lambda = lambda, learning_samples = learning_samples, 
              observed_values = observed_values, X_tilde = X_tilde, gam = gam, 
              alpha = alpha, loc = loc, sca = sca)$total
        }
      }
    } else {
      if (isTRUE(reparam)) {
        function(theta) {
          phi_variant(lambda = exp(theta), learning_samples = learning_samples, 
                      observed_values = observed_values, gam = gam, 
                      alpha = alpha, loc = loc, sca = sca)$total
        }
      } else {
        function(lambda) {
          phi_variant(lambda = lambda, learning_samples = learning_samples, 
                      observed_values = observed_values, gam = gam, 
                      alpha = alpha, loc = loc, sca = sca)$total
        }
      }
    }
  } else NULL
  
  fit <- if (optimizer == "nlminb") {
    if (isTRUE(reparam)) {
      nlminb(start = theta00, objective = obj_box)
    } else {
      nlminb(start = lambda00, objective = obj_box, lower = rep(0, len))
    }
  } else if (optimizer == "optim_bfgs") {
    if (isTRUE(reparam)) {
      optim(par = theta00, fn = obj_box, method = "BFGS")
    } else {
      optim(par = lambda00, fn = obj_box, method = "L-BFGS-B", lower = rep(0, len))
    }
  } else {
    gradient_descent(lambda_start = lambda00, 
                     learning_samples = learning_samples, 
                     observed_values = observed_values, 
                     gam = gam, 
                     max_iter = max_iter, 
                     stag_pat = stag_pat, 
                     bootstrap = bootstrap, 
                     X_tilde = X_tilde, 
                     plot_conv = plot_conv, 
                     alpha = alpha, loc = loc, sca = sca, 
                     bs = bs, eta = eta, 
                     optimizer = optimizer, 
                     reparam = reparam, 
                     eval_every = eval_every, 
                     full_grad = full_grad)
  }
  
  lambda_end <- if (optimizer == "nlminb" || optimizer == "optim_bfgs") {
    if (isTRUE(reparam)) exp(fit$par) else fit$par
  } else fit$best_lambda
  
  best_phi <- if (optimizer == "nlminb") {
    fit$objective
  } else if (optimizer == "optim_bfgs") {
    fit$value
  } else {
    fit$best_phi
  }
  
  best_em <- if (optimizer == "sgd" || optimizer == "adam") {
    fit$best_em
  } else {
    if (isTRUE(bootstrap)) {
      phi(lambda = lambda_end, learning_samples = learning_samples, 
          observed_values = observed_values, X_tilde = X_tilde, gam = gam, 
          alpha = alpha, loc = loc, sca = sca)$em
    } else {
      phi_variant(lambda = lambda_end, learning_samples = learning_samples, 
                  observed_values = observed_values, gam = gam, alpha = alpha, 
                  loc = loc, sca = sca)$em
    }
  }
  
  best_ws <- if (optimizer == "sgd" || optimizer == "adam") {
    fit$best_ws
  } else {
    if (isTRUE(bootstrap)) {
      phi(lambda = lambda_end, learning_samples = learning_samples, 
          observed_values = observed_values, X_tilde = X_tilde, gam = gam, 
          alpha = alpha, loc = loc, sca = sca)$ws
    } else {
      phi_variant(lambda = lambda_end, learning_samples = learning_samples, 
                  observed_values = observed_values, gam = gam, alpha = alpha, 
                  loc = loc, sca = sca)$ws
    }
  }
  
  X_hat <- max(forecast_sample * lambda_end)
  
  return(list(X_hat = X_hat, lambda_end = lambda_end, best_phi = best_phi, 
              best_em = best_em, best_ws = best_ws))
  
}


run_extrapolation <- function(Z, t_end, FH_end, len, overlap, gam, plot_ = TRUE, 
                              max_iter = 10000L, stag_pat = 1000L, bs = 1L, 
                              eta = 1e-2, plot_conv = FALSE, bootstrap = FALSE, 
                              optimizer = "adam", alpha = 1, loc = 0, sca = 1, 
                              reparam = FALSE, eval_every = 1L, full_grad = FALSE) {
  preds <- numeric(FH_end)
  lambdas <- matrix(0, FH_end, len)
  best_phis <- numeric(FH_end)
  best_ems <- numeric(FH_end)
  best_wss <- numeric(FH_end)
  
  for (h in seq_len(FH_end)) {
    pred_list <- extrapolation(Z = Z, t_end = t_end, len = len, 
                               overlap = overlap, gam = gam, FH = h, 
                               max_iter = max_iter, bootstrap = bootstrap, 
                               eta = eta, bs = bs, stag_pat = stag_pat, 
                               alpha = alpha, loc = loc, sca = sca, 
                               plot_conv = plot_conv, reparam = reparam, 
                               eval_every = eval_every, optimizer = optimizer, 
                               full_grad = full_grad)
    preds[h] <- pred_list$X_hat
    lambdas[h, ] <- pred_list$lambda_end
    best_phis[h] <- pred_list$best_phi
    best_ems[h] <- pred_list$best_em
    best_wss[h] <- pred_list$best_ws
  }
  
  if (isTRUE(plot_)) {
    plot_prediction(Z = Z, preds = preds, t_end = t_end, FH_end = FH_end, 
                    plot_start = t_end - 2*FH_end, plot_end = t_end + FH_end, 
                    y_label = "", main_title = "", em = NULL)
  }
  
  return(list(preds = preds, lambdas = lambdas, best_phis = best_phis, 
              best_ems = best_ems, best_wss = best_wss))
}


excursion_metric <- function(X, X_hat, alpha = 1, loc = 0, sca = 1) {
  stopifnot(length(X) == length(X_hat))
  n <- length(X)
  
  E_max <- 0
  E_hat <- 0
  
  for (i in seq_len(n)) {
    E_max <- E_max + H_cdf(max(X[i], X_hat[i]), alpha = alpha, loc = loc, sca = sca)
    E_hat <- E_hat + H_cdf(X_hat[i], alpha = alpha, loc = loc, sca = sca)
  }
  
  Em <- (1/n) * (2*E_max - E_hat) - 0.5
  
  return(Em)
}

# functions for plotting
plot_convergence <- function(iter_hist, phi_hist, last_improve_iter, best_phi) {
  
  op <- par(no.readonly = TRUE); on.exit(par(op))
  
  if (all(phi_hist > 0)) {
    plot(iter_hist, phi_hist, type = "l", lwd = 1.5, xlab = "", ylab = "", log = "y")
  } else {
    plot(iter_hist, phi_hist, type = "l", lwd = 1.5, xlab = "", ylab = "")
  }
  ytex <- latex2exp::TeX("$\\Phi\\left(\\lambda^{(k)}\\right)$")
  mtext(ytex, side = 2, line = 2)
  mtext("Iterations", side = 1, line = 2)
  
  phi_start <- phi_hist[1L]
  
  points(iter_hist[1L], phi_start, pch = 16, cex = 1.5, col = "orange")
  points(iter_hist[last_improve_iter], best_phi, pch = 16, cex = 1.5, col = "red")
}


plot_prediction <- function(Z, preds, t_end, FH_end, plot_start, plot_end, 
                            y_label = "Z", main_title = NULL, em = NULL, em_star = NULL) {
  plot_range <- plot_start:plot_end
  truth_idx <- (t_end + 1L):(t_end + FH_end)
  Z_segment <- Z[plot_range]
  
  plot(plot_range, Z_segment, type = "n", xlab = "t", ylab = y_label, 
       ylim = range(c(Z_segment, preds)), main = main_title, log = "y")
  
  v10 <- plot_range[plot_range %% 10L == 0L]
  abline(v = v10, col = "grey70", lwd = 0.5)
  
  lines(plot_range, Z_segment, lwd = 1.5, col = "cornflowerblue")
  
  lines(c(t_end, t_end + 1L), c(Z[t_end], preds[1L]), lwd = 1.5, 
        col = "cornflowerblue", lty = "dashed")
  
  lines(truth_idx, preds, lwd = 1.5, col = "cornflowerblue", lty = "dashed")
  
  if (!is.null(em) || !is.null(em_star)) {
    ylim_right <- range(c(em, em_star), finite = TRUE)
    
    if (!is.null(em)) {
      par(new = TRUE)
      plot(truth_idx, em, type = "l", lwd = 3, col = "red", axes = FALSE, lty = 3, 
           xlab = "", ylab = "", ylim = ylim_right, xlim = c(plot_start, plot_end))
    }
    
    if (!is.null(em_star)) {
      par(new = TRUE)
      plot(truth_idx, em_star, type = "l", lwd = 1.5, col = "darkgreen", axes = FALSE, 
           xlab = "", ylab = "", ylim = ylim_right, xlim = c(plot_start, plot_end))
    }
    
    axis(4, col.axis = "black", col = "black")
    mtext("em (right axis)", side = 4, line = 3, col = "black")
  }
}


################################################################################
################################################################################
################################################################################

# function to run the "main experiments" (for Figures 4, 5, 6)
run_main_experiment <- function(read_in_path, write_path, gam, t_end, FH_end, len, overlap = 0) {
  M_test <- readRDS(read_in_path)
  numb <- nrow(M_test)
  X_hat_matrix <- matrix(0, numb, FH_end)
  dir.create(write_path, recursive = TRUE, showWarnings = FALSE)
  for (i in seq_len(numb)) {
    result <- run_extrapolation(Z = M_test[i, ], t_end = t_end, FH_end = FH_end, 
                                len = len, overlap = overlap, gam = gam, plot_ = TRUE, 
                                max_iter = 10000L, stag_pat = 1000L, bs = 1L, 
                                eta = 0.1, plot_conv = FALSE, bootstrap = FALSE, 
                                optimizer = "adam", reparam = TRUE, 
                                eval_every = 1L, full_grad = FALSE)
    X_hat_matrix[i, ] <- result$preds
    write.table(matrix(result$preds, nrow = 1), 
                file = paste0(write_path, "/X_hat_matrix_gamma=", as.integer(gam), ".txt"), 
                append = TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE)
  }
  
  em <- numeric(FH_end)
  
  for (l in seq_len(FH_end)) {
    em[l] <- excursion_metric(X = M_test[, (t_end+1):(t_end+FH_end)][, l], X_hat = X_hat_matrix[, l])
    cat(em[l], if (l < FH_end) " " else "\n", append = TRUE, 
        file = paste0(write_path, "/em_gamma=", as.integer(gam), ".txt"))
  }
  
  return(list(X_hat_matrix = X_hat_matrix, em = em))
}

################################################################################

# adam bs=1

run_main_experiment(read_in_path = "Data/20step_prediction/M_br_1.7.rds", 
                    write_path = "Results/results20step_prediction/results_br_07/new", 
                    gam = 1.0, t_end = 2121, FH_end = 20, len = 21)

run_main_experiment(read_in_path = "Data/20step_prediction/M_sm_1.7.rds", 
                    write_path = "Results/results20step_prediction/results_sm_07/new", 
                    gam = 1.0, t_end = 2121, FH_end = 20, len = 21)

run_main_experiment(read_in_path = "Data/20step_prediction/M_sm_1.7.rds", 
                    write_path = "Results/results20step_prediction/results_sm_07/new", 
                    gam = 2.0, t_end = 2121, FH_end = 20, len = 21)

run_main_experiment(read_in_path = "Data/20step_prediction/M_eg_1.7.rds", 
                    write_path = "Results/results20step_prediction/results_eg_07/new", 
                    gam = 1.0, t_end = 2121, FH_end = 20, len = 21)

################################################################################

# adam bs=1

run_main_experiment(read_in_path = "Data/20step_prediction/M_br_1.3.rds", 
                    write_path = "Results/results20step_prediction/results_br_03/new", 
                    gam = 0.0, t_end = 2121, FH_end = 20, len = 21)

run_main_experiment(read_in_path = "Data/20step_prediction/M_br_1.3.rds", 
                    write_path = "Results/results20step_prediction/results_br_03/new", 
                    gam = 3.0, t_end = 2121, FH_end = 20, len = 21)

run_main_experiment(read_in_path = "Data/20step_prediction/M_sm_1.3.rds", 
                    write_path = "Results/results20step_prediction/results_sm_03/new", 
                    gam = 0.0, t_end = 2121, FH_end = 20, len = 21)

run_main_experiment(read_in_path = "Data/20step_prediction/M_eg_1.3.rds", 
                    write_path = "Results/results20step_prediction/results_eg_03/new", 
                    gam = 0.0, t_end = 2121, FH_end = 20, len = 21)

################################################################################

sample_indices <- function(M_test, n) {
  numb <- nrow(M_test)
  sample.int(numb, n)
}

# function that creates the images 4, 5, 6
create_final_image <- function(read_in_path_data, 
                               read_in_path_preds, 
                               read_in_path_em, 
                               t_end = 2121, FH_end = 20, y_label = "X(t)", 
                               seed = 100, em_star = NULL) {
  M_test <- readRDS(read_in_path_data)
  numb <- nrow(M_test)
  set.seed(seed)
  j0 <- sample_indices(M_test = M_test, n = 1L)
  Z <- as.numeric(M_test[j0, ])
  X_hats <- as.matrix(read.table(read_in_path_preds, header = FALSE, sep = "", 
                                 fill = TRUE))
  preds <- as.numeric(X_hats[j0, ])
  em <- scan(read_in_path_em, quiet = TRUE)
  
  plot_prediction(Z = Z, preds = preds, t_end = t_end, FH_end = FH_end, 
                  plot_start = t_end - FH_end, plot_end = t_end + FH_end, 
                  em = em, y_label = y_label, main_title = NULL, em_star = em_star)
  return(j0)
}


# create final images

t0_vec <- c(2122:2141)

em_star_br <- 1 - 2 / (2*pnorm(0.5*br_pars["1.7"]*sqrt(t0_vec-2121)) + 1)
em_star_sm <- 1 - 2 / (2*pnorm((t0_vec-2121)/(2*sm_pars["1.7"])) + 1)
em_star_eg <- 1 - 2 / (2 + (1/sqrt(2))*sqrt(1-exp(-(t0_vec-2121)/eg_pars["1.7"])))

create_final_image(read_in_path_data = "Data/20step_prediction/M_br_1.7.rds", 
                   read_in_path_preds = "Results/results20step_prediction/results_br_07/new/X_hat_matrix_gamma=1.txt", 
                   read_in_path_em = "Results/results20step_prediction/results_br_07/new/em_gamma=1.txt", 
                   t_end = 2121, FH_end = 20, y_label = "B(t)", em_star = em_star_br)

create_final_image(read_in_path_data = "Data/20step_prediction/M_sm_1.7.rds", 
                   read_in_path_preds = "Results/results20step_prediction/results_sm_07/new/X_hat_matrix_gamma=1.txt", 
                   read_in_path_em = "Results/results20step_prediction/results_sm_07/new/em_gamma=1.txt", 
                   t_end = 2121, FH_end = 20, y_label = "S(t)", em_star = em_star_sm)

create_final_image(read_in_path_data = "Data/20step_prediction/M_sm_1.7.rds", 
                   read_in_path_preds = "Results/results20step_prediction/results_sm_07/new/X_hat_matrix_gamma=2.txt", 
                   read_in_path_em = "Results/results20step_prediction/results_sm_07/new/em_gamma=2.txt", 
                   t_end = 2121, FH_end = 20, y_label = "S(t)", em_star = em_star_sm)

create_final_image(read_in_path_data = "Data/20step_prediction/M_eg_1.7.rds", 
                   read_in_path_preds = "Results/results20step_prediction/results_eg_07/new/X_hat_matrix_gamma=1.txt", 
                   read_in_path_em = "Results/results20step_prediction/results_eg_07/new/em_gamma=1.txt", 
                   t_end = 2121, FH_end = 20, y_label = "G(t)", em_star = em_star_eg)

em_br <- as.numeric(read.table("Results/results20step_prediction/results_br_07/new/em_gamma=1.txt", header = FALSE))
em_sm <- as.numeric(read.table("Results/results20step_prediction/results_sm_07/new/em_gamma=1.txt", header = FALSE))
em_eg <- as.numeric(read.table("Results/results20step_prediction/results_eg_07/new/em_gamma=1.txt", header = FALSE))


plot(1:20, em_br, type = "l", col = "red", ylim = range(c(em_br, em_sm, em_eg)))
lines(1:20, em_sm, type = "l", col = "blue")
lines(1:20, em_eg, type = "l", col = "green")

################################################################################

t0_vec <- c(2122:2141)

em_star_br <- 1 - 2 / (2*pnorm(0.5*br_pars["1.3"]*sqrt(t0_vec-2121)) + 1)
em_star_sm <- 1 - 2 / (2*pnorm((t0_vec-2121)/(2*sm_pars["1.3"])) + 1)
em_star_eg <- 1 - 2 / (2 + (1/sqrt(2))*sqrt(1-exp(-(t0_vec-2121)/eg_pars["1.3"])))

create_final_image(read_in_path_data = "Data/20step_prediction/M_br_1.3.rds", 
                   read_in_path_preds = "Results/results20step_prediction/results_br_03/new/X_hat_matrix_gamma=0.txt", 
                   read_in_path_em = "Results/results20step_prediction/results_br_03/new/em_gamma=0.txt", 
                   t_end = 2121, FH_end = 20, y_label = "B(t)", em_star = em_star_br)

create_final_image(read_in_path_data = "Data/20step_prediction/M_br_1.3.rds", 
                   read_in_path_preds = "Results/results20step_prediction/results_br_03/new/X_hat_matrix_gamma=3.txt", 
                   read_in_path_em = "Results/results20step_prediction/results_br_03/new/em_gamma=3.txt", 
                   t_end = 2121, FH_end = 20, y_label = "B(t)", em_star = em_star_br)

create_final_image(read_in_path_data = "Data/20step_prediction/M_sm_1.3.rds", 
                   read_in_path_preds = "Results/results20step_prediction/results_sm_03/new/X_hat_matrix_gamma=0.txt", 
                   read_in_path_em = "Results/results20step_prediction/results_sm_03/new/em_gamma=0.txt", 
                   t_end = 2121, FH_end = 20, y_label = "S(t)", em_star = em_star_sm)

create_final_image(read_in_path_data = "Data/20step_prediction/M_eg_1.3.rds", 
                   read_in_path_preds = "Results/results20step_prediction/results_eg_03/new/X_hat_matrix_gamma=0.txt", 
                   read_in_path_em = "Results/results20step_prediction/results_eg_03/new/em_gamma=0.txt", 
                   t_end = 2121, FH_end = 20, y_label = "G(t)", em_star = em_star_eg)

em_br <- as.numeric(read.table("Results/results20step_prediction/results_br_03/new/em_gamma=0.txt", header = FALSE))
em_sm <- as.numeric(read.table("Results/results20step_prediction/results_sm_03/new/em_gamma=0.txt", header = FALSE))
em_eg <- as.numeric(read.table("Results/results20step_prediction/results_eg_03/new/em_gamma=0.txt", header = FALSE))


plot(1:20, em_br, type = "l", col = "red", ylim = range(c(em_br, em_sm, em_eg)))
lines(1:20, em_sm, type = "l", col = "blue")
lines(1:20, em_eg, type = "l", col = "green")




