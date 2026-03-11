frechet_loglik <- function(theta, x) {
  # x: data vector X_1,...,X_M
  x_min <- min(x)
  alpha <- exp(theta[1])
  sig <- exp(theta[2])
  mu <- x_min - exp(theta[3])
  y <- (x - mu) / sig
  M <- length(x)
  ll <- M * log(alpha / sig) - sum(y^(-alpha)) - (alpha + 1) * sum(log(y))
  return(ll)
}


frechet_negloglik <- function(theta, x) {
  -frechet_loglik(theta = theta, x = x)
}


fit_frechet_mle <- function(x) {
  x_min <- min(x)
  # simple starting values
  alpha_start <- 1
  sig_start <- sd(x)
  mu_start <- x_min - 0.1
  theta_start <- c(log(alpha_start), log(sig_start), log(x_min - mu_start))
  opt <- optim(par = theta_start, fn = frechet_negloglik, x = x, 
               method = "BFGS", control = list(maxit = 10000L))
  alpha_hat <- exp(opt$par[1])
  sig_hat <- exp(opt$par[2])
  mu_hat <- x_min - exp(opt$par[3])
  return(
    list(
      par = c(alpha = alpha_hat, mu = mu_hat, sig = sig_hat), 
      logLik = -opt$value, 
      convergence = opt$convergence
    )
  )
}


fill_missing_years <- function(x) {
  yrs_have <- as.integer(names(x))
  yrs_full <- seq(min(yrs_have), max(yrs_have))
  out <- rep(NA_real_, length(yrs_full))
  names(out) <- yrs_full
  out[names(x)] <- x
  return(out)
}


fill_missing_months <- function(x) {
  nms <- names(x)
  d0 <- as.Date(paste0(nms, "-01"))
  
  # extend to full years: Jan of min year ... Dec of max year
  d_min <- as.Date(sprintf("%d-01-01", as.integer(format(min(d0), "%Y"))))
  d_max <- as.Date(sprintf("%d-12-01", as.integer(format(max(d0), "%Y"))))
  
  d <- seq(d_min, d_max, by = "month")
  nm <- format(d, "%Y-%m")
  out <- rep(NA_real_, length(nm))
  storage.mode(out) <- storage.mode(x)
  names(out) <- nm
  out[nms] <- x
  return(out)
}


monthly_to_annual_max <- function(x) {
  yrs <- substr(names(x), 1, 4)
  out <- tapply(x, yrs, function(z) if (all(is.na(z))) NA_real_ else max(z, na.rm = TRUE))
  return(out[order(as.integer(names(out)))])
}


read_in_annual_data <- function(file_name) {
  df_annual <- read.csv(file_name)
  emxp_annual <- df_annual$EMXP
  names(emxp_annual) <- df_annual$DATE
  emxp_annual <- fill_missing_years(emxp_annual)
  return(emxp_annual)
}


read_in_monthly_data <- function(file_name) {
  df_monthly <- read.csv(file_name)
  emxp_monthly <- df_monthly$EMXP
  names(emxp_monthly) <- df_monthly$DATE
  emxp_monthly <- fill_missing_months(emxp_monthly)
  return(emxp_monthly)
}


################################################################################
################################################################################

emxp_annual <- read_in_annual_data(file_name = "Data/rainfall_data/GM000004199.csv")
emxp_annual_2 <- read_in_annual_data(file_name = "Data/rainfall_data/GME00126262.csv")


# which years are missing?
which(is.na(emxp_annual))
which(is.na(emxp_annual_2))


emxp_monthly <- read_in_monthly_data(file_name = "Data/rainfall_data/GM000004199 Kopie.csv")
emxp_monthly_2 <- read_in_monthly_data(file_name = "Data/rainfall_data/GME00126262 Kopie.csv")


# which months are missing?
which(is.na(emxp_monthly))
which(is.na(emxp_monthly_2))




# get yearly maximums from monthly maximums
emxp_annual_from_months <- monthly_to_annual_max(emxp_monthly)
emxp_annual_from_months_2 <- monthly_to_annual_max(emxp_monthly_2)



# which years are still missing?
which(is.na(emxp_annual_from_months))
which(is.na(emxp_annual_from_months_2))

################################################################################

# fill the two missing years (1946, 1947) with observations from other weather station
names(emxp_annual_from_months)[is.na(emxp_annual_from_months)]

emxp_annual_from_months["1946"] <- emxp_annual_from_months_2["1946"]
emxp_annual_from_months["1947"] <- emxp_annual_from_months_2["1947"]

################################################################################

# Figure (8) left

yrs <- as.integer(names(emxp_annual_from_months))

plot(yrs, emxp_annual_from_months, type = "l", col = "cornflowerblue", lwd = 1.5, 
     xlab = "Year", ylab = "Annual Daily Maximum Rainfall (mm)", main = "", 
     xaxt = "n")

v20 <- seq(ceiling(min(yrs)/20)*20, floor(max(yrs)/20)*20, by = 20)
abline(v = v20, col = "lightgrey")

v5 <- seq(ceiling(min(yrs)/5)*5, floor(max(yrs)/5)*5, by = 5)
axis(1, at = v5, labels = FALSE, tck = -0.02)
axis(1, at = v20, labels = TRUE, tck = -0.05, las = 2)

################################################################################

# Figure (8) right

fit <- fit_frechet_mle(emxp_annual_from_months)
fit$convergence

alpha_hat <- fit$par["alpha"]
mu_hat <- fit$par["mu"]
sigma_hat <- fit$par["sig"]

alpha_hat
mu_hat
sigma_hat

# first plot: empirical CDF
plot(ecdf(emxp_annual_from_months), main = "", xlab = "x", ylab = "CDF", 
     do.points = FALSE, verticals = TRUE, col = "cornflowerblue", lwd = 2)

# x-range over which to plot the fitted CDF
x_grid <- seq(min(emxp_annual_from_months), max(emxp_annual_from_months), length.out = 2000)

vec <- H_cdf(x_grid, alpha = alpha_hat, loc = mu_hat, sca = sigma_hat)

# add fitted Fréchet CDF as a smooth line
lines(x_grid, vec, col = "red1")

legend("bottomright", legend = c("Empirical CDF", "Fitted Fréchet CDF"), 
       col = c("cornflowerblue", "red1"), lwd = c(2, 1))

################################################################################
################################################################################

# test bootstrap with specific gamma
n_bootstrap <- 100
forecast_length <- 3L
preds_100 <- matrix(0, n_bootstrap, forecast_length)
gam <- 2.0

set.seed(1002003)
for (i in seq_len(n_bootstrap)) {
  result_fig9 <- run_extrapolation(Z = emxp_annual_from_months, t_end = 144, 
                                   FH_end = forecast_length, len = 4, overlap = 0, 
                                   gam = gam, plot_ = TRUE, max_iter = 10000L, 
                                   stag_pat = 1000L, bs = 1L, eta = 0.1, 
                                   plot_conv = FALSE, bootstrap = TRUE, 
                                   optimizer = "adam", alpha = alpha_hat, 
                                   loc = mu_hat, sca = sigma_hat, 
                                   reparam = TRUE, eval_every = 1L)
  preds_100[i, ] <- result_fig9$preds
}


################################################################################
################################################################################

# function to create figure 9
create_fig9 <- function(gam, preds_100) {
  
  # test non-bootstrap with specific gamma
  result_fig9 <- run_extrapolation(Z = emxp_annual_from_months, t_end = 144, 
                                   FH_end = 3, len = 4, overlap = 0, gam = gam, 
                                   plot_ = TRUE, max_iter = 10000L, stag_pat = 1000L, 
                                   bs = 1L, eta = 0.1, plot_conv = FALSE, 
                                   bootstrap = FALSE, optimizer = "adam", 
                                   alpha = alpha_hat, loc = mu_hat, 
                                   sca = sigma_hat, reparam = TRUE, 
                                   eval_every = 1L)
  
  specific_prediction <- result_fig9$preds
  
  # plot
  env_max <- apply(preds_100, 2, max)
  env_min <- apply(preds_100, 2, min)
  years_pred <- 2023:2025
  first_pred_year <- years_pred[1]
  idx_prev <- which(names(emxp_annual_from_months) == first_pred_year) - 1
  year_prev <- names(emxp_annual_from_months)[idx_prev]
  x_prev <- emxp_annual_from_months[idx_prev]
  years_env <- c(year_prev, years_pred)
  env_max_e <- c(x_prev, env_max)
  env_min_e <- c(x_prev, env_min)
  spec_pred <- as.numeric(result_fig9$preds)
  spec_pred_e <- c(x_prev, spec_pred)
  x_min <- years_pred[1] - 10
  x_max <- tail(years_pred, 1)
  keep_obs <- names(emxp_annual_from_months) >= x_min & names(emxp_annual_from_months) <= x_max
  ylim <- range(emxp_annual_from_months[keep_obs], env_min_e, env_max_e, spec_pred_e)
  plot(names(emxp_annual_from_months)[keep_obs], emxp_annual_from_months[keep_obs], 
       type = "l", lwd = 2, col = "cornflowerblue", xlab = "Year", main = "", 
       ylab = "Annual Daily Maximum Rainfall (mm)", xlim = c(x_min, x_max), ylim = ylim)
  polygon(x = c(years_env, rev(years_env)), y = c(env_max_e, rev(env_min_e)), 
          col = adjustcolor("orange", alpha.f = 0.3), border = NA)
  lines(years_env, env_max_e, lwd = 2, col = "orange")
  lines(years_env, env_min_e, lwd = 2, col = "orange")
  lines(years_env, spec_pred_e, lwd = 2, col = "cornflowerblue", lty = 2)
}

set.seed(10032026)
create_fig9(gam = 2.0, preds_100 = preds_100)





