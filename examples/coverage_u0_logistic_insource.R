# 250-rep u0 coverage using the IN-SOURCE analytic noise-propagation
# covariance now returned by estimate_IC (cov_method = "noise_propagation").
# Same logistic setup / seeds as the earlier studies (n=128, nr=0.10, 95% CI).
# We read cov_u0 (new) and cov_u0_resid (old) straight off res$boundary_state.

suppressMessages({library(deSolve); library(parallel); library(devtools)})
invisible(devtools::load_all(quiet = TRUE))

f <- function(u, p, t) c(p[1] * u[1] - p[2] * u[1]^2)
p_star  <- c(1, 1/10); u0_true <- c(0.1); npoints <- 128L
t_eval  <- seq(0, 10, length.out = npoints)
modelODE <- function(tv, s, pp) list(as.vector(f(s, pp, tv)))
sol   <- deSolve::ode(y = u0_true, times = t_eval, func = modelODE,
                      parms = p_star, rtol = 1e-12, atol = 1e-14)
U_clean <- sol[, 2, drop = FALSE]; tt <- sol[, 1, drop = FALSE]
D <- ncol(U_clean); nr <- 0.10
noise_sd <- nr * sqrt(mean(as.vector(U_clean)^2))
n_reps <- 250L; z <- qnorm(0.975)

one_rep <- function(rep_id) {
  set.seed(10000L + rep_id)
  U <- U_clean + rnorm(npoints * D, 0, noise_sd)
  out <- tryCatch({
    res <- solveWendy(f, U, tt, method = "IRLS",
                      control = list(estimate_IC = TRUE,
                                     estimate_trajectory = FALSE,
                                     test_fun_type = "SSL"))
    bs <- res$boundary_state
    list(u0hat = as.numeric(bs$u0hat),
         se_new = sqrt(diag(bs$cov_u0)),
         se_old = if (!is.null(bs$cov_u0_resid)) sqrt(diag(bs$cov_u0_resid)) else NA_real_,
         method = bs$cov_method, err = NA_character_)
  }, error = function(e) list(u0hat = NA, se_new = NA, se_old = NA,
                              method = NA, err = conditionMessage(e)))
  data.frame(rep_id = rep_id, u0hat = out$u0hat[1],
             se_new = out$se_new[1], se_old = out$se_old[1],
             method = out$method,
             cov_new = abs(out$u0hat[1] - u0_true) <= z * out$se_new[1],
             cov_old = abs(out$u0hat[1] - u0_true) <= z * out$se_old[1],
             err = out$err)
}

ncores <- max(1L, parallel::detectCores() - 1L)
cat(sprintf("Running %d reps (in-source noiseprop) on %d cores...\n", n_reps, ncores))
t0 <- Sys.time()
df <- do.call(rbind, parallel::mclapply(seq_len(n_reps), one_rep, mc.cores = ncores))
cat(sprintf("Done in %.1fs\n", as.numeric(difftime(Sys.time(), t0, units = "secs"))))

ok <- df[is.finite(df$se_new), ]
sd_emp <- sd(ok$u0hat); mc <- function(p) sqrt(p*(1-p)/nrow(ok))
cat("\n===== u0 coverage with IN-SOURCE noise-propagation cov_u0 (logistic) =====\n")
cat(sprintf("reps usable           : %d   (cov_method: %s)\n",
            nrow(ok), paste(unique(ok$method), collapse = ",")))
cat(sprintf("empirical SD(u0hat)   : %.5f\n", sd_emp))
cat(sprintf("%-24s %8s %8s %10s\n", "estimator", "mean SE", "calib", "95% cov"))
cat(sprintf("%-24s %8.4f %8.3f %9.1f%%\n", "cov_u0_resid (old)",
            mean(ok$se_old), mean(ok$se_old)/sd_emp, 100*mean(ok$cov_old)))
cat(sprintf("%-24s %8.4f %8.3f %9.1f%%   (+/- %.1f%% MC)\n", "cov_u0 (noise prop)",
            mean(ok$se_new), mean(ok$se_new)/sd_emp, 100*mean(ok$cov_new), 100*mc(mean(ok$cov_new))))
cat("nominal-vs-empirical (noise prop):\n")
for (lvl in c(0.50, 0.80, 0.90, 0.95)) {
  zc <- qnorm(0.5 + lvl/2)
  cat(sprintf("  %2.0f%% -> %5.1f%%\n", 100*lvl,
              100*mean(abs(ok$u0hat - u0_true) <= zc * ok$se_new)))
}
cat("==========================================================================\n")
saveRDS(df, "examples/coverage_u0_logistic_insource.rds")
