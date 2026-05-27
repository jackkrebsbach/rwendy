# Coverage study for the closed-form LS covariance of u0 (estimate_IC).
# 250 Monte-Carlo reps on the logistic system, full solveWendy pipeline
# (p estimated each rep), additive Gaussian noise at nr = 0.10, nominal 95% CI.
#
# For each rep we form the Wald interval  u0hat +/- z * sqrt(cov_u0[d,d])
# and check whether the true u0 falls inside. Empirical coverage is the
# fraction of reps that cover. We also report the calibration ratio
# mean(model SE) / sd(u0hat): ~1 means the formula's SE matches the actual
# sampling spread of the estimator.

suppressMessages({
  library(deSolve)
  library(parallel)
  library(devtools)
})
invisible(devtools::load_all(quiet = TRUE))

## ---- problem setup (mirrors example_logistic.R) ----------------------------
f <- function(u, p, t) c(p[1] * u[1] - p[2] * u[1]^2)

p_star  <- c(1, 1/10)
u0_true <- c(0.1)
npoints <- 128L
t_span  <- c(0.0, 10)
t_eval  <- seq(t_span[1], t_span[2], length.out = npoints)

modelODE <- function(tvec, state, parameters) list(as.vector(f(state, parameters, tvec)))
sol   <- deSolve::ode(y = u0_true, times = t_eval, func = modelODE,
                      parms = p_star, rtol = 1e-12, atol = 1e-14)
U_clean <- sol[, 2, drop = FALSE]
tt      <- sol[, 1, drop = FALSE]

nr       <- 0.10
noise_sd <- nr * sqrt(mean(as.vector(U_clean)^2))

n_reps <- 250L
z      <- qnorm(0.975)          # 95% two-sided
D      <- ncol(U_clean)

## ---- one replicate ---------------------------------------------------------
one_rep <- function(rep_id) {
  set.seed(10000L + rep_id)     # per-rep seed -> reproducible under parallelism
  U <- U_clean + rnorm(npoints * D, mean = 0, sd = noise_sd)

  out <- tryCatch({
    res <- solveWendy(f, U, tt, method = "IRLS",
                      control = list(
                        estimate_IC         = TRUE,
                        estimate_trajectory = FALSE,
                        test_fun_type       = "SSL"))
    bs    <- res$boundary_state
    u0hat <- as.numeric(res$u0hat)
    cv    <- bs$cov_u0
    list(u0hat = u0hat, cov_u0 = cv,
         converged = isTRUE(bs$converged), diverged = isTRUE(bs$diverged),
         phat = res$phat, err = NA_character_)
  }, error = function(e) list(u0hat = NA, cov_u0 = NULL, converged = NA,
                              diverged = NA, phat = NA, err = conditionMessage(e)))

  se        <- if (!is.null(out$cov_u0)) sqrt(diag(out$cov_u0)) else rep(NA_real_, D)
  lo        <- out$u0hat - z * se
  hi        <- out$u0hat + z * se
  covered   <- u0_true >= lo & u0_true <= hi

  data.frame(rep_id    = rep_id,
             u0hat     = out$u0hat[1],
             se        = se[1],
             lo        = lo[1],
             hi        = hi[1],
             covered   = covered[1],
             abs_err   = abs(out$u0hat[1] - u0_true[1]),
             converged = out$converged,
             diverged  = out$diverged,
             null_cov  = is.null(out$cov_u0),
             err       = out$err)
}

## ---- run -------------------------------------------------------------------
ncores <- max(1L, parallel::detectCores() - 1L)
cat(sprintf("Running %d reps on %d cores (sigma = %.4f, nr = %.2f)...\n",
            n_reps, ncores, noise_sd, nr))

t0  <- Sys.time()
rows <- parallel::mclapply(seq_len(n_reps), one_rep, mc.cores = ncores)
df   <- do.call(rbind, rows)
cat(sprintf("Done in %.1fs\n", as.numeric(difftime(Sys.time(), t0, units = "secs"))))

## ---- summarise -------------------------------------------------------------
ok      <- df[!df$null_cov & is.finite(df$se), ]
n_ok    <- nrow(ok)
cov_hat <- mean(ok$covered)
mc_se   <- sqrt(cov_hat * (1 - cov_hat) / n_ok)     # Monte-Carlo SE of coverage

cat("\n================ u0 coverage summary (logistic) ================\n")
cat(sprintf("reps requested        : %d\n", n_reps))
cat(sprintf("reps with usable cov  : %d  (NULL cov: %d, diverged: %d)\n",
            n_ok, sum(df$null_cov), sum(df$diverged, na.rm = TRUE)))
cat(sprintf("true u0               : %.4f\n", u0_true))
cat(sprintf("mean u0hat            : %.5f   (bias = %+.5f)\n",
            mean(ok$u0hat), mean(ok$u0hat) - u0_true))
cat(sprintf("empirical SD(u0hat)   : %.5f   <- actual sampling spread\n", sd(ok$u0hat)))
cat(sprintf("mean model SE         : %.5f   <- sqrt(cov_u0) from the formula\n", mean(ok$se)))
cat(sprintf("calibration ratio     : %.3f   (mean SE / SD; 1.0 = well calibrated)\n",
            mean(ok$se) / sd(ok$u0hat)))
cat(sprintf("RMSE(u0hat)           : %.5f\n", sqrt(mean(ok$abs_err^2))))
cat(sprintf("\nEMPIRICAL 95%% COVERAGE : %.1f%%  (+/- %.1f%% MC-SE)\n",
            100 * cov_hat, 100 * mc_se))
cat("================================================================\n")

saveRDS(df, "examples/coverage_u0_logistic.rds")
write.csv(df, "examples/coverage_u0_logistic.csv", row.names = FALSE)
cat("Per-rep results saved to examples/coverage_u0_logistic.{rds,csv}\n")
