# Coverage study for cov_u0, CONDITIONING ON THE TRUE PARAMETERS.
#
# Identical to coverage_u0_logistic.R (logistic, n=128, nr=0.10, 95% CI,
# 250 reps) except estimate_IC is called with p = p_star instead of the
# pipeline's p-hat. We still run solveWendy each rep so the symbolic
# evaluators (f_, dF_dt_, ...) and the SSL window r_c are byte-for-byte the
# same as the p-hat study; only the parameter vector fed to estimate_IC
# changes. This isolates whether the closed-form LS variance s^2/B^TB is
# itself calibrated, free of p-hat sampling noise.

suppressMessages({
  library(deSolve)
  library(parallel)
  library(devtools)
})
invisible(devtools::load_all(quiet = TRUE))

## ---- problem setup (identical to the p-hat study) --------------------------
f <- function(u, p, t) c(p[1] * u[1] - p[2] * u[1]^2)

p_star  <- c(1, 1/10)
u0_true <- c(0.1)
npoints <- 128L
t_span  <- c(0.0, 10)
t_eval  <- seq(t_span[1], t_span[2], length.out = npoints)

# Build the symbolic evaluators (same as solveWendy lines 88-129). The LLVM
# visitor closures don't survive fork(), so this is called inside each worker.
build_evaluators <- function(f, D, J) {
  u_expr <- do.call(c, lapply(1:D, function(i) symengine::S(paste0("u", i))))
  p_expr <- do.call(c, lapply(1:J, function(i) symengine::S(paste0("p", i))))
  t_expr <- symengine::S("t")
  f_expr <- f(u_expr, p_expr, t_expr)
  vars   <- c(p_expr, u_expr, t_expr)
  Ju_sym  <- compute_symbolic_jacobian(f_expr, u_expr)
  dF_sym  <- compute_symbolic_total_time_deriv(f_expr,  u_expr, f_expr, t_expr)
  d2F_sym <- compute_symbolic_total_time_deriv(dF_sym,  u_expr, f_expr, t_expr)
  d3F_sym <- compute_symbolic_total_time_deriv(d2F_sym, u_expr, f_expr, t_expr)
  list(f_      = build_fn(f_expr,  vars),
       J_u     = build_fn(Ju_sym,  vars),
       dF_dt_  = build_fn(dF_sym,  vars),
       d2F_dt2_= build_fn(d2F_sym, vars),
       d3F_dt3_= build_fn(d3F_sym, vars))
}

modelODE <- function(tvec, state, parameters) list(as.vector(f(state, parameters, tvec)))
sol   <- deSolve::ode(y = u0_true, times = t_eval, func = modelODE,
                      parms = p_star, rtol = 1e-12, atol = 1e-14)
U_clean <- sol[, 2, drop = FALSE]
tt      <- sol[, 1, drop = FALSE]

nr       <- 0.10
noise_sd <- nr * sqrt(mean(as.vector(U_clean)^2))

n_reps <- 250L
z      <- qnorm(0.975)
D      <- ncol(U_clean)

## ---- one replicate ---------------------------------------------------------
one_rep <- function(rep_id) {
  set.seed(10000L + rep_id)              # SAME seeds as the p-hat study
  U <- U_clean + rnorm(npoints * D, mean = 0, sd = noise_sd)

  out <- tryCatch({
    # Run the pipeline only to obtain the SSL window r_c on this dataset
    # (identical selection to the p-hat study); rebuild evaluators locally
    # and re-run estimate_IC at p_star instead of res$phat.
    res <- solveWendy(f, U, tt, method = "IRLS",
                      control = list(
                        estimate_IC         = TRUE,
                        estimate_trajectory = FALSE,
                        test_fun_type       = "SSL"))
    ev <- build_evaluators(f, D, length(p_star))
    bs <- estimate_IC(U, ev$f_, ev$dF_dt_, ev$d2F_dt2_, ev$d3F_dt3_,
                      tt, p_star, res$rc, J_u = ev$J_u, sigma = noise_sd,
                      update_trap_u0 = FALSE)
    list(u0hat = as.numeric(bs$u0hat), cov_u0 = bs$cov_u0,
         converged = isTRUE(bs$converged), diverged = isTRUE(bs$diverged),
         phat = res$phat, err = NA_character_)
  }, error = function(e) list(u0hat = NA, cov_u0 = NULL, converged = NA,
                              diverged = NA, phat = NA, err = conditionMessage(e)))

  se      <- if (!is.null(out$cov_u0)) sqrt(diag(out$cov_u0)) else rep(NA_real_, D)
  lo      <- out$u0hat - z * se
  hi      <- out$u0hat + z * se
  covered <- u0_true >= lo & u0_true <= hi

  data.frame(rep_id = rep_id, u0hat = out$u0hat[1], se = se[1],
             lo = lo[1], hi = hi[1], covered = covered[1],
             abs_err = abs(out$u0hat[1] - u0_true[1]),
             converged = out$converged, diverged = out$diverged,
             null_cov = is.null(out$cov_u0), err = out$err)
}

## ---- run -------------------------------------------------------------------
ncores <- max(1L, parallel::detectCores() - 1L)
cat(sprintf("Running %d reps at TRUE p on %d cores (sigma = %.4f, nr = %.2f)...\n",
            n_reps, ncores, noise_sd, nr))

t0   <- Sys.time()
rows <- parallel::mclapply(seq_len(n_reps), one_rep, mc.cores = ncores)
df   <- do.call(rbind, rows)
cat(sprintf("Done in %.1fs\n", as.numeric(difftime(Sys.time(), t0, units = "secs"))))

## ---- summarise -------------------------------------------------------------
ok      <- df[!df$null_cov & is.finite(df$se), ]
n_ok    <- nrow(ok)
cov_hat <- mean(ok$covered)
mc_se   <- sqrt(cov_hat * (1 - cov_hat) / n_ok)
zsc     <- (ok$u0hat - u0_true) / ok$se

cat("\n========== u0 coverage at TRUE p (logistic) ==========\n")
cat(sprintf("reps usable           : %d  (NULL cov: %d, diverged: %d)\n",
            n_ok, sum(df$null_cov), sum(df$diverged, na.rm = TRUE)))
cat(sprintf("mean u0hat            : %.5f   (bias = %+.5f)\n",
            mean(ok$u0hat), mean(ok$u0hat) - u0_true))
cat(sprintf("empirical SD(u0hat)   : %.5f\n", sd(ok$u0hat)))
cat(sprintf("mean model SE         : %.5f\n", mean(ok$se)))
cat(sprintf("calibration ratio     : %.3f   (mean SE / SD; 1.0 = calibrated)\n",
            mean(ok$se) / sd(ok$u0hat)))
cat(sprintf("std-resid sd(z)       : %.3f   (1.0 = calibrated; var overstated %.2fx)\n",
            sd(zsc), 1 / var(zsc)))
cat(sprintf("RMSE(u0hat)           : %.5f\n", sqrt(mean(ok$abs_err^2))))
cat(sprintf("\nEMPIRICAL 95%% COVERAGE : %.1f%%  (+/- %.1f%% MC-SE)\n",
            100 * cov_hat, 100 * mc_se))
for (lvl in c(0.50, 0.80, 0.90, 0.95)) {
  zc <- qnorm(0.5 + lvl / 2)
  cat(sprintf("  nominal %2.0f%% -> empirical %5.1f%%\n",
              100 * lvl, 100 * mean(abs(ok$u0hat - u0_true) <= zc * ok$se)))
}
cat("======================================================\n")

saveRDS(df, "examples/coverage_u0_logistic_pstar.rds")
write.csv(df, "examples/coverage_u0_logistic_pstar.csv", row.names = FALSE)
