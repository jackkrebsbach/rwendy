# Coverage study: NOISE-PROPAGATION conditional variance + parameter term.
#
# The estimator u0 = B^T rhs / B^T B is a deterministic function of the data
# U, so propagate the data noise through it by the delta method:
#
#     Cov_noise(u0 | p) = sigma^2 * (d u0 / d vec(U)) (d u0 / d vec(U))^T
#
# (numeric Jacobian w.r.t. U; captures the nonlinear f(p,U,t) term, the
# linear Vp_BL U term, and the implicit EM(u0) coupling). This uses ONLY the
# propagated noise -- unlike the shipped cov_LS = ||B u0 - rhs||^2/(K_bl-1)
# B^TB, whose residual norm also absorbs deterministic EM/trapezoidal
# truncation error and is therefore ~1.8x too wide in SE.
#
# Total reported covariance adds the Fisher parameter term:
#     Cov_total = Cov_noise(u0|p) + J Cov(phat) J^T,   J = d u0 / d p.
#
# sigma is the library's own estimate (estimate_std, k=6). Same logistic
# setup / seeds as the other studies (n=128, nr=0.10, 95% CI, 250 reps).

suppressMessages({
  library(deSolve); library(parallel); library(devtools); library(numDeriv)
})
invisible(devtools::load_all(quiet = TRUE))

f <- function(u, p, t) c(p[1] * u[1] - p[2] * u[1]^2)

p_star  <- c(1, 1/10)
u0_true <- c(0.1)
npoints <- 128L
t_eval  <- seq(0, 10, length.out = npoints)

build_evaluators <- function(f, D, J) {
  u_expr <- do.call(c, lapply(1:D, function(i) symengine::S(paste0("u", i))))
  p_expr <- do.call(c, lapply(1:J, function(i) symengine::S(paste0("p", i))))
  t_expr <- symengine::S("t"); f_expr <- f(u_expr, p_expr, t_expr)
  vars   <- c(p_expr, u_expr, t_expr)
  Ju <- compute_symbolic_jacobian(f_expr, u_expr)
  dF <- compute_symbolic_total_time_deriv(f_expr, u_expr, f_expr, t_expr)
  d2 <- compute_symbolic_total_time_deriv(dF, u_expr, f_expr, t_expr)
  d3 <- compute_symbolic_total_time_deriv(d2, u_expr, f_expr, t_expr)
  list(f_ = build_fn(f_expr, vars), J_u = build_fn(Ju, vars),
       dF_dt_ = build_fn(dF, vars),
       d2F_dt2_ = build_fn(d2, vars), d3F_dt3_ = build_fn(d3, vars))
}

modelODE <- function(tv, s, pp) list(as.vector(f(s, pp, tv)))
sol     <- deSolve::ode(y = u0_true, times = t_eval, func = modelODE,
                        parms = p_star, rtol = 1e-12, atol = 1e-14)
U_clean <- sol[, 2, drop = FALSE]
tt      <- sol[, 1, drop = FALSE]
M       <- npoints; D <- ncol(U_clean)
nr      <- 0.10
noise_sd<- nr * sqrt(mean(as.vector(U_clean)^2))
n_reps  <- 250L
z       <- qnorm(0.975)

one_rep <- function(rep_id) {
  set.seed(10000L + rep_id)
  U <- U_clean + rnorm(M * D, 0, noise_sd)
  out <- tryCatch({
    res    <- solveWendy(f, U, tt, method = "IRLS",
                         control = list(estimate_IC = TRUE,
                                        estimate_trajectory = FALSE,
                                        test_fun_type = "SSL"))
    phat   <- res$phat
    pcov   <- summary(res)$param_cov
    sig    <- estimate_std(U, k = 6)
    ev     <- build_evaluators(f, D, length(p_star))
    ic <- function(p, Umat)
      estimate_IC(Umat, ev$f_, ev$dF_dt_, ev$d2F_dt2_, ev$d3F_dt3_,
                  tt, p, res$rc, J_u = ev$J_u, sigma = sig,
                  update_trap_u0 = FALSE)
    bs   <- ic(phat, U)

    # noise propagation: d u0 / d vec(U)  (D x M*D), forward differences
    JU   <- numDeriv::jacobian(function(uv) ic(phat, matrix(uv, M, D))$u0hat,
                               as.vector(U), method = "simple")
    cov_noise <- sig^2 * (JU %*% t(JU))                    # D x D

    # parameter term
    Jp   <- numDeriv::jacobian(function(p) ic(p, U)$u0hat, phat)   # D x J
    cov_param <- Jp %*% pcov %*% t(Jp)

    list(u0hat = as.numeric(bs$u0hat),
         se_old   = sqrt(diag(bs$cov_u0_resid)),
         se_noise = sqrt(diag(cov_noise)),
         se_total = sqrt(diag(cov_noise + cov_param)),
         se_par   = sqrt(diag(cov_param)),
         sig = sig, err = NA_character_)
  }, error = function(e) list(u0hat = NA, se_old = NA, se_noise = NA,
                              se_total = NA, se_par = NA, sig = NA,
                              err = conditionMessage(e)))

  data.frame(rep_id = rep_id, u0hat = out$u0hat[1], sig = out$sig[1],
             se_old = out$se_old[1], se_noise = out$se_noise[1],
             se_par = out$se_par[1], se_total = out$se_total[1],
             cov_old   = abs(out$u0hat[1] - u0_true) <= z * out$se_old[1],
             cov_noise = abs(out$u0hat[1] - u0_true) <= z * out$se_noise[1],
             cov_total = abs(out$u0hat[1] - u0_true) <= z * out$se_total[1],
             err = out$err)
}

ncores <- max(1L, parallel::detectCores() - 1L)
cat(sprintf("Running %d reps (noise propagation) on %d cores...\n", n_reps, ncores))
t0   <- Sys.time()
rows <- parallel::mclapply(seq_len(n_reps), one_rep, mc.cores = ncores)
df   <- do.call(rbind, rows)
cat(sprintf("Done in %.1fs\n", as.numeric(difftime(Sys.time(), t0, units = "secs"))))

ok <- df[is.finite(df$se_total), ]
sd_emp <- sd(ok$u0hat)
mc <- function(p) sqrt(p * (1 - p) / nrow(ok))
cat("\n=========== u0 coverage: noise-propagation variance (logistic) ===========\n")
cat(sprintf("reps usable           : %d   (mean sigma_hat = %.4f, true = %.4f)\n",
            nrow(ok), mean(ok$sig), noise_sd))
cat(sprintf("empirical SD(u0hat)   : %.5f   <- target for the SE\n", sd_emp))
cat(sprintf("%-22s %8s %10s %12s\n", "estimator", "mean SE", "calib", "95% cov"))
cat(sprintf("%-22s %8.4f %10.3f %10.1f%%   (shipped, too wide)\n",
            "cov_LS (old)",   mean(ok$se_old),   mean(ok$se_old)/sd_emp,   100*mean(ok$cov_old)))
cat(sprintf("%-22s %8.4f %10.3f %10.1f%%\n",
            "noise-only",     mean(ok$se_noise), mean(ok$se_noise)/sd_emp, 100*mean(ok$cov_noise)))
cat(sprintf("%-22s %8.4f %10.3f %10.1f%%   (+/- %.1f%% MC)\n",
            "noise + param",  mean(ok$se_total), mean(ok$se_total)/sd_emp, 100*mean(ok$cov_total),
            100*mc(mean(ok$cov_total))))
cat(sprintf("\nmean param SE = %.5f  (param share of total var = %.2f%%)\n",
            mean(ok$se_par), 100 * mean(ok$se_par^2 / ok$se_total^2)))
cat("nominal-vs-empirical (noise+param):\n")
for (lvl in c(0.50, 0.80, 0.90, 0.95)) {
  zc <- qnorm(0.5 + lvl/2)
  cat(sprintf("  %2.0f%% -> %5.1f%%\n", 100*lvl,
              100*mean(abs(ok$u0hat - u0_true) <= zc * ok$se_total)))
}
cat("==========================================================================\n")

saveRDS(df, "examples/coverage_u0_logistic_noiseprop.rds")
write.csv(df, "examples/coverage_u0_logistic_noiseprop.csv", row.names = FALSE)
