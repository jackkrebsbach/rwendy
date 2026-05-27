# Coverage study for cov_u0 with PARAMETER COVARIANCE PROPAGATED IN.
#
# Law of total variance:  Cov(u0hat) ~= Cov_LS(u0|p) + J Cov(phat) J^T,
# where J = d u0hat / d p (numeric Jacobian) and Cov(phat) is the Fisher
# param covariance (G^T S^-1 G)^-1 reported by summary(res)$param_cov.
#
# Same logistic setup / seeds as the other two studies (n=128, nr=0.10,
# 95% CI, 250 reps). We report coverage under BOTH the conditional LS
# covariance (cov_LS) and the augmented covariance (cov_total) so the
# effect of the parameter term is directly visible.

suppressMessages({
  library(deSolve)
  library(parallel)
  library(devtools)
  library(numDeriv)
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
  t_expr <- symengine::S("t")
  f_expr <- f(u_expr, p_expr, t_expr)
  vars   <- c(p_expr, u_expr, t_expr)
  Ju  <- compute_symbolic_jacobian(f_expr, u_expr)
  dF  <- compute_symbolic_total_time_deriv(f_expr, u_expr, f_expr, t_expr)
  d2  <- compute_symbolic_total_time_deriv(dF, u_expr, f_expr, t_expr)
  d3  <- compute_symbolic_total_time_deriv(d2, u_expr, f_expr, t_expr)
  list(f_ = build_fn(f_expr, vars), J_u = build_fn(Ju, vars),
       dF_dt_ = build_fn(dF, vars),
       d2F_dt2_ = build_fn(d2, vars), d3F_dt3_ = build_fn(d3, vars))
}

modelODE <- function(tv, s, pp) list(as.vector(f(s, pp, tv)))
sol     <- deSolve::ode(y = u0_true, times = t_eval, func = modelODE,
                        parms = p_star, rtol = 1e-12, atol = 1e-14)
U_clean <- sol[, 2, drop = FALSE]
tt      <- sol[, 1, drop = FALSE]
nr      <- 0.10
noise_sd<- nr * sqrt(mean(as.vector(U_clean)^2))
n_reps  <- 250L
z       <- qnorm(0.975)
D       <- ncol(U_clean)

one_rep <- function(rep_id) {
  set.seed(10000L + rep_id)
  U <- U_clean + rnorm(npoints * D, 0, noise_sd)
  out <- tryCatch({
    res  <- solveWendy(f, U, tt, method = "IRLS",
                       control = list(estimate_IC = TRUE,
                                      estimate_trajectory = FALSE,
                                      test_fun_type = "SSL"))
    pcov <- summary(res)$param_cov
    ev   <- build_evaluators(f, D, length(p_star))
    u0_of_p <- function(p)
      estimate_IC(U, ev$f_, ev$dF_dt_, ev$d2F_dt2_, ev$d3F_dt3_,
                  tt, p, res$rc, J_u = ev$J_u, sigma = noise_sd,
                  update_trap_u0 = FALSE)$u0hat
    bs   <- estimate_IC(U, ev$f_, ev$dF_dt_, ev$d2F_dt2_, ev$d3F_dt3_,
                        tt, res$phat, res$rc, J_u = ev$J_u, sigma = noise_sd,
                        update_trap_u0 = FALSE)
    Jp   <- numDeriv::jacobian(u0_of_p, res$phat)   # D x J
    cov_LS    <- bs$cov_u0_resid
    cov_param <- Jp %*% pcov %*% t(Jp)              # D x D
    cov_total <- cov_LS + cov_param
    list(u0hat = as.numeric(bs$u0hat),
         se_ls    = sqrt(diag(cov_LS)),
         se_param = sqrt(diag(cov_param)),
         se_total = sqrt(diag(cov_total)),
         err = NA_character_)
  }, error = function(e) list(u0hat = NA, se_ls = NA, se_param = NA,
                              se_total = NA, err = conditionMessage(e)))

  data.frame(rep_id = rep_id, u0hat = out$u0hat[1],
             se_ls = out$se_ls[1], se_param = out$se_param[1],
             se_total = out$se_total[1],
             cov_ls    = abs(out$u0hat[1] - u0_true) <= z * out$se_ls[1],
             cov_total = abs(out$u0hat[1] - u0_true) <= z * out$se_total[1],
             err = out$err)
}

ncores <- max(1L, parallel::detectCores() - 1L)
cat(sprintf("Running %d reps (param-propagated) on %d cores...\n", n_reps, ncores))
t0   <- Sys.time()
rows <- parallel::mclapply(seq_len(n_reps), one_rep, mc.cores = ncores)
df   <- do.call(rbind, rows)
cat(sprintf("Done in %.1fs\n", as.numeric(difftime(Sys.time(), t0, units = "secs"))))

ok <- df[is.finite(df$se_total), ]
sd_emp <- sd(ok$u0hat)
cat("\n===== u0 coverage with parameter covariance propagated (logistic) =====\n")
cat(sprintf("reps usable           : %d\n", nrow(ok)))
cat(sprintf("empirical SD(u0hat)   : %.5f\n", sd_emp))
cat(sprintf("mean SE (LS only)     : %.5f   (calib ratio %.3f)\n",
            mean(ok$se_ls), mean(ok$se_ls) / sd_emp))
cat(sprintf("mean SE (param only)  : %.5f\n", mean(ok$se_param)))
cat(sprintf("mean SE (total)       : %.5f   (calib ratio %.3f)\n",
            mean(ok$se_total), mean(ok$se_total) / sd_emp))
cat(sprintf("param share of total var (mean) : %.2f%%\n",
            100 * mean(ok$se_param^2 / ok$se_total^2)))
cat(sprintf("\n95%% COVERAGE  LS only : %.1f%%\n", 100 * mean(ok$cov_ls)))
cat(sprintf("95%% COVERAGE  total   : %.1f%%\n", 100 * mean(ok$cov_total)))
cat("=======================================================================\n")

saveRDS(df, "examples/coverage_u0_logistic_propagate.rds")
write.csv(df, "examples/coverage_u0_logistic_propagate.csv", row.names = FALSE)
