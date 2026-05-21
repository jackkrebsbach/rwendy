# Sweep the number of boundary-layer test functions per side (n_bl) used inside
# estimate_u0() on the logistic problem. Holds r_c fixed (SSL change-point) so
# the BL window width is unchanged; only the count (and the implicit placement
# step inside build_boundary_layer_block) varies.
#
# Reports per-(noise, n_bl) median absolute error on u0 / uM, the LM-derived
# corner covariance (a proxy for prior tightness), and the effective placement
# step that build_boundary_layer_block chooses.

library(deSolve)

invisible({devtools::load_all(quiet = TRUE)})

# ---- problem -----------------------------------------------------------------

f <- function(u, p, t) c(p[1] * u[1] - p[2] * u[1]^2)
p_star  <- c(1, 1/10)
u0_true <- c(0.1)
npoints <- 128
t_span  <- c(0.0, 10)
t_eval  <- seq(t_span[1], t_span[2], length.out = npoints)

modelODE <- function(tvec, state, parameters) { list(as.vector(f(state, parameters, tvec))) }
sol     <- deSolve::ode(y = u0_true, times = t_eval, func = modelODE,
                        parms = p_star, rtol = 1e-12, atol = 1e-14)
U_clean <- sol[, 2, drop = FALSE]
uM_true <- as.numeric(U_clean[nrow(U_clean), ])
tt      <- sol[, 1, drop = FALSE]

# ---- one shot: run estimate_u0 over several n_bl values with the same phat ----

run_sweep_one_dataset <- function(U, tt, n_bl_grid) {
  res <- solveWendy(
    f, U, tt,
    method  = "IRLS",
    control = list(estimate_u0 = FALSE, test_fun_type = "SSL")
  )
  ctx <- res$opt_ctx
  r_c <- res$rc

  out <- lapply(n_bl_grid, function(n_bl) {
    if (n_bl > r_c) return(NULL)  # n_bl > r_c is nonsensical (no peaks available)
    v <- estimate_u0(U, ctx$f_, ctx$dF_dt_, ctx$d2F_dt2_, ctx$d3F_dt3_,
                     tt, res$phat, r_c,
                     corners = "both",
                     n_bl_per_side = n_bl)
    # Placement step used by build_boundary_layer_block (for reporting).
    step_used <- max(1L, floor(r_c / n_bl) - 2L)
    list(n_bl = n_bl, step = step_used,
         u0hat = v$u0hat, uMhat = v$uMhat,
         var_u0 = as.numeric(v$cov_u0), var_uM = as.numeric(v$cov_uM))
  })
  list(r_c = r_c, phat = res$phat, sweeps = out)
}

# ---- sweep over noise ratios and seeds ---------------------------------------

n_bl_grid    <- c(3, 4, 5, 6, 7, 8, 10, 14)
noise_ratios <- c(0.01, 0.02, 0.05, 0.10, 0.20)
n_seeds      <- 20
seeds        <- 1000 + seq_len(n_seeds)

results <- data.frame()

cat("=== BL count sweep (logistic, SSL, npoints=128) ===\n")
for (nr in noise_ratios) {
  for (seed in seeds) {
    set.seed(seed)
    U_vec    <- as.vector(U_clean)
    noise_sd <- nr * sqrt(mean(U_vec^2))
    U        <- U_clean + rnorm(npoints, mean = 0, sd = noise_sd)

    sw <- tryCatch(run_sweep_one_dataset(U, tt, n_bl_grid),
                   error = function(e) NULL)
    if (is.null(sw)) next

    for (s in sw$sweeps) {
      if (is.null(s)) next
      results <- rbind(results, data.frame(
        noise_ratio = nr,
        seed        = seed,
        r_c         = sw$r_c,
        n_bl        = s$n_bl,
        step        = s$step,
        u0_err      = abs(s$u0hat - u0_true),
        uM_err      = abs(s$uMhat - uM_true),
        var_u0      = s$var_u0,
        var_uM      = s$var_uM
      ))
    }
  }
  cat(sprintf("  nr=%.2f done\n", nr))
}

cat(sprintf("\nFitted SSL r_c values: %s (median %d)\n",
            paste(sort(unique(results$r_c)), collapse=", "),
            median(results$r_c)))

# ---- summary by (noise, n_bl) -----------------------------------------------

agg <- aggregate(cbind(u0_err, uM_err, var_u0) ~ noise_ratio + n_bl + step,
                 data = results, FUN = median)
agg <- agg[order(agg$noise_ratio, agg$n_bl), ]
cat("\n=== Median errors and prior variance ===\n")
print(agg, row.names = FALSE, digits = 4)

# Identify the best n_bl per noise ratio (lowest median u0_err).
best <- do.call(rbind, lapply(noise_ratios, function(nr) {
  sub <- agg[agg$noise_ratio == nr, ]
  if (nrow(sub) == 0) return(NULL)
  sub[which.min(sub$u0_err), ]
}))
cat("\n=== Best n_bl per noise ratio (min median u0_err) ===\n")
print(best, row.names = FALSE, digits = 4)

# Across noise ratios, is there a globally good choice?
global_rank <- aggregate(u0_err ~ n_bl, data = agg, FUN = mean)
global_rank <- global_rank[order(global_rank$u0_err), ]
cat("\n=== n_bl ranked by mean-over-noise median u0_err ===\n")
print(global_rank, row.names = FALSE, digits = 4)

out_csv <- file.path("examples", "bl_count_compare_results.csv")
write.csv(results, out_csv, row.names = FALSE)
cat(sprintf("\nWrote %d rows to %s\n", nrow(results), out_csv))
