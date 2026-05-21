# Compare estimate_u0 variants on the logistic problem.
#
# The augmented BL weak residual in estimate_u0() can be solved over
#   - "both": both U[1, ] and U[M, ] (default, 2*D unknowns),
#   - "u0":   only U[1, ]            (D unknowns, U[M, ] held at observation),
#   - "uM":   only U[M, ]            (D unknowns, U[1, ] held at observation).
#
# This script sweeps a grid of noise ratios and seeds, runs all three variants
# against the same fitted parameter vector (so phat is identical across
# variants), and reports absolute error against the *true* corner states. A
# noiseless sanity check runs first.

library(deSolve)

invisible({devtools::load_all()})

# ---- problem -----------------------------------------------------------------

f <- function(u, p, t) c(p[1] * u[1] - p[2] * u[1]^2)
p_star  <- c(1, 1/10)
u0_true <- c(0.1)
npoints <- 128
t_span  <- c(0.0, 10)
t_eval  <- seq(t_span[1], t_span[2], length.out = npoints)

modelODE <- function(tvec, state, parameters) { list(as.vector(f(state, parameters, tvec))) }
sol      <- deSolve::ode(y = u0_true, times = t_eval, func = modelODE,
                         parms = p_star, rtol = 1e-12, atol = 1e-14)
U_clean  <- sol[, 2, drop = FALSE]
uM_true  <- as.numeric(U_clean[nrow(U_clean), ])
tt       <- sol[, 1, drop = FALSE]

# Helper: run solveWendy without doing the boundary state internally, then call
# estimate_u0 once per `corners` variant with the same phat / r_c. This isolates
# the comparison to the corner solve itself.
run_variants <- function(U, tt, corners_list = c("both", "u0", "uM")) {
  res <- solveWendy(
    f, U, tt,
    method  = "IRLS",
    control = list(estimate_u0 = FALSE, test_fun_type = "SSL")
  )
  ctx <- res$opt_ctx
  r_c_bl <- if (identical(ctx$control$test_fun_type, "SSL")) res$rc else res$min_radius

  variants <- lapply(corners_list, function(cc) {
    estimate_u0(U, ctx$f_, ctx$dF_dt_, ctx$d2F_dt2_, ctx$d3F_dt3_,
                tt, res$phat, r_c_bl, corners = cc)
  })
  names(variants) <- corners_list
  list(phat = res$phat, r_c = r_c_bl, variants = variants)
}

# ---- 1. sanity check on noiseless data ---------------------------------------

cat("=== Sanity check: noiseless logistic ===\n")
clean <- run_variants(U_clean, tt)
for (cc in names(clean$variants)) {
  v <- clean$variants[[cc]]
  cat(sprintf("  corners=%-4s  |u0hat - u0|=%.3e  |uMhat - uM|=%.3e\n",
              cc,
              abs(v$u0hat - u0_true),
              abs(v$uMhat - uM_true)))
}
stopifnot(abs(clean$variants$both$u0hat - u0_true) < 1e-2,
          abs(clean$variants$u0$u0hat   - u0_true) < 1e-2,
          abs(clean$variants$uM$uMhat   - uM_true) < 1e-2)

# ---- 2. noise sweep ----------------------------------------------------------

noise_ratios <- c(0.01, 0.02, 0.05, 0.10, 0.20)
n_seeds      <- 20
seeds        <- 1000 + seq_len(n_seeds)

results <- data.frame()

cat("\n=== Noise sweep ===\n")
for (nr in noise_ratios) {
  for (seed in seeds) {
    set.seed(seed)
    U_vec    <- as.vector(U_clean)
    noise_sd <- nr * sqrt(mean(U_vec^2))
    noise    <- rnorm(npoints, mean = 0, sd = noise_sd)
    U        <- U_clean + noise

    out <- tryCatch(run_variants(U, tt),
                    error = function(e) { message("failure: ", e$message); NULL })
    if (is.null(out)) next

    p_rel <- rel_err(out$phat, p_star)
    u0_obs_err <- abs(U[1, 1]              - u0_true)
    uM_obs_err <- abs(U[nrow(U), 1]        - uM_true)
    for (cc in names(out$variants)) {
      v <- out$variants[[cc]]
      results <- rbind(results, data.frame(
        noise_ratio = nr,
        seed        = seed,
        corners     = cc,
        r_c         = out$r_c,
        phat_relerr = p_rel,
        u0_err      = abs(v$u0hat - u0_true),
        uM_err      = abs(v$uMhat - uM_true),
        u0_obs_err  = u0_obs_err,
        uM_obs_err  = uM_obs_err
      ))
    }
  }
  cat(sprintf("  nr=%.2f done (%d seeds)\n", nr, n_seeds))
}

# ---- 3. summary --------------------------------------------------------------

agg <- aggregate(cbind(u0_err, uM_err, phat_relerr) ~ noise_ratio + corners,
                 data = results, FUN = median)
agg <- agg[order(agg$noise_ratio, agg$corners), ]
cat("\n=== Median errors by (noise_ratio, corners) ===\n")
print(agg, row.names = FALSE, digits = 4)

# Pairwise: how often does u0-only beat both for u0, and uM-only beat both for uM?
wide <- reshape(
  results[, c("noise_ratio", "seed", "corners", "u0_err", "uM_err")],
  idvar     = c("noise_ratio", "seed"),
  timevar   = "corners",
  direction = "wide"
)
win_table <- do.call(rbind, lapply(noise_ratios, function(nr) {
  sub <- wide[wide$noise_ratio == nr, ]
  data.frame(
    noise_ratio        = nr,
    n                  = nrow(sub),
    u0only_beats_both  = mean(sub$u0_err.u0 < sub$u0_err.both),
    uMonly_beats_both  = mean(sub$uM_err.uM < sub$uM_err.both),
    median_u0_both     = median(sub$u0_err.both),
    median_u0_u0only   = median(sub$u0_err.u0),
    median_uM_both     = median(sub$uM_err.both),
    median_uM_uMonly   = median(sub$uM_err.uM)
  )
}))
cat("\n=== Win rate of single-corner variants vs both ===\n")
print(win_table, row.names = FALSE, digits = 4)

cat("\n=== Also: corner error vs raw observation ===\n")
obs_table <- do.call(rbind, lapply(noise_ratios, function(nr) {
  sub <- wide[wide$noise_ratio == nr, ]
  obs_u0 <- abs(results$u0_obs_err[results$noise_ratio == nr & results$corners == "both"])
  obs_uM <- abs(results$uM_obs_err[results$noise_ratio == nr & results$corners == "both"])
  data.frame(
    noise_ratio     = nr,
    median_u0_obs   = median(obs_u0),
    median_u0_both  = median(sub$u0_err.both),
    median_u0_u0only = median(sub$u0_err.u0),
    median_uM_obs   = median(obs_uM),
    median_uM_both  = median(sub$uM_err.both),
    median_uM_uMonly = median(sub$uM_err.uM)
  )
}))
print(obs_table, row.names = FALSE, digits = 4)

# Save raw runs for later analysis.
out_csv <- file.path("examples", "estimate_u0_compare_results.csv")
write.csv(results, out_csv, row.names = FALSE)
cat(sprintf("\nWrote %d rows to %s\n", nrow(results), out_csv))
