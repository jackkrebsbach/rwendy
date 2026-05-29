# Hyperparameter sweep + cross-method comparison: JSTATE vs JOINT.
#
# Grid:
#   systems x noise x mp1 = {logistic, lv, lorenz} x {0.05, 0.15, 0.35} x {64, 256, 512}
#   seeds   = 5
#   JSTATE  = 36 configs   (jstate_lambda x jstate_pen x jstate_data_pen)
#   JOINT   = 24 configs   (joint_lambda x joint_basis_rank x joint_deriv_pen x joint_state_warmstart)
# Total runs: 27 cells x 5 seeds x 60 configs = 8100 fits.
#
# Outputs (under /tmp/):
#   bench_jstate_joint.csv     long-format results, one row per fit
#   bench_best_configs.csv     best config (median p_rel) per (system, nr, mp1, method)
#   bench_p_rel_grid.png       3x3 grid of best p_rel vs mp1, JSTATE vs JOINT
#   bench_state_rmse_grid.png  same for state RMSE
#
# Re-run analysis only (post-sweep): set RECOMPUTE <- FALSE below.

suppressMessages(invisible(devtools::load_all(".", quiet = TRUE)))
suppressMessages({
  library(deSolve)
  library(parallel)
})

RECOMPUTE   <- TRUE
OUT_DIR     <- "/tmp"
OUT_CSV     <- file.path(OUT_DIR, "bench_jstate_joint.csv")
OUT_BEST    <- file.path(OUT_DIR, "bench_best_configs.csv")
OUT_PLOT_P  <- file.path(OUT_DIR, "bench_p_rel_grid.png")
OUT_PLOT_S  <- file.path(OUT_DIR, "bench_state_rmse_grid.png")

# ----- Systems --------------------------------------------------------------
SYSTEMS <- list(
  logistic = list(
    f = function(u, p, t) c(p[1] * u[1] - p[2] * u[1]^2),
    p_star = c(1, 1/10), u0 = c(0.1), tspan = c(0, 10)
  ),
  lv = list(
    f = function(u, p, t) c(p[1]*u[1] + p[2]*u[1]*u[2],
                            p[3]*u[2] + p[4]*u[1]*u[2]),
    p_star = c(1, -0.1, -1.5, 0.075), u0 = c(10, 5), tspan = c(0, 10)
  ),
  # Lorenz: tspan kept short (3 ~ 1.5 Lyapunov times) so mp1=64 isn't catastrophic.
  lorenz = list(
    f = function(u, p, t) c(p[1] * (u[2] - u[1]),
                            u[1] * (p[2] - u[3]) - u[2],
                            u[1] * u[2] - p[3] * u[3]),
    p_star = c(10, 28, 8/3), u0 = c(-8, 10, 27), tspan = c(0, 3)
  )
)

# ----- Hyperparameter grids -------------------------------------------------
JSTATE_CONFIGS <- (function() {
  out <- list()
  for (lam in list("auto", 1e3, 1e4))
    for (pen in list(0, 1, 100, 1e4))
      for (dpn in list(0.1, 1, 10))
        out[[length(out) + 1L]] <- list(jstate_lambda   = lam,
                                        jstate_pen      = pen,
                                        jstate_data_pen = dpn)
  out
})()

JOINT_CONFIGS <- (function() {
  out <- list()
  for (lam in list("auto", 300, 3000))
    for (rk  in list("auto", "full"))
      for (dp  in list(0, 1e-4))
        for (ws  in list(FALSE, TRUE))
          out[[length(out) + 1L]] <- list(joint_lambda          = lam,
                                          joint_basis_rank      = rk,
                                          joint_deriv_pen       = dp,
                                          joint_state_warmstart = ws)
  out
})()

# ----- Sweep dimensions -----------------------------------------------------
NOISE_LEVELS <- c(0.05, 0.15, 0.35)
MP1_LEVELS   <- c(64, 256, 512)
N_SEEDS      <- 5L
SYS_NAMES    <- names(SYSTEMS)

# ----- Helpers --------------------------------------------------------------
gen_data <- function(sys, mp1, nr, seed) {
  t_eval <- seq(sys$tspan[1], sys$tspan[2], length.out = mp1)
  modelODE <- function(tv, st, pa) list(as.vector(sys$f(st, pa, tv)))
  sol <- deSolve::ode(y = sys$u0, times = t_eval, func = modelODE,
                      parms = sys$p_star, rtol = 1e-10, atol = 1e-12)
  Uc <- sol[, -1, drop = FALSE]
  set.seed(seed)
  noise_sd <- nr * sqrt(mean(as.vector(Uc)^2))
  U <- Uc + matrix(rnorm(length(Uc), sd = noise_sd), nrow = nrow(Uc))
  list(U = U, tt = matrix(sol[, 1], ncol = 1), Uc = Uc)
}

config_to_str <- function(cfg) {
  vals <- vapply(cfg, function(v) {
    if (is.logical(v))      as.character(v)
    else if (is.numeric(v)) formatC(v, format = "g", digits = 4)
    else                    as.character(v)
  }, character(1))
  paste(names(cfg), vals, sep = "=", collapse = ", ")
}

NA_ROW <- function(spec, elapsed) {
  data.frame(system = spec$system, nr = spec$nr, mp1 = spec$mp1, seed = spec$seed,
             method = spec$method, config_id = spec$config_id,
             config_str = config_to_str(spec$config),
             p_rel = NA_real_, state_rmse = NA_real_, elapsed = elapsed,
             converged = FALSE, iter = NA_integer_, stringsAsFactors = FALSE)
}

run_one <- function(spec) {
  sys <- SYSTEMS[[spec$system]]
  d <- tryCatch(gen_data(sys, spec$mp1, spec$nr, spec$seed),
                error = function(e) NULL)
  if (is.null(d)) return(NA_ROW(spec, 0))
  t1 <- proc.time()[3]
  res <- tryCatch(
    suppressWarnings(suppressMessages(
      solveWendy(sys$f, d$U, d$tt, method = spec$method, control = spec$config)
    )),
    error = function(e) NULL
  )
  elapsed <- proc.time()[3] - t1
  if (is.null(res) || is.null(res$phat) || anyNA(res$phat)) return(NA_ROW(spec, elapsed))
  uhat <- res$data$uhat
  state_rmse <- if (!is.null(uhat) && all(dim(uhat) == dim(d$Uc)))
                  sqrt(mean((uhat - d$Uc)^2)) else NA_real_
  data.frame(system = spec$system, nr = spec$nr, mp1 = spec$mp1, seed = spec$seed,
             method = spec$method, config_id = spec$config_id,
             config_str = config_to_str(spec$config),
             p_rel = rel_err(res$phat, sys$p_star),
             state_rmse = state_rmse,
             elapsed = elapsed,
             converged = isTRUE(res$data$converged %||% TRUE),
             iter = as.integer(res$data$iterations %||% NA),
             stringsAsFactors = FALSE)
}

# ----- Build spec list ------------------------------------------------------
build_specs <- function() {
  specs <- list()
  for (snm in SYS_NAMES) for (nr in NOISE_LEVELS) for (mp1 in MP1_LEVELS)
    for (seed in seq_len(N_SEEDS)) {
      for (i in seq_along(JSTATE_CONFIGS))
        specs[[length(specs) + 1L]] <- list(
          system = snm, nr = nr, mp1 = mp1, seed = seed,
          method = "JSTATE", config = JSTATE_CONFIGS[[i]],
          config_id = sprintf("JSTATE_%02d", i))
      for (i in seq_along(JOINT_CONFIGS))
        specs[[length(specs) + 1L]] <- list(
          system = snm, nr = nr, mp1 = mp1, seed = seed,
          method = "JOINT", config = JOINT_CONFIGS[[i]],
          config_id = sprintf("JOINT_%02d", i))
    }
  specs
}

# ----- Sweep (chunked by (system, mp1) for memory safety + intermediate saves) -
N_CORES <- 4L   # mp1=512 jacobians are ~40 MB/worker -> keep cores modest

build_chunk_specs <- function(snm, mp1) {
  out <- list()
  for (nr in NOISE_LEVELS) for (seed in seq_len(N_SEEDS)) {
    for (i in seq_along(JSTATE_CONFIGS))
      out[[length(out) + 1L]] <- list(
        system = snm, nr = nr, mp1 = mp1, seed = seed,
        method = "JSTATE", config = JSTATE_CONFIGS[[i]],
        config_id = sprintf("JSTATE_%02d", i))
    for (i in seq_along(JOINT_CONFIGS))
      out[[length(out) + 1L]] <- list(
        system = snm, nr = nr, mp1 = mp1, seed = seed,
        method = "JOINT", config = JOINT_CONFIGS[[i]],
        config_id = sprintf("JOINT_%02d", i))
  }
  out
}

if (RECOMPUTE) {
  chunks <- expand.grid(system = SYS_NAMES, mp1 = MP1_LEVELS, stringsAsFactors = FALSE)
  per_chunk <- length(NOISE_LEVELS) * N_SEEDS * (length(JSTATE_CONFIGS) + length(JOINT_CONFIGS))
  cat(sprintf("Total runs: %d  (%d chunks x %d specs/chunk, %d cores)\n",
              nrow(chunks) * per_chunk, nrow(chunks), per_chunk, N_CORES))
  all_rows <- list()
  t_global <- proc.time()[3]
  for (ci in seq_len(nrow(chunks))) {
    snm <- chunks$system[ci]; mp1 <- chunks$mp1[ci]
    specs <- build_chunk_specs(snm, mp1)
    t0 <- proc.time()[3]
    rows <- parallel::mclapply(specs, run_one, mc.cores = N_CORES,
                                mc.preschedule = FALSE)
    el <- proc.time()[3] - t0
    rows_df <- rows[vapply(rows, is.data.frame, logical(1))]
    df_chunk <- do.call(rbind, rows_df)
    n_err <- length(rows) - length(rows_df)
    n_fail <- if (nrow(df_chunk) > 0) sum(is.na(df_chunk$p_rel)) else 0
    cat(sprintf("[%d/%d] %-8s mp1=%-4d  %d/%d ok, %d NA, %d crash, %.1fs\n",
                ci, nrow(chunks), snm, mp1,
                nrow(df_chunk) - n_fail, length(specs), n_fail, n_err, el))
    all_rows[[ci]] <- df_chunk
    # Save partial after each chunk
    df_partial <- do.call(rbind, all_rows)
    write.csv(df_partial, OUT_CSV, row.names = FALSE)
  }
  cat(sprintf("Full sweep done in %.1f min\n", (proc.time()[3] - t_global) / 60))
  df <- do.call(rbind, all_rows)
  write.csv(df, OUT_CSV, row.names = FALSE)
  cat(sprintf("Wrote %s (%d rows, %d NA p_rel)\n",
              OUT_CSV, nrow(df), sum(is.na(df$p_rel))))
} else {
  df <- read.csv(OUT_CSV, stringsAsFactors = FALSE)
}

# ----- Best config per cell -------------------------------------------------
# Use median p_rel across seeds (NA-safe). Cells where every seed failed -> NA.
median_safe <- function(x) if (all(is.na(x))) NA_real_ else median(x, na.rm = TRUE)
agg <- aggregate(cbind(p_rel, state_rmse) ~ system + nr + mp1 + method + config_id + config_str,
                 data = df, FUN = median_safe, na.action = na.pass)
best <- do.call(rbind, lapply(
  split(agg, list(agg$system, agg$nr, agg$mp1, agg$method), drop = TRUE),
  function(g) if (all(is.na(g$p_rel))) g[1, ] else g[which.min(g$p_rel), ]))
best <- best[order(best$system, best$nr, best$mp1, best$method), ]
rownames(best) <- NULL
write.csv(best, OUT_BEST, row.names = FALSE)
cat(sprintf("Wrote %s (%d rows)\n", OUT_BEST, nrow(best)))

# ----- Cross-method comparison ---------------------------------------------
wide <- reshape(best[, c("system", "nr", "mp1", "method", "p_rel", "state_rmse")],
                idvar = c("system", "nr", "mp1"), timevar = "method", direction = "wide")
wide$winner_p     <- ifelse(wide$p_rel.JSTATE     < wide$p_rel.JOINT,     "JSTATE", "JOINT")
wide$winner_state <- ifelse(wide$state_rmse.JSTATE < wide$state_rmse.JOINT, "JSTATE", "JOINT")
wide$p_ratio      <- wide$p_rel.JSTATE     / wide$p_rel.JOINT
wide$state_ratio  <- wide$state_rmse.JSTATE / wide$state_rmse.JOINT
cat("\n========== Best median p_rel per cell ==========\n")
print(wide[, c("system","nr","mp1","p_rel.JSTATE","p_rel.JOINT","winner_p","p_ratio")],
      row.names = FALSE, digits = 4)
cat("\n========== Best median state_rmse per cell ==========\n")
print(wide[, c("system","nr","mp1","state_rmse.JSTATE","state_rmse.JOINT","winner_state","state_ratio")],
      row.names = FALSE, digits = 4)
cat(sprintf("\nJSTATE wins on p_rel: %d/%d cells   on state_rmse: %d/%d cells\n",
            sum(wide$winner_p     == "JSTATE", na.rm=TRUE), nrow(wide),
            sum(wide$winner_state == "JSTATE", na.rm=TRUE), nrow(wide)))

# ----- Plots ----------------------------------------------------------------
plot_grid <- function(file, metric, title) {
  png(file, width = 1200, height = 900, res = 110)
  par(mfrow = c(3, 3), mar = c(4, 4, 2.5, 1), oma = c(0, 0, 2, 0))
  for (snm in SYS_NAMES) for (nr in NOISE_LEVELS) {
    sub <- best[best$system == snm & best$nr == nr, ]
    ys <- sub[[metric]]
    if (all(is.na(ys))) {
      plot.new(); title(main = sprintf("%s, nr=%.2f", snm, nr)); next
    }
    ylim <- range(ys, finite = TRUE)
    plot(NA, xlim = range(MP1_LEVELS), ylim = ylim, log = "y", xaxt = "n",
         xlab = "mp1", ylab = paste("best", metric),
         main = sprintf("%s, nr=%.2f", snm, nr))
    axis(1, at = MP1_LEVELS)
    for (m in c("JSTATE", "JOINT")) {
      ss <- sub[sub$method == m, ]
      ss <- ss[order(ss$mp1), ]
      lines(ss$mp1, ss[[metric]], type = "b",
            col = ifelse(m == "JSTATE", "#1f77b4", "#d62728"),
            pch = 19, lwd = 2)
    }
    legend("topright", legend = c("JSTATE", "JOINT"),
           col = c("#1f77b4", "#d62728"), pch = 19, lwd = 2, bty = "n", cex = 0.8)
  }
  mtext(title, outer = TRUE, cex = 1.2, font = 2)
  dev.off()
}
plot_grid(OUT_PLOT_P, "p_rel",      "Best parameter rel. error per cell (lower is better)")
plot_grid(OUT_PLOT_S, "state_rmse", "Best state RMSE per cell (lower is better)")
cat(sprintf("Wrote %s\n", OUT_PLOT_P))
cat(sprintf("Wrote %s\n", OUT_PLOT_S))
