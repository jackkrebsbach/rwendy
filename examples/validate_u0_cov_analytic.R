# Validate the ANALYTIC noise-propagation covariance for u0 against the
# numeric-Jacobian propagation (which routes through the real estimate_IC).
# Nothing in R/ is modified; this reconstructs estimate_IC's internal setup
# faithfully and computes the analytic Cov_noise, then compares.

suppressMessages({library(deSolve); library(devtools); library(numDeriv)})
invisible(devtools::load_all(quiet = TRUE))

f <- function(u, p, t) c(p[1] * u[1] - p[2] * u[1]^2)
build_evaluators <- function(f, D, J) {
  u_expr <- do.call(c, lapply(1:D, function(i) symengine::S(paste0("u", i))))
  p_expr <- do.call(c, lapply(1:J, function(i) symengine::S(paste0("p", i))))
  t_expr <- symengine::S("t"); f_expr <- f(u_expr, p_expr, t_expr)
  vars   <- c(p_expr, u_expr, t_expr)
  Ju_sym <- compute_symbolic_jacobian(f_expr, u_expr)
  dF <- compute_symbolic_total_time_deriv(f_expr, u_expr, f_expr, t_expr)
  d2 <- compute_symbolic_total_time_deriv(dF, u_expr, f_expr, t_expr)
  d3 <- compute_symbolic_total_time_deriv(d2, u_expr, f_expr, t_expr)
  list(f_ = build_fn(f_expr, vars), J_u = build_fn(Ju_sym, vars),
       dF_dt_ = build_fn(dF, vars), d2F_dt2_ = build_fn(d2, vars),
       d3F_dt3_ = build_fn(d3, vars))
}

# Analytic Cov_noise(u0|p), reconstructing estimate_IC's setup (lines 127-191).
analytic_cov_noise <- function(U, ev, tt, p, r_c, sigma,
                               n_bl = NULL, em_order = 4L) {
  M <- nrow(U); D <- ncol(U); J <- length(p)
  tt_vec <- as.vector(tt); dt <- mean(diff(tt_vec))
  r_c  <- min(r_c, floor((M - 1L) / 2L))
  n_bl <- if (is.null(n_bl)) max(3L, as.integer(ceiling(r_c / 4))) else max(1L, as.integer(n_bl))
  step <- max(1L, as.integer(floor((r_c - 1L) / max(1L, n_bl - 1L))))
  bl_left <- lapply(0:4, function(ord)
    build_boundary_layer_block(psi, tt_vec, r_c, order = ord, n_bl = n_bl, step = step))
  K_bl <- nrow(bl_left[[1]])
  apply_trap <- function(Bm) { Bm[, 1] <- Bm[, 1] * 0.5; Bm[, M] <- Bm[, M] * 0.5; Bm }
  V_BL  <- apply_trap(bl_left[[1]]); Vp_BL <- apply_trap(bl_left[[2]])
  bl_phi_t1 <- matrix(0, K_bl, 5); for (ord in 0:4) bl_phi_t1[, ord + 1] <- bl_left[[ord + 1]][, 1]
  B <- bl_phi_t1[, 1]; BtB <- sum(B * B)
  c2 <- dt^2 / 12; c4 <- if (em_order >= 4L) dt^4 / 720 else 0
  em_correction <- function(u0_curr) {
    u_t1 <- as.vector(u0_curr); inp_t1 <- matrix(c(p, u_t1, tt_vec[1]), ncol = 1L)
    fd_t1 <- list(as.vector(ev$f_(inp_t1)), as.vector(ev$dF_dt_(inp_t1)),
                  as.vector(ev$d2F_dt2_(inp_t1)), as.vector(ev$d3F_dt3_(inp_t1)))
    EM <- matrix(0, K_bl, D)
    for (k in seq_len(K_bl)) {
      phi_t1_k <- as.list(bl_phi_t1[k, ])
      g1 <- g_deriv_at_endpoint(phi_t1_k, fd_t1, u_t1, order = 1L)
      g3 <- if (em_order >= 4L) g_deriv_at_endpoint(phi_t1_k, fd_t1, u_t1, order = 3L) else rep(0, D)
      EM[k, ] <- c2 * g1 - c4 * g3
    }
    EM
  }
  # u0 from the REAL estimate_IC (so the linearization point is identical)
  u0 <- estimate_IC(U, ev$f_, ev$dF_dt_, ev$d2F_dt2_, ev$d3F_dt3_, tt, p, r_c,
                    J_u = ev$J_u, sigma = sigma, update_trap_u0 = FALSE)$u0hat
  I_D <- diag(D)
  h <- 1e-6 * max(1, sqrt(sum(u0^2))); A <- matrix(0, D, D)
  for (e_i in seq_len(D)) {
    up <- u0; up[e_i] <- up[e_i] + h; dn <- u0; dn[e_i] <- dn[e_i] - h
    dEM_e <- (em_correction(up) - em_correction(dn)) / (2 * h)
    A[, e_i] <- as.numeric(crossprod(B, dEM_e))
  }
  Minv <- solve(BtB * I_D + A)
  bV <- as.numeric(crossprod(B, V_BL)); bVp <- as.numeric(crossprod(B, Vp_BL))
  acc <- matrix(0, D, D)
  for (m in seq_len(M)) {
    if (bV[m] == 0 && bVp[m] == 0) next
    Ju_m <- matrix(as.vector(ev$J_u(c(p, U[m, ], tt_vec[m]))), D, D)
    C_m  <- -dt * bV[m] * Ju_m - dt * bVp[m] * I_D
    DU_m <- Minv %*% C_m
    acc  <- acc + tcrossprod(DU_m)
  }
  list(cov = sigma^2 * acc, u0 = u0)
}

## ---- compare on several seeds ----------------------------------------------
p_star <- c(1, 1/10); u0_true <- 0.1; np <- 128L; M <- np; D <- 1L
t_eval <- seq(0, 10, length.out = np)
modelODE <- function(tv, s, pp) list(as.vector(f(s, pp, tv)))
sol <- deSolve::ode(y = u0_true, times = t_eval, func = modelODE, parms = p_star, rtol = 1e-12, atol = 1e-14)
U_clean <- sol[, 2, drop = FALSE]; tt <- sol[, 1, drop = FALSE]
sigma_true <- 0.10 * sqrt(mean(as.vector(U_clean)^2))
ev <- build_evaluators(f, D, 2L)

cat(sprintf("%-6s %12s %12s %10s\n", "seed", "SE_analytic", "SE_numeric", "rel_diff"))
for (s in 1:8) {
  set.seed(10000L + s); U <- U_clean + rnorm(np, 0, sigma_true)
  res <- solveWendy(f, U, tt, method = "IRLS",
                    control = list(estimate_IC = TRUE, estimate_trajectory = FALSE, test_fun_type = "SSL"))
  phat <- res$phat; rc <- res$rc
  an <- analytic_cov_noise(U, ev, tt, phat, rc, sigma_true)
  ic <- function(uv) estimate_IC(matrix(uv, M, D), ev$f_, ev$dF_dt_, ev$d2F_dt2_, ev$d3F_dt3_,
                                 tt, phat, rc, J_u = ev$J_u, sigma = sigma_true,
                                 update_trap_u0 = FALSE)$u0hat
  JU <- numDeriv::jacobian(ic, as.vector(U), method = "simple")
  cov_num <- sigma_true^2 * (JU %*% t(JU))
  se_an <- sqrt(diag(an$cov)); se_num <- sqrt(diag(cov_num))
  cat(sprintf("%-6d %12.6f %12.6f %9.2f%%\n", 10000L + s, se_an, se_num,
              100 * abs(se_an - se_num) / se_num))
}
