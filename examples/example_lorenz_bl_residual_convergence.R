# %%
# Validate O(h^6) convergence of the augmented boundary-layer residual.
#
# The math: for the true solution to the ODE and a smooth test function phi(t)
# supported on [t_1, t_1 + r_time] with phi(t_1) != 0 and phi vanishing at the
# right endpoint, the augmented residual in integral form is
#
#   r_aug = h*trap(phi*F) + h*trap(phi'*u) + phi(t_1)*u(t_1) - phi(t_M)*u(t_M)
#           - h^2/12  [g'(t_M)   - g'(t_1)]
#   h^4/720 [g'''(t_M) - g'''(t_1)]
#
# with g(t) = phi(t) F(p,u(t),t) + phi'(t) u(t). At (p*, U*) and any smooth phi,
# IBP says r_aug = 0 exactly in the continuous limit. The trapezoid + 2-term
# Euler-Maclaurin discretization picks up O(h^6) truncation error, so:
#
#   ||r_aug||_2  ~  O(h^6)   as h -> 0
#
# This script tests that empirically. It bypasses the SSL/MSG plumbing entirely
# and assembles r_aug from the helpers (build_boundary_layer_block,
# g_deriv_at_endpoint) using a single left-BL test function shape fixed in TIME.

library(deSolve)
library(ggplot2)
library(devtools)
invisible({ devtools::load_all("/Users/jackkrebsbach/Documents/ml/rwendy") })

f <- function(u, p, t) c(
  p[1] * (u[2] - u[1]),
  u[1] * (p[2] - u[3]) - u[2],
  u[1] * u[2] - p[3] * u[3]
)

p_star  <- c(10.0, 28.0, 8.0 / 3.0)
u0_true <- c(-8, 10, 27)
t_span  <- c(0, 10)
D       <- 3L
J_n     <- 3L

modelODE <- function(tvec, state, parms) list(as.vector(f(state, parms, tvec)))

u_expr <- do.call(c, lapply(seq_len(D),  function(i) symengine::S(paste0("u", i))))
p_expr <- do.call(c, lapply(seq_len(J_n),function(i) symengine::S(paste0("p", i))))
t_expr <- symengine::S("t")
vars   <- c(p_expr, u_expr, t_expr)

f_sym   <- f(u_expr, p_expr, t_expr)
dF_sym  <- compute_symbolic_total_time_deriv(f_sym,  u_expr, f_sym, t_expr)
d2F_sym <- compute_symbolic_total_time_deriv(dF_sym, u_expr, f_sym, t_expr)
d3F_sym <- compute_symbolic_total_time_deriv(d2F_sym,u_expr, f_sym, t_expr)
f_       <- build_fn(f_sym,   vars)
dF_dt_   <- build_fn(dF_sym,  vars)
d2F_dt2_ <- build_fn(d2F_sym, vars)
d3F_dt3_ <- build_fn(d3F_sym, vars)

r_time <- 1.0
p_psi  <- 16L
tau    <- symengine::S("tau")
psi_sym <- (1 - (tau / r_time)^2)^p_psi
psi_d <- function(order) {
  expr <- psi_sym
  for (k in seq_len(order)) expr <- symengine::D(expr, "tau")
  symengine::lambdify(expr, tau)
}
psi_fns <- lapply(0:4, psi_d)   # phi^(0..4) as functions of tau = t - t_1

trap_weights <- function(M) c(0.5, rep(1, M - 2L), 0.5)

eval_fd <- function(u_pt, t_pt) {
  input <- matrix(c(p_star, u_pt, t_pt), ncol = 1L)
  list(as.vector(f_(input)),
       as.vector(dF_dt_(input)),
       as.vector(d2F_dt2_(input)),
       as.vector(d3F_dt3_(input)))
}

# IMPORTANT: dt is chosen so r_time / dt is a positive integer. Otherwise the
# right edge of phi's support t_1 + r_time falls between two grid points, the
# `in_supp` mask drops/adds a partial sample as we sweep h, and the residual
# error gets dominated by that grid-alignment artifact instead of the EM
# truncation. With r_time = 1, dt = 1/N for integer N gives clean convergence.
n_per_rtime <- c(10L, 20L, 40L, 80L, 160L, 320L)        # samples per r_time
dt_vec      <- r_time / n_per_rtime
npoints_vec <- as.integer((t_span[2] - t_span[1]) / dt_vec + 1)
err <- numeric(length(dt_vec))
dts <- dt_vec

for (i in seq_along(dt_vec)) {
  npoints <- npoints_vec[i]
  t_eval  <- seq(t_span[1], t_span[2], length.out = npoints)
  dt      <- dt_vec[i]

  sol <- deSolve::ode(y = u0_true, times = t_eval, func = modelODE,
                      parms = p_star, rtol = 1e-15, atol = 1e-15)
  U  <- sol[, -1]                       # M x D
  tt <- sol[, 1]                        # M-vector
  M  <- length(tt)
  t1 <- tt[1]
  tM <- tt[M]

  in_supp <- tt >= t1 & tt <= t1 + r_time
  tau_in  <- tt[in_supp] - t1
  phi_grid  <- rep(0, M); phi_grid [in_supp] <- psi_fns[[1]](tau_in)
  phi1_grid <- rep(0, M); phi1_grid[in_supp] <- psi_fns[[2]](tau_in)

  input_full <- rbind(matrix(rep(p_star, M), nrow = J_n),
                     t(U),
                     matrix(tt, nrow = 1L))
  F_eval <- f_(input_full)              # M x D

  # Integral terms (trapezoid sums times h)
  w <- trap_weights(M)
  trap_phiF  <- dt * as.vector(t(F_eval) %*% (w * phi_grid))   # D-vector
  trap_phipU <- dt * as.vector(t(U)      %*% (w * phi1_grid))  # D-vector

  # Boundary terms: phi(t_1)*U[1,] - phi(t_M)*U[M,] (raw phi values)
  bdry <- phi_grid[1] * U[1, ] - phi_grid[M] * U[M, ]

  # EM correction: -h^2/12 * (g'(tM)-g'(t1)) + h^4/720 * (g'''(tM)-g'''(t1))
  # g = phi*F + phi'*u   —>   needs phi^(0..4) at the endpoints
  phi_at <- function(tau_pt) lapply(psi_fns, function(fn) fn(tau_pt))
  phi_t1 <- phi_at(0)                                  # tau = 0  -> peak
  phi_tM <- if (tM > t1 + r_time) {
    list(0, 0, 0, 0, 0)                               # outside support
  } else {
    phi_at(tM - t1)
  }
  fd_t1 <- eval_fd(U[1, ], t1)
  fd_tM <- eval_fd(U[M, ], tM)

  g1_t1 <- g_deriv_at_endpoint(phi_t1, fd_t1, U[1, ], order = 1L)
  g1_tM <- g_deriv_at_endpoint(phi_tM, fd_tM, U[M, ], order = 1L)
  g3_t1 <- g_deriv_at_endpoint(phi_t1, fd_t1, U[1, ], order = 3L)
  g3_tM <- g_deriv_at_endpoint(phi_tM, fd_tM, U[M, ], order = 3L)

  EM <- -(dt^2 / 12) * (g1_tM - g1_t1) + (dt^4 / 720) * (g3_tM - g3_t1)

  r_aug   <- trap_phiF + trap_phipU + bdry + EM      # D-vector
  err[i]  <- sqrt(sum(r_aug^2))

  cat(sprintf("dt = %.6f   ||trap_phiF|| = %.3e   ||trap_phipU|| = %.3e   ||bdry|| = %.3e   ||EM|| = %.3e   ||r_aug||_2 = %.3e\n",
              dt, sqrt(sum(trap_phiF^2)), sqrt(sum(trap_phipU^2)),
              sqrt(sum(bdry^2)), sqrt(sum(EM^2)), err[i]))
}

# Empirical slope (least-squares on log-log)
slope <- unname(coef(lm(log(err) ~ log(dts)))[2])
cat(sprintf("\nEmpirical slope: %.3f  (expected 6)\n", slope))

# Reference slope lines anchored at the first point
h2_ref <- err[1] * (dts / dts[1])^2
h4_ref <- err[1] * (dts / dts[1])^4
h6_ref <- err[1] * (dts / dts[1])^6

df     <- data.frame(dt = dts, error = err)
df_ref <- data.frame(
  dt    = rep(dts, 3),
  value = c(h2_ref, h4_ref, h6_ref),
  rate  = factor(rep(c("O(h^2)", "O(h^4)", "O(h^6)"), each = length(dts)),
                 levels = c("O(h^2)", "O(h^4)", "O(h^6)"))
)

ggplot() +
  geom_line(data = df_ref, aes(x = dt, y = value, linetype = rate),
            color = "gray40", alpha = 0.8, linewidth = 0.6) +
  geom_line(data = df, aes(x = dt, y = error), color = "firebrick", linewidth = 0.9) +
  geom_point(data = df, aes(x = dt, y = error),
             color = "firebrick", fill = "white", shape = 21, size = 3, stroke = 1.2) +
  scale_x_log10(breaks = dts, labels = scales::label_scientific()(dts)) +
  scale_y_log10(labels = scales::label_scientific()) +
  scale_linetype_manual(values = c("O(h^2)" = "dashed",
                                   "O(h^4)" = "dotted",
                                   "O(h^6)" = "longdash")) +
  annotation_logticks(sides = "bl") +
  labs(
    title    = sprintf("Augmented residual (integral form):  empirical slope %.2f", slope),
    subtitle = sprintf("Lorenz, p = p*, U = U_true (clean), r_time = %.2f", r_time),
    x        = "Step size  h",
    y        = expression("||" * r[aug] * "||"[2]),
    linetype = "reference"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom",
        panel.grid.minor = element_blank())
