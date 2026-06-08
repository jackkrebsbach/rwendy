
# %%
# WSINDy model discovery: recover the ODE right-hand side from data alone,
# then refine the discovered model through the normal WENDy pipeline.
# Benchmarks follow MathBioCU/WSINDy_ODE (Messenger & Bortz, SIAM MMS 2021).
library(deSolve)

invisible({devtools::load_all()})

set.seed(8675309)

simulate <- function(f, p_star, u0, t_eval, nr) {
  modelODE <- function(tvec, state, parameters) { list(as.vector(f(state, parameters, tvec))) }
  sol <- deSolve::ode(y = u0, times = t_eval, func = modelODE, parms = p_star)
  U_clean <- sol[, -1, drop = FALSE]
  noise_sd <- nr * sqrt(mean(U_clean^2))
  U <- U_clean + matrix(rnorm(length(U_clean), sd = noise_sd), nrow = nrow(U_clean))
  list(U = U, tt = matrix(sol[, 1], ncol = 1), U_clean = U_clean)
}

# Map a wsindy fit onto the true (term label, state) -> coefficient support,
# returning the discovered coefficients in the order of `truth`.
support_coefs <- function(ws, truth) {
  vapply(seq_len(nrow(truth)), function(i) {
    ws$W[match(truth$label[i], rownames(ws$W)), truth$state[i]]
  }, numeric(1))
}

# 1. Lotka-Volterra (reference setup): u1' = 3*u1 - u1*u2, u2' = -6*u2 + u1*u2
f_lv <- function(u, p, t) {
  c(p[1] * u[1] + p[2] * u[1] * u[2],
    p[3] * u[2] + p[4] * u[1] * u[2])
}
p_lv <- c(3, -1, -6, 1)
truth_lv <- data.frame(label = c("u1", "u1*u2", "u2", "u1*u2"),
                       state = c(1, 1, 2, 2))

set.seed(8675309)
dat <- simulate(f_lv, p_lv, u0 = c(1, 1), t_eval = seq(0, 10, by = 0.001), nr = 0.1)

ws_lv <- solveWSINDy(dat$U, dat$tt)
print(ws_lv)

stopifnot(identical(unname(colSums(ws_lv$W != 0)), c(2, 2)))           # exact support
cat(sprintf("\nLV discovery rel_err = %.4f (noise 5%%)\n",
            rel_err(support_coefs(ws_lv, truth_lv), p_lv)))

# One-call discovery + WENDy refinement, vs the known-f fit
res_disc  <- solveWendy(U = dat$U, tt = dat$tt, method = "IRLS")       # f = NULL -> discovery
res_known <- solveWendy(f_lv, dat$U, dat$tt, method = "IRLS")

cat(sprintf("p̂ (discovered f, IRLS) = [%s]  rel_err = %.4f\n",
            paste(sprintf("%.4f", res_disc$phat), collapse = ", "),
            rel_err(res_disc$phat, p_lv)))
cat(sprintf("p̂ (known f,      IRLS) = [%s]  rel_err = %.4f\n",
            paste(sprintf("%.4f", res_known$phat), collapse = ", "),
            rel_err(res_known$phat, p_lv)))
stopifnot(!is.null(res_disc$wsindy))

# Consistency: generated closure reproduces Theta %*% W on data rows, and
# detect_n_params sees every p[k]
i <- c(1, 500, 1001)
recon <- t(vapply(i, function(m) ws_lv$f(dat$U[m, ], ws_lv$p0, dat$tt[m]), numeric(2)))
stopifnot(max(abs(recon - (ws_lv$Theta %*% ws_lv$W)[i, ])) < 1e-8)
stopifnot(length(ws_lv$p0) == sum(ws_lv$W != 0))

# 2. Lorenz (reference setup, x0 = (-8, 7, 27)): cross-term recovery
f_lorenz <- function(u, p, t) {
  c(p[1] * (u[2] - u[1]),
    u[1] * (p[2] - u[3]) - u[2],
    u[1] * u[2] - p[3] * u[3])
}
p_lorenz <- c(10, 28, 8 / 3)
# Full weight vector on the discovered support (7 terms)
truth_lorenz <- data.frame(
  label = c("u1", "u2", "u1", "u2", "u1*u3", "u1*u2", "u3"),
  state = c(1, 1, 2, 2, 2, 3, 3))
w_lorenz <- c(-10, 10, 28, -1, -1, 1, -8 / 3)

set.seed(8675309)
dat <- simulate(f_lorenz, p_lorenz, u0 = c(-8, 7, 27),
                t_eval = seq(0, 10, by = 0.002), nr = 0.05)

ws_lor <- solveWSINDy(dat$U, dat$tt)
print(ws_lor)

stopifnot(identical(unname(colSums(ws_lor$W != 0)), c(2, 3, 2)))       # 7-term support
cat(sprintf("\nLorenz discovery rel_err = %.4f (noise 5%%)\n",
            rel_err(support_coefs(ws_lor, truth_lorenz), w_lorenz)))

# 3. Logistic (D = 1): u' = u - 0.1*u^2
f_log <- function(u, p, t) c(p[1] * u[1] - p[2] * u[1]^2)
truth_log <- data.frame(label = c("u1", "u1^2"), state = c(1, 1))

set.seed(8675309)
dat <- simulate(f_log, c(1, 0.1), u0 = c(0.1), t_eval = seq(0, 10, by = 0.01), nr = 0.05)

ws_log <- solveWSINDy(dat$U, dat$tt)
print(ws_log)
stopifnot(identical(unname(colSums(ws_log$W != 0)), c(2)))
cat(sprintf("\nLogistic discovery rel_err = %.4f (noise 5%%)\n",
            rel_err(support_coefs(ws_log, truth_log), c(1, -0.1))))

# Rescaling off still recovers here
ws_log_raw <- solveWSINDy(dat$U, dat$tt, control = list(wsindy_rescale = FALSE))
stopifnot(identical(which(ws_log_raw$W != 0), which(ws_log$W != 0)))

# 4. MSTLS sanity on a synthetic system with known sparse solution
K <- 60; J <- 12
G_syn <- matrix(rnorm(K * J), K, J)
w_true <- numeric(J); w_true[c(2, 7)] <- c(1.5, -0.8)
fit <- wsindy_mstls(G_syn, as.vector(G_syn %*% w_true), rep(1, J),
                    10^seq(-4, 0, length.out = 100))
stopifnot(max(abs(fit$w - w_true)) < 1e-10)

# Duplicated extra term exercises the rank-deficient ginv fallback
ws_dup <- solveWSINDy(dat$U, dat$tt,
                      control = list(wsindy_extra_terms = c("u[1]^2", "u[1]^2")))
stopifnot(inherits(ws_dup, "wsindy"))

cat("\nAll WSINDy example checks passed.\n")
