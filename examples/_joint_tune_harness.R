# Shared harness for JOINT hyperparameter tuning across systems.
suppressMessages(invisible(devtools::load_all(".", quiet = TRUE)))
suppressMessages(library(deSolve))

# ---- system definitions --------------------------------------------------
sys_logistic <- list(
  name   = "logistic",
  f      = function(u, p, t) c(p[1] * u[1] - p[2] * u[1]^2),
  p_star = c(1, 1/10),
  u0     = c(0.1),
  tspan  = c(0, 10)
)
sys_lv <- list(
  name   = "lotka_volterra",
  f      = function(u, p, t) c(p[1]*u[1] + p[2]*u[1]*u[2], p[3]*u[2] + p[4]*u[1]*u[2]),
  p_star = c(1, -0.1, -1.5, 0.075),
  u0     = c(10, 5),
  tspan  = c(0, 20)
)
sys_lorenz <- list(
  name   = "lorenz",
  f      = function(u, p, t) c(p[1]*(u[2]-u[1]), u[1]*(p[2]-u[3]) - u[2], u[1]*u[2] - p[3]*u[3]),
  p_star = c(10, 28, 8/3),
  u0     = c(-8, 10, 27),
  tspan  = c(0, 10)
)
SYSTEMS <- list(logistic = sys_logistic, lotka_volterra = sys_lv, lorenz = sys_lorenz)

# ---- data generation -----------------------------------------------------
gen_data <- function(sys, mp1, nr, seed = 1) {
  t_eval <- seq(sys$tspan[1], sys$tspan[2], length.out = mp1)
  modelODE <- function(tv, st, pa) list(as.vector(sys$f(st, pa, tv)))
  sol <- deSolve::ode(y = sys$u0, times = t_eval, func = modelODE, parms = sys$p_star,
                      rtol = 1e-10, atol = 1e-10)
  Uc  <- sol[, -1, drop = FALSE]
  noise_sd <- nr * sqrt(mean(as.vector(Uc)^2))
  set.seed(seed)
  U <- Uc + matrix(rnorm(length(Uc), sd = noise_sd), nrow = nrow(Uc))
  list(U = U, tt = matrix(sol[, 1], ncol = 1), Uc = Uc, noise_sd = noise_sd, t_eval = t_eval)
}

# relative Frobenius state error (scale-invariant, comparable across systems)
rel_state <- function(Uhat, Uc) norm(Uhat - Uc, "F") / norm(Uc, "F")

# Build the JOINT system once (optimize = FALSE) and collect scaling diagnostics.
build_joint <- function(sys, mp1, nr, seed = 1, n_bl = 20) {
  d <- gen_data(sys, mp1, nr, seed)
  wobj <- solveWendy(sys$f, d$U, d$tt, method = "JOINT",
                     control = list(optimize = FALSE, n_bl = n_bl))
  # MSG-IRLS warm start + weak-residual variance scale
  ctx <- wobj$opt_ctx
  p_irls <- wendy:::.run_irls(ctx, NULL, ctx$system$G, ctx$system$b)$p
  Sd <- tryCatch(mean(diag(wobj$S(p_irls))), error = function(e) NA)
  list(wobj = wobj, d = d, p_irls = p_irls,
       diag = list(system = sys$name, mp1 = mp1, nr = nr, D = ncol(d$U),
                   K = nrow(wobj$V), sigma = as.numeric(ctx$sig)[1],
                   mean_diagS = Sd, sig_scale = sqrt(Sd)))
}

# Evaluate one (lambda, mu) via optimizeWendy (reuses the built system).
eval_joint <- function(bj, sys, lambda, mu = 0, n_bl = 20, thresh = 1e-3) {
  ctl <- list(joint_lambda = lambda, joint_l1 = mu, n_bl = n_bl, joint_l1_thresh = thresh)
  r <- tryCatch(optimizeWendy(bj$wobj, control = ctl),
                error = function(e) NULL)
  if (is.null(r)) return(NULL)
  list(lambda = lambda, mu = mu,
       rel_state = rel_state(r$data$uhat, bj$d$Uc),
       p_relerr  = rel_err(r$phat, sys$p_star),
       n_active  = r$data$n_active, n_coef = r$data$n_coef,
       outer = r$data$outer_iters)
}
