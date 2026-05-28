devtools::load_all(".", quiet = TRUE)

## Logistic: u' = p1*u - p2*u^2,  p_star = (1, 1/10)
f <- function(u, p, t) c(p[1] * u[1] - p[2] * u[1]^2)
p_star  <- c(1, 1/10)
u0      <- c(0.1)
npoints <- 128
t_eval  <- seq(0, 10, length.out = npoints)

modelODE <- function(tvec, state, parameters) list(as.vector(f(state, parameters, tvec)))
sol <- deSolve::ode(y = u0, times = t_eval, func = modelODE, parms = p_star)

set.seed(8675309 + 1)
U_vec    <- as.vector(sol[, -1])
noise_sd <- 0.05 * sqrt(mean(U_vec^2))
U  <- sol[, 2, drop = FALSE] + rnorm(npoints, sd = noise_sd)
tt <- sol[, 1, drop = FALSE]

# MSG residual + MSG-IRLS warm start (default test functions);
# SSL + BL basis with n_bl = 20.
ctrl <- list(n_bl = 20)

## --- 1. Build only, verify analytic Jacobian against numDeriv -------------
## (check both rho=0 and rho>0 so the derivative-penalty block is exercised)
res0 <- solveWendy(f, U, tt, method = "JOINT",
                   control = modifyList(ctrl, list(optimize = FALSE, joint_lambda = 1)))
ctx  <- res0$opt_ctx
bb   <- wendy:::.joint_bases(ctx)
B    <- bb$B; Bd <- bb$Bpp   # second-derivative (curvature) basis
p0_irls <- wendy:::.run_irls(ctx, NULL, ctx$system$G, ctx$system$b)$p
theta0  <- c(as.vector(B %*% ctx$U_sys), p0_irls)
cat(sprintf("Residual dim: K(MSG)=%d  Kb(SSL+BL basis)=%d  mp1=%d\n",
            nrow(ctx$V), nrow(B), nrow(ctx$U_sys)))
for (rho in c(0, 5)) {
  ctx$control$joint_deriv_pen <- rho
  fns  <- wendy:::.joint_make_fns(ctx, B, Bd)
  Jana <- fns$jacobian(theta0)
  Jnum <- numDeriv::jacobian(fns$residual, theta0)
  cat(sprintf("Jacobian (rho=%g): dims=%s  max|ana-num|=%.3e  rel=%.3e\n",
              rho, paste(dim(Jana), collapse="x"),
              max(abs(Jana - Jnum)), max(abs(Jana - Jnum)) / max(abs(Jnum))))
}
ctx$control$joint_deriv_pen <- 0

## --- 2. JOINT fit, sweep joint_lambda -------------------------------------
uclean <- sol[, 2]
cat(sprintf("\nMSG-IRLS p0 = (%.4f, %.4f)  rel_err=%.4f\n",
            p0_irls[1], p0_irls[2], rel_err(p0_irls, p_star)))
cat(sprintf("noisy-data state RMSE = %.4f\n\n", sqrt(mean((U[,1] - uclean)^2))))

for (lam in c(1, 1e2, 1e4, 1e6)) {
  res <- solveWendy(f, U, tt, method = "JOINT",
                    control = modifyList(ctrl, list(joint_lambda = lam)))
  cat(sprintf("lambda=%.0e  phat=(%.4f, %.4f)  rel_err=%.4f  uhat_RMSE=%.4f  conv=%s it=%s\n",
              lam, res$phat[1], res$phat[2], rel_err(res$phat, p_star),
              sqrt(mean((res$data$uhat[,1] - uclean)^2)),
              res$data$converged, res$data$iterations))
}

## --- 3. Hard-sparse L1 on coefficients (debiased, exact zeros) ------------
cat("\nHard-sparse L1 sweep (lambda = 1e4, debias on):\n")
for (mu in c(0, 1e-2, 1e-1, 1, 10)) {
  res <- solveWendy(f, U, tt, method = "JOINT",
                    control = modifyList(ctrl, list(joint_lambda = 1e4, joint_l1 = mu)))
  cc       <- as.vector(res$data$c)
  n_exact0 <- sum(cc == 0)
  cat(sprintf("mu=%-6g phat=(%.4f, %.4f)  rel_err=%.4f  uhat_RMSE=%.4f  active=%d/%d  exact_zeros=%d\n",
              mu, res$phat[1], res$phat[2], rel_err(res$phat, p_star),
              sqrt(mean((res$data$uhat[,1] - uclean)^2)),
              res$data$n_active, res$data$n_coef, n_exact0))
}

## debias on vs off at a fixed mu
cat("\ndebias on vs off (lambda = 1e4, mu = 1):\n")
for (db in c(TRUE, FALSE)) {
  res <- solveWendy(f, U, tt, method = "JOINT",
                    control = modifyList(ctrl, list(joint_lambda = 1e4, joint_l1 = 1,
                                                    joint_l1_debias = db)))
  cc <- as.vector(res$data$c)
  cat(sprintf("debias=%-5s active=%d/%d  exact_zeros=%d  uhat_RMSE=%.4f  rel_err=%.4f\n",
              db, res$data$n_active, res$data$n_coef, sum(cc == 0),
              sqrt(mean((res$data$uhat[,1] - uclean)^2)), rel_err(res$phat, p_star)))
}
