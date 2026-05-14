# %%
library(deSolve)
invisible({ devtools::load_all() })


f <- function(u, p, t) {
  c(
    p[1] * (u[2] - u[1]),
    u[1] * (p[2] - u[3]) - u[2],
    u[1] * u[2] - p[3] * u[3]
  )
}

p_star <- c(10.0, 28.0, 8.0 / 3.0)
u0     <- c(-8, 10, 27)
t_eval <- seq(0, 10, length.out = 256)

modelODE <- function(tvec, state, parms) list(as.vector(f(state, parms, tvec)))

sol <- deSolve::ode(y = u0, times = t_eval, func = modelODE, parms = p_star,
                    rtol = 1e-12, atol = 1e-12)

nr       <- 0.15
noise_sd <- nr * sqrt(mean(as.array(sol[, -1])^2))

set.seed(8675309 + 1)
noise <- matrix(rnorm(nrow(sol) * 3, sd = noise_sd), nrow = nrow(sol))
U  <- sol[, -1] + noise
tt <- matrix(sol[, 1], ncol = 1)


res <- solveWendy(f, U, tt, method = "IRLS",
                  control = list(optimize = FALSE, estimate_u0 = FALSE))

sys       <- res$opt_ctx$system        # original:  closures take (p)
sys_joint <- res$system_joint          # joint:     closures take (U, p, tt) / (U, p, tt, b)

U_fit  <- res$U                        # interpolated U used inside the problems
tt_fit <- res$tt

b_joint <- sys_joint$build_compute_b(U_fit, tt_fit)
sig     <- res$sig

p_test <- p_star


val_orig  <- res$wnll(p_test)
val_joint <- sys_joint$wnll(U_fit, p_test, tt_fit, b_joint, sig)

cat("=== wnll ===\n")
cat(sprintf("  original : %.10f\n", val_orig))
cat(sprintf("  joint    : %.10f\n", val_joint))
cat(sprintf("  abs diff : %.2e\n\n", abs(val_orig - val_joint)))


grad_orig  <- res$J_wnll(p_test)
grad_joint <- sys_joint$J_wnll(U_fit, p_test, tt_fit, b_joint, sig)

cat("=== J_wnll ===\n")
cat(sprintf("  original : %s\n", paste(sprintf("%.6f", grad_orig),  collapse = ", ")))
cat(sprintf("  joint    : %s\n", paste(sprintf("%.6f", grad_joint), collapse = ", ")))
cat(sprintf("  max diff : %.2e\n\n", max(abs(grad_orig - grad_joint))))


hess_orig  <- res$H_wnll(p_test)
hess_joint <- sys_joint$H_wnll(U_fit, p_test, tt_fit, b_joint, sig)

cat("=== H_wnll ===\n")
cat("  original:\n");  print(round(hess_orig,  6))
cat("  joint:\n");     print(round(hess_joint, 6))
cat(sprintf("  max diff : %.2e\n\n", max(abs(hess_orig - hess_joint))))


eps <- 1e-5

# Numerical gradient of the original wnll
grad_fd <- numDeriv::grad(res$wnll, p_test)

cat("=== J_wnll vs finite-difference (original) ===\n")
cat(sprintf("  analytic : %s\n", paste(sprintf("%.6f", grad_orig), collapse = ", ")))
cat(sprintf("  fd       : %s\n", paste(sprintf("%.6f", grad_fd),   collapse = ", ")))
cat(sprintf("  max diff : %.2e\n\n", max(abs(grad_orig - grad_fd))))

# Numerical Hessian of the original wnll
hess_fd <- numDeriv::hessian(res$wnll, p_test)

cat("=== H_wnll vs finite-difference (original) ===\n")
cat("  analytic:\n"); print(round(hess_orig, 4))
cat("  fd:\n");       print(round(hess_fd,   4))
cat(sprintf("  max diff : %.2e\n\n", max(abs(hess_orig - hess_fd))))

# Numerical gradient of the joint wnll (wrapping out U/tt/b/sig)
wnll_joint_p <- function(p) sys_joint$wnll(U_fit, p, tt_fit, b_joint, sig)
grad_fd_joint <- numDeriv::grad(wnll_joint_p, p_test)

cat("=== J_wnll vs finite-difference (joint) ===\n")
cat(sprintf("  analytic : %s\n", paste(sprintf("%.6f", grad_joint),    collapse = ", ")))
cat(sprintf("  fd       : %s\n", paste(sprintf("%.6f", grad_fd_joint), collapse = ", ")))
cat(sprintf("  max diff : %.2e\n\n", max(abs(grad_joint - grad_fd_joint))))


# Validate J_u_wnll (analytic gradient w.r.t. U) against finite differences
grad_U_analytic <- sys_joint$J_u_wnll(U_fit, p_test, tt_fit, b_joint, sig)

wnll_joint_U <- function(U_vec) {
  U_mat <- matrix(U_vec, nrow = nrow(U_fit), ncol = ncol(U_fit))
  b_    <- sys_joint$build_compute_b(U_mat, tt_fit)
  sys_joint$wnll(U_mat, p_test, tt_fit, b_, sig)
}
grad_U_fd <- matrix(numDeriv::grad(wnll_joint_U, c(U_fit)),
                    nrow = nrow(U_fit), ncol = ncol(U_fit))

cat("=== J_u_wnll vs finite-difference ===\n")
cat(sprintf("  max abs diff : %.2e\n", max(abs(grad_U_analytic - grad_U_fd))))
cat(sprintf("  rel error    : %.2e\n\n", max(abs(grad_U_analytic - grad_U_fd)) /
                                          max(abs(grad_U_fd))))


# Validate J_sig_wnll (analytic gradient w.r.t. sig) against finite differences
sig_vec     <- as.numeric(sig$cpu())
grad_sig_analytic <- sys_joint$J_sig_wnll(U_fit, p_test, tt_fit, b_joint, sig)

wnll_joint_sig <- function(sig_r) {
  sig_ <- torch::torch_tensor(sig_r, dtype = torch::torch_float64())
  sys_joint$wnll(U_fit, p_test, tt_fit, b_joint, sig_)
}
grad_sig_fd <- numDeriv::grad(wnll_joint_sig, sig_vec)

cat("=== J_sig_wnll vs finite-difference ===\n")
cat(sprintf("  analytic : %s\n", paste(sprintf("%.6f", grad_sig_analytic), collapse = ", ")))
cat(sprintf("  fd       : %s\n", paste(sprintf("%.6f", grad_sig_fd),       collapse = ", ")))
cat(sprintf("  max diff : %.2e\n\n", max(abs(grad_sig_analytic - grad_sig_fd))))
