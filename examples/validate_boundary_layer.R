# %%
# Validate the boundary-layer augmentation for JWENDy.
# Runs the joint system with control$include_boundary_layer = FALSE and TRUE
# on a Lorenz dataset and inspects the resulting matrices and residuals.

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
sol      <- deSolve::ode(y = u0, times = t_eval, func = modelODE, parms = p_star,
                         rtol = 1e-12, atol = 1e-12)

nr       <- 0.15
noise_sd <- nr * sqrt(mean(as.array(sol[, -1])^2))

set.seed(8675309 + 1)
noise <- matrix(rnorm(nrow(sol) * 3, sd = noise_sd), nrow = nrow(sol))
U  <- sol[, -1] + noise
tt <- matrix(sol[, 1], ncol = 1)

res_off <- solveWendy(f, U, tt, method = "IRLS",
                     control = list(optimize = FALSE, estimate_u0 = FALSE,
                                    test_fun_type = "SSL",
                                    include_boundary_layer = FALSE))

res_on  <- solveWendy(f, U, tt, method = "IRLS",
                     control = list(optimize = FALSE, estimate_u0 = FALSE,
                                    test_fun_type = "SSL",
                                    include_boundary_layer = TRUE))

sys_off <- res_off$system_joint
sys_on  <- res_on$system_joint

U_fit   <- res_off$U
tt_vec  <- as.vector(res_off$tt)
sig     <- res_off$sig

cat("=== Matrix shapes ===\n")
cat(sprintf("  BL off : K = %d  (K_interior=%d, K_bl=%d)\n",
            sys_off$K, sys_off$K, 0L))
prob_on <- res_on$wendy_problems_joint[[1]]
cat(sprintf("  BL on  : K = %d  (K_interior=%d, K_bl=%d)\n",
            sys_on$K, prob_on$bl$K_interior, prob_on$bl$K_bl))

b_off <- sys_off$build_compute_b(U_fit, tt_vec)
b_on  <- sys_on$build_compute_b(U_fit, tt_vec)

g_off <- sys_off$g(U_fit, p_star, tt_vec)
g_on  <- sys_on$g(U_fit, p_star, tt_vec)

r_off <- as.numeric((g_off - b_off)$cpu())
r_on  <- as.numeric((g_on  - b_on )$cpu())

K_off <- length(r_off) / 3L
K_on  <- length(r_on)  / 3L

cat("\n=== Residual norms at p_star ===\n")
cat(sprintf("  ||r_interior||      = %.4f   (K=%d)\n", sqrt(sum(r_off^2)), K_off))
cat(sprintf("  ||r_augmented||     = %.4f   (K=%d)\n", sqrt(sum(r_on^2)),  K_on))
cat(sprintf("  ||r_aug[interior]|| = %.4f   (first K_int rows of augmented r)\n",
            sqrt(sum(r_on[seq_len(prob_on$bl$K_interior * 3L)]^2))))
cat(sprintf("  ||r_aug[BL]||       = %.4f   (last K_bl rows of augmented r)\n",
            sqrt(sum(tail(r_on, prob_on$bl$K_bl * 3L)^2))))

cat("\n=== BL off identity check ===\n")

g_off_v <- as.numeric(g_off$cpu())
g_on_v  <- as.numeric(g_on$cpu())
g_on_interior <- g_on_v[seq_len(prob_on$bl$K_interior * 3L)]
cat(sprintf("  ||g_off - g_on[interior]|| = %.4e   (should NOT be ~0: SVDs differ)\n",
            sqrt(sum((g_off_v - g_on_interior)^2))))
cat("  (BL rows are appended after SVD so K_interior coincides between modes.)\n")

wnll_off  <- sys_off$wnll(U_fit, p_star, tt_vec, b_off, sig)
wnll_on   <- sys_on$wnll (U_fit, p_star, tt_vec, b_on,  sig)
J_wnll_off <- sys_off$J_wnll(U_fit, p_star, tt_vec, b_off, sig)
J_wnll_on  <- sys_on$J_wnll (U_fit, p_star, tt_vec, b_on,  sig)

cat("\n=== Objective and gradient ===\n")
cat(sprintf("  wnll      :  off = %.4f    on = %.4f\n", wnll_off, wnll_on))
cat(sprintf("  J_wnll off:  [%s]\n", paste(sprintf("%.4f", J_wnll_off), collapse = ", ")))
cat(sprintf("  J_wnll on :  [%s]\n", paste(sprintf("%.4f", J_wnll_on),  collapse = ", ")))

gU_off <- sys_off$J_u_wnll(U_fit, p_star, tt_vec, b_off, sig)
gU_on  <- sys_on$J_u_wnll (U_fit, p_star, tt_vec, b_on,  sig)

cat("\n=== J_u_wnll ===\n")
cat(sprintf("  shape off = %d x %d   shape on = %d x %d\n",
            nrow(gU_off), ncol(gU_off), nrow(gU_on), ncol(gU_on)))
cat(sprintf("  ||J_u_wnll off|| = %.4f   ||J_u_wnll on|| = %.4f\n",
            sqrt(sum(gU_off^2)), sqrt(sum(gU_on^2))))
cat(sprintf("  max |dU at t_1|  off = %.4f   on = %.4f\n",
            max(abs(gU_off[1, ])), max(abs(gU_on[1, ]))))
cat(sprintf("  max |dU at t_M|  off = %.4f   on = %.4f\n",
            max(abs(gU_off[nrow(gU_off), ])), max(abs(gU_on[nrow(gU_on), ]))))

wnll_p_on <- function(p) {
  b <- sys_on$build_compute_b(U_fit, tt_vec)
  sys_on$wnll(U_fit, p, tt_vec, b, sig)
}
J_wnll_fd <- numDeriv::grad(wnll_p_on, p_star)

cat("\n=== J_wnll vs finite-difference (BL on) ===\n")
cat(sprintf("  analytic : %s\n", paste(sprintf("%.6f", J_wnll_on), collapse = ", ")))
cat(sprintf("  fd       : %s\n", paste(sprintf("%.6f", J_wnll_fd), collapse = ", ")))
cat(sprintf("  max diff : %.2e   rel err : %.2e\n",
            max(abs(J_wnll_on - J_wnll_fd)),
            max(abs(J_wnll_on - J_wnll_fd)) / max(abs(J_wnll_fd))))
cat("  NOTE: analytic Jp_r captures d/dp of the integral term over the augmented\n")
cat("  V, but treats the EM correction's p-dependence (via f, dF_dt, ...) as\n")
cat("  constant — symmetric to the J_u_wnll boundary-term U-dependence drop.\n")
cat("  Expect a few-percent FD mismatch as long as the EM block is non-negligible.\n")

# Plot the boundary-layer test functions
prob   <- res_on$wendy_problems_joint[[1]]
V_mat  <- as.array(prob$V$cpu())
K_int  <- prob$bl$K_interior
K_bl   <- prob$bl$K_bl
bl_idx <- (K_int + 1):(K_int + K_bl)

# First left-BL row: peak at t_1, support sweeping right
plot(tt_vec, V_mat[K_int + 1, ], type = "l",
     xlab = "t", ylab = expression(phi[BL]),
     main = "First left-BL test function (peak at t_1)")

# All BL rows on one figure
matplot(tt_vec, t(V_mat[bl_idx, ]),
        type = "l", lty = 1,
        xlab = "t", ylab = expression(phi[BL]),
        main = sprintf("All %d BL test functions", K_bl))