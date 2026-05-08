
# %%
# library(wendy)
library(deSolve)
library(MASS)
library(uGMAR)
library(numDeriv)

invisible({devtools::load_all()})

f <- function(u, p, t) {
  u1 <- p[1] * u[1] + p[2] * u[1] * u[2]
  u2 <- p[3] * u[2] + p[4] * u[1] * u[2]
  c(u1, u2)
}

p_star <- c(1, -0.1, -1.5, 0.075)
u0 <- c(10,5)
p0 <- c(2, -0.1, -1, 0.25)
npoints <- 256
t_span <- c(0, 20)
t_eval <- seq(t_span[1], t_span[2], length.out = npoints);

modelODE <- function(tvec, state, parameters) { list(as.vector(f(state, parameters, tvec))) }
sol <- deSolve::ode(y = u0, times = t_eval, func = modelODE, parms = p_star)

# Additive Gaussian Noise
nr <- 0.5
U_vec <- as.vector(sol[,-1])
noise_sd <- nr * sqrt(mean(U_vec^2))
noise <- matrix(
  rnorm(nrow(sol) * (ncol(sol) - 1), mean = 0, sd = noise_sd),
  nrow = nrow(sol)
)
U <- sol[, -1] + noise
# U[,1] <- mean(U[,2])
tt <- matrix(sol[, 1], ncol = 1)

res <- solveWendy(f, U, tt, method = "ROOT")

Vp <- as.array(res$V_prime$contiguous())
V <- as.array(res$V$contiguous())

F_ <- function(U, p){
    mp1 <- nrow(U)
    J <- length(p)
    tU <- t(U)
    ttt   <- matrix(tt, nrow = 1L)
    p_mat <- matrix(rep(p, mp1), ncol = mp1, nrow = J)
    input <- rbind(p_mat, tU, ttt)
    Fp <- res$f(input)
    return(Fp)
}

# Root on just state
root_target <- function(U, p){
  U <- matrix(U, ncol = 2)
  as.vector(-Vp %*% U -  V %*% F_(U, p))
}

# U <- res$state$U_star

# u0_vec <- as.vector(U)
# max_iter <- 100
# tol <- 1e-8
# lambda <- 1e-8

# u <- u0_vec
# pn <- res$phat

# G <- root_target(u, pn)
# J <- jacobian(\(u) root_target(u, pn), u)

# for (iter in seq_len(max_iter)) {
#   cat("n=", iter, "||G||:", sqrt(sum(G^2)), sprintf("err(p̂): %.6f", rel_err(pn, p_star)), "\n")

#   JtJ  <- t(J) %*% J
#   delta <- solve(JtJ + lambda * diag(length(u)), t(J) %*% G)
#   u_new <- u - delta

#   pn    <- solveWendy(f, matrix(u_new, ncol = 2), tt, method = "IRLS", control = list(test_fun_type = "MSG", optimize = TRUE))$phat
#   G_new <- root_target(u_new, pn)

#   # Broyden rank-1 Jacobian update
#   s <- -delta
#   y <- G_new - G
#   J <- J + ((y - J %*% s) %*% t(s)) / as.numeric(t(s) %*% s)

#   u <- u_new
#   G <- G_new

#   if (max(abs(delta)) < tol) break
# }

# Un <- matrix(u, ncol = 2)
# Un <- sweep(Un, 2, colMeans(U - Un), "+")

# frob_rel_err <- function(A, B) norm(A - B, "F") / norm(B, "F")
# sol_true <- sol[, -1]

# cat(sprintf("\nIRLS p̂:  %s  rel_err=%.6f\n", paste(round(res$phat, 4), collapse=" "), rel_err(res$phat, p_star)))
# cat(sprintf("Broyden p̂: %s  rel_err=%.6f\n", paste(round(pn, 4), collapse=" "), rel_err(pn, p_star)))
# cat(sprintf("U_star  state frob_rel_err=%.6f\n", frob_rel_err(as.matrix(U), sol_true)))
# cat(sprintf("Un      state frob_rel_err=%.6f\n", frob_rel_err(Un, sol_true)))

# plot(tt, sol_true[,1], col = "blue", ylim = range(c(sol_true[,1], as.matrix(U)[,1], Un[,1])),
#      main = "State estimation", ylab = "u1")
# points(tt, sol[,2] + noise[,2],   col = "green",  pch = 1)
# points(tt, as.matrix(U)[,1], col = "black",  pch = 1)
# points(tt, Un[,1],           col = "purple", pch = 1)

# --- Joint Gauss-Newton on (u, p) simultaneously ---
# U <- res$state$U_star

# u0_vec <- as.vector(U)
# max_iter <- 100
# tol <- 1e-8
# lambda <- 1e-8

# u <- u0_vec
# pn <- res$phat

# nu <- length(u0_vec)
# np <- length(res$phat)

# joint_target <- function(theta) {
#   root_target(theta[seq_len(nu)], theta[nu + seq_len(np)])
# }

# theta   <- c(u0_vec, res$phat)
# G_j     <- joint_target(theta)
# J_joint <- jacobian(joint_target, theta)

# for (iter in seq_len(max_iter)) {
#   cat(sprintf("unweighted iter=%d ||G||=%.6f err(p̂)=%.6f\n",
#               iter, sqrt(sum(G_j^2)), rel_err(theta[nu + seq_len(np)], p_star)))

#   JtJ_j <- t(J_joint) %*% J_joint
#   delta  <- solve(JtJ_j + lambda * diag(length(theta)), t(J_joint) %*% G_j)
#   theta_new <- theta - delta

#   G_new   <- joint_target(theta_new)
#   s <- -delta; y <- G_new - G_j
#   J_joint <- J_joint + ((y - J_joint %*% s) %*% t(s)) / as.numeric(t(s) %*% s)

#   theta <- theta_new
#   G_j   <- G_new

#   if (max(abs(delta)) < tol) break
# }

# u_joint  <- theta[seq_len(nu)]
# p_joint  <- theta[nu + seq_len(np)]
# Un_joint <- matrix(u_joint, ncol = 2)

# frob_rel_err <- function(A, B) norm(A - B, "F") / norm(B, "F")
# U_star   <- as.array(res$state$U_star)
# sol_true <- sol[, -1]

# cat(sprintf("\np̂_IRLS:           %s  rel_err=%.6f\n", paste(round(res$phat,   4), collapse=" "), rel_err(res$phat,   p_star)))
# cat(sprintf("p̂_ROOT:    %s  rel_err=%.6f\n", paste(round(p_joint,   4), collapse=" "), rel_err(p_joint,   p_star)))
# cat(sprintf("RKS state err=%.6f\n", frob_rel_err(U_star,    sol_true)))
# cat(sprintf("WENDy ROOT err=%.6f\n",  frob_rel_err(Un_joint,  sol_true)))

# plot(tt, sol[,2], col = "blue", ylim = range(c(sol[,2], U[,1], Un_joint[,1])), main = "State estimation", ylab = "u1")
# points(tt, sol[,2] + noise[,2],   col = "green",  pch = 1)
# points(tt, Un_joint[,1],  col = "purple", pch = 1)
# points(tt, U[,1],  col = "black", pch = 1)

# -- (GLS Gauss-Newton using S(p)) --
U <- res$state$U_star

u0_vec <- as.vector(U)
max_iter <- 100
tol <- 1e-8
lambda <- 1e-8

u <- u0_vec
pn <- res$phat

nu <- length(u0_vec)
np <- length(res$phat)

joint_target <- function(theta) {
  root_target(theta[seq_len(nu)], theta[nu + seq_len(np)])
}

theta_w   <- c(u0_vec, res$phat)
G_w <- joint_target(theta_w)
J_w <- jacobian(joint_target, theta_w)

S_inv  <- solve(as.array(res$S(res$phat)$contiguous()))

for (iter in seq_len(max_iter)) {
  p_cur  <- theta_w[nu + seq_len(np)]

  cat(sprintf("iter=%d ||G||=%.6f err(p̂)=%.6f\n",
              iter, sqrt(sum(G_w^2)), rel_err(p_cur, p_star)))

  JtSiJ <- t(J_w) %*% S_inv %*% J_w
  delta  <- solve(JtSiJ + lambda * diag(length(theta_w)), t(J_w) %*% S_inv %*% G_w)
  theta_new <- theta_w - delta

  G_new <- joint_target(theta_new)
  s <- -delta; y <- G_new - G_w
  J_w <- J_w + ((y - J_w %*% s) %*% t(s)) / as.numeric(t(s) %*% s)

  theta_w <- theta_new
  G_w     <- G_new

  if (max(abs(delta)) < tol) break
}

u_w  <- theta_w[seq_len(nu)]
p_w  <- theta_w[nu + seq_len(np)]
Un_w <- matrix(u_w, ncol = 2)
sol_true <- sol[, -1]

cat(sprintf("\np̂_IRLS: %s  rel_err=%.8f\n", paste(round(res$phat, 4), collapse=" "), rel_err(res$phat, p_star)))
cat(sprintf("p̂_ROOT: %s   rel_err=%.8f\n", paste(round(p_w, 4), collapse=" "), rel_err(p_w, p_star)))
# cat(sprintf("MEAN(p̂_ROOT, p̂_IRLS):                   rel_err=%.8f\n", rel_err((res$phat + p_w) / 2, p_star)))

cat(sprintf("Û_RKS state err=%.6f\n", frob_rel_err(U,    sol_true)))
cat(sprintf("Û_ROOT err=%.6f\n", frob_rel_err(Un_w, sol_true)))

plot(tt, sol[,2], col = "blue", ylim = range(c(sol[,2], U[,1], Un_joint[,1])), main = "State estimation", ylab = "u1", cex=0.5)
points(tt, sol[,2] + noise[,2],   col = "green",  pch = 1, cex=0.5)
points(tt, Un_w[,1],  col = "purple", pch = 1, cex=0.5)
points(tt, U[,1],  col = "black", pch = 1)

# Vpinv <- ginv(Vp)
# Un <- U
# pn <- res$phat
# for(i in seq_along(20)){
#   Un_est <- -Vpinv %*% V %*% F_(Un, p = pn)
#   Un[,1] <- Un_est[,1]
#   res2 <- solveWendy(f, matrix(Un, ncol = 2), tt, method = "IRLS", control = list(test_fun_type = "MSG", optimize = TRUE))
#   pn <- res2$phat
# }

# Un <- Un + mean(Un - U)

# plot(tt, sol[,2], col = "blue", ylim = range(c(sol[,2], U[,1])))
# points(tt, U[,1], col = "green")
# points(tt, -Un[,1] + 10, col = "purple")

# sol_hat <- deSolve::ode(u0, t_eval, modelODE, res$phat)

# plot(U[,1],U[,2], cex = 0.5)
# points(sol_hat[,2], sol_hat[,3], cex = 0.5, col = "red")

# plot(res$wendy_problems[[1]]$min_radius_radii, res$wendy_problems[[1]]$min_radius_errors)
# abline(v = res$wendy_problems[[1]]$min_radius, col = "red")

# print(res$phat)