{
library(deSolve)

source("./R/symbolics.R")
source("./R/test_functions.R")
source("./R/noise.R")
source("./R/optim.R")
source("./R/weak_residual.R")
source("./R/wendy.R")

f <- function(u, p, t) {
  c(p[1] * u[1] - p[2] * u[1]^2)
}

noise_sd <- 0.1
p_star <- c(1, 1);
u0 <- c(0.01);
p0 <- c(0.5, 0.5);
npoints <- 256;
t_span <- c(0.005, 10);
t_eval <- seq(t_span[1], t_span[2], length.out = npoints);

modelODE <- function(tvec, state, parameters) { list(as.vector(f(state, parameters, tvec))) }
sol <- deSolve::ode(y = u0, times = t_eval, func = modelODE, parms = p_star)

# Additive Gaussian Noise
U <- matrix(c(sol[, 2] + rnorm(npoints, mean = 0, sd = noise_sd)), ncol = 1)
tt <- matrix(sol[, 1], ncol = 1)

# Log Normal Noise
#noisy <- sol[, 2] * exp(rnorm(npoints, mean = 0, sd = noise_sd))
#U <- matrix(noisy, ncol = 1)
#tt <- matrix(sol[, 1], ncol = 1)
}

res <- solveWendy(f, p0, U, tt, method = "MLE", optimize = T)
#sol_hat <- deSolve::ode(u0, t_eval, modelODE, res$phat)

finite_J <- calc_gradient(p_star, res$wnll)
analytic_J <- res$J_wnll(p_star)

cat(finite_J)
cat(analytic_J)

finite_H <- calc_hessian(p_star, res$wnll)
analytic_H <- res$H_wnll(p_star)

cat(finite_H)
cat(analytic_H)

# Vp <- as.matrix(res$V_prime)
# U_state <- U
# svd_result <- svd(Vp)
# k <- length(svd_result$d)
# u_k <- svd_result$u[, 1:k]
# d_k <- svd_result$d[1:k]
# v_k <- svd_result$v[, 1:k]
# V_prime_pinv <- v_k %*% diag(1/d_k) %*% t(u_k)
#
# mp1 <- nrow(U)
# J <- length(p0)
#
# for(i in seq(1:1)){
#   p_mat <- matrix(rep(res$phat, mp1), ncol = mp1, nrow = J)
#   input <- rbind(p_mat, t(U_state), t(tt))
#   FU_ <- res$f(input)
#   U_state <- -1 * V_prime_pinv %*% (as.matrix(res$V) %*% FU_)
# }
#
# plot(U_state, cex = 0.5, col = "blue", ylim = c(-1,1.5))
# points(U, cex = 0.5)
# points(sol[, 2], cex = 0.5, col = "green")
#points(U, cex = 0.5, col = "red")

 # {
 # finite_g <- calc_gradient(p_star, res$wnll)
 # #finite_g <- gradient_4th_order(res$wnll, p_star)
 # J_p <- res$J_wnll(p_star)
 #
 # finite_h <- calc_hessian(p_star, res$wnll)
 # H_p <- res$H_wnll(p_star)
 #
 # cat("\n")
 # print(finite_g)
 # print(J_p)

#cat("\n")
#print(finite_h)
#print(H_p)

#cat("\n")
#print(res$phat)
#}




