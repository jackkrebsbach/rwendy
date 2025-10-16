# Logistic Example
{
set.seed(8675309)
library(RColorBrewer)
library(tidyverse)
library(deSolve)
library(symengine)
library(trust)
library(uGMAR)
library(wendy)
library(trustOptim)

source("./R/noise.R")
source("./R/symbolics.R")
source("./R/test_functions.R")
source("./R/weak_residual.R")
source("./R/wendy.R")

f <- function(u, p, t) { c(p[1] * u[1] - p[2] * u[1]^2) }

logistic <- function(u, p, t) {
  list(p[[1]] * u[[1]] - p[[2]] * u[[1]]^2)
}

noise_sd <- 0.05;
p_star <- c(1, 1);
u0 <- c(0.01);
p0 <- c(0.5, 0.5);
npoints <- 400;
t_span <- c(0, 10);
t_eval <- seq(t_span[1], t_span[2], length.out = npoints);

modelODE <- function(tvec, state, parameters) { list(as.vector(f(state, parameters, tvec))) }
sol <- deSolve::ode(y = u0, times = t_eval, func = modelODE, parms = p_star)
U <- matrix(c(sol[, 2] + rnorm(npoints, mean = 0, sd = noise_sd)), ncol = 1)
tt <- matrix(sol[, 1], ncol = 1)
#test_funs <- WendySolver(logistic, U, p0, tt, log_level = "info", compute_svd_ = T, optimize_ = FALSE)

res <- solveWendy(f, p0, U, tt, lip = T, optimize = F, compute_svd = T)

g <- res$g
G <- res$G
b <- res$b
b1 <- res$b1

p_hat <- solve(t(G) %*% G, t(G) %*% b1)
#p_hat <- p_star

covR <- res$S(p_hat)
eig <- eigen(covR, symmetric = TRUE)
r_w <- eig$vectors %*% diag(1/sqrt(eig$values)) %*% t(eig$vectors) %*% r

r <- (G %*% p_hat) - b1

covR<- res$S(p_hat)
cfact <- chol(covR)


dens_r   <- density(r)
dens_rw  <- density(r_w)

cols <- brewer.pal(4, "Set1")

plot(dens_r, main="OLS vs IRLS Residuals", xlab="Residual value", ylab="Density", lwd=2, col=cols[1], xlim=c(-5,5))
lines(dens_rw, col=cols[2], lwd=2, lty=4)
curve(dnorm(x, mean=0, sd=1), add=TRUE, col=cols[3], lwd=2)
legend("topright", legend=c("OLS Residuals", "IRLS Residuals", "Standard Normal"), col=cols[1:3], lwd=2)

}

{
# Lorenz Example
#set.seed(8675309)
library(deSolve)
library(symengine)
library(RColorBrewer)
library(tidyverse)
library(trust)
library(uGMAR)
library(wendy)
library(trustOptim)

source("./R/noise.R")
source("./R/symbolics.R")
source("./R/test_functions.R")
source("./R/weak_residual.R")
source("./R/wendy.R")

f <- function(u, p, t) {
  du1 <- p[1] * (u[2] - u[1])
  du2 <- u[1] * (p[2] - u[3]) - u[2]
  du3 <- u[1] * u[2] - p[3] * u[3]
  c(du1, du2, du3)
}

lorenz <- function(u, p, t) {
  du1 <- p[[1]] * (u[[2]] - u[[1]])
  du2 <- u[[1]] * (p[[2]] - u[[3]]) - u[[2]]
  du3 <- u[[1]] * u[[2]] - p[[3]] * u[[3]]
  list(du1, du2, du3)
}

noise_sd <- 0.05
p_star <- c(10.0, 28.0, 8.0 / 3.0)
p0 <- c(12.0, 21, 4.0)
u0 <- c(2, 1, 1)
npoints <- 100
t_span <- c(0, 10)
t_eval <- seq(t_span[1], t_span[2], length.out = npoints)

modelODE <- function(tvec, state, parameters) { list(as.vector(f(state, parameters, tvec))) }

sol <- deSolve::ode(y = u0, times = t_eval, func = modelODE, parms = p_star)

noise <- matrix(
  rnorm(nrow(sol) * (ncol(sol) - 1), mean = 0, sd = noise_sd),
  nrow = nrow(sol)
)

U <- sol[, -1] + noise
tt <- matrix(sol[, 1], ncol = 1)

res <- solveWendy(f, p0, U, tt, lip = T, optimize = F, compute_svd = T)

g <- res$g
G <- res$G
b <- res$b
b1 <- res$b1

p_hat <- solve(t(G) %*% G, t(G) %*% b1)
#p_hat <- p_star
r <- (G %*% p_hat) - b1

covR <- res$S(p_hat)
eig <- eigen(covR, symmetric = TRUE)
r_w <- eig$vectors %*% diag(1/sqrt(eig$values)) %*% t(eig$vectors) %*% r

r <- r / sd(r)
r_w <- r_w / sd(r_w)

dens_r   <- density(r)
dens_rw  <- density(r_w)

cols <- brewer.pal(4, "Set1")

plot(dens_r, main="OLS vs IRLS Residuals", xlab="Residual value", ylab="Density", lwd=2, col=cols[1], xlim=c(-5,5))
lines(dens_rw, col=cols[2], lwd=2, lty=4)
curve(dnorm(x, mean=0, sd=1), add=TRUE, col=cols[3], lwd=2)
legend("topright", legend=c("OLS Residuals", "IRLS Residuals", "Standard Normal"), col=cols[1:3], lwd=2)

}


#res$J_wnll(p0)
#calc_gradient(p0, res$wnll)

#res$H_wnll(p0)
#calc_hessian(p0, res$wnll)

#test_funs <- WendySolver(lorenz, U, p0, tt, log_level = "info", compute_svd_ = F, optimize_ = FALSE)
#V_r <- res$V

#V <- test_funs$V_prime






