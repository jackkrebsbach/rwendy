
# %%
library(wendy)
library(deSolve)
library(devtools)
library(ggplot2)

# u = c(S, I);  p = c(nu, mu, V)
# S0 = initial susceptible larvae (density?)

N0 <- 300  # cohort size: larvae in the population
S0 <- 250  # Number initial Suceptible

f  <- function(u, p, t) {
  S <- u[1]
  I <- u[2]
  nu <- p[1]
  mu <- p[2]
  V <- p[3]
  tr <- nu * S * I * (S / S0)^V
  c(-tr, tr - mu * I)
}

# true parameters 
p_star <- c(nu = 2.5/N0, mu = 0.4, V = 1.5) 

# initial condition
u0 <- c(S = S0, I = 1)                  

# ~15 weeks
n <- 40
Tmax <- 15                         

sol    <- pracma::rk4sys(function(t, u) f(u, p_star, t), 0, Tmax, u0, n - 1)
u_true <- sol$y                            
tt     <- matrix(sol$x, ncol = 1)

set.seed(1)
sig <- 0.05 * sqrt(colMeans(u_true^2))     # 5% additive Gaussian noise
U   <- u_true + cbind(rnorm(n, sd = sig[1]), rnorm(n, sd = sig[2]))

par(mar = c(4,4,1,1), col.axis = "#1B1D1E",
    col.lab = "#1B1D1E", fg = "#A2A4A3", cex = 1.2)
plot(tt, U[,1], pch = 16, col = adjustcolor("#CFB87C", .9),
     xlab = "t (weeks)", ylab = "larvae", ylim = range(U))
points(tt, U[,2], pch = 16, col = adjustcolor("#1D70A1", .9))
lines(tt, u_true[,1], col = "#CFB87C", lwd = 5)
lines(tt, u_true[,2], col = "#1D70A1", lwd = 5)
legend("right", bty = "n",
       legend = c("S (susceptible)", "I (infected)"),
       col = c("#CFB87C", "#1D70A1"), lwd = 5)


res <- solveWendy(f, U, tt, p0 = c(0.01, 1, 1), method = "MLE")  
phat <- res$phat
names(phat) <- paste0(names(p_star), "\u0302")

print(phat)
print(p_star)                         # truth

cat("\np̂ relative error: ")
cat(rel_err(res$phat, p_star))