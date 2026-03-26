library(devtools); library(deSolve)
invisible(devtools::load_all())

# Lorenz system (lip=FALSE, 3 params)
f <- function(u, p, t) {
  c(p[1]*(u[2]-u[1]),
    u[1]*(p[2]-u[3])-u[2],
    u[1]*u[2]-p[3]*u[3])
}

set.seed(42)
p_star <- c(10.0, 28.0, 8/3)
u0 <- c(-8, 10, 27)
t_eval <- seq(0, 2, length.out = 64)   # short, fast
modelODE <- function(t, y, p) list(as.vector(f(y, p, t)))
sol <- deSolve::ode(y=u0, times=t_eval, func=modelODE, parms=p_star)
U  <- sol[,-1] + matrix(rnorm(nrow(sol)*3, sd=0.5), ncol=3)
tt <- sol[,1]

res <- solveWendy(f, p_star, U, tt, lip=FALSE, method="MLE",
                  control=list(optimize=FALSE))

L    <- res$L
Jp_L <- res$wendy_problems[[1]]$Jp_L
S    <- res$S
Jp_S <- res$Jp_S

h <- 1e-5
J <- 3; D <- 3

cat("=== Check Jp_L vs numerical ∂L/∂p_j ===\n")
L0 <- L(p_star)
for (j in seq_len(J)) {
  ej <- rep(0,J); ej[j] <- h
  dL_num <- (L(p_star + ej) - L(p_star - ej)) / (2*h)
  dL_ana <- Jp_L(p_star)[,,j]
  err <- max(abs(dL_ana - dL_num))
  cat(sprintf("  j=%d  max|∂L/∂p_j error|= %.3e\n", j, err))
}

cat("\n=== Check Jp_S vs numerical ∂S/∂p_j ===\n")
S0 <- S(p_star)
for (j in seq_len(J)) {
  ej <- rep(0,J); ej[j] <- h
  dS_num <- (S(p_star + ej) - S(p_star - ej)) / (2*h)
  dS_ana <- Jp_S(p_star)[,,j]
  err <- max(abs(dS_ana - dS_num))
  cat(sprintf("  j=%d  max|∂S/∂p_j error|= %.3e\n", j, err))
}

cat("\n=== Check J_wnll components at p_star ===\n")
Sp   <- S(p_star)
cholF <- chol(Sp)
r    <- res$g(p_star) - res$b
S_inv_r <- backsolve(cholF, forwardsolve(t(cholF), r))
J_rp <- res$Jp_r(p_star)   # K*D x J
J_Sp <- Jp_S(p_star)        # K*D x K*D x J

for (j in seq_len(J)) {
  J_S_j <- J_Sp[,,j]
  J_r_j <- J_rp[,j]
  tmp        <- J_S_j %*% S_inv_r
  prt0       <- 2.0 * sum(J_r_j * S_inv_r)
  prt1       <- -1.0 * sum(S_inv_r * tmp)
  fact       <- backsolve(cholF, forwardsolve(t(cholF), J_S_j))
  logDetPart <- sum(diag(fact))
  cat(sprintf("  j=%d  prt0=%.4f  prt1=%.4f  logDet=%.4f  total=%.4f\n",
              j, prt0, prt1, logDetPart, 0.5*(prt0+prt1+logDetPart)))
}

cat("\n=== Analytical vs numerical gradient ===\n")
cat("J_wnll(p_star):", res$J_wnll(p_star), "\n")
grad_num <- sapply(seq_len(J), function(j) {
  ej <- rep(0,J); ej[j] <- h
  (res$wnll(p_star+ej) - res$wnll(p_star-ej))/(2*h)
})
cat("numerical grad:", grad_num, "\n")
