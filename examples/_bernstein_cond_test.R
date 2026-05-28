# Conditioning wall of the GLOBAL Bernstein basis vs the (orthonormal) SVD basis.
# A degree-n Bernstein design matrix has n+1 columns; as n grows the columns go
# nearly linearly dependent, so the smallest singular value -> 0, the condition
# number blows up, and the numerical rank falls below n+1 (-> rank-deficient fit,
# the NA we saw). The SVD basis is orthonormal by construction: every singular
# value is 1, so it supplies as many well-conditioned, full-domain modes as the
# data supports.
source("examples/_joint_tune_harness.R")

bernstein_basis <- function(tt, degree) {
  x  <- as.vector(tt); x <- (x - min(x)) / (max(x) - min(x))
  Bb <- vapply(0:degree, function(k) choose(degree, k) * x^k * (1 - x)^(degree - k),
               numeric(length(x)))
  t(Bb)                                              # (degree+1) x mp1
}

tt <- seq(0, 10, length.out = 256)
cat("Global Bernstein design matrix (mp1 = 256):\n")
cat(sprintf("%-5s %-8s %-12s %-10s %-12s\n",
            "deg", "n+1", "cond(kappa)", "num_rank", "min_sing"))
for (d in c(8, 16, 24, 32, 40, 48, 56, 64)) {
  Bt  <- t(bernstein_basis(tt, d))                   # mp1 x (d+1)
  sv  <- svd(Bt)$d
  tol <- sv[1] * .Machine$double.eps * max(dim(Bt))
  cat(sprintf("%-5d %-8d %-12.2e %-10d %-12.2e\n",
              d, d + 1L, sv[1] / sv[length(sv)], sum(sv > tol), sv[length(sv)]))
}

cat("\nSVD basis (current default) on the same grids:\n")
for (nm in c("logistic", "lotka_volterra", "lorenz")) {
  sys <- SYSTEMS[[nm]]; mp1 <- if (nm == "lorenz") 256 else 128
  bj  <- build_joint(sys, mp1, nr = 0.05, seed = 1, n_bl = 20)
  B   <- wendy:::.joint_bases(bj$wobj$opt_ctx)$B     # Kb x mp1, orthonormal rows
  sv  <- svd(t(B))$d                                 # design-matrix singular values
  cat(sprintf("%-15s Kb=%-3d  cond=%.3e  min_sing=%.4f  max_sing=%.4f\n",
              nm, nrow(B), sv[1] / sv[length(sv)], min(sv), max(sv)))
}
