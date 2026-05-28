# Can a global Bernstein polynomial basis represent the clean trajectories?
# Bernstein basis of degree n on [0,1]: B_{k,n}(x) = C(n,k) x^k (1-x)^(n-k),
# k = 0..n  (n+1 functions). Map t -> x = (t - t0)/T. We fit uhat = Bb^T c to the
# CLEAN data by least squares and report the relative Frobenius error.
source("examples/_joint_tune_harness.R")

# Bb: (degree+1) x mp1, row k+1 = B_{k,degree}(x_m).
bernstein_basis <- function(tt, degree) {
  x  <- as.vector(tt); x <- (x - min(x)) / (max(x) - min(x))   # -> [0,1]
  Bb <- vapply(0:degree, function(k) choose(degree, k) * x^k * (1 - x)^(degree - k),
               numeric(length(x)))
  t(Bb)                                                        # (degree+1) x mp1
}

# Relative Frobenius error of the least-squares Bernstein fit to clean data Uc.
bern_fit_rel <- function(Uc, Bb) {
  Bt     <- t(Bb)                                              # mp1 x (degree+1)
  fitted <- apply(Uc, 2, function(u) Bt %*% lm.fit(Bt, u)$coefficients)
  norm(fitted - Uc, "F") / norm(Uc, "F")
}

# Piecewise (local-support) Bernstein: n_seg contiguous segments, degree d each.
# Each of the (d+1) Bernstein functions per segment is zero outside its segment.
# Total basis size = n_seg * (d+1).
bernstein_basis_pw <- function(tt, degree, n_seg) {
  tt <- as.vector(tt); mp1 <- length(tt)
  edges <- unique(round(seq(1, mp1 + 1, length.out = n_seg + 1)))
  rows  <- vector("list", 0)
  for (s in seq_len(length(edges) - 1L)) {
    idx <- edges[s]:(edges[s + 1L] - 1L)
    xs  <- tt[idx]; xs <- (xs - min(xs)) / (max(xs) - min(xs) + 1e-12)
    Bseg <- vapply(0:degree, function(k) choose(degree, k) * xs^k * (1 - xs)^(degree - k),
                   numeric(length(idx)))
    for (k in seq_len(degree + 1L)) {
      r <- numeric(mp1); r[idx] <- Bseg[, k]; rows[[length(rows) + 1L]] <- r
    }
  }
  do.call(rbind, rows)
}

degrees <- c(8, 12, 16, 24, 32, 48)
cat("== Global Bernstein (one polynomial over [0,T]) ==\n")
cat(sprintf("%-15s %-5s | %s\n", "system", "mp1",
            paste(sprintf("deg=%-2d", degrees), collapse = " ")))
for (nm in c("logistic", "lotka_volterra", "lorenz")) {
  sys <- SYSTEMS[[nm]]; mp1 <- if (nm == "lorenz") 256 else 128
  Uc  <- gen_data(sys, mp1, nr = 0.05, seed = 1)$Uc            # clean trajectory
  tt  <- seq(sys$tspan[1], sys$tspan[2], length.out = mp1)
  errs <- sapply(degrees, function(d) bern_fit_rel(Uc, bernstein_basis(tt, d)))
  cat(sprintf("%-15s %-5d | %s\n", nm, mp1,
              paste(sprintf("%.4f", errs), collapse = "  ")))
}

cat("\n== Piecewise Bernstein (degree d per segment, n_seg segments) ==\n")
cat(sprintf("%-15s | %s\n", "system(deg,nseg->Kb)",
            "rel_err"))
configs <- list(c(4, 8), c(4, 16), c(4, 32), c(8, 8), c(8, 16), c(16, 8))
for (nm in c("logistic", "lotka_volterra", "lorenz")) {
  sys <- SYSTEMS[[nm]]; mp1 <- if (nm == "lorenz") 256 else 128
  Uc  <- gen_data(sys, mp1, nr = 0.05, seed = 1)$Uc
  tt  <- seq(sys$tspan[1], sys$tspan[2], length.out = mp1)
  for (cfg in configs) {
    d <- cfg[1]; ns <- cfg[2]
    Bb <- bernstein_basis_pw(tt, d, ns)
    cat(sprintf("%-15s d=%-2d nseg=%-2d Kb=%-3d  rel_err=%.4f\n",
                nm, d, ns, nrow(Bb), bern_fit_rel(Uc, Bb)))
  }
}
