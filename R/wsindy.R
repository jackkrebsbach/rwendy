# WSINDy: Weak-form Sparse Identification of Nonlinear Dynamics for ODEs.
#
# Discovers the right-hand side u' = f(u) from data by sparse regression on a
# candidate library in the weak form, following Messenger & Bortz, "Weak SINDy:
# Galerkin-Based Data-Driven Model Selection" (SIAM MMS, 2021) and the
# reference implementation MathBioCU/WSINDy_ODE. Test functions are built from
# the package's psi primitive (test_functions.R).
#
# Per equation d (each state gets its own test-function support mt_d, pt_d):
#   Theta(U) : mp1 x J   candidate terms evaluated at the (scaled) data
#   G_d = V_d %*% Theta  : K_d x J
#   b_d = -(Vp_d %*% U[,d])
# MSTLS (modified sequential-thresholding least squares) then selects a sparse
# w_d with G_d w_d ~ b_d, with thresholds applied to physical-unit coefficients.

# Library generation ===========================================================

# All monomial powers beta in N^D with sum(beta) <= deg, constant row first,
# ordered by total degree then numeric lexicographic. J = choose(deg + D, D).
wsindy_poly_powers <- function(D, deg) {
  B <- as.matrix(expand.grid(rep(list(0:deg), D)))
  B <- B[rowSums(B) <= deg, , drop = FALSE]
  dimnames(B) <- NULL
  ord <- do.call(order, c(list(rowSums(B)), lapply(seq_len(D), function(d) B[, d])))
  B[ord, , drop = FALSE]
}

# Display name for a monomial powers vector, e.g. "1", "u1", "u1^2*u2".
wsindy_mono_label <- function(b) {
  if (all(b == 0)) return("1")
  parts <- character(0)
  for (d in seq_along(b)) {
    if (b[d] == 1) parts <- c(parts, paste0("u", d))
    else if (b[d] > 1) parts <- c(parts, paste0("u", d, "^", b[d]))
  }
  paste(parts, collapse = "*")
}

# Code name for a monomial powers vector, e.g. "1", "u[1]", "u[1]^2*u[2]".
wsindy_mono_code <- function(b) {
  if (all(b == 0)) return("1")
  parts <- character(0)
  for (d in seq_along(b)) {
    if (b[d] == 1) parts <- c(parts, paste0("u[", d, "]"))
    else if (b[d] > 1) parts <- c(parts, paste0("u[", d, "]^", b[d]))
  }
  paste(parts, collapse = "*")
}

# Build the candidate library: monomials up to total degree wsindy_poly_deg,
# optional trig terms sin(k*u_d)/cos(k*u_d) for k in wsindy_trig_freqs, and
# user-supplied wsindy_extra_terms (code strings like "u[1]*u[3]^2").
# Returns parallel structures indexed by term j = 1..J:
#   B     : J x D integer powers (trig/extra rows are all-zero placeholders)
#   kind  : "poly" | "trig" | "extra"
#   trig  : list, NULL or list(fun, k, d)
#   label : display name ("u1^2*u2", "sin(2*u1)")
#   code  : closure code ("u[1]^2*u[2]", "sin(2*u[1])")
wsindy_library <- function(D, control) {
  deg <- as.integer(control$wsindy_poly_deg)
  if (is.na(deg) || deg < 1) stop("control$wsindy_poly_deg must be an integer >= 1.")

  B <- wsindy_poly_powers(D, deg)
  J_poly <- nrow(B)
  kind  <- rep("poly", J_poly)
  trig  <- vector("list", J_poly)
  label <- apply(B, 1, wsindy_mono_label)
  code  <- apply(B, 1, wsindy_mono_code)

  for (k in control$wsindy_trig_freqs) {
    for (d in seq_len(D)) {
      for (fun in c("sin", "cos")) {
        arg_lab <- if (k == 1) paste0("u", d) else sprintf("%g*u%d", k, d)
        arg_cod <- if (k == 1) paste0("u[", d, "]") else sprintf("%g*u[%d]", k, d)
        B     <- rbind(B, rep(0L, D))
        kind  <- c(kind, "trig")
        trig  <- c(trig, list(list(fun = fun, k = k, d = d)))
        label <- c(label, sprintf("%s(%s)", fun, arg_lab))
        code  <- c(code, sprintf("%s(%s)", fun, arg_cod))
      }
    }
  }

  for (term in control$wsindy_extra_terms) {
    B     <- rbind(B, rep(0L, D))
    kind  <- c(kind, "extra")
    trig  <- c(trig, list(NULL))
    label <- c(label, gsub("u\\[(\\d+)\\]", "u\\1", term))
    code  <- c(code, term)
  }

  list(B = B, kind = kind, trig = trig, label = label, code = code,
       J = length(kind), D = D, deg = deg)
}

# Evaluate the library: Theta is mp1 x J, column j = term j at every grid
# point. Polynomial terms are evaluated on the SCALED data U/scale_x (the
# reference's build_theta); trig and extra terms on the raw data.
wsindy_theta <- function(U, tt, lib, scale_x = NULL) {
  mp1 <- nrow(U); D <- ncol(U); J <- lib$J
  Us <- if (is.null(scale_x)) U else sweep(U, 2, scale_x, "/")
  Theta <- matrix(1, mp1, J)
  for (j in seq_len(J)) {
    if (lib$kind[j] == "poly") {
      b <- lib$B[j, ]
      col <- rep(1, mp1)
      for (d in seq_len(D)) if (b[d] != 0) col <- col * Us[, d]^b[d]
      Theta[, j] <- col
    } else if (lib$kind[j] == "trig") {
      tg <- lib$trig[[j]]
      Theta[, j] <- if (tg$fun == "sin") sin(tg$k * U[, tg$d]) else cos(tg$k * U[, tg$d])
    } else {
      # extra: rewrite u[d] -> U[, d] and evaluate over the grid (t available)
      expr <- parse(text = gsub("u\\[(\\d+)\\]", "U[, \\1]", lib$code[j]))[[1]]
      Theta[, j] <- rep_len(eval(expr, list(U = U, t = as.vector(tt))), mp1)
    }
  }
  if (any(!is.finite(Theta))) {
    bad <- unique(lib$label[colSums(!is.finite(Theta)) > 0])
    stop("Non-finite values in the WSINDy library for term(s): ",
         paste(bad, collapse = ", "),
         ". Lower wsindy_poly_deg or rescale the data.")
  }
  Theta
}

# Test function support selection ==============================================
#
# Port of findcornerpts/findcorners (MathBioCU/WSINDy_ODE) for the piecewise
# polynomial test function psi(t; r, p) = (1 - (t/r)^2)^p.

# Critical wavenumber k* of one state: corner of the cumulative magnitude
# spectrum and of the log spectrum (averaged), via the two-piece-linear corner
# finder get_corner_index (itself a port of the reference's getcorner).
wsindy_corner_k <- function(x, tt) {
  T_ <- length(tt)
  half <- seq_len(ceiling(T_ / 2))
  mag <- Mod(fft(x))
  mag <- c(mag[(floor(T_ / 2) + 1):T_], mag[seq_len(floor(T_ / 2))])  # fftshift
  Ufft <- pmax(mag[half] / sqrt(2 * length(half)), .Machine$double.xmin)
  Umax <- which.max(Ufft)
  if (Umax < 3) return(1L)
  wn <- ((seq_len(T_) - 1) - floor(T_ / 2)) * (2 * pi) / (max(tt) - min(tt))
  xx <- wn[half]
  t1 <- get_corner_index(cumsum(Ufft[1:Umax]), xx[1:Umax])
  t2 <- get_corner_index(log(Ufft[1:Umax]), xx[1:Umax])
  tstarind <- floor((t1 + t2) / 2)
  max(Umax - tstarind + 1, 1)
}

# Support half-width mt and degree pt for one state (findcorners, phi_class 1):
# mt solves l(m) = log((2m-1)/m^2) (4 pi^2 k^2 m^2 - 3 T^2 tauhat^2)
#                  - 2 T^2 tauhat^2 log(tau) = 0,
# pt = max(2, floor(log(tau)/log((2 mt - 1)/mt^2))).
wsindy_tf_params <- function(x, tt, tau = 1e-16, tauhat = 1) {
  T_ <- length(tt)
  k <- wsindy_corner_k(x, tt)
  l_fun <- function(m) {
    log((2 * m - 1) / m^2) * (4 * pi^2 * k^2 * m^2 - 3 * T_^2 * tauhat^2) -
      2 * T_^2 * tauhat^2 * log(tau)
  }
  mnew <- tryCatch(stats::uniroot(l_fun, c(1 + 1e-8, 2 / sqrt(tau)))$root,
                   error = function(e) T_ / 2 / k)
  if (mnew > T_ / 2 - 1) mnew <- T_ / 2 / k
  mt <- max(2L, min(floor((T_ - 1) / 2), as.integer(ceiling(mnew))))
  pt <- max(2, floor(log(tau) / log((2 * mt - 1) / mt^2)))
  list(mt = mt, pt = pt, k = k)
}

# Test function rows: same row construction as build_test_function_matrix
# (psi values / derivative on the support grid, L2-normalized), but placing at
# most max_K uniformly strided centers so large datasets stay tractable.
# Returns V (order 0) and Vp (order 1) together, K x mp1 each.
wsindy_tf_pair <- function(test_function, tt, radius, max_K = Inf) {
  len_tt <- length(tt)
  dt <- mean(diff(tt))
  if ((2 * radius + 1) > len_tt - 2) radius <- floor((len_tt - 3) / 2)
  diameter <- 2L * as.integer(radius) + 1L

  r <- dt * radius
  lin <- seq(-r, r, length.out = diameter)
  xx <- lin[2:(diameter - 1)]
  norm_factor <- sqrt(sum(test_function(xx, r)^2))
  row0 <- c(0, test_function(xx, r) / norm_factor, 0)
  ph1  <- test_function_derivative(test_function, radius, dt, 1L)
  row1 <- c(0, ph1(xx) / norm_factor, 0)

  starts_all <- 2:(len_tt - diameter)
  stride <- max(1L, as.integer(ceiling(length(starts_all) / max_K)))
  starts <- starts_all[seq(1, length(starts_all), by = stride)]
  K <- length(starts)

  V  <- matrix(0, K, len_tt)
  Vp <- matrix(0, K, len_tt)
  for (i in seq_len(K)) {
    cols <- starts[i]:(starts[i] + diameter - 1L)
    V[i, cols]  <- row0
    Vp[i, cols] <- row1
  }
  list(V = V, Vp = Vp, radius = as.integer(radius), K = K)
}

# MSTLS sparse regression ======================================================

# Least squares with rank-deficiency fallback (min-norm via ginv). With
# gamma > 0 this solves the l2-regularized (ridge) problem
#   argmin_w ||A w - y||^2 + gamma^2 ||w||^2  =  (A'A + gamma^2 I)^-1 A'y,
# the gamma -> 0 limit of the paper's GLS objective (Algorithm step 4/5). The
# ridge is what makes near-collinear libraries identifiable (sin vs polynomial
# terms for the pendulum; the 252-term Linear-5D library) under noise.
wsindy_lstsq <- function(A, y, gamma = 0) {
  if (gamma > 0) {
    AtA <- crossprod(A)
    diag(AtA) <- diag(AtA) + gamma^2
    out <- tryCatch(solve(AtA, crossprod(A, y)), error = function(e) NULL)
    if (is.null(out) || any(!is.finite(out))) out <- as.vector(MASS::ginv(A) %*% y)
    return(as.vector(out))
  }
  out <- tryCatch(qr.solve(A, y, tol = 1e-12), error = function(e) NULL)
  if (is.null(out) || any(!is.finite(out))) out <- as.vector(MASS::ginv(A) %*% y)
  out
}

# Modified sequential-thresholding least squares for one equation, following
# sparsifyDynamics + wsindy_pde_RGLS_seq (MathBioCU/WSINDy_ODE):
#  - G has columns in SCALED units; M_diag converts scaled coefficients to
#    physical units (w_phys = M_diag * w_scaled)
#  - thresholds act on PHYSICAL coefficients with
#    bnds_j = ||b|| / ||g_j|| * M_j, keep lambda*max(1,bnds) <= |w| <= min(1,bnds)/lambda
#  - loss(lambda) = 2*alpha*||G(w_l - w_ls)||/||G w_ls|| + 2*(1-alpha)*nnz/J,
#    minimized over the lambda grid (smallest minimizer); alpha = 0.5 is the
#    published loss ||G(w - w_ls)||/||G w_ls|| + nnz/J
wsindy_mstls <- function(G, b, M_diag, lambdas, alpha = 0.5, gamma = 0) {
  J <- ncol(G)
  degenerate <- list(w = numeric(J), lambda_star = NA_real_,
                     losses = rep(NA_real_, length(lambdas)),
                     w_ls = numeric(J), iterations = 0L, degenerate = TRUE)

  bnorm <- sqrt(sum(b^2))
  if (bnorm < .Machine$double.eps) return(degenerate)

  wls_s <- wsindy_lstsq(G, b, gamma)           # scaled-unit (ridge) LS solution
  Gwls <- sqrt(sum((G %*% wls_s)^2))
  if (Gwls < .Machine$double.eps) { degenerate$w_ls <- M_diag * wls_s; return(degenerate) }

  gnorm <- sqrt(colSums(G^2))
  gnorm[gnorm < .Machine$double.eps] <- Inf    # dead columns can never enter
  bnds <- bnorm / gnorm * M_diag

  iterate <- function(lambda) {
    LB <- lambda * pmax(1, bnds)
    UB <- (1 / lambda) * pmin(1, bnds)
    w <- M_diag * wls_s                        # physical units throughout
    small <- rep(FALSE, J)
    its <- 0L
    for (j in seq_len(J)) {
      small_new <- (abs(w) < LB) | (abs(w) > UB)
      its <- j
      if (identical(small_new, small)) break
      small <- small_new
      w[small] <- 0
      keep <- which(!small)
      if (length(keep) == 0L) break
      w[keep] <- M_diag[keep] * wsindy_lstsq(G[, keep, drop = FALSE], b, gamma)
    }
    list(w = w, its = its)
  }

  results <- lapply(lambdas, iterate)
  losses <- vapply(results, function(r) {
    2 * alpha * sqrt(sum((G %*% (r$w / M_diag - wls_s))^2)) / Gwls +
      2 * (1 - alpha) * sum(r$w != 0) / J
  }, numeric(1))
  istar <- which.min(losses)                   # first = smallest minimizing lambda

  list(w = results[[istar]]$w, lambda_star = lambdas[istar], losses = losses,
       w_ls = M_diag * wls_s, iterations = results[[istar]]$its,
       degenerate = FALSE)
}

# Discovered model -> f(u, p, t) ==============================================

# Generate the closure, p0, and the parameter map in ONE state-major loop so
# the p[k] indices, p0 order, and param_map cannot drift apart. Every selected
# term is emitted as p[k]*theta_k (bare p[k] for the constant), so F(0) = 0 and
# the linear-in-p WENDy path (b = b_raw - g0 with g0 = 0) stays consistent.
wsindy_make_f <- function(W, lib) {
  D <- ncol(W)
  k <- 0L
  components <- character(D)
  p0 <- numeric(0)
  pmap <- list()

  for (d in seq_len(D)) {
    nz <- which(W[, d] != 0)
    if (length(nz) == 0L) {
      components[d] <- "0"
      next
    }
    pieces <- character(length(nz))
    for (i in seq_along(nz)) {
      j <- nz[i]
      k <- k + 1L
      pieces[i] <- if (identical(lib$code[j], "1")) {
        sprintf("p[%d]", k)
      } else {
        sprintf("p[%d]*%s", k, lib$code[j])
      }
      p0 <- c(p0, W[j, d])
      pmap[[k]] <- data.frame(param = k, state = d, term = j,
                              label = lib$label[j], coef = W[j, d],
                              stringsAsFactors = FALSE)
    }
    components[d] <- paste(pieces, collapse = " + ")
  }

  fn_txt <- sprintf("function(u, p, t) c(%s)", paste(components, collapse = ", "))
  f <- eval(parse(text = fn_txt), envir = asNamespace("wendy"))

  list(f = f, p0 = p0,
       param_map = if (length(pmap)) do.call(rbind, pmap) else NULL)
}

# Symbolic RHS with the numeric coefficients embedded (full precision), for
# format_ode_system / summary and the WENDy hand-off display. Built by simply
# evaluating the generated closure on symbolic states with numeric p0.
wsindy_make_f_sym <- function(f, p0, D) {
  u_expr <- sym_build(lapply(seq_len(D), function(i) sym_symbol(paste0("u", i))))
  t_expr <- sym_symbol("t")
  f(u_expr, p0, t_expr)
}

# Entry point ==================================================================

#' Discover ODE structure with Weak-form SINDy (WSINDy)
#'
#' Identifies the right-hand side of an ODE system \eqn{\dot{u} = f(u)} from
#' noisy time-series data using the weak-form Sparse Identification of
#' Nonlinear Dynamics algorithm of Messenger and Bortz (SIAM Multiscale Model.
#' Simul., 2021), following the reference implementation
#' \code{MathBioCU/WSINDy_ODE}. A library of candidate terms (monomials up to a
#' total degree, optionally trigonometric and custom terms) is integrated
#' against compactly supported test functions \eqn{\psi(t) = (1-(t/r)^2)^p}
#' whose support radius and degree are selected per state from the data's
#' Fourier spectrum; per-equation modified sequential-thresholding least
#' squares (MSTLS) then selects a sparse coefficient vector. The discovered
#' terms are assembled into an \code{f(u, p, t)} closure and starting values
#' \code{p0} ready for [solveWendy()].
#'
#' Assumes additive Gaussian measurement noise on a (close to) uniform time
#' grid; discovery quality improves with sampling density (the reference
#' examples use \eqn{\Delta t} of 0.01 to 0.001). For calibrated estimates and
#' a parameter covariance, call \code{solveWendy(f = NULL, U, tt, ...)}, which
#' runs this discovery and the WENDy refinement in one step.
#'
#' @param U Numeric matrix. Rows are observed states at the time points in
#'   \code{tt}; columns are state variables.
#' @param tt Numeric vector (or one-column matrix). Time points corresponding
#'   to the rows of \code{U}.
#' @param control Named list of control parameters, as produced by
#'   [wendy_control()]; the \code{wsindy_*} entries configure the candidate
#'   library, test-function selection, and the MSTLS sweep. See
#'   [wendy_control()].
#'
#' @return An object of class \code{wsindy} with elements including \code{W}
#'   (term-by-state coefficient matrix with term labels as row names), \code{f}
#'   (generated \code{f(u, p, t)} closure referencing \code{p[1..J]}),
#'   \code{p0} (nonzero coefficients in closure order), \code{f_sym} (symbolic
#'   RHS with numeric coefficients), \code{param_map}, per-equation
#'   test-function parameters (\code{mt}, \code{pt}, \code{K}) and
#'   \code{lambda_star}/\code{losses} over the lambda grid, the per-equation
#'   weak systems \code{G}/\code{b}, the scales \code{scale_x}/\code{M_diag},
#'   and residual diagnostics.
#'
#' @examples
#' # Lotka-Volterra: u1' = u1 - 0.1*u1*u2, u2' = -1.5*u2 + 0.075*u1*u2
#' f <- function(u, p, t) {
#'   c(p[1] * u[1] + p[2] * u[1] * u[2],
#'     p[3] * u[2] + p[4] * u[1] * u[2])
#' }
#' p_star <- c(1, -0.1, -1.5, 0.075)
#' u0     <- c(10, 5)
#' tt     <- seq(0, 10, by = 0.01)
#'
#' modelODE <- function(tvec, state, parameters) {
#'   list(as.vector(f(state, parameters, tvec)))
#' }
#' sol <- deSolve::ode(y = u0, times = tt, func = modelODE, parms = p_star)
#'
#' set.seed(8675309)
#' U_clean  <- sol[, -1]
#' noise_sd <- 0.05 * sqrt(mean(U_clean^2))
#' U <- U_clean + matrix(rnorm(length(U_clean), sd = noise_sd), nrow = length(tt))
#'
#' ws <- solveWSINDy(U, tt)
#' ws          # discovered equations
#' ws$p0       # coefficients, ready as p0 for solveWendy(ws$f, U, tt, ...)
#'
#' @export
solveWSINDy <- function(U, tt, control = NULL) {
  if (is.null(control)) {
    control <- default_control
  } else if (!all(names(default_control) %in% names(control))) {
    # Raw user list (standalone path): validate and merge. A fully-merged list
    # (handed over by solveWendy) is used as-is to avoid double-merging.
    unknown <- setdiff(names(control), names(default_control))
    if (length(unknown) > 0) {
      stop("Unknown control parameter(s): ", paste(unknown, collapse = ", "),
           ". See ?wendy_control for valid options.")
    }
    control <- modifyList(default_control, control)
  }

  U  <- as.matrix(U)
  tt <- as.vector(tt)
  if (nrow(U) != length(tt)) stop("nrow(U) must equal length(tt).")
  D   <- ncol(U)
  mp1 <- nrow(U)
  if (mp1 < 10) stop("WSINDy needs at least 10 time points.")

  lambdas <- control$wsindy_lambdas
  if (is.null(lambdas)) lambdas <- 10^seq(-4, 0, length.out = 100)
  alpha <- control$wsindy_alpha_loss

  # Ridge strength gamma. "auto" estimates the noise ratio sigma_NR from the
  # data (high-order FD filter, estimate_std) and sets gamma = sqrt(sigma_NR),
  # the paper's prescription gamma = sqrt(sigma_NR) (sec 3.3). On (near-)clean
  # data the estimate is ~0 and gamma is floored to 0, preserving the
  # machine-precision noise-free recovery; under noise it regularizes the
  # near-collinear / large libraries that otherwise collapse (pendulum, Lorenz,
  # Linear-5D). A numeric gamma is used verbatim.
  gamma <- control$wsindy_gamma
  sigma_NR_hat <- NA_real_
  if (is.character(gamma) && identical(gamma, "auto")) {
    sigma_hat <- tryCatch(estimate_std(U), error = function(e) rep(0, D))
    rms <- sqrt(colMeans(U^2)); rms[rms < .Machine$double.eps] <- Inf
    sigma_NR_hat <- mean(sigma_hat / rms)
    gamma <- if (is.finite(sigma_NR_hat) && sigma_NR_hat > 1e-3) sqrt(sigma_NR_hat) else 0
  } else if (is.null(gamma) || !is.numeric(gamma) || !is.finite(gamma) || gamma < 0) {
    gamma <- 0
  }

  lib <- wsindy_library(D, control)

  # Test-function parameters per state from the data spectrum
  tfp <- lapply(seq_len(D), function(d) {
    wsindy_tf_params(U[, d], tt, tau = control$wsindy_tau, tauhat = control$wsindy_tauhat)
  })

  # State scales from the raw data: scale_x_n = ||u_n^deg||^(1/deg)
  # (scale_Theta = 2 in the reference); M_diag maps scaled -> physical coefs.
  scale_x <- rep(1, D)
  if (isTRUE(control$wsindy_rescale)) {
    for (n in seq_len(D)) {
      s <- sum(U[, n]^(2 * lib$deg))^(1 / (2 * lib$deg))
      if (is.finite(s) && s > .Machine$double.eps) scale_x[n] <- s
    }
  }
  M_diag <- vapply(seq_len(lib$J), function(j) {
    if (lib$kind[j] == "poly") 1 / prod(scale_x^lib$B[j, ]) else 1
  }, numeric(1))

  Theta_s <- wsindy_theta(U, tt, lib, scale_x = scale_x)

  # Per-equation weak system + MSTLS
  W <- matrix(0, lib$J, D, dimnames = list(lib$label, paste0("u", seq_len(D))))
  per_eq <- vector("list", D)
  G_list <- vector("list", D)
  b_list <- vector("list", D)
  mt <- integer(D); pt <- numeric(D); k_corner <- integer(D); K <- integer(D)
  rel_residual <- numeric(D); r2 <- numeric(D); cond_G <- numeric(D)
  warnings <- character(0)

  for (d in seq_len(D)) {
    pt_d <- tfp[[d]]$pt
    psi_p <- function(t, r) psi(t, r, p = pt_d)
    tfr <- wsindy_tf_pair(psi_p, tt, tfp[[d]]$mt, max_K = control$wsindy_K)
    G_d <- tfr$V %*% Theta_s
    b_d <- as.vector(-(tfr$Vp %*% U[, d]))

    # Optional GLS whitening with Cov = (1-a) Vp Vp' + a diag(Vp Vp')
    # (the ODE paper's leading-order weak-residual covariance)
    a_gls <- control$wsindy_gls
    if (is.numeric(a_gls) && a_gls > 0) {
      Cov <- tcrossprod(tfr$Vp)
      Cov <- (1 - a_gls) * Cov + a_gls * diag(diag(Cov))
      R <- tryCatch(chol(Cov), error = function(e) NULL)
      if (!is.null(R)) {
        G_d <- forwardsolve(t(R), G_d)
        b_d <- forwardsolve(t(R), b_d)
      } else {
        warnings <- c(warnings, sprintf("u%d: GLS covariance not positive definite; using OLS", d))
      }
    }

    res <- wsindy_mstls(G_d, b_d, M_diag, lambdas, alpha = alpha, gamma = gamma)
    W[, d] <- res$w
    per_eq[[d]] <- res
    G_list[[d]] <- G_d
    b_list[[d]] <- b_d
    mt[d] <- tfr$radius; pt[d] <- pt_d; k_corner[d] <- tfp[[d]]$k; K[d] <- tfr$K

    # Diagnostics on the system MSTLS saw (scaled coefficients w/M)
    r_d <- G_d %*% (res$w / M_diag) - b_d
    bn <- sqrt(sum(b_d^2))
    rel_residual[d] <- if (bn > 0) sqrt(sum(r_d^2)) / bn else NA_real_
    ss <- sum((b_d - mean(b_d))^2)
    r2[d] <- if (ss > 0) 1 - sum(r_d^2) / ss else NA_real_
    sv <- svd(G_d, nu = 0, nv = 0)$d
    cond_G[d] <- if (min(sv) > 0) max(sv) / min(sv) else Inf

    nnz <- sum(res$w != 0)
    if (res$degenerate) {
      warnings <- c(warnings, sprintf("u%d: weak derivative ~ 0 (state near equilibrium); equation set to 0", d))
    } else if (nnz == 0) {
      warnings <- c(warnings, sprintf("u%d: no terms selected", d))
    } else {
      imin <- which.min(res$losses)
      if (imin %in% c(1L, length(res$losses))) {
        warnings <- c(warnings, sprintf("u%d: lambda* at grid boundary; consider widening wsindy_lambdas", d))
      }
      if (is.finite(rel_residual[d]) && rel_residual[d] > 0.5) {
        warnings <- c(warnings, sprintf("u%d: poor weak-form fit (rel residual %.2f)", d, rel_residual[d]))
      }
    }
    if (K[d] < lib$J) {
      warnings <- c(warnings, sprintf(
        "u%d: K = %d test functions < J = %d library terms; identification underdetermined", d, K[d], lib$J))
    }
  }

  if (all(W == 0)) {
    warning("WSINDy selected no terms in any equation; the discovered model is empty. ",
            "Inspect $losses, or adjust wsindy_poly_deg / wsindy_lambdas.")
  }

  # Discovered model as a WENDy-ready closure + symbolic form
  model <- wsindy_make_f(W, lib)
  f_sym <- wsindy_make_f_sym(model$f, model$p0, D)

  res <- structure(
    list(
      W = W,
      terms = data.frame(index = seq_len(lib$J), label = lib$label,
                         code = lib$code, kind = lib$kind,
                         stringsAsFactors = FALSE),
      f = model$f,
      p0 = model$p0,
      f_sym = f_sym,
      param_map = model$param_map,
      lambda_star = vapply(per_eq, `[[`, numeric(1), "lambda_star"),
      losses = lapply(per_eq, `[[`, "losses"),
      lambdas = lambdas,
      w_ls = lapply(per_eq, `[[`, "w_ls"),
      G = G_list, b = b_list,
      Theta = sweep(Theta_s, 2, M_diag, "/"),    # raw-unit library evaluation
      mt = mt, pt = pt, k_corner = k_corner, K = K,
      scale_x = scale_x, M_diag = M_diag,
      gamma = gamma, sigma_NR_hat = sigma_NR_hat,
      rescale = isTRUE(control$wsindy_rescale),
      condition_number = cond_G,
      rel_residual = rel_residual, r2 = r2,
      warnings = warnings,
      control = control
    ),
    class = "wsindy"
  )
  attr(res, "call") <- match.call()
  attr(res, "n_obs") <- mp1
  attr(res, "n_states") <- D
  attr(res, "n_terms") <- lib$J
  res
}

# S3 methods ===================================================================

#' Print method for wsindy objects
#'
#' @param x A wsindy object returned by [solveWSINDy()]
#' @param digits Number of significant digits for the displayed coefficients
#' @param ... Additional arguments (ignored)
#' @return Invisibly returns the input object
#' @method print wsindy
#' @export
print.wsindy <- function(x, digits = 4, ...) {
  D <- attr(x, "n_states")
  cat("WSINDy Model Discovery\n")
  cat("======================\n\n")
  cat("Observations:", attr(x, "n_obs"),
      "  States:", D,
      "  Library terms:", attr(x, "n_terms"), "\n")
  cat(sprintf("Test functions (per eq): mt = [%s], pt = [%s], K = [%s]\n\n",
              paste(x$mt, collapse = ", "),
              paste(x$pt, collapse = ", "),
              paste(x$K, collapse = ", ")))

  cat("Discovered system:\n")
  for (d in seq_len(D)) {
    nz <- which(x$W[, d] != 0)
    rhs <- if (length(nz) == 0L) {
      "0"
    } else {
      pieces <- vapply(nz, function(j) {
        lab <- rownames(x$W)[j]
        if (identical(lab, "1")) sprintf("%+.*g", digits, x$W[j, d])
        else sprintf("%+.*g*%s", digits, x$W[j, d], lab)
      }, character(1))
      sub("^\\+", "", paste(pieces, collapse = " "))
    }
    cat(sprintf("  du%d/dt = %s\n", d, rhs))
  }

  cat(sprintf("\nNonzero terms per equation: [%s]\n",
              paste(colSums(x$W != 0), collapse = ", ")))
  cat(sprintf("lambda*: [%s]\n",
              paste(sprintf("%.2e", x$lambda_star), collapse = ", ")))
  cat(sprintf("Relative weak residual: [%s]\n",
              paste(sprintf("%.3f", x$rel_residual), collapse = ", ")))
  if (length(x$warnings)) {
    cat("\nWarnings:\n")
    for (w in x$warnings) cat("  -", w, "\n")
  }
  invisible(x)
}

#' Summary method for wsindy objects
#'
#' @param object A wsindy object returned by [solveWSINDy()]
#' @param ... Additional arguments (ignored)
#' @return An object of class \code{summary.wsindy}
#' @method summary wsindy
#' @export
summary.wsindy <- function(object, ...) {
  D <- attr(object, "n_states")

  # Loss-landscape flatness: number of lambdas within 5% of the minimum loss.
  # A wide flat region indicates a robustly identified model (ODE paper).
  flatness <- vapply(object$losses, function(l) {
    if (all(is.na(l))) return(NA_integer_)
    sum(l <= min(l, na.rm = TRUE) * 1.05, na.rm = TRUE)
  }, integer(1))

  summ <- list(
    call = attr(object, "call"),
    n_obs = attr(object, "n_obs"),
    n_states = D,
    n_terms = attr(object, "n_terms"),
    ode_system = format_ode_system(object$f_sym, format = "ascii"),
    equations = data.frame(
      Equation     = paste0("du", seq_len(D), "/dt"),
      Nonzero      = colSums(object$W != 0),
      mt           = object$mt,
      pt           = object$pt,
      K            = object$K,
      Lambda       = object$lambda_star,
      RelResidual  = object$rel_residual,
      R2           = object$r2,
      CondG        = object$condition_number,
      LossFlatness = flatness
    ),
    param_map = object$param_map,
    scale_x = object$scale_x,
    rescale = object$rescale,
    warnings = object$warnings
  )
  class(summ) <- "summary.wsindy"
  summ
}

#' Print method for summary.wsindy objects
#'
#' @param x A summary.wsindy object
#' @param digits Number of significant digits to display
#' @param ... Additional arguments (ignored)
#' @return Invisibly returns the input object
#' @method print summary.wsindy
#' @export
print.summary.wsindy <- function(x, digits = 4, ...) {
  cat("WSINDy Discovery Summary\n")
  cat("========================\n\n")
  cat("Call:\n")
  print(x$call)
  cat("\n")
  cat("Observations:", x$n_obs, "  States:", x$n_states,
      "  Library terms:", x$n_terms, "\n")
  if (x$rescale) {
    cat(sprintf("State scales: [%s]\n",
                paste(sprintf("%.4g", x$scale_x), collapse = ", ")))
  } else {
    cat("Rescaling disabled.\n")
  }

  cat("\nDiscovered system:\n")
  cat(x$ode_system, "\n\n")

  cat("Per-equation diagnostics:\n")
  print(x$equations, digits = digits, row.names = FALSE)

  if (!is.null(x$param_map)) {
    cat("\nParameters:\n")
    print(x$param_map, digits = digits, row.names = FALSE)
  }
  if (length(x$warnings)) {
    cat("\nWarnings:\n")
    for (w in x$warnings) cat("  -", w, "\n")
  }
  invisible(x)
}

#' Extract the coefficient matrix from a wsindy object
#'
#' @param object A wsindy object
#' @param ... Additional arguments (ignored)
#' @return Term-by-state numeric coefficient matrix \code{W}
#' @method coef wsindy
#' @export
coef.wsindy <- function(object, ...) {
  object$W
}
