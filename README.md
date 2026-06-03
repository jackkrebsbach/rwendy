# WENDy

[![Binder](https://mybinder.org/badge_logo.svg)](
https://mybinder.org/v2/gh/jackkrebsbach/binder/HEAD?urlpath=notebooks/wendy_demo.ipynb
)

**Weak-form Estimation of Nonlinear Dynamics (WENDy)** is an algorithm for
estimating the parameters of a system of ordinary differential equations (ODEs)
from noisy time-series data.

Instead of numerically differentiating noisy observations — or repeatedly
forward-solving the ODE inside an optimizer — WENDy converts the ODE to its
**weak form** by integrating against compactly supported test functions. This
sidesteps numerical differentiation entirely and turns parameter estimation
into a regression problem. The result is an estimator that is robust to large
measurement noise and, for stiff or high-dimensional systems, often *faster*
and *more accurate* than forward-solver-based nonlinear least squares.

## Installation

```r
# install.packages("remotes")
remotes::install_github("jackkrebsbach/rwendy")
```

## Quick start

```r
library(wendy)
library(deSolve)

# Logistic growth
f <- function(u, p, t) {
  c(p[1] * u[1] - p[2] * u[1]^2)
}
p_star <- c(1, 1/10)
tt     <- seq(0, 10, length.out = 60)

modelODE <- function(tvec, state, parameters) {
  list(as.vector(f(state, parameters, tvec)))
}
sol <- deSolve::ode(y = 0.5, times = tt, func = modelODE, parms = p_star)

set.seed(1)
# Additive Gaussian noise
U <- sol[, -1, drop = FALSE] + matrix(rnorm(length(tt), sd = 0.1), length(tt))

res <- solveWendy(f, U, matrix(tt, ncol = 1), method = "IRLS")
res$phat            # estimated (r, K)
summary(res)        # estimates, standard errors, covariance, diagnostics
```

## Defining the ODE

The right-hand side `f(u, p, t)` takes the state vector `u`, the parameter
vector `p`, and time `t`, and **returns a vector** with one entry per state
variable. For example, the Lorenz system:

```r
f <- function(u, p, t) {
  du1 <- p[1] * (u[2] - u[1])
  du2 <- u[1] * (p[2] - u[3]) - u[2]
  du3 <- u[1] * u[2] - p[3] * u[3]
  c(du1, du2, du3)
}
```

`U` is then a matrix whose columns are the state variables and whose rows
correspond to the time points in `tt`.

## Estimation methods

Choose with the `method` argument of `solveWendy()`:

| Method     | Description                                              | Notes |
|------------|----------------------------------------------------------|-------|
| `"IRLS"`   | Iteratively reweighted weak-form least squares (default) | General-purpose; robust to noise. |
| `"MLE"`    | Weak-form maximum likelihood via trust-region optimization | Best statistical efficiency; handles nonlinear-in-parameters ODEs. |
| `"OLS"`    | Unweighted weak-form least squares                       | Fast; best on well-resolved / low-noise data. |
| `"OE"`     | Output error (forward-solve nonlinear least squares)     | Requires an initial guess `p0`. |
| `"HYBRID"` | WENDy estimate used to initialize an output-error refine | Useful for stiff / hard-to-fit systems. |

WENDy also supports two noise models via `noise_dist`: `"addgaussian"`
(additive Gaussian, the default) and `"lognormal"` (multiplicative log-normal).

```r
res <- solveWendy(f, U, tt, method = "MLE", noise_dist = "addgaussian")
```

## Working with the result

`solveWendy()` returns an object of class `wendy`. The usual extractor methods
are available:

```r
res$phat              # point estimate (also: coef(res))
coef(res)             # estimated parameter vector
summary(res)          # full summary 
residuals(res)        # weak-form residuals
```

`summary(res)` reports the fitted ODE system, parameter estimates with
**standard errors** and the **parameter covariance** matrix (inverse Fisher
information), a weak negative log-likelihood value, weak-residual quantiles,
and convergence diagnostics.

## Beyond point estimates

- **Initial-condition estimation** — set `estimate_IC = TRUE` to recover the
  initial state; the estimate is returned in `res$u0hat` (with covariance).
- **State trajectory smoothing** — set `estimate_trajectory = TRUE` to obtain a
  denoised trajectory in `res$state` (`res$state$U_star`), using either an
  extended RTS smoother (`smoother = "erts"`, default) or a Matérn-5/2 Gaussian
  process (`smoother = "gp"`).

```r
res <- solveWendy(
  f, U, tt,
  method  = "IRLS",
  control = list(estimate_IC = TRUE, estimate_trajectory = TRUE)
)
res$u0hat              # estimated initial condition
res$state$U_star       # smoothed state trajectory
```

## Control options

Fine-grained behavior is configured through the `control` argument, which
accepts a plain named list or a call to `wendy_control()`; entries you supply
override the defaults. See `?wendy_control` for the full set of options.
Some commonly used ones:

| Option                  | Default | Description |
|-------------------------|---------|-------------|
| `optimize`              | `TRUE`  | Run the IRLS/MLE optimization (vs. setup only). |
| `estimate_IC`           | `FALSE` | Estimate the initial condition. |
| `estimate_trajectory`   | `FALSE` | Estimate the smoothed state trajectory. |
| `noise_sd`              | `NA`    | Known noise standard deviation, if available. |
| `test_fun_type`         | `"MSG"` | Test functions: multi-scale global (`"MSG"`) or single-scale local (`"SSL"`). |
| `test_fun`              | `"phi"` | C-infinity bump (`"phi"`) or piecewise polynomial (`"psi"`) test function. |
| `k_max`                 | `200`   | Maximum number of test functions. |
| `radius_params`         | `2^(0:3)` | Radius scale factors for MSG test functions. |
| `interpolation_method`  | `NULL`  | `"gp"`, `"spline"`, `"linear"`, `"cubic"`, `"loess"`, `"kernel"`, or `"poly_ls_N"`. |
| `smoother`              | `"erts"`| State smoother for `estimate_trajectory`: `"erts"` or `"gp"`. |
| `apply_fn`              | `NULL`  | Custom apply for multi-start (e.g. `parallel::mclapply`); `NULL` uses `lapply`. |


## Examples

The [`examples/`](examples/) directory contains worked end-to-end scripts for a
range of benchmark systems:

- Logistic growth, Lotka–Volterra, Duffing oscillator
- FitzHugh–Nagumo, Hindmarsh–Rose (neuroscience)
- Goodwin (2D/3D), protein transduction (biochemistry)
- SIR with time-dependent infectivity, HIV dynamics
- Lorenz and Lorenz-96 (chaotic systems)

## References

- Bortz, D. M., Messenger, D. A., & Dukic, V. (2023). *Direct Estimation of
  Parameters in ODE Models Using WENDy: Weak-form Estimation of Nonlinear
  Dynamics.* **Bulletin of Mathematical Biology**, 85(110).
  doi:[10.1007/s11538-023-01208-6](https://doi.org/10.1007/S11538-023-01208-6)
- Rummel, N., Messenger, D. A., Becker, S., Dukic, V., & Bortz, D. M. (2025).
  *WENDy for Nonlinear-in-Parameters ODEs.* arXiv:2502.08881.
  doi:[10.48550/arXiv.2502.08881](https://doi.org/10.48550/arXiv.2502.08881)
- Tran, A., & Bortz, D. M. (2026). *Weak Form Scientific Machine Learning:
  Test Function Construction for System Identification.* **SIAM Journal on
  Scientific Computing** (accepted).
