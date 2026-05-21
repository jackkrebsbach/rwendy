# WENDy

## System dependencies

`symengine` requires native C/C++ libraries. Install them once via your system
package manager:

```bash
zypper install cmake gmp-devel mpfr-devel mpc-devel    ## openSUSE
dnf    install cmake gmp-devel mpfr-devel libmpc-devel ## Fedora
apt    install cmake libgmp-dev libmpfr-dev libmpc-dev ## Debian / Ubuntu
brew   install cmake gmp mpfr libmpc                   ## macOS
```

Windows users installing from source need [Rtools](https://cran.r-project.org/bin/windows/Rtools/).

## Install WENDy

Once CRAN-released:

```r
install.packages("wendy")
```

Development version from GitHub:

```r
# install.packages("remotes")
remotes::install_github("jackkrebsbach/rwendy")
```

All R-level dependencies (`symengine`, `trust`, `minpack.lm`, `numDeriv`,
`numbers`, `MASS`, `deSolve`, `FME`) are declared in `DESCRIPTION` and pulled in
automatically by `install.packages` / `remotes::install_github`.

## Quick start

```r
library(wendy)
library(deSolve)

# Logistic growth: u' = r * u * (1 - u / K), p = (r, K)
f <- function(u, p, t) p[1] * u[1] * (1 - u[1] / p[2])
p_star <- c(0.5, 10)
tt <- seq(0, 10, length.out = 60)

modelODE <- function(tvec, state, parameters) {
  list(as.vector(f(state, parameters, tvec)))
}
sol <- deSolve::ode(y = 1, times = tt, func = modelODE, parms = p_star)
set.seed(1)
U <- sol[, -1, drop = FALSE] + matrix(rnorm(length(tt), sd = 0.1), length(tt))

res <- solveWendy(f, U, matrix(tt, ncol = 1), method = "IRLS")
res$phat
```
