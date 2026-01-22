# WENDy 

## System level dependencies

Symengine needs system level dependencies. From the the terminal run one of the following according to your OS.

```bash
zypper install cmake gmp-devel mpfr-devel mpc-devel    ## openSUSE
dnf    install cmake gmp-devel mpfr-devel libmpc-devel  ## Fedora
apt    install cmake libgmp-dev libmpfr-dev libmpc-dev  ## Debian
brew   install cmake gmp mpfr libmpc                    ## Mac OS
```

## R dependencies 

```r
install.packages(c("symengine", "trust", "stats", "numbers", "torch", "minpack.lm", "deSolve"))
```

If you are having trouble installing the library *torch* you can follow [this guide.](https://cran.r-project.org/web/packages/torch/vignettes/installation.html)


To install WENDy run the following. Note that if the github repository is private you will have to set your github credentials.

```r
library(remotes)
remotes::install_github("jackkrebsbach/rwendy")
```