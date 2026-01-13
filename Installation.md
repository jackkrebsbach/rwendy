# WENDy 

WENDy depends on other packages. To install them run 

```r
install.packages(c("symengine", "trust", "stats", "numbers", "torch", "minpack.lm"))
```
Symengine needs system level dependencies. From the command line run one of the following:

```bash
zypper install cmake gmp-devel mpfr-devel mpc-devel    ## openSUSE
dnf    install cmake gmp-devel mpfr-devel libmpc-devel  ## Fedora
apt    install cmake libgmp-dev libmpfr-dev libmpc-dev  ## Debian
brew   install cmake gmp mpfr libmpc                    ## Mac OS
```

To simulate data using the examples install **deSolve**

```r
install.packages("deSolve")
```

To install WENDy run. Note that if the github repository is private you will have to set your github credentials.

```r
library(devtools)
devtools::install_github("jackkrebsbach/rwendy")
```