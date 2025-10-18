# WENDy 

WENDy depends on other packages. To install them run 

```r
install.packages(c("symengine", "trust", "stats"))
```
Symengine needs system level dependencies. From the command line run one of the following:

```bash
## openSUSE
zypper install cmake gmp-devel mpfr-devel mpc-devel    
## Fedora
dnf    install cmake gmp-devel mpfr-devel libmpc-devel 
## Debian
apt    install cmake libgmp-dev libmpfr-dev libmpc-dev 
## Mac OS
brew   install cmake gmp mpfr libmpc                   

```

To simulate data using the examples install **deSolve**

```r
install.packages("deSolve")
```

To install WENDy run

```r
devtools::install_github("jackkrebsbach/rwendy")
```
