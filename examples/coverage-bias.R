library(wendy)
library(deSolve)
library(symengine)
library(trust)
library(tidyverse)
library(dplyr)
library(purrr)
library(tibble)
library(dplyr)
library(purrr)

# Method 0 is the coverage and bias, Method 1 is MLE Paper
addAdditiveGaussian <- function(noise_ratio, U, strategy) {
  state_sd <- apply(U, 2, sd)
  D <- ncol(U)
  mp1 <- nrow(U)
  switch(strategy,
         "0" = {
            x <- rnorm(mp1 * length(state_sd), mean = 0, sd = rep(state_sd, each = mp1))
            noise <- matrix(x, nrow = mp1, ncol = D)
            U + noise
          },{
            noise <- sqrt(noise_ratio * (norm(U, "F")^2 / mp1))
            U + noise
         }
  )
}

# Method 0 is the coverage and bias, Method 1 is MLE Paper
addMultiplicativeLogNormal <- function(noise_ratio, U, strategy) {
  state_sd <- apply(U, 2, sd)
  D <- ncol(U)
  mp1 <- nrow(U)
  switch(strategy,
         "0" = {
           col_range <- apply(U, 2, function(col) max(col) - min(col))
           sigma <- sqrt((log(col_range * noise_ratio) + 1) / 2)
           x <- rnorm(mp1 * D, mean = 0, sd = rep(state_sd, each = mp1))
           noise <- matrix(x, nrow = mp1, ncol = D)
           U * exp(noise)
         },{
           noise <- rnorm(mp1 * D, mean = 0, sd = sqrt(noise_ratio))
           U * exp(noise)
         }
  )
}

# Method 0 is the coverage and bias, Method 1 is MLE Paper
add_noise <- function(U, noise_ratio = 0.05, distribution = "AdditiveGaussian", strategy = 1){
  strategy <- as.character(strategy)
  switch(distribution, "AdditiveGaussian" = addAdditiveGaussian(noise_ratio, U, strategy),
                       "MultiplicativeLogNormal" = addMultiplicativeLogNormal(noise_ratio, U, strategy)
         )
}

f <- function(u, p, t) {
  c(p[1] * u[1] - p[2] * u[1]^2)
}

# Experiment params
distributions <- c("AdditiveGaussian", "MultiplicativeLogNormal")
noise_strategies <- c(0, 1)
methods <- c("IRLS", "MLE")
noise_ratios <- seq(0, 0.1, by = 0.05)
nreps <- 10

# ODE params
p_star <- c(1, 1)
u0 <- c(0.01)
p0_range <- list(c(0, 10), c(0, 10))
npoints <- 103
t_span <- c(0.001, 10)
t_eval <- seq(t_span[1], t_span[2], length.out = npoints)

modelODE <- function(tvec, state, parameters) {
  list(as.vector(f(state, parameters, tvec)))
}
sol <- deSolve::ode(y = u0, times = t_eval, func = modelODE, parms = p_star)
U <- sol[, 2, drop = F]
tt <- sol[, 1, drop = F]

params <- expand.grid(
  distribution = distributions,
  strategy = noise_strategies,
  method = methods,
  noise_ratio = noise_ratios,
  stringsAsFactors = FALSE
)

run_single <- function(dist, strat, meth, nr, p0_range) {
  dist <- "AdditiveGaussian"
  strat <- 0
  meth <- "IRLS"
  nr <- 0.05

  p0 <- mapply(\(r) runif(1, r[1], r[2]), p0_range)
  U_noisy <- add_noise(U, noise_ratio = nr, distribution = dist, strategy = strat)
  res <- solveWendy(f, p0, U_noisy, tt, method = meth)

  tibble(
    distribution = dist,
    strategy = strat,
    method = meth,
    noise_ratio = nr,
    converged = res$converged,
    iterations = res$iterations,
    final_sw_pvalue = if(!is.null(res$final_sw_pvalue)) res$final_sw_pvalue else NA,
    result = list(res)
  )
}

results <- params %>%
  crossing(rep = 1:nreps) %>%
  pmap_df(\(distribution, strategy, method, noise_ratio, rep) {
    run_single(distribution, strategy, method, noise_ratio, p0_range) %>%
      mutate(replicate = rep)
  })

head(results)

results_summary <- results %>%
  group_by(distribution, strategy, method, noise_ratio) %>%
  summarise(
    converged_rate = mean(converged),
    mean_iterations = mean(iterations),
    mean_sw_pvalue = mean(final_sw_pvalue, na.rm = TRUE),
    n_replicates = n(),
    .groups = "drop"
  )

print(results_summary)

saveRDS(results, "full_results.rds")
write_csv(results_summary, "results_summary.csv")
