# %%
# Three-state SIR with known population and only infected prevalence measured.
# S(0) and R(0) are NOT supplied to the fit. The three-state truth below is used
# only to simulate infected measurements and to score/plot the reconstruction.
# Closed population, constant positive rates, absolute prevalence (not incidence
# or an unknown reporting fraction), and a window covering growth and decline.
# See examples/example_latent_sir.md for the identifiability derivation.
library(deSolve)
library(plotly)
invisible(devtools::load_all())

population <- 10000
p_star <- c(beta = 0.4, gamma = 0.12)
u0 <- c(S = 8500, I = 500, R = 1000)
npoints <- 20L
t_span <- c(0, 40)
tt <- seq(t_span[1], t_span[2], length.out = npoints)

# Begin with all three physical compartments, in people.
sir_full <- function(t, state, p) {
  infection <- p[1] * state[1] * state[2] / population
  recovery <- p[2] * state[2]
  list(c(-infection, infection - recovery, recovery))
}
sol <- deSolve::ode(u0, tt, sir_full, p_star, rtol = 1e-11, atol = 1e-11)
truth <- sol[, c("S", "I", "R")]
stopifnot(abs(sum(u0) - population) < 1e-8)

# set.seed(8675309)
nr <- 0.1
noise_sd <- nr * sqrt(mean(truth[, "I"]^2))
observations <- matrix(NA_real_, npoints, 3L,
                       dimnames = list(NULL, c("S", "I", "R")))
observations[, "I"] <- truth[, "I"] + rnorm(npoints, sd = noise_sd)

# Use i=I/N and v=S/(S+R), the susceptible share of noninfected people.
# Then S/N=(1-i)*v and R/N=(1-i)*(1-v). Box bounds 0<i,v<1
f <- function(u, p, t) {
  i <- u[1]; v <- u[2]
  c(p[1]*i*(1-i)*v-p[2]*i, -p[1]*i*v*(1-v)-p[2]*i*v/(1-i))
}

U <- cbind(I = observations[, "I"] / population, susceptible_share = NA_real_)

elapsed <- system.time({
  fit <- solveWendyGP(f, U, tt,
    parameter_lower = c(0, 0),
    parameter_upper = c(Inf,Inf),
    state_lower = c(0, 0),
    state_upper = c(1, 1))
})

i_hat <- fit$U_hat[, 1]
v_hat <- fit$U_hat[, 2]

state_hat <- population * cbind(S = (1-i_hat)*v_hat, I = i_hat,
                                R = (1-i_hat)*(1-v_hat))
initial_hat <- state_hat[1, ]
conservation_error <- max(abs(rowSums(state_hat) - population))

infection_information <- function(times, p, s0, i0, noise_fraction) {
  beta <- p[1]
  gamma <- p[2]

  rhs <- function(t, state, unused) {
    i <- state[1]
    s <- state[2]
    Z <- matrix(state[-(1:2)], 2L, 4L)
    A <- matrix(c(beta*s-gamma, -beta*s, beta*i, -beta*i), 2L, 2L)
    B <- cbind(c(s*i, -s*i), c(-i, 0), c(0, 0), c(0, 0))
    list(c(beta*s*i-gamma*i, -beta*s*i, as.vector(A %*% Z + B)))
  }

  Z0 <- cbind(c(0, 0), c(0, 0), c(0, 1), c(1, 0))

  trajectory <- deSolve::ode(c(i0, s0, as.vector(Z0)), times, rhs, NULL, rtol = 1e-10, atol = 1e-12)

  J <- trajectory[, c(4, 6, 8, 10), drop = FALSE] / noise_fraction
  colnames(J) <- c("beta", "gamma", "s0", "i0")

  scales <- sqrt(colSums(J^2))
  empty <- list(rank=0L,condition=Inf,local_noise_sd=setNames(rep(NA_real_,4L),colnames(J)),
    r0_noise_sd=NA_real_,jacobian=J)
  if (any(!is.finite(scales) | scales <= 0)) return(empty)
  spectrum <- svd(sweep(J, 2, scales, "/"))
  rank <- sum(spectrum$d > max(spectrum$d) * 1e-8)
  if (rank < 4L) {
    empty$rank <- rank; empty$singular_values <- spectrum$d
    return(empty)
  }
  covariance <- sweep(sweep(spectrum$v %*% diag(1 / spectrum$d^2) %*% t(spectrum$v), 1, scales, "/"), 2, scales, "/")
  list(rank = rank, singular_values = spectrum$d,
       condition = max(spectrum$d) / min(spectrum$d),
       local_noise_sd = setNames(sqrt(diag(covariance)), colnames(J)),
       r0_noise_sd = sqrt(sum(covariance[3:4, 3:4])), jacobian = J)
}
info <- infection_information(tt, fit$phat, initial_hat["S"] / population, initial_hat["I"] / population, noise_sd / population)

print(fit)
cat(sprintf("\nElapsed: %.3f s", elapsed[["elapsed"]] ))
print(data.frame(parameter = names(p_star), truth = unname(p_star), estimate = unname(fit$phat) ))

ref <- deSolve::ode(u0, fit$tt, sir_full, p_star, rtol = 1e-11, atol = 1e-11)[, c("S", "I", "R")]
component_colors <- c(S = "#0072B2", I = "#D55E00", R = "#009E73")
fig <- plot_ly()
for (component in colnames(state_hat)) {
  color <- component_colors[[component]]
  fig <- fig |>
    add_trace(x = fit$tt, y = ref[, component] / population, type = "scatter", mode = "lines",
              name = paste(component, "truth"), legendgroup = component,
              line = list(color = color, width = 1.5), opacity = 0.65) |>
    add_trace(x = fit$tt, y = state_hat[, component] / population, type = "scatter", mode = "lines",
              name = paste(component, if (component == "I") "fitted" else "inferred"),
              legendgroup = component, line = list(color = color, width = 2.5, dash = "dash"))
}
fig <- fig |>
  add_trace(x = tt, y = observations[, "I"] / population, type = "scatter", mode = "markers",
            name = "I observations", legendgroup = "I",
            marker = list(color = component_colors[["I"]], size = 8), opacity = 0.6) |>
  layout(title = list(text = paste0("SIR: infected observed",
    "<br><sup>Estimated beta = ", sprintf("%.3f", fit$phat[1]),
    ", gamma = ", sprintf("%.3f", fit$phat[2]),
    " · True: ", paste(sprintf("%.3f", p_star), collapse = ", "), "</sup>")),
    xaxis = list(title = "Time"), yaxis = list(title = "Proportion of population"),
    hovermode = "x unified", legend = list(groupclick = "toggleitem"))
print(fig)