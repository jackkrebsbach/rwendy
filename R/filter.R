# 𝐂^{∞} Bump test function
# 𝜱(t; r, n) = exp(η/ [(1 - (t/r)²)]₊)
phi <- function(t, r, eta = 9) {
  exp(-eta / (1 - (t / r)^2))
}

# Piecewise polynomial test function -
# Advantageous because we have access to analytic form of Fourier coefficients
# 𝛹(t; r, p) = [(1- (t/r)²)]ᵖ₊
psi <- function(t, r, p = 16){
  (1 - ( t / r)^2)^p
} 

wendy_filter <- function(U, tt, p, f, radius, test_function = phi){
    dt <- mean(diff(tt))
    tf_prime <- test_function_derivative(test_function, radius, dt, order = 1)
    V <- build_full_test_function_matrix(test_function, tt, radius, order = 0)
    Vp <- build_full_test_function_matrix(test_function, tt, radius, order = 1)
}

