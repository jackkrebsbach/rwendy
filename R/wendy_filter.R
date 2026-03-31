
# For the Euler-Maclaurin correction we need:
# ∇ₜf(p,u,t) and ∇ᵤf(p,u,t)

# Retuns the correction for one integral rₖ = ∫₀ᵀ g(t) dt
# g(t) = φ(t)f(p,u,t) + φ′(t)u(t)
# g′(t) = φ′(t)f(p,u,t) + φ(t)fₜ(p,u,t) + φ′′(t)u(t) + φ′(t)u′(t)
#       = 2φ′(t)f(p,u,t) + φ(t)fₜ(p,u,t) + φ′′(t)u(t) 
# Where fₜ(p,u,t) = ∂f/∂u f(p,u,t) +  ∂f/∂t
em_correction <- function(
  test_fun_0, test_fun_der_0, J_u, J_t, f_, u_0
){
  
  # ∂f/∂t 

}