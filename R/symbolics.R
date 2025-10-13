library(symengine)

compute_symbolic_jacobian <- function(f_expr, vars) {
  dims <- dim(f_expr)

  if (is.null(dims)) dims <- length(f_expr)

  n_vars <- length(vars)

  f_flat <- as.vector(f_expr)

  deriv_list <- lapply(f_flat, function(f_i) {
    lapply(vars, function(v) D(f_i, v))
  })

  deriv_flat <- unlist(deriv_list, recursive = FALSE)

  out_dims <- c(dims, n_vars)

  J <- array(deriv_flat, dim = out_dims)

  return(J)
}

build_fn <- function(expr_array, vars) {
  dims <- dim(expr_array)
  if (is.null(dims)) {
    expr_vec <- symengine::Vector(expr_array)
  } else {
    expr_flat <- array(expr_array)
    expr_vec <- symengine::Vector(expr_flat)
  }
  visitor <- DoubleVisitor(expr_vec,
                           args = vars,
                           perform_cse = TRUE,
                           llvm_opt_level = if (symengine_have_component("llvm")) 3L else -1L)
  function(input) {
    vals <- visitor_call(visitor, input, do_transpose = FALSE)
    if (is.null(dims)) return(vals)
    fn <- array(vals, dim = rev(dims)) # To match column wise vectorization
    fn <- aperm(fn, length(dims):1)
  }
}

lognormal_transform <- function(f_sym){
  # Dimension of system
  D <- length(f_sym)

  # Get symbolic variables
  u <- do.call(c, lapply(1:D, \(i) S(paste0("u", i))))

  # Build substitution arguments: u1, exp(u1), u2, exp(u2), ...
  sub_args <- unlist(Map(list, u, lapply(u, exp)))

  # Perform: f_sym[i] / uᵢ, substitute uᵢ ↦ exp(uᵢ)
  logu <- do.call(c, lapply(1:D, \(i) {
    fud <- f_sym[i] / u[i]
    do.call(symengine::subs, c(list(fud), sub_args)) # subs is a symengine command
  }))
}






