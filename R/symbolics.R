# Detect the number of parameters by scanning f's body for param[N] references,
# where param is the name of the second argument of f (the parameter vector).
# Returns the maximum index found, i.e. the number of parameters.
# @keywords internal
detect_n_params <- function(f) {
  param_name <- names(formals(f))[[2]]
  body_str   <- paste(deparse(body(f)), collapse = " ")
  pattern    <- paste0(param_name, "\\[\\d+\\]")
  hits       <- regmatches(body_str, gregexpr(pattern, body_str))[[1]]
  if (length(hits) == 0L) return(0L)
  max(as.integer(regmatches(hits, gregexpr("\\d+", hits))))
}

# Detect the maximum polynomial degree of state variables in the symbolic RHS.
# Parses the string representation of f_expr and looks for u_i^N patterns.
# Returns an integer >= 1.
# @keywords internal
detect_max_state_order <- function(f_expr, u_expr) {
  expr_str  <- paste(as.character(f_expr), collapse = " ")
  u_names   <- as.character(u_expr)
  max_order <- 1L
  for (u_name in u_names) {
    escaped <- gsub("([.|()\\^{}+?$*])", "\\\\\\1", u_name)
    pattern <- paste0(escaped, "\\^(\\d+)")
    hits    <- regmatches(expr_str, gregexpr(pattern, expr_str))[[1]]
    for (h in hits) {
      pwr       <- as.integer(sub(paste0(escaped, "\\^"), "", h))
      max_order <- max(max_order, pwr)
    }
  }
  max_order
}

compute_symbolic_jacobian <- function(f_expr, vars) {
  dims <- dim(f_expr)

  if (is.null(dims)) dims <- length(f_expr)

  n_vars <- length(vars)
  f_flat <- as.vector(f_expr)

  deriv_list <- lapply(f_flat, function(f_i) {
    lapply(vars, function(v) symengine::D(f_i, v))
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
  visitor <- symengine::DoubleVisitor(expr_vec,
                           args = vars,
                           perform_cse = TRUE,
                           llvm_opt_level = -1L)
  function(input) {
    symengine::visitor_call(visitor, input, do_transpose = TRUE)
  }
}

# Total time derivative of a symbolic vector g(u(t), t) along the trajectory u'=f.
# Returns a vector of symbolic expressions of the same length as g_sym.
compute_symbolic_total_time_deriv <- function(g_sym, u_expr, f_expr, t_expr) {
  g_flat <- as.vector(g_sym)
  f_flat <- as.vector(f_expr)
  u_flat <- as.vector(u_expr)
  do.call(c, lapply(g_flat, function(g_i) {
    total <- symengine::D(g_i, t_expr)
    for (k in seq_along(f_flat))
      total <- total + symengine::D(g_i, u_flat[[k]]) * f_flat[[k]]
    total
  }))
}

lognormal_transform <- function(f_sym){
  D <- length(f_sym)
  u <- do.call(c, lapply(1:D, function(i) symengine::S(paste0("u", i))))
  sub_args <- unlist(Map(list, u, lapply(u, exp)))

  logu <- do.call(c, lapply(1:D, function(i) {
    fud <- f_sym[i] / u[i]
    do.call(symengine::subs, c(list(fud), sub_args))
  }))
}