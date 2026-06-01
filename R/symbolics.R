# Symbolic engine ==============================================================
#
# WENDy needs a small set of symbolic operations on the ODE right-hand side
# f(u, p, t): create symbols, differentiate, substitute, stringify, render LaTeX,
# and compile to a fast numeric evaluator. These are provided by one of two
# interchangeable backends:
#
#   * "symengine" - fast, LLVM-compiled, but a heavy external dependency.
#   * "native"    - pure base-R fallback (the `rsym` class below), slower in
#                    principle but dependency-free.
#
# The file is organised in four sections:
#   1. Backend selection            - sym_backend()
#   2. Primitive dispatch           - sym_*() generics routing to a backend
#   3. Native backend               - the `rsym` S3 class + native_*() impls
#   4. Backend-agnostic algorithms  - Jacobians, total time derivatives, etc.
#
# Sections 1-3 keep all representation details behind the sym_*() primitives, so
# the algorithms in section 4 (and the callers in wendy.R) are written once and
# work unchanged on either backend.
#
# Representation conventions:
#   native    : scalars and vectors/arrays are `rsym` objects (shape via dim_sym).
#   symengine : scalars are Basic; 1-D vectors are symengine Vectors; arrays are
#               base-R list-arrays of Basic scalars (column-major).


# 1. Backend selection ---------------------------------------------------------

.wendy_sym <- new.env(parent = emptyenv())

# Resolve the active backend ("symengine" or "native").
# Override with options(wendy.symbolic_backend = "native" | "symengine");
# otherwise symengine is used if installed, else the native fallback.
sym_backend <- function() {
  ov <- getOption("wendy.symbolic_backend", NULL)
  if (!is.null(ov)) {
    ov <- match.arg(ov, c("symengine", "native"))
    if (ov == "symengine" && !requireNamespace("symengine", quietly = TRUE))
      stop("options(wendy.symbolic_backend='symengine') set but symengine is not installed.")
    return(ov)
  }
  if (is.null(.wendy_sym$backend))
    .wendy_sym$backend <- if (requireNamespace("symengine", quietly = TRUE))
      "symengine" else "native"
  .wendy_sym$backend
}


# 2. Primitive dispatch --------------------------------------------------------

# Create a scalar symbol from a name.
sym_symbol <- function(name) {
  if (sym_backend() == "native") native_symbol(name) else symengine::S(name)
}

# Number of elements in a symbolic vector/array.
sym_length <- function(x) length(x)

# Shape of a symbolic vector/array (NULL for a flat vector).
sym_dim <- function(x) {
  if (sym_backend() == "native") native_dim(x) else dim(x)
}

# i-th element (column-major) as a scalar.
sym_elt <- function(x, i) {
  if (sym_backend() == "native") native_build(list(x[[i]])) else as.vector(x)[[i]]
}

# Assemble a symbolic vector/array from a list of scalars (column-major).
sym_build <- function(scalar_list, dims = NULL) {
  if (sym_backend() == "native") return(native_build(scalar_list, dims = dims))
  if (is.null(dims) || length(dims) <= 1L) do.call(c, scalar_list)
  else array(scalar_list, dim = dims)
}

# Partial derivative of a scalar expression w.r.t. a scalar variable.
sym_diff <- function(scalar, var) {
  if (sym_backend() == "native") native_diff(scalar, var) else symengine::D(scalar, var)
}

# Substitute scalars `from` with scalars `to` in a scalar expression.
sym_subs <- function(scalar, from, to) {
  if (sym_backend() == "native") return(native_subs(scalar, from, to))
  do.call(symengine::subs, c(list(scalar), unlist(Map(list, from, to))))
}

# Character (ascii) representation, one string per element.
sym_strings <- function(x) {
  if (sym_backend() == "native") return(native_strings(x))
  as.character(symengine::Vector(if (is.null(dim(x))) as.vector(x) else array(x)))
}

# LaTeX string for a scalar expression.
sym_latex <- function(scalar) {
  if (sym_backend() == "native") native_latex(scalar) else symengine::codegen(scalar, "latex")
}

# DoubleVisitor opt level: a non-negative value JIT-compiles via LLVM (much
# faster per eval), but symengine silently falls back to a tree-walking visitor
# for negative values OR when it was built without LLVM. So use the JIT only
# when LLVM is actually available, else -1. The (cheap) capability check is
# cached; override the JIT level with options(wendy.symengine_llvm_opt_level=).
sym_llvm_opt_level <- function() {
  if (is.null(.wendy_sym$has_llvm))
    .wendy_sym$has_llvm <- tryCatch(
      isTRUE(symengine::symengine_have_component("llvm")),
      error = function(e) FALSE)
  if (.wendy_sym$has_llvm)
    as.integer(getOption("wendy.symengine_llvm_opt_level", 3L)) else -1L
}

# Compile a symbolic vector/array into a vectorised numeric evaluator.
# Contract: input (n_vars x n_pts) -> output (n_pts x n_out), n_out = #elements.
sym_compile <- function(x, vars) {
  if (sym_backend() == "native") return(native_compile(x, vars))
  dims <- dim(x)
  expr_vec <- if (is.null(dims)) symengine::Vector(x) else symengine::Vector(array(x))
  visitor <- symengine::DoubleVisitor(expr_vec, args = vars,
                                      perform_cse = TRUE,
                                      llvm_opt_level = sym_llvm_opt_level())
  function(input) symengine::visitor_call(visitor, input, do_transpose = TRUE)
}

# Return an R function of one scalar variable.
sym_lambdify <- function(scalar, var) {
  if (sym_backend() == "native") native_lambdify(scalar, var) else symengine::lambdify(scalar, var)
}

# 3. Native backend ------------------------------------------------------------
#
# The `rsym` S3 class is a list of base-R language objects (calls / names /
# numeric constants) carrying an optional `dim_sym` attribute. A scalar is an
# `rsym` of length 1; a vector/array stores its elements in column-major order
# with `dim_sym` recording the shape.
#
# Operator overloading (Ops/Math group generics) lets a user's f(u, p, t) be
# evaluated symbolically with the exact source code it uses numerically.
# Differentiation uses stats::D(), substitution uses substitute(), and
# compilation builds a vectorised closure via eval(). Supported functions are
# those in stats::D()'s derivative table (+ - * / ^, exp, log, sqrt, sin, cos,
# tan, ...), which covers the ODE models WENDy targets.

# Constructor: wrap a language object (or list of them) as an rsym.
new_rsym <- function(exprs, dims = NULL) {
  if (!is.list(exprs)) exprs <- list(exprs)
  structure(exprs, dim_sym = dims, class = "rsym")
}

# Strip class/dim and return the bare list of language objects.
rsym_exprs <- function(x) {
  attr(x, "dim_sym") <- NULL
  class(x) <- NULL
  x
}

# Coerce an operand (rsym or numeric) to a plain list of elements.
.rsym_operand <- function(x) if (inherits(x, "rsym")) rsym_exprs(x) else as.list(x)

#' @export
Ops.rsym <- function(e1, e2) {
  op <- .Generic
  if (missing(e2)) {                       # unary + / -
    a <- .rsym_operand(e1)
    return(new_rsym(lapply(a, function(ai) call(op, ai))))
  }
  a <- .rsym_operand(e1)
  b <- .rsym_operand(e2)
  n <- max(length(a), length(b))
  res <- vector("list", n)
  for (k in seq_len(n)) {
    ak <- a[[(k - 1L) %% length(a) + 1L]]
    bk <- b[[(k - 1L) %% length(b) + 1L]]
    res[[k]] <- call(op, ak, bk)
  }
  new_rsym(res)
}

#' @export
Math.rsym <- function(x, ...) {
  op <- .Generic
  new_rsym(lapply(rsym_exprs(x), function(xi) call(op, xi)))
}

#' @export
`[.rsym` <- function(x, i) new_rsym(rsym_exprs(x)[i])

#' @export
`[[.rsym` <- function(x, i) new_rsym(rsym_exprs(x)[[i]])

#' @export
c.rsym <- function(...) {
  parts <- lapply(list(...), function(z)
    if (inherits(z, "rsym")) rsym_exprs(z) else as.list(z))
  new_rsym(do.call(c, parts))
}

#' @export
length.rsym <- function(x) length(rsym_exprs(x))

#' @export
dim.rsym <- function(x) attr(x, "dim_sym")

#' @export
as.character.rsym <- function(x, ...)
  vapply(rsym_exprs(x), function(e) paste(deparse(e), collapse = ""), character(1))

#' @export
print.rsym <- function(x, ...) {
  cat("<rsym>", if (!is.null(dim(x))) paste0("[", paste(dim(x), collapse = "x"), "]"), "\n")
  cat(as.character(x), sep = "\n")
  invisible(x)
}

native_symbol <- function(name) new_rsym(as.name(name))

native_dim <- function(x) attr(x, "dim_sym")

native_build <- function(scalar_list, dims = NULL)
  new_rsym(lapply(scalar_list, function(s) rsym_exprs(s)[[1]]), dims = dims)

native_diff <- function(scalar, var) {
  e <- rsym_exprs(scalar)[[1]]
  if (is.numeric(e)) return(new_rsym(0))
  vn <- as.character(rsym_exprs(var)[[1]])
  new_rsym(stats::D(e, vn))
}

native_subs <- function(scalar, from, to) {
  env <- list()
  for (i in seq_along(from))
    env[[as.character(rsym_exprs(from[[i]])[[1]])]] <- rsym_exprs(to[[i]])[[1]]
  new_rsym(do.call("substitute", list(rsym_exprs(scalar)[[1]], env)))
}

native_strings <- function(x) as.character(x)

# Best-effort LaTeX: base R has no symbolic LaTeX printer, so drop the explicit
# multiplication stars from the deparse. Used only for human-readable summaries.
native_latex <- function(scalar)
  gsub("\\*", " ", paste(deparse(rsym_exprs(scalar)[[1]]), collapse = ""))

# Build a vectorised numeric evaluator. Contract matches symengine's
# visitor_call(do_transpose = TRUE): input is (n_vars x n_pts) [a bare vector is
# treated as a single point]; output is (n_pts x n_out), n_out = #elements.
# Vectorisation is over points: each variable is bound to a whole row vector, so
# one eval() per output expression drives R's C-level arithmetic over all points.
native_compile <- function(x, vars) {
  var_names <- vapply(rsym_exprs(vars), as.character, character(1))
  exprs <- rsym_exprs(x)
  n_out <- length(exprs)
  function(input) {
    if (is.null(dim(input))) input <- matrix(input, ncol = 1L)
    n_pts <- ncol(input)
    env <- new.env(parent = baseenv())
    for (i in seq_along(var_names)) assign(var_names[i], input[i, ], envir = env)
    out <- matrix(0, n_pts, n_out)
    for (j in seq_len(n_out))
      out[, j] <- rep(eval(exprs[[j]], env), length.out = n_pts)
    out
  }
}

# Return an R function of a single (scalar symbol) variable, vectorised.
native_lambdify <- function(scalar, var) {
  e <- rsym_exprs(scalar)[[1]]
  vn <- as.character(rsym_exprs(var)[[1]])
  function(val) {
    env <- new.env(parent = baseenv())
    assign(vn, val, envir = env)
    eval(e, env)
  }
}


# 4. Backend-agnostic algorithms -----------------------------------------------

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
  expr_str  <- paste(sym_strings(f_expr), collapse = " ")
  u_names   <- sym_strings(u_expr)
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
  # Build derivatives in column-major flat order: state index fastest, var index slowest.
  # The resulting array has dims c(dims, n_vars) so that arr[i, ..., v] = ∂f_i / ∂vars[v]
  # under R's native column-major reshape.
  dims <- sym_dim(f_expr)
  if (is.null(dims)) dims <- sym_length(f_expr)
  n_vars <- sym_length(vars)
  n_f    <- sym_length(f_expr)

  deriv_list <- vector("list", n_vars * n_f)
  idx <- 1L
  for (v in seq_len(n_vars)) {
    var_v <- sym_elt(vars, v)
    for (i in seq_len(n_f)) {
      deriv_list[[idx]] <- sym_diff(sym_elt(f_expr, i), var_v)
      idx <- idx + 1L
    }
  }
  sym_build(deriv_list, dims = c(dims, n_vars))
}

build_fn <- function(expr_array, vars) {
  sym_compile(expr_array, vars)
}

# Total time derivative of a symbolic vector g(u(t), t) along the trajectory u'=f.
# Returns a vector of symbolic expressions of the same length as g_sym.
compute_symbolic_total_time_deriv <- function(g_sym, u_expr, f_expr, t_expr) {
  n_g <- sym_length(g_sym)
  n_f <- sym_length(f_expr)
  terms <- lapply(seq_len(n_g), function(gi) {
    g_i   <- sym_elt(g_sym, gi)
    total <- sym_diff(g_i, t_expr)
    for (k in seq_len(n_f))
      total <- total + sym_diff(g_i, sym_elt(u_expr, k)) * sym_elt(f_expr, k)
    total
  })
  sym_build(terms)
}

lognormal_transform <- function(f_sym){
  D <- sym_length(f_sym)
  u    <- lapply(seq_len(D), function(i) sym_symbol(paste0("u", i)))
  expu <- lapply(u, exp)
  terms <- lapply(seq_len(D), function(i) {
    fud <- sym_elt(f_sym, i) / u[[i]]
    sym_subs(fud, u, expu)
  })
  sym_build(terms)
}
