#' Build a WENDy problem from a dataset.
#'
#' WENDyProblem is constructed per dataset and they are later stacked by
#' build_wendy_system().
#'
#' @param wendy_data  List with elements U (mp1 x D), tt (mp1-vector or column
#'   matrix), var (mp1 x D variance matrix).
#' @param f_,J_u,J_up,J_p,J_pp,J_upp  Callable symbolic derivatives.
#' @param J       Number of parameters.
#' @param lip     Logical; TRUE when f is linear in parameters.
#' @param sig     Torch scalar: estimated noise SD.
#' @param device  Torch device.
#' @param control Control list (uses $test_fun_type, $scale_by_var, etc.).
#' @return A list of class "WENDyProblem".
#' @keywords internal
build_wendy_problem <- function(wendy_data, f_, J_u, J_up, J_p, J_pp, J_upp, J, lip, sig, device, control) {
  U   <- wendy_data$U
  tt  <- as.vector(wendy_data$tt)
  var <- wendy_data$var
  D   <- ncol(U)
  mp1 <- nrow(U)

  # Ensure sig is a D-vector for einsum 'a' indices; expand scalar if needed
  if (sig$numel() != D) {
    sig <- sig * torch::torch_ones(D, dtype = torch::torch_float64(), device = device)
  }

  tf <- if (control$test_fun_type == "SSL") {
    build_full_test_function_matrices_ssl(U, tt, control)
  } else {
    build_full_test_function_matrices_msg(U, tt, control, control$compute_svd)
  }
  K   <- nrow(tf$V)

  V  <- torch::torch_tensor(tf$V,       dtype = torch::torch_float64(), device = device)
  Vp <- torch::torch_tensor(tf$V_prime, dtype = torch::torch_float64(), device = device)

  F_ <- build_F(U, tt, f_, J, device)

  G  <- build_G_matrix(V, U, tt, F_, J, device)
  g0 <- torch::torch_mm(V, F_(rep(0, J)))$reshape(c(-1))

  g <- if (!lip) {
    local({ Vt <- V; Fl <- F_
      function(p) torch::torch_matmul(Vt, Fl(p))$reshape(c(-1))
    })
  } else {
    build_g_linear(G, device)
  }

  U_tensor <- torch::torch_tensor(U, dtype = torch::torch_float64(), device = device)
  b_raw    <- -1 * torch::torch_mm(Vp, U_tensor)$reshape(c(-1))
  b        <- if (!lip) b_raw else b_raw - g0

  Jp_r <- if (!lip) build_Jp_r(J_p,  K, D, J, mp1, V, U, tt, device)
          else      build_Jp_r_linear(G)

  Hp_r <- build_Hp_r(J_pp, K, D, J, mp1, V, U, tt, device)

  L0 <- build_L0(K, D, mp1, Vp, sig, device)

  L <- if (!lip) build_L(      U, tt, J_u,  K, V, L0, sig, J, device)
       else      build_L_linear(U, tt, J_u,  K, V, L0, sig, J, device)

  Jp_L <- if (!lip) build_Jp_L(      U, tt, J_up, K, J, D, V,  sig, device)
          else      build_Jp_L_linear(U, tt, J_u,  K, V, L0, sig, J, device)

  Hp_L <- build_Hp_L(U, tt, J_upp, K, J, D, V, sig, device)

  W <- if (!is.null(control$scale_by_var) && !is.na(control$scale_by_var)) {
    var_vec <- c(t(var))
    torch::torch_diag(torch::torch_tensor(var_vec, dtype = torch::torch_float64(), device = device))
  } else {
    NULL
  }

  structure(
    list(
      # Data
      U = U, tt = tt, var = var,
      # Test-function tensors and dimensions
      V = V, Vp = Vp, K = K, D = D, mp1 = mp1, J = J, min_radius = tf$min_radius,
      # Residual components
      F_ = F_, g = g, g0 = g0, b = b, G = G,
      Jp_r = Jp_r, Hp_r = Hp_r,
      # Covariance factor components
      L0 = L0, L = L, Jp_L = Jp_L, Hp_L = Hp_L,
      # Variance weights
      W = W
    ),
    class = "WENDyProblem"
  )
}


#' Assemble a list of WENDy problems into system
#'
#' Returns a plain list with everything the optimizers need:
#'   g, b, G, Jp_r, Hp_r, L, Jp_L, Hp_L, W,
#'   S, Jp_S, wnll, J_wnll, H_wnll, K, D, J.
#'
#' @param wendy_problems  List of WENDyProblem objects.
#' @param lip             Logical; TRUE when f is linear in parameters.
#' @param diag_reg        Diagonal regularisation added to S.
#' @param scale_by_var    Logical; if TRUE build block-diagonal W from per-problem variances.
#' @param device          Torch device.
#' @keywords internal
build_wendy_system <- function(wendy_problems, lip, diag_reg, scale_by_var, device) {
  J      <- wendy_problems[[1]]$J
  D      <- wendy_problems[[1]]$D
  mp1    <- wendy_problems[[1]]$mp1
  K_list <- sapply(wendy_problems, `[[`, "K")
  K      <- sum(K_list)
  M      <- length(wendy_problems)

  if (M == 1L) {
    prob <- wendy_problems[[1]]
    g    <- prob$g
    b    <- prob$b
    G    <- prob$G
    Jp_r <- prob$Jp_r
    Hp_r <- prob$Hp_r
    L    <- prob$L
    Jp_L <- prob$Jp_L
    Hp_L <- prob$Hp_L
    W    <- prob$W
  } else {
    g_fns <- lapply(wendy_problems, `[[`, "g")
    g     <- local({ gs <- g_fns
      function(p) torch::torch_cat(lapply(gs, function(gi) gi(p)))
    })

    b <- torch::torch_cat(lapply(wendy_problems, `[[`, "b"))
    G <- torch::torch_cat(lapply(wendy_problems, `[[`, "G"), dim = 1L)

    Jp_r_fns <- lapply(wendy_problems, `[[`, "Jp_r")
    Jp_r     <- local({ 
      fns <- Jp_r_fns
      function(p){ 
        torch::torch_cat(lapply(fns, function(fn) fn(p)), dim = 1L)
      }
    })

    Hp_r_fns <- lapply(wendy_problems, `[[`, "Hp_r")
    Hp_r     <- local({ 
      fns <- Hp_r_fns
      function(p){ 
        torch::torch_cat(lapply(fns, function(fn) fn(p)), dim = 1L)
      }
    })

    L    <- build_L_block(lapply(wendy_problems, `[[`, "L"),    K_list, D, mp1,    device)
    Jp_L <- build_Jp_L_block(lapply(wendy_problems, `[[`, "Jp_L"), K_list, D, mp1, J, device)
    Hp_L <- build_Hp_L_block(lapply(wendy_problems, `[[`, "Hp_L"), K_list, D, mp1, J, device)

    W <- if (!is.null(scale_by_var) && !is.na(scale_by_var)) {
      var_vec <- unlist(lapply(wendy_problems, function(pr) c(t(pr$var))))
      torch::torch_diag(torch::torch_tensor(var_vec, dtype = torch::torch_float64(), device = device))
    } else {
      NULL
    }
  }

  S      <- build_S(L, W, diag_reg = diag_reg)
  Jp_S   <- build_J_S(L, Jp_L, J, K, D, W, diag_reg = diag_reg)
  wnll   <- build_wnll(S, g, b, K, D)
  J_wnll <- build_J_wnll(S, Jp_S, Jp_r, g, b, J)
  H_wnll <- if (!lip) {
    build_H_wnll(S, Jp_S, L, Jp_L, Hp_L, Jp_r, Hp_r, g, b, J, W, diag_reg = diag_reg)
  } else {
    build_H_wnll_linear(S, Jp_S, L, Jp_L, Jp_r, g, b, J, W, diag_reg = diag_reg)
  }

  list(
    g = g, b = b, G = G,
    Jp_r = Jp_r, Hp_r = Hp_r,
    L = L, Jp_L = Jp_L, Hp_L = Hp_L, W = W,
    S = S, Jp_S = Jp_S,
    wnll = wnll, J_wnll = J_wnll, H_wnll = H_wnll,
    K = K, D = D, J = J
  )
}
