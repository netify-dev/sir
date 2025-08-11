
#' @useDynLib sir, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom stats dpois dnorm dbinom qnorm rnorm rgamma
NULL

# ---- Helper functions for data handling ----

#' Helper to prepare Z (4D array) for C++ consumption
#'
#' C++ (RcppArmadillo) handles 4D arrays best when passed as a list of 3D cubes.
#' We convert Z (m x m x q x T) into a list of q (m x m x T) cubes.
#'
#' @param Z The 4D array of exogenous covariates.
#' @return A list of 3D arrays.
prepare_Z_list <- function(Z) {
  if (is.null(Z) || length(dim(Z)) == 0) {
    return(list())
  }

  dims <- dim(Z)
  ndims <- length(dims)
  m <- dims[1]

  if (ndims == 3) {
    # Assume (m x m x T) and q=1
    T_len <- dims[3]
    # Ensure it remains 3D even if T=1 or m=1
    if (is.null(dim(Z))) {
        Z <- array(Z, dim=c(m, m, T_len))
    }
    return(list(Z))
  }

  if (ndims != 4) {
    stop("Z must be a 3D or 4D array.")
  }

  # Z is (m x m x q x T)
  q <- dims[3]
  T_len <- dims[4]

  Z_list <- lapply(1:q, function(k) {
      # Extract the (m x m x T) cube for the k-th covariate
      Zk <- Z[,,k,]
      # Ensure it remains a 3D array even if T=1 or m=1
      if (is.null(dim(Zk)) || length(dim(Zk)) < 3) {
          Zk <- array(Zk, dim=c(m, m, T_len))
      }
      return(Zk)
  })

  return(Z_list)
}

#' Flatten Y for GLM input
#' @param Y (m x m x T) array.
#' @return Vectorized Y (column-major).
flatten_Y <- function(Y) {
    return(c(Y))
}

#' Flatten Z for GLM input
#'
#' Converts Z (m x m x q x T) into a (m*m*T) x q matrix, matching the flattening of Y.
#'
#' @param Z (m x m x q x T) array.
#' @return (m*m*T) x q matrix.
flatten_Z <- function(Z) {
    if (is.null(Z)) return(NULL)

    dims <- dim(Z)
    ndims <- length(dims)

    if (ndims == 3) {
        # (m x m x T), q=1
        Z_flat <- matrix(c(Z), ncol=1)
        cname <- if (!is.null(dimnames(Z)[[3]])) dimnames(Z)[[3]][1] else "Z1"
        colnames(Z_flat) <- cname
        return(Z_flat)
    }

    # (m x m x q x T)
    q <- dims[3]
    # Permute Z to (m, m, T, q) for correct flattening order when using c() or matrix()
    Z_perm <- aperm(Z, c(1, 2, 4, 3))
    Z_flat <- matrix(Z_perm, ncol = q)

    cnames <- if (!is.null(dimnames(Z)[[3]])) dimnames(Z)[[3]] else paste0("Z", 1:q)
    colnames(Z_flat) <- cnames
    return(Z_flat)
}


# ---- Core Model Components (R wrappers using Cpp) ----

#' Calculate the bilinear predictor (eta)
#'
#' Computes eta_{i,j,t} = Z_{i,j,t}^T theta + alpha^T X_{i,j,t} beta.
#' Uses optimized C++ functions for tensor products.
#'
#' @param tab Parameter vector (theta, alpha[-1], beta).
#' @param W (m x m x p) array of influence covariates.
#' @param X (m x m x T) array for the bilinear part.
#' @param Z (m x m x q x T) array of exogenous covariates.
#' @return An (m x m x T) array of linear predictors.
#' @export
eta_tab <- function(tab, W, X, Z) {
  p <- if (is.null(W)) 0 else dim(W)[3]
  q <- if (is.null(Z)) 0 else dim(Z)[3]
  m <- dim(X)[1]
  T_len <- dim(X)[3]

  # Parse parameters
  if (q > 0) {
      theta <- tab[1:q]
  } else {
      theta <- numeric(0)
  }

  if (p > 0) {
      if (p > 1) {
          alpha_start = q + 1
          alpha_end = q + p - 1
          alpha <- c(1, tab[alpha_start:alpha_end])
      } else {
          alpha <- c(1)
      }
      beta_start = q + p
      beta <- tab[beta_start:length(tab)]
  }


  # Bilinear part: AXB
  if (p > 0) {
      # Build A, B using optimized C++ amprod wrapper
      A <- cpp_amprod_W_v(W, alpha)
      B <- cpp_amprod_W_v(W, beta)

      # Calculate AXB using optimized C++ tprod wrapper
      AXB <- cpp_tprod_A_X_Bt(X, A, B)
  } else {
      AXB <- array(0, dim=c(m, m, T_len))
  }

  # Exogenous part: ZT
  if (q > 0) {
    # We use the R implementation of amprod for the 4D Z contraction,
    # as the C++ optimization focused on the 3D cases.
    # Ensure Z is 4D before calling amprod on mode 3
    if (length(dim(Z)) == 3) {
        Z <- array(Z, dim=c(m, m, 1, T_len))
    }
    ZT  <- amprod(Z, matrix(theta, nrow=1), 3)
    # Squeeze the result dimension (m x m x 1 x T) -> (m x m x T)
    ZT <- array(ZT, dim=c(m, m, T_len))

  } else {
    ZT <- array(0, dim=c(m, m, T_len))
  }

  ZT + AXB
}

#' Calculate Negative Log-Likelihood
#'
#' @param tab Parameter vector.
#' @param Y Outcome array.
#' @param W Influence covariates.
#' @param X Bilinear covariates.
#' @param Z Exogenous covariates.
#' @param family Distribution family.
#' @return The negative log-likelihood value.
#' @export
mll_sir <- function(tab, Y, W, X, Z, family) {
  ETA <- eta_tab(tab, W, X, Z)

  if (family == "poisson") {
    # Use dpois for numerical stability
    # Ensure lambda > 0
    lambda <- exp(ETA)
    # Stabilization for very large lambda
    lambda[lambda > 1e300] <- 1e300
    nll <- -sum(dpois(Y, lambda = lambda, log = TRUE), na.rm = TRUE)
  } else if (family == "normal") {
    # Use dnorm (assuming sd=1, as the fitting procedure estimates sigma^2 separately if needed)
    nll <- -sum(dnorm(Y, mean = ETA, sd = 1, log = TRUE), na.rm = TRUE)
  } else if (family == "binomial") {
    # Use dbinom for stability
    prob <- 1/(1+exp(-ETA))
    # Stabilization for probabilities near 0 or 1
    prob[prob < 1e-15] <- 1e-15
    prob[prob > 1 - 1e-15] <- 1 - 1e-15
    nll <- -sum(dbinom(Y, size = 1, prob = prob, log = TRUE), na.rm = TRUE)
  } else {
    stop("Unsupported family in mll_sir.")
  }

  return(nll)
}


# ---- Functions from tfunctions.r (R Fallbacks) ----
# Including R implementations of tensor functions for use cases not covered by the optimized Cpp (e.g. 4D ZT calculation)

#' Matricization (R implementation)
#' @keywords internal
mat<-function(A,k)
{
  # Handle vector case
  if (is.null(dim(A))) {
      if (k==1) return(matrix(A, ncol=1))
      else stop("Invalid mode for vector.")
  }

  Ak<-t(apply(A,k,"c"))
  # Ensure Ak is a matrix and handle dimension issues during transposition
  if(!is.matrix(Ak)) Ak <- matrix(Ak, nrow=dim(A)[k])
  if(is.matrix(Ak) && nrow(Ak)!=dim(A)[k])  { Ak<-t(Ak) }
  Ak
}

#' Array-matrix product (R implementation)
#'
#' This R implementation is used for operations not covered by the specialized Cpp optimizations (like 4D Z).
#' @keywords internal
amprod<-function(A,M,k)
{
  if(is.vector(M)) { M<-matrix(M,nrow=1) } # Treat vectors as row matrices by default for this implementation

  K<-length(dim(A))
  if(is.null(K)) { # A is a vector
      if (k!=1) stop("Invalid mode for vector.")
      AM <- M %*% A
      return(c(AM)) # Return vector
  }

  A_mat <- mat(A,k)

  if (ncol(M) != nrow(A_mat)) {
      # Handle potential transposition issues if dim(A)[k] is 1
      if (ncol(M) == ncol(A_mat) && nrow(A_mat) == 1 && dim(A)[k] == 1) {
           A_mat <- t(A_mat)
      } else {
        # Provide detailed dimension info for debugging
        stop(paste("Dimension mismatch in R amprod. k=", k, "dim(M) =", paste(dim(M), collapse="x"), ", dim(A)[k] =", dim(A)[k], ", nrow(A_mat)=", nrow(A_mat)))
      }
  }

  AM<-M %*% A_mat

  # Determine new dimensions
  dims_A <- dim(A)
  new_dims <- c(dim(M)[1], dims_A[-k])

  if (length(new_dims) <= 1) return(c(AM)) # Return vector if result is 1D or scalar

  AMA<-array(AM, dim=new_dims)

  # Handle permutation
  if (K > 1) {
    # Calculate permutation vector
    perm_order <- c(k, (1:K)[-k])
    # We need the inverse permutation to restore the original order relative to the new dimensions
    inv_perm <- match(1:K, perm_order)

    if (any(is.na(inv_perm))) {
        return(AMA) # Should not happen
    }
    AMA <- aperm(AMA, inv_perm)
  }
  AMA
}

#' Tucker product (R implementation)
#' @keywords internal
tprod<-function(A,B,modes=1:length(B))
{
  X<-A
  for(i in seq_along(modes)) {
    k <- modes[i]
    M <- B[[i]]
    X<-amprod(X, M, k)
  }
  X
}

# ---- Functions from mltrHelpers.R (Data Prep) ----

#' Cast directed dyadic variable into array
#' @importFrom reshape2 acast
#' @export
castArray = function(dyadData, var, monadic=FALSE, row=FALSE){
	# Assumes dyadData has columns i, j, t, and the value.var
	arr = reshape2::acast(dyadData, i ~ j ~ t, value.var=var)
	arr[is.na(arr)] = 0
	if(monadic){
	    # Logic for handling monadic variables (placing them on the diagonal)
		if(row){
			for(t in 1:dim(arr)[3]){
			    # Handle case where setdiff might return empty vector
				diag_vals = apply(arr[,,t], 1, function(x){
				    val <- setdiff(unique(x),0)
				    if (length(val) > 0) return(val[1]) else return(0)
				})
				diag(arr[,,t]) = diag_vals
			}
		} else {
			for(t in 1:dim(arr)[3]){
				diag_vals = apply(arr[,,t], 2, function(x){
				    val <- setdiff(unique(x),0)
				    if (length(val) > 0) return(val[1]) else return(0)
				})
				diag(arr[,,t]) = diag_vals
			}
		}
	}
	return(arr)
}

# (createRelCovar and prepMLTR can be added here from the original scripts if desired by the user)
