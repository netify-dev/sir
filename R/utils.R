
#' @useDynLib sir, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom stats dpois dnorm dbinom qnorm rnorm rgamma rpois rbinom vcov model.matrix
NULL

# ---- Helper functions for data handling ----

#' Prepare Z Array for C++ Consumption
#'
#' @description
#' Converts a 4D array of exogenous covariates into a format optimized for C++ processing.
#' RcppArmadillo handles 4D arrays most efficiently when decomposed into a list of 3D cubes.
#'
#' @details
#' This function handles the dimensional restructuring needed for computation
#' in the C++ backend. The transformation preserves the covariate structure while
#' enabling vectorized operations in Armadillo.
#' 
#' Input dimensions:
#' \itemize{
#'   \item 3D array (m x m x T): Interpreted as single covariate (q=1)
#'   \item 4D array (m x m x q x T): Multiple covariates
#' }
#' 
#' Output structure:
#' \itemize{
#'   \item List of length q
#'   \item Each element is an (m x m x T) array for one covariate
#' }
#'
#' @param Z Array of exogenous covariates. Can be:
#'   \itemize{
#'     \item NULL: Returns empty list
#'     \item 3D array (m x m x T): Single time-varying covariate
#'     \item 4D array (m x m x q x T): Multiple covariates
#'   }
#'   
#' @return List of 3D arrays, where each element corresponds to one covariate
#'   across all time periods. Length equals q (number of covariates).
#'   
#' @examples
#' \dontrun{
#' # Single covariate
#' Z_single <- array(rnorm(10*10*5), dim=c(10,10,5))
#' Z_list <- prepare_Z_list(Z_single)
#' length(Z_list)  # Returns 1
#' 
#' # Multiple covariates  
#' Z_multi <- array(rnorm(10*10*3*5), dim=c(10,10,3,5))
#' Z_list <- prepare_Z_list(Z_multi)
#' length(Z_list)  # Returns 3
#' }
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

#' Flatten Y Array for GLM Input
#' 
#' @description
#' Converts a 3D network outcome array into a vector suitable for GLM estimation.
#' Uses column-major ordering (R's default) to ensure consistency with design matrices.
#' 
#' @details
#' The flattening follows R's column-major convention:
#' \itemize{
#'   \item Elements are ordered: Y[1,1,1], Y[2,1,1], ..., Y[m,1,1], Y[1,2,1], ...
#'   \item This matches how design matrices are constructed
#'   \item Preserves the correspondence between outcomes and covariates
#' }
#' 
#' @param Y Three-dimensional array (m x m x T) of network outcomes.
#'   
#' @return Numeric vector of length m*m*T containing flattened outcomes.
flatten_Y <- function(Y) {
    return(c(Y))
}

#' Flatten Z Array for GLM Input
#'
#' @description
#' Converts a 4D covariate array into a design matrix compatible with flattened Y.
#' Ensures proper alignment between outcomes and covariates for GLM estimation.
#'
#' @details
#' The flattening process maintains the correspondence between each Y[i,j,t] and
#' its associated covariates Z[i,j,:,t]. The resulting matrix has:
#' \itemize{
#'   \item Rows: One for each outcome observation (m*m*T total)
#'   \item Columns: One for each covariate (q total)
#'   \item Order: Matches the flattening of Y (column-major)
#' }
#' 
#' Special handling for 3D input (single covariate):
#' \itemize{
#'   \item Automatically detected and converted to column matrix
#'   \item Preserves any dimension names for interpretability
#' }
#'
#' @param Z Covariate array. Can be:
#'   \itemize{
#'     \item NULL: Returns NULL
#'     \item 3D array (m x m x T): Converted to (m*m*T) x 1 matrix
#'     \item 4D array (m x m x q x T): Converted to (m*m*T) x q matrix
#'   }
#'   
#' @return Design matrix with dimensions (m*m*T) x q, or NULL if Z is NULL.
#'   Column names are preserved from dimension names or auto-generated.
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

#' Calculate Linear Predictor (eta) for SIR Model
#'
#' @description
#' Computes the linear predictor for the Social Influence Regression model,
#' combining exogenous effects and network influence through bilinear terms.
#'
#' @details
#' The linear predictor has two components:
#' 
#' \deqn{\eta_{i,j,t} = \theta^T z_{i,j,t} + \sum_{k,l} X_{k,l,t} A_{i,k} B_{j,l}}
#' 
#' Where:
#' \itemize{
#'   \item First term: Linear effect of exogenous covariates
#'   \item Second term: Bilinear network influence effect
#' }
#' 
#' The influence matrices are parameterized as:
#' \itemize{
#'   \item \eqn{A = I + \sum_{r=1}^{p-1} \alpha_r W_r} (sender effects)
#'   \item \eqn{B = \beta_0 I + \sum_{r=1}^{p} \beta_r W_r} (receiver effects)
#' }
#' 
#' Note: The first alpha is fixed at 1 for identifiability.
#' 
#' \strong{Computational Strategy:}
#' \itemize{
#'   \item Uses optimized C++ routines for matrix products
#'   \item Exploits sparsity when present
#'   \item Minimizes memory allocation through in-place operations
#' }
#'
#' @param tab Numeric vector of parameters. For the standard model, ordered as
#'   [theta, alpha_2:p, beta_1:p] (length q + 2p - 1). For fix_receiver mode,
#'   ordered as [theta, alpha_1:p] (length q + p).
#'
#' @param W Three-dimensional array (m x m x p) of influence covariates.
#'   Each slice W[,,r] parameterizes the influence structure.
#'   If NULL or p=0, no network influence is included.
#'
#' @param X Three-dimensional array (m x m x T) carrying network influence.
#'   Typically lagged outcomes. Required even if W is NULL (can be zeros).
#'
#' @param Z Array of exogenous covariates. Can be:
#'   \itemize{
#'     \item NULL or q=0: No exogenous effects
#'     \item 3D array (m x m x T): Single covariate
#'     \item 4D array (m x m x q x T): Multiple covariates
#'   }
#'
#' @param fix_receiver Logical. If TRUE, B is fixed to the identity matrix and
#'   tab is parsed as [theta, alpha_1:p] with all alpha free. Default FALSE.
#'
#' @return Three-dimensional array (m x m x T) of linear predictors.
#'   Each element eta[i,j,t] is the linear predictor for outcome Y[i,j,t]
#'   before applying the link function.
#'   
#' @examples
#' \dontrun{
#' # Setup
#' m <- 10; T <- 5; p <- 2; q <- 1
#' W <- array(rnorm(m*m*p), dim=c(m,m,p))
#' X <- array(rnorm(m*m*T), dim=c(m,m,T))
#' Z <- array(rnorm(m*m*q*T), dim=c(m,m,q,T))
#' 
#' # Parameter vector
#' tab <- c(0.5,      # theta (q=1)
#'          0.2,      # alpha_2 (alpha_1=1 fixed)
#'          0.3, 0.4) # beta_1, beta_2
#'          
#' # Compute linear predictor
#' eta <- eta_tab(tab, W, X, Z)
#' dim(eta)  # Returns c(10, 10, 5)
#' }
#' @export
eta_tab <- function(tab, W, X, Z, fix_receiver=FALSE) {
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

  if (fix_receiver && p > 0) {
      # fix_receiver: tab = [theta, alpha_1:p], B = I
      alpha <- tab[(q+1):(q+p)]
  } else if (p > 0) {
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
  if (fix_receiver && p > 0) {
      A <- cpp_amprod_W_v(W, alpha)
      B <- diag(m)
      AXB <- cpp_tprod_A_X_Bt(X, A, B)
  } else if (p > 0) {
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

#' Calculate Negative Log-Likelihood for SIR Model
#'
#' @description
#' Computes the negative log-likelihood for the Social Influence Regression model
#' under the specified distributional family. Used as the objective function
#' for parameter estimation.
#'
#' @details
#' The likelihood depends on the distributional family:
#' 
#' \strong{Poisson:} For count outcomes
#' \deqn{L = \prod_{i,j,t} \frac{e^{-\lambda_{ijt}} \lambda_{ijt}^{y_{ijt}}}{y_{ijt}!}}
#' where \eqn{\lambda_{ijt} = \exp(\eta_{ijt})}
#' 
#' \strong{Normal:} For continuous outcomes
#' \deqn{L = \prod_{i,j,t} \frac{1}{\sqrt{2\pi\sigma^2}} \exp\left(-\frac{(y_{ijt} - \mu_{ijt})^2}{2\sigma^2}\right)}
#' where \eqn{\mu_{ijt} = \eta_{ijt}} (identity link)
#' 
#' \strong{Binomial:} For binary outcomes
#' \deqn{L = \prod_{i,j,t} p_{ijt}^{y_{ijt}} (1-p_{ijt})^{1-y_{ijt}}}
#' where \eqn{p_{ijt} = \frac{1}{1 + \exp(-\eta_{ijt})}} (logit link)
#' 
#' \strong{Numerical Stability:}
#' \itemize{
#'   \item Uses log-space computations to avoid underflow
#'   \item Bounds probabilities away from 0 and 1 for binomial
#'   \item Caps extreme values of lambda for Poisson
#'   \item Handles NA values by exclusion from likelihood
#' }
#'
#' @param tab Numeric vector of parameters [theta, alpha_2:p, beta_1:p].
#'   
#' @param Y Three-dimensional array (m x m x T) of observed outcomes.
#'   Can contain NA values which are automatically excluded.
#'   
#' @param W Three-dimensional array (m x m x p) of influence covariates,
#'   or NULL for no network influence.
#'   
#' @param X Three-dimensional array (m x m x T) carrying network influence.
#'   
#' @param Z Array of exogenous covariates (3D or 4D), or NULL.
#'   
#' @param family Character string specifying the distribution:
#'   \itemize{
#'     \item "poisson": Count data with log link
#'     \item "normal": Continuous data with identity link (assumes sigma=1)
#'     \item "binomial": Binary data with logit link
#'   }
#'
#' @param fix_receiver Logical. If TRUE, B is fixed to identity and tab is
#'   parsed as [theta, alpha_1:p]. Default FALSE.
#'
#' @return Numeric scalar giving the negative log-likelihood.
#'   Lower values indicate better fit. Used for optimization.
#'
#' @note The normal family assumes unit variance (sigma=1) for simplicity.
#'   The actual variance is estimated separately if needed.
#'   
#' @examples
#' \dontrun{
#' # Poisson example
#' Y <- array(rpois(1000, 2), dim=c(10,10,10))
#' nll <- mll_sir(tab, Y, W, X, Z, "poisson")
#' 
#' # Binomial example with missing data
#' Y_binary <- array(rbinom(1000, 1, 0.3), dim=c(10,10,10))
#' Y_binary[1,1,1] <- NA  # Missing value
#' nll <- mll_sir(tab, Y_binary, W, X, Z, "binomial")
#' }
#' @export
mll_sir <- function(tab, Y, W, X, Z, family, fix_receiver=FALSE) {
  ETA <- eta_tab(tab, W, X, Z, fix_receiver=fix_receiver)

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

#' Cast Directed Dyadic Data into Array Format
#' 
#' @description
#' Transforms long-format dyadic data (edge list) into a 3D array suitable for
#' SIR model analysis. Handles both dyadic and monadic variables with proper
#' placement in the network structure.
#' 
#' @details
#' This function converts network data from "long" format (one row per edge/time)
#' to "array" format (3D array indexed by sender, receiver, time).
#' 
#' \strong{Expected Input Format:}
#' Data frame with columns:
#' \itemize{
#'   \item i: Sender node identifier
#'   \item j: Receiver node identifier  
#'   \item t: Time period
#'   \item var: The variable value for this edge
#' }
#' 
#' \strong{Monadic Variables:}
#' When monadic=TRUE, the function treats the variable as node-level rather than
#' edge-level. The values are placed on the diagonal of each time slice:
#' \itemize{
#'   \item row=TRUE: Uses sender (i) attributes
#'   \item row=FALSE: Uses receiver (j) attributes
#' }
#' 
#' \strong{Missing Data:}
#' Missing edges in the input are filled with zeros in the output array.
#' This assumes that absence of an edge means zero interaction.
#' 
#' @param dyadData Data frame in long format with columns i, j, t, and the
#'   value variable specified by var parameter.
#'   
#' @param var Character string naming the column containing the values to
#'   be placed in the array.
#'   
#' @param monadic Logical indicating whether the variable is monadic (node-level)
#'   rather than dyadic (edge-level). Default is FALSE.
#'   
#' @param row Logical, only used when monadic=TRUE. If TRUE, uses sender
#'   attributes; if FALSE, uses receiver attributes. Default is FALSE.
#'   
#' @return Three-dimensional array (m x m x T) where:
#'   \itemize{
#'     \item First dimension: Sender nodes (i)
#'     \item Second dimension: Receiver nodes (j)
#'     \item Third dimension: Time periods (t)
#'   }
#'   Missing edges are filled with zeros.
#'   
#' @examples
#' \dontrun{
#' # Create example dyadic data
#' dyad_data <- data.frame(
#'   i = c(1,1,2,2,3,3),
#'   j = c(2,3,1,3,1,2),
#'   t = c(1,1,1,1,1,1),
#'   trade = c(100,150,80,120,90,110)
#' )
#' 
#' # Convert to array
#' trade_array <- castArray(dyad_data, "trade")
#' dim(trade_array)  # 3 x 3 x 1
#' 
#' # Monadic example (node GDP)
#' node_data <- data.frame(
#'   i = rep(1:3, each=3),
#'   j = rep(1:3, 3),
#'   t = 1,
#'   gdp = rep(c(1000,2000,1500), each=3)
#' )
#' 
#' gdp_array <- castArray(node_data, "gdp", monadic=TRUE, row=TRUE)
#' diag(gdp_array[,,1])  # Shows node GDPs on diagonal
#' }
#' 
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
