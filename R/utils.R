
#' @useDynLib sir, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom stats dpois dnorm dbinom qnorm rnorm rgamma rpois rbinom runif vcov model.matrix
NULL

# helper functions for data handling

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
	# (m x m x T) with q=1
	T_len <- dims[3]
	# keep 3D even if T=1 or m=1
	if (is.null(dim(Z))) {
		Z <- array(Z, dim=c(m, m, T_len))
	}
	return(list(Z))
  }

  if (ndims != 4) {
	cli::cli_abort("Z must be a 3D or 4D array.")
  }

  # Z is (m x m x q x T)
  q <- dims[3]
  T_len <- dims[4]

  Z_list <- lapply(1:q, function(k) {
	  # extract (m x m x T) cube for covariate k
	  Zk <- Z[,,k,]
	  # keep 3D even if T=1 or m=1
	  if (is.null(dim(Zk)) || length(dim(Zk)) < 3) {
		  Zk <- array(Zk, dim=c(m, m, T_len))
	  }
	  return(Zk)
  })

  return(Z_list)
}

# convert 4D W array (m x m x p x T) to a list of T cubes for C++
prepare_W_field <- function(W) {
	T_len <- dim(W)[4]
	m <- dim(W)[1]
	p <- dim(W)[3]
	lapply(seq_len(T_len), function(t) {
		Wt <- W[,,,t]
		if (is.null(dim(Wt)) || length(dim(Wt)) < 3) {
			Wt <- array(Wt, dim = c(m, m, p))
		}
		Wt
	})
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
	# permute to (m, m, T, q) for correct flattening order
	Z_perm <- aperm(Z, c(1, 2, 4, 3))
	Z_flat <- matrix(Z_perm, ncol = q)

	cnames <- if (!is.null(dimnames(Z)[[3]])) dimnames(Z)[[3]] else paste0("Z", 1:q)
	colnames(Z_flat) <- cnames
	return(Z_flat)
}


# core model components (R wrappers around C++)

#' Calculate Linear Predictor (eta) for SIR Model
#'
#' @description
#' Computes the linear predictor \eqn{\eta_{i,j,t} = \theta^T z_{i,j,t} +
#' \sum_{k,l} X_{k,l,t} A_{i,k} B_{j,l}}, combining exogenous covariate
#' effects with the bilinear network influence term. Handles both static
#' (3D) and dynamic (4D) W arrays.
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
  dynamic_W <- !is.null(W) && length(dim(W)) == 4
  p <- if (is.null(W)) 0 else dim(W)[3]
  q <- if (is.null(Z)) 0 else if (length(dim(Z)) == 3) 1 else dim(Z)[3]
  n1 <- dim(X)[1]
  n2 <- dim(X)[2]
  m <- n1
  T_len <- dim(X)[3]

  # parse parameters
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
		  alpha_start <- q + 1
		  alpha_end <- q + p - 1
		  alpha <- c(1, tab[alpha_start:alpha_end])
	  } else {
		  alpha <- c(1)
	  }
	  beta_start <- q + p
	  beta <- tab[beta_start:length(tab)]
  }

  # bilinear part: AXB
  if (dynamic_W && p > 0) {
	  # dynamic W: compute A_t * X_t * B_t' per period
	  AXB <- array(0, dim = c(n1, n2, T_len))
	  for (t in seq_len(T_len)) {
		  W_t <- W[,,,t]
		  A_t <- cpp_amprod_W_v(W_t, alpha)
		  if (fix_receiver) {
			  AXB[,,t] <- A_t %*% X[,,t]
		  } else {
			  B_t <- cpp_amprod_W_v(W_t, beta)
			  AXB[,,t] <- A_t %*% X[,,t] %*% t(B_t)
		  }
	  }
  } else if (fix_receiver && p > 0) {
	  A <- cpp_amprod_W_v(W, alpha)
	  B <- diag(n2)
	  AXB <- cpp_tprod_A_X_Bt(X, A, B)
  } else if (p > 0) {
	  A <- cpp_amprod_W_v(W, alpha)
	  B <- cpp_amprod_W_v(W, beta)
	  AXB <- cpp_tprod_A_X_Bt(X, A, B)
  } else {
	  AXB <- array(0, dim=c(n1, n2, T_len))
  }

  # exogenous part: ZT
  if (q > 0) {
	if (length(dim(Z)) == 3) {
		Z <- array(Z, dim=c(n1, n2, 1, T_len))
	}
	ZT  <- amprod(Z, matrix(theta, nrow=1), 3)
	ZT <- array(ZT, dim=c(n1, n2, T_len))
  } else {
	ZT <- array(0, dim=c(n1, n2, T_len))
  }

  ZT + AXB
}

#' Calculate Negative Log-Likelihood for SIR Model
#'
#' @description
#' Computes the negative log-likelihood for the SIR model under the specified
#' distributional family. Diagonal entries (self-ties) and NA values are
#' excluded. Used as the objective function for parameter estimation.
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

  # exclude diagonal (self-ties) to match C++ gradient/hessian which skips i == j
  n1 <- dim(Y)[1]
  n2 <- dim(Y)[2]
  T_len <- dim(Y)[3]
  if (n1 == n2) {
	for (t in seq_len(T_len)) {
	  diag(Y[,,t]) <- NA
	  diag(ETA[,,t]) <- NA
	}
  }

  if (family == "poisson") {
	lambda <- exp(ETA)
	lambda[lambda > 1e300] <- 1e300
	# guard against lambda underflowing to 0
	lambda[lambda < 1e-300] <- 1e-300
	nll <- -sum(dpois(Y, lambda = lambda, log = TRUE), na.rm = TRUE)
  } else if (family == "normal") {
	nll <- -sum(dnorm(Y, mean = ETA, sd = 1, log = TRUE), na.rm = TRUE)
  } else if (family == "binomial") {
	prob <- 1 / (1 + exp(-ETA))
	prob[prob < 1e-15] <- 1e-15
	prob[prob > 1 - 1e-15] <- 1 - 1e-15
	nll <- -sum(dbinom(Y, size = 1, prob = prob, log = TRUE), na.rm = TRUE)
  } else {
	cli::cli_abort("Unsupported family in {.fn mll_sir}: {.val {family}}.")
  }

  return(nll)
}


# r fallback tensor functions (used for 4D Z calculations not covered by C++)

#' Matricization (R implementation)
#' @keywords internal
mat<-function(A,k)
{
  # handle vector case
  if (is.null(dim(A))) {
	  if (k==1) return(matrix(A, ncol=1))
	  else cli::cli_abort("Invalid mode for vector in {.fn mat}.")
  }

  Ak<-t(apply(A,k,"c"))
  # ensure Ak is a matrix, handle dimension issues
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
  if(is.vector(M)) { M<-matrix(M,nrow=1) } # treat vectors as row matrices

  K<-length(dim(A))
  if(is.null(K)) { # a is a vector
	  if (k!=1) cli::cli_abort("Invalid mode for vector in {.fn amprod}.")
	  AM <- M %*% A
	  return(c(AM)) # Return vector
  }

  A_mat <- mat(A,k)

  if (ncol(M) != nrow(A_mat)) {
	  # handle transposition issues if dim(A)[k] is 1
	  if (ncol(M) == ncol(A_mat) && nrow(A_mat) == 1 && dim(A)[k] == 1) {
		   A_mat <- t(A_mat)
	  } else {
		cli::cli_abort("Dimension mismatch in {.fn amprod}: k={.val {k}}, dim(M)={.val {paste(dim(M), collapse='x')}}, dim(A)[k]={.val {dim(A)[k]}}, nrow(A_mat)={.val {nrow(A_mat)}}.")
	  }
  }

  AM<-M %*% A_mat

  # determine new dimensions
  dims_A <- dim(A)
  new_dims <- c(dim(M)[1], dims_A[-k])

  if (length(new_dims) <= 1) return(c(AM)) # Return vector if result is 1D or scalar

  AMA<-array(AM, dim=new_dims)

  # handle permutation
  if (K > 1) {
	# calculate permutation vector
	perm_order <- c(k, (1:K)[-k])
	# inverse permutation to restore original dimension order
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

# data prep helpers

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
#' @param dyad_data Data frame in long format with columns i, j, t, and the
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
#' trade_array <- cast_array(dyad_data, "trade")
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
#' gdp_array <- cast_array(node_data, "gdp", monadic=TRUE, row=TRUE)
#' diag(gdp_array[,,1])  # Shows node GDPs on diagonal
#' }
#' 
#' @export
cast_array <- function(dyad_data, var, monadic=FALSE, row=FALSE){
	# reshape long-format dyadic data into 3D array
	i_levels <- sort(unique(as.character(dyad_data$i)))
	j_levels <- sort(unique(as.character(dyad_data$j)))
	t_levels <- sort(unique(as.character(dyad_data$t)))
	arr <- array(0, dim = c(length(i_levels), length(j_levels), length(t_levels)),
				 dimnames = list(i_levels, j_levels, t_levels))
	i_idx <- match(as.character(dyad_data$i), i_levels)
	j_idx <- match(as.character(dyad_data$j), j_levels)
	t_idx <- match(as.character(dyad_data$t), t_levels)
	for (row_n in seq_len(nrow(dyad_data))) {
		arr[i_idx[row_n], j_idx[row_n], t_idx[row_n]] <- dyad_data[[var]][row_n]
	}
	arr[is.na(arr)] <- 0
	if(monadic){
		# place monadic variable on diagonal
		if(row){
			for(t in 1:dim(arr)[3]){
				diag_vals <- apply(arr[,,t], 1, function(x){
					val <- setdiff(unique(x),0)
					if (length(val) > 0) return(val[1]) else return(0)
				})
				diag(arr[,,t]) <- diag_vals
			}
		} else {
			for(t in 1:dim(arr)[3]){
				diag_vals <- apply(arr[,,t], 2, function(x){
					val <- setdiff(unique(x),0)
					if (length(val) > 0) return(val[1]) else return(0)
				})
				diag(arr[,,t]) <- diag_vals
			}
		}
	}
	return(arr)
}

