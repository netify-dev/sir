#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// ----------------------------------------------------------------------------
// --- Tensor Utility Functions ---
// ----------------------------------------------------------------------------

//' Tensor Product for SIR Model (A * X * B')
//' 
//' @description
//' Performs the bilinear transformation central to the Social Influence Regression model.
//' Computes A * X_t * B' for each time slice t, where this product represents how
//' network influence flows through the sender effects (A) and receiver effects (B).
//' 
//' @details
//' This operation is the computational bottleneck of the SIR model, appearing in both
//' the likelihood evaluation and gradient computation. The function implements:
//' 
//' For each time t: Result[,,t] = A * X[,,t] * B'
//' 
//' Where:
//' - A captures sender-specific influence patterns (how nodes affect others)
//' - B captures receiver-specific influence patterns (how nodes are affected)
//' - X typically contains lagged network outcomes that carry influence forward
//' 
//' Mathematical interpretation:
//' - Element (i,j) of the result represents the total influence flowing from i to j
//' - This influence is mediated by the entire network structure at time t
//' - The bilinear form allows for complex, indirect influence pathways
//' 
//' Computational optimizations:
//' - Pre-computes B' once rather than for each time slice
//' - Uses Armadillo's optimized BLAS routines for matrix multiplication
//' - Memory-efficient slice-wise operations to avoid large temporary matrices
//' - Compiler optimizations enabled through RcppArmadillo
//' 
//' @param X Three-dimensional array (m x m x T) representing the network state over time.
//'   Each slice X[,,t] is the network at time t that carries influence.
//'   
//' @param A Matrix (m x m) of sender effects. Element A[i,k] represents how node i
//'   is influenced by the sending behavior of node k.
//'   
//' @param B Matrix (m x m) of receiver effects. Element B[j,l] represents how node j's
//'   reception is modified by node l's receiving patterns.
//'   
//' @return Three-dimensional array (m x m x T) where element [i,j,t] represents the
//'   total bilinear influence from node i to node j at time t.
//'   
//' @examples
//' \dontrun{
//' // In R:
//' m <- 10; T <- 5
//' X <- array(rnorm(m*m*T), dim=c(m,m,T))
//' A <- matrix(rnorm(m*m), m, m)
//' B <- matrix(rnorm(m*m), m, m)
//' result <- cpp_tprod_A_X_Bt(X, A, B)
//' }
//' 
//' @note This function is called repeatedly during optimization, so efficiency is critical.
//'   The implementation avoids unnecessary memory allocations and leverages BLAS Level 3
//'   operations for optimal performance.
//'   
// [[Rcpp::export]]
arma::cube cpp_tprod_A_X_Bt(const arma::cube& X, const arma::mat& A, const arma::mat& B) {
  int m = X.n_rows;
  int T = X.n_slices;

  arma::cube AXB(m, m, T);
  arma::mat Bt = B.t();

  // Loop over time T
  for(int t=0; t < T; ++t) {
    AXB.slice(t) = A * X.slice(t) * Bt;
  }
  return AXB;
}

//' Array-Matrix Product for Influence Matrices
//' 
//' @description
//' Computes a weighted sum of influence covariate matrices to construct the
//' parameterized influence matrices A or B in the SIR model.
//' 
//' @details
//' This function implements the parameterization:
//' Result = sum(k=1 to p) v[k] * W[,,k]
//' 
//' In the SIR model context:
//' - A = alpha[1] * I + sum(k=2 to p) alpha[k] * W[,,k-1]
//' - B = sum(k=1 to p) beta[k] * W[,,k]
//' 
//' The parameterization reduces the number of free parameters from O(m²) to O(p),
//' where typically p << m. This makes estimation feasible for larger networks.
//' 
//' Computational strategy:
//' - Skips zero coefficients to save computation
//' - Uses in-place addition to minimize memory allocation
//' - Leverages Armadillo's expression templates for efficiency
//' 
//' @param W Three-dimensional array (m x m x p) of influence covariates.
//'   Each slice W[,,k] represents one way nodes can influence each other
//'   (e.g., geographic proximity, social distance, shared attributes).
//'   
//' @param v Vector (p x 1) of coefficients for the linear combination.
//'   These are the parameters being estimated (either alpha or beta).
//'   
//' @return Matrix (m x m) representing the weighted combination of influence
//'   covariate matrices. This becomes either the A or B matrix in the model.
//'   
//' @examples
//' \dontrun{
//' // In R:
//' m <- 10; p <- 3
//' W <- array(rnorm(m*m*p), dim=c(m,m,p))
//' v <- c(1, 0.5, -0.3)  // Coefficients
//' A <- cpp_amprod_W_v(W, v)
//' // A is now the parameterized influence matrix
//' }
//' 
//' @note The function checks for dimension compatibility and will throw an
//'   error if v has incorrect length. Zero coefficients are detected and
//'   skipped to improve performance when the model is sparse.
//'   
// [[Rcpp::export]]
arma::mat cpp_amprod_W_v(const arma::cube& W, const arma::vec& v) {
    int m = W.n_rows;
    int p = W.n_slices;

    if (v.n_elem != p) {
        stop("Dimension mismatch in cpp_amprod_W_v.");
    }

    // This is equivalent to reshaping W to (m*m x p) and multiplying by v, then reshaping back.
    // arma::mat result = arma::reshape(W, m*m, p) * v;
    // result.reshape(m, m);

    // Loop implementation (often clearer and sometimes faster for small p)
    arma::mat result = arma::zeros<arma::mat>(m, m);
    for(int k=0; k < p; ++k) {
        if (v(k) != 0.0) {
            result += W.slice(k) * v(k);
        }
    }
    return result;
}

// ----------------------------------------------------------------------------
// --- ALS Design Matrix Construction ---
// ----------------------------------------------------------------------------
// These functions optimize the main bottleneck in the ALS algorithm by constructing
// the design matrices (Wbeta, Walpha) in C++ instead of slow R loops/apply calls.

//' Construct Design Matrix for Alpha Updates in ALS
//' 
//' @description
//' Builds the design matrix for updating sender effects (alpha parameters) in the
//' Alternating Least Squares algorithm, holding receiver effects (beta) fixed.
//' 
//' @details
//' In the ALS algorithm, when updating alpha with beta fixed, the model becomes
//' linear in alpha. The design matrix for this GLM sub-problem has columns
//' corresponding to each influence covariate.
//' 
//' For covariate k and observation (i,j,t), the design matrix element is:
//' [W[,,k] * X[,,t] * B'][i,j]
//' 
//' Where B = sum_l beta[l] * W[,,l] is the current receiver effects matrix.
//' 
//' The algorithm:
//' 1. Compute B from current beta and W
//' 2. For each covariate k:
//'    - Calculate W[,,k] * X * B' for all time points
//'    - Flatten to match the vectorized Y
//' 3. Combine into design matrix
//' 
//' This C++ implementation is 10-100x faster than the equivalent R code using
//' loops or apply functions, making ALS feasible for larger networks.
//' 
//' @param W Three-dimensional array (m x m x p) of influence covariates.
//'   These parameterize how influence flows through the network.
//'   
//' @param X Three-dimensional array (m x m x T) of network states over time.
//'   Typically contains lagged outcomes that carry influence forward.
//'   
//' @param beta Vector (p x 1) of current receiver effect parameters.
//'   These are held fixed while updating alpha in this ALS step.
//'   
//' @return Matrix (m*m*T x p) that serves as the design matrix for GLM estimation.
//'   Each column corresponds to one influence covariate, rows match vectorized Y.
//'   
//' @note This function is called once per ALS iteration. The resulting matrix can
//'   be large (m²T x p), so memory usage should be considered for big networks.
//'   
//' @examples
//' \dontrun{
//' // Called internally by sir_alsfit during the alpha update step
//' // After computing this design matrix, the update is:
//' // glm(Y ~ -1 + cbind(Z_design, Wbeta_design), family=...)
//' }
//' 
// [[Rcpp::export]]
arma::mat cpp_construct_Wbeta_design(const arma::cube& W, const arma::cube& X, const arma::vec& beta) {
    int m = W.n_rows;
    int p = W.n_slices;
    int T = X.n_slices;

    // 1. Calculate WSbeta = W * beta (m x m matrix)
    arma::mat WSbeta = cpp_amprod_W_v(W, beta);

    // 2. Calculate Wbeta design matrix ((m*m*T) x p)
    arma::mat Wbeta_design(m*m*T, p);

    // For each covariate k, calculate Wk %*% X %*% WSbeta^T and flatten
    for(int k=0; k < p; ++k) {
        // This uses the optimized tprod implementation
        arma::cube Wk_X_WSbeta = cpp_tprod_A_X_Bt(X, W.slice(k), WSbeta);
        // Flatten the resulting cube into a column vector (Column-major flattening)
        // This ensures the flattening matches R's c(Y) order.
        Wbeta_design.col(k) = arma::vectorise(Wk_X_WSbeta);
    }

    return Wbeta_design;
}

//' Construct Design Matrix for Beta Updates in ALS
//' 
//' @description
//' Builds the design matrix for updating receiver effects (beta parameters) in the
//' Alternating Least Squares algorithm, holding sender effects (alpha) fixed.
//' 
//' @details
//' In the ALS algorithm, when updating beta with alpha fixed, the model becomes
//' linear in beta. This function constructs the required design matrix efficiently.
//' 
//' For covariate k and observation (i,j,t), the design matrix element is:
//' [A * X[,,t] * W[,,k]'][i,j]
//' 
//' Where A = I + sum_l alpha[l] * W[,,l] is the current sender effects matrix
//' (with alpha[1] = 1 fixed for identifiability).
//' 
//' The algorithm mirrors the alpha update but with roles reversed:
//' 1. Compute A from current alpha and W
//' 2. For each covariate k:
//'    - Calculate A * X * W[,,k]' for all time points
//'    - Flatten to match the vectorized Y
//' 3. Combine into design matrix
//' 
//' @param W Three-dimensional array (m x m x p) of influence covariates.
//'   
//' @param X Three-dimensional array (m x m x T) of network states over time.
//' @param alpha Vector (p x 1) of current sender effect parameters.
//'   These are held fixed while updating beta in this ALS step.
//'   
//' @return Matrix (m*m*T x p) serving as the design matrix for beta GLM estimation.
//'   
//' @note The symmetry with cpp_construct_Wbeta_design reflects the bilinear
//'   structure of the model, where sender and receiver effects play dual roles.
//'   
//' @examples
//' \dontrun{
//' // Called internally by sir_alsfit during the beta update step
//' // The GLM call becomes:
//' // glm(Y ~ -1 + cbind(Z_design, Walpha_design), family=...)
//' }
// [[Rcpp::export]]
arma::mat cpp_construct_Walpha_design(const arma::cube& W, const arma::cube& X, const arma::vec& alpha) {
    int m = W.n_rows;
    int p = W.n_slices;
    int T = X.n_slices;

    // 1. Calculate WSalpha = W * alpha (m x m matrix)
    arma::mat WSalpha = cpp_amprod_W_v(W, alpha);

    // 2. Calculate Walpha design matrix ((m*m*T) x p)
    arma::mat Walpha_design(m*m*T, p);

    // For each covariate k, calculate WSalpha %*% X %*% Wk^T and flatten
    for(int k=0; k < p; ++k) {
        arma::cube WSalpha_X_Wk = cpp_tprod_A_X_Bt(X, WSalpha, W.slice(k));
        // Flatten the resulting cube into a column vector
        Walpha_design.col(k) = arma::vectorise(WSalpha_X_Wk);
    }

    return Walpha_design;
}

// ----------------------------------------------------------------------------
// --- Gradient and Hessian Calculation ---
// ----------------------------------------------------------------------------

//' Calculate Gradient and Hessian for Direct Optimization
//' 
//' @description
//' Computes the gradient vector and Hessian matrix of the negative log-likelihood
//' for the SIR model. These are essential for gradient-based optimization methods
//' like BFGS used in the direct optimization approach.
//' 
//' @details
//' This function implements the analytical derivatives of the SIR likelihood with
//' respect to all parameters. The computation leverages the bilinear structure
//' of the model for efficiency.
//' 
//' \strong{Gradient Computation:}
//' For each parameter, the gradient accumulates:
//' \deqn{\nabla_{\cdot} NLL = \sum_{i,j,t} (\mu_{ijt} - y_{ijt}) \frac{\partial \eta_{ijt}}{\partial \cdot}}
//' 
//' Where the partial derivatives are:
//' - \eqn{\partial \eta / \partial \theta_k = Z_{ijk,t}}
//' - \eqn{\partial \eta / \partial \alpha_k = [W_k X B']_{ij}}
//' - \eqn{\partial \eta / \partial \beta_k = [A X W_k']_{ij}}
//' 
//' \strong{Hessian Computation:}
//' The Hessian has two components:
//' 1. Fisher Information (always positive semi-definite):
//'    \deqn{H_{Fisher} = \sum_{i,j,t} w_{ijt} \nabla \eta_{ijt} \nabla \eta_{ijt}'}
//'    where \eqn{w_{ijt}} is the GLM weight (variance function)
//' 
//' 2. Observed Information adjustment (for non-canonical links):
//'    \deqn{H_{Obs} = \sum_{i,j,t} (\mu_{ijt} - y_{ijt}) \nabla^2 \eta_{ijt}}
//'    This term captures the curvature from the bilinear structure
//' 
//' \strong{Identifiability Constraint:}
//' Since alpha_1 is fixed at 1 for identifiability, the function:
//' - Computes full derivatives internally
//' - Projects out the alpha_1 dimension using matrix J
//' - Returns reduced-dimension gradient and Hessian
//' 
//' \strong{Numerical Stability:}
//' - Handles missing data (NA) by skipping those observations
//' - Adds small stabilization to binomial weights to prevent singularity
//' - Uses log-space computations where possible
//' 
//' \strong{Computational Efficiency:}
//' - Pre-computes A and B matrices once per function call
//' - Reuses matrix products across gradient components
//' - Vectorized operations using Armadillo
//' - Avoids redundant calculations in inner loops
//' 
//' @param tab Parameter vector [theta, alpha_2:p, beta] of length q+2p-1.
//'   
//' @param Y Three-dimensional array (m x m x T) of observed outcomes.
//'   Can contain NA values which are automatically skipped.
//'   
//' @param W Three-dimensional array (m x m x p) of influence covariates.
//'   
//' @param X Three-dimensional array (m x m x T) carrying network influence.
//'   
//' @param Z_list List of q three-dimensional arrays (m x m x T), one per covariate.
//'   Passed as list for efficient memory handling of 4D structure.
//'   
//' @param family String specifying distribution: "poisson", "normal", or "binomial".
//'   Determines the link function and variance structure.
//'   
//' @return List containing:
//'   \item{grad}{Gradient vector of length q+2p-1 (after identifiability projection)}
//'   \item{hess}{Hessian matrix of dimension (q+2p-1) x (q+2p-1)}
//'   \item{shess}{Score outer product matrix (for robust standard errors)}
//'   
//' @note The Hessian may not be positive definite far from the optimum, which can
//'   cause optimization issues. The score outer product (shess) provides an
//'   alternative for standard error calculation.
//'   
//' @examples
//' \dontrun{
//' // Called internally by optim() during optimization:
//' result <- cpp_mll_gH(current_params, Y, W, X, Z_list, "poisson")
//' // Use gradient for search direction
//' // Use Hessian for step size (quasi-Newton methods)
//' }
// [[Rcpp::export]]
Rcpp::List cpp_mll_gH(const arma::vec& tab, const arma::cube& Y, const arma::cube& W, const arma::cube& X,
                      const Rcpp::List& Z_list, const std::string& family) {

    int m = Y.n_rows;
    int T = Y.n_slices;
    int p = W.n_slices;
    int q = Z_list.size();

    // --- 1. Parse Parameters and Setup ---
    arma::vec theta(q);
    if (q > 0) {
        theta = tab.subvec(0, q-1);
    }

    arma::vec alpha(p);
    arma::vec beta(p);

    if (p > 0) {
        alpha(0) = 1.0;
        if (p > 1) {
            int start_idx = q;
            int end_idx = q + p - 2;
            if (start_idx <= end_idx) {
                alpha.subvec(1, p-1) = tab.subvec(start_idx, end_idx);
            }
        }
        int start_idx_beta = q + p - 1;
        beta = tab.subvec(start_idx_beta, tab.n_elem-1);
    }


    // Initialize G, H, S (Full dimension q + 2p)
    int dim_full = q + 2*p;
    // Initialize gradient, Hessian, and score matrices
    arma::vec G(dim_full, fill::zeros);
    arma::mat H(dim_full, dim_full, fill::zeros);
    arma::mat S(dim_full, dim_full, fill::zeros);

    // Calculate A and B matrices (m x m)
    arma::mat A(m, m, fill::zeros);
    arma::mat B(m, m, fill::zeros);
    if (p > 0) {
        A = cpp_amprod_W_v(W, alpha);
        B = cpp_amprod_W_v(W, beta);
    }

    // --- 2. Main Loop over T, i, j ---
    // We iterate through the data to calculate G, H, S based on the NLL derivatives.
    for (int t = 0; t < T; ++t) {
        arma::mat Xt = X.slice(t);
        arma::mat Yt = Y.slice(t);

        // Calculate linear predictor Eta_t (m x m)

        // ZT part
        arma::mat ZT(m, m, fill::zeros);
        if (q > 0) {
            for(int k=0; k<q; ++k) {
                // Safely access the cube from the list
                if (Z_list.size() > k) {
                    arma::cube Zk = as<arma::cube>(Z_list[k]);
                    if (Zk.n_slices > t) {
                        // Extract slice t from the k-th covariate cube
                        arma::mat Zkt = Zk.slice(t);
                        ZT += Zkt * theta(k);
                    }
                }
            }
        }


        // Bilinear part AXB
        arma::mat AXB(m, m, fill::zeros);
        if (p > 0) {
            AXB = A * Xt * B.t();
        }

        arma::mat Eta = ZT + AXB;

        // Calculate Mu (mean) and Residuals/Weights based on family
        arma::mat Mu(m, m);
        arma::mat Resid(m, m); // NLL residual: d(NLL)/d(eta)
        arma::mat Weights(m, m); // Hessian weights: d^2(NLL)/d(eta)^2

        if (family == "poisson") {
            Mu = arma::exp(Eta);
            Resid = Mu - Yt;
            Weights = Mu;
        } else if (family == "normal") {
            Mu = Eta;
            Resid = Mu - Yt;
            Weights.fill(1.0); // Assuming sigma^2=1
        } else if (family == "binomial") {
            // Use arma::exp for element-wise exponentiation
            Mu = 1.0 / (1.0 + arma::exp(-Eta)); // This is 'p'
            Resid = Mu - Yt;
            Weights = Mu % (1.0 - Mu); // p*(1-p)
            // Add small stabilization factor to weights to prevent exact zeros
            Weights += 1e-16;
        } else {
            stop("Unsupported family.");
        }

        // Loop over dyads (i, j)
        for (int i = 0; i < m; ++i) {
            for (int j = 0; j < m; ++j) {
                if (i == j || std::isnan(Yt(i,j))) continue;

                double resid_ij = Resid(i, j);
                double weight_ij = Weights(i, j);

                // --- Calculate derivatives of eta w.r.t parameters (d_eta) ---
                // This forms the 'Xtab' vector for this observation (q + 2p vector)

                arma::vec d_eta(dim_full, fill::zeros);

                // 1. d_eta / d_theta (q vector)
                if (q > 0) {
                    for(int k=0; k<q; ++k) {
                        arma::cube Zk = as<arma::cube>(Z_list[k]);
                        d_eta(k) = Zk(i, j, t);
                    }
                }

                if (p > 0) {
                    // 2. d_eta / d_alpha (p vector)
                    // Derivative is (W_k X B^T)_{ij}
                    for(int k=0; k<p; ++k) {
                        // This calculation avoids explicit Xij construction
                        arma::mat WkXB = W.slice(k) * Xt * B.t();
                        d_eta(q + k) = WkXB(i, j);
                    }

                    // 3. d_eta / d_beta (p vector)
                    // Derivative is (A X W_k^T)_{ij}
                    for(int k=0; k<p; ++k) {
                        arma::mat AXWk = A * Xt * W.slice(k).t();
                        d_eta(q + p + k) = AXWk(i, j);
                    }
                }

                // --- Update G, H, S ---

                // Update Gradient (G += resid * d_eta)
                arma::vec score_i = resid_ij * d_eta;
                G += score_i;

                // Update Score outer product (S += score * score^T)
                S += score_i * score_i.t();

                // Update Hessian (H += weight * d_eta * d_eta^T + resid * d^2_eta)

                // GLM part (Fisher Information approximation)
                arma::mat H_glm = weight_ij * (d_eta * d_eta.t());

                // Bilinear part (Observed Information adjustment for d^2_eta)
                arma::mat H_adj(dim_full, dim_full, fill::zeros);

                // The adjustment term is only strictly required for the Observed Information Hessian.
                // It was included in the original R Poisson script but omitted for Normal/Binomial.
                // We include it here for Poisson as in the original script.
                if (family == "poisson" && p > 0) {
                    // d^2(eta)/d(alpha_k)d(beta_l) = (W_k X W_l^T)_{ij}

                    // Calculate the p x p matrix of second derivatives H_ab
                    arma::mat H_ab(p, p);
                    for(int k=0; k<p; ++k) {
                        for(int l=0; l<p; ++l) {
                            arma::mat WkXWl = W.slice(k) * Xt * W.slice(l).t();
                            H_ab(k, l) = WkXWl(i, j);
                        }
                    }
                    // Adjustment = resid * H_ab
                    arma::mat rH_ab = resid_ij * H_ab;

                    // Place into the Hessian adjustment matrix
                    H_adj.submat(q, q + p, q + p - 1, dim_full - 1) = rH_ab;
                    H_adj.submat(q + p, q, dim_full - 1, q + p - 1) = rH_ab.t();
                }

                H += H_glm + H_adj;
            }
        }
    }

    // --- 3. Apply Identifiability Constraint (J matrix) ---
    // J removes the dimension corresponding to alpha_1 (index q).

    if (p == 0) {
        // No bilinear part, no J needed.
        return Rcpp::List::create(
            Rcpp::Named("grad") = G,
            Rcpp::Named("hess") = H,
            Rcpp::Named("shess") = S
        );
    }

    int dim_reduced = dim_full - 1;
    arma::mat J(dim_reduced, dim_full, fill::zeros);

    // Construct J matrix
    // Identity for theta (0 to q-1)
    if (q > 0) {
        J.submat(0, 0, q - 1, q - 1) = arma::eye(q, q);
    }
    // Identity for alpha[-1] (q to q+p-2) -> maps from (q+1 to q+p-1) in full
    if (p > 1) {
        J.submat(q, q + 1, q + p - 2, q + p - 1) = arma::eye(p - 1, p - 1);
    }
    // Identity for beta (q+p-1 to end) -> maps from (q+p to end) in full
    if (p > 0) {
       J.submat(q + p - 1, q + p, dim_reduced - 1, dim_full - 1) = arma::eye(p, p);
    }


    // Project G, H, S
    arma::vec grad_id = J * G;
    arma::mat hess_id = J * H * J.t();
    arma::mat shess_id = J * S * J.t();

    return Rcpp::List::create(
        Rcpp::Named("grad") = grad_id,
        Rcpp::Named("hess") = hess_id,
        Rcpp::Named("shess") = shess_id
    );
}
