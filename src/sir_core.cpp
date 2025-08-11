#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

using namespace Rcpp;
using namespace arma;

// ----------------------------------------------------------------------------
// --- Tensor Utility Functions ---
// ----------------------------------------------------------------------------

// Efficient calculation of tprod(X, list(A, B)) specialized for SIR model.
// Calculates A %*% X_t %*% B^T for each slice t of X.
// X is (m x m x T), A and B are (m x m). Result is (m x m x T).
// [[Rcpp::export]]
arma::cube cpp_tprod_A_X_Bt(const arma::cube& X, const arma::mat& A, const arma::mat& B) {
  int m = X.n_rows;
  int T = X.n_slices;

  arma::cube AXB(m, m, T);
  arma::mat Bt = B.t();

  // Parallelize the loop over time T
  #pragma omp parallel for schedule(static)
  for(int t=0; t < T; ++t) {
    AXB.slice(t) = A * X.slice(t) * Bt;
  }
  return AXB;
}

// Efficient calculation of amprod(W, v, 3) specialized for SIR model.
// Calculates the linear combination of slices of W using vector v.
// W is (m x m x p), v is (p x 1). Result is (m x m).
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

//' Construct Wbeta design matrix for ALS update (theta, alpha)
//'
//' Efficiently calculates and flattens the Wbeta component of the design matrix.
//' Wbeta[,,k,t] corresponds to W[,,k] %*% X[,,t] %*% (W*beta)^T.
//'
//' @param W (m x m x p) array of influence covariates.
//' @param X (m x m x T) array, typically lagged outcomes.
//' @param beta (p x 1) vector of current beta estimates.
//' @return (m*m*T) x p matrix for the GLM design matrix.
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
    // Parallelize loop over p features
    #pragma omp parallel for schedule(static)
    for(int k=0; k < p; ++k) {
        // This uses the optimized tprod implementation
        arma::cube Wk_X_WSbeta = cpp_tprod_A_X_Bt(X, W.slice(k), WSbeta);
        // Flatten the resulting cube into a column vector (Column-major flattening)
        // This ensures the flattening matches R's c(Y) order.
        Wbeta_design.col(k) = arma::vectorise(Wk_X_WSbeta);
    }

    return Wbeta_design;
}

//' Construct Walpha design matrix for ALS update (theta, beta)
//'
//' Efficiently calculates and flattens the Walpha component of the design matrix.
//' Walpha[,,k,t] corresponds to (W*alpha) %*% X[,,t] %*% W[,,k]^T.
//'
//' @param W (m x m x p) array of influence covariates.
//' @param X (m x m x T) array, typically lagged outcomes.
//' @param alpha (p x 1) vector of current alpha estimates.
//' @return (m*m*T) x p matrix for the GLM design matrix.
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
    #pragma omp parallel for schedule(static)
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

// This implementation uses the mathematically robust approach of pre-calculating A and B
// and deriving the gradients based on that structure.
// eta = Z*theta + A*X*B^T, where A=W*alpha, B=W*beta.

//' Calculate Gradient and Hessian for SIR models (NLL)
//'
//' Computes the gradient, Hessian, and score cross-product for the Negative Log-Likelihood (NLL)
//' of the SIR model. This implementation is optimized by avoiding explicit calculation of Xij
//' inside the loops, relying instead on the (A, B) structure.
//'
//' @param tab Parameter vector (theta, alpha[-1], beta).
//' @param Y (m x m x T) array of outcomes.
//' @param W (m x m x p) array of influence covariates.
//' @param X (m x m x T) array for the bilinear part.
//' @param Z_list List of q (m x m x T) arrays for exogenous covariates (passed this way for 4D handling).
//' @param family Distribution family ("poisson", "normal", "binomial").
//' @return A list containing the gradient (grad), Hessian (hess), and score cross-product (shess) of the NLL.
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
    // We use temporary storage for parallel reduction if OpenMP is used for the main loop
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
    // Parallelizing over T is often effective.

    // Note: Parallelizing this loop requires thread-safe updates to G, H, S (e.g., using reduction or private copies).
    // Given the complexity of H and S updates, we keep the loop serial here for correctness,
    // as the major bottlenecks (tensor products) are already optimized/parallelized elsewhere.
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
