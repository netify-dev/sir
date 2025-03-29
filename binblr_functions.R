#### ---- functions for fitting bilinear logistic (binary) regression models
#### ---- thoughts:
#### ---- 1. Maintain the same rank-1 structure used in the Poisson code (alpha_1 = 1 for identifiability).
#### ---- 2. Use an iterative block coordinate descent or a direct `optim`-based approach.
#### ---- 3. For a binary outcome, we assume a Bernoulli likelihood with logistic link:
#### ----    Y_{i,j,t} ~ Bernoulli( p_{i,j,t} ),  where p_{i,j,t} = 1 / (1 + exp(-eta_{i,j,t})).
#### ---- 4. Negative log-likelihood is:  - sum_{i,j,t}[ Y_{i,j,t} * eta - log(1 + exp(eta)) ].


#### ---- gradient and Hessian of minus log likelihood
mll_gH_binom <- function(tab, Y, W, X, Z)
{
  # This function calculates:
  #   1. The gradient of the negative log-likelihood (NLL) under a logistic bilinear model
  #   2. The Hessian of the NLL
  #   3. The cross-product of score vectors (S matrix) for possible robust SEs
  #
  # Model:
  #   Y_{i,j,t} ~ Bernoulli( p_{i,j,t} ),  where p_{i,j,t} = logistic( eta_{i,j,t} )
  #   eta_{i,j,t} = theta^T Z_{i,j,t} + alpha^T (some transform of X_{i,j,t}) beta
  #
  # Identifiability:
  #   We fix alpha_1 = 1, then store alpha[-1], so total dim for tab = q + (p-1) + p.
  #
  # Args:
  #   tab : parameter vector in the order (theta, alpha[-1], beta).
  #   Y   : m x m x T array of {0,1} binary outcomes
  #   W   : m x m x p array of nodal/dyadic influence covariates
  #   X   : an object that encodes how Y_{i,j,t-1} (or other structures) feed into the bilinear term
  #   Z   : m x m x q x T array of exogenous dyadic covariates for theta
  #
  # Returns a list with:
  #   grad : the gradient vector of length (q + 2p - 1)
  #   hess : the Hessian matrix (q + 2p - 1) x (q + 2p - 1)
  #   shess: the cross-product of individual score vectors (for robust/sandwich variance)

  m <- dim(Y)[1]
  p <- dim(W)[3]   # number of influence covariates
  q <- dim(Z)[3]   # number of exogenous dyad-level covariates

  # parse parameter vector
  theta <- tab[1:q]
  alpha <- c(1, tab[(q+1):(q + p - 1)])
  beta  <- tab[-(1:(q + p - 1))]

  # We'll track: gll (gradient), Hll (Hessian), Sll (score outer products)
  gll   <- rep(0, q + 2*p)  # dimension pre-J
  Hll   <- matrix(0, q + 2*p, q + 2*p)
  Sll   <- matrix(0, q + 2*p, q + 2*p)

  for(i in 1:m){
    for(j in setdiff(1:m, i))
    {
      # outcome: Y_{i,j, ] is length T
      y_ij <- Y[i,j, ]

      # covariates for theta: Z_{i,j,,] => T x q
      Zij <- t(Z[i,j,,])   # dimension T x q

      # build bilinear X_{i,j,t}
      Xij <- tprod(X, list(t(W[i,,]), t(W[j,,])))

      # linear predictor: eta_{t} = Z_ij[t, ] %*% theta + alpha^T Xij[t] beta
      # dimension of Xij is (p, p, T)
      # we'll flatten that with alpha, beta
      eta_vec <- Zij %*% theta + c( tprod(Xij, list(matrix(alpha,1,p), matrix(beta,1,p))) )

      # logistic transform:
      p_vec  <- 1 / (1 + exp(-eta_vec))  # length T

      # negative log-likelihood piece for each t is - [y_{t} * eta - log(1 + exp(eta))]
      # partial derivative wrt eta_t: (p_t - y_t)
      resid <- (p_vec - y_ij)  # length T

      # partial derivatives w.r.t (theta, alpha, beta)
      #   * w.r.t theta_k => Z_{ij,k}
      #   * w.r.t alpha   => X_{ij} * beta
      #   * w.r.t beta    => X_{ij}^T alpha
      Xb <- t( amprod(Xij, matrix(beta,1,p), 2 )[ ,1, ] )  # T x p
      Xa <- t( amprod(Xij, matrix(alpha,1,p), 1 )[1,, ] )  # T x p

      # combine into big design matrix: (theta part, alpha part, beta part)
      # dimension T x (q + 2p). We'll remove alpha_1 row with J later.
      Xtab <- cbind(Zij, Xb, Xa)  # T x (q + 2*p)

      # gradient accumulates: sum_t [ resid_t * partial_t ]
      # do row-wise multiplication
      eX  <- sweep(Xtab, 1, resid, "*")
      gll <- gll + colSums(eX)

      # crossprod for robust:
      Sll <- Sll + crossprod(eX)

      # Hessian: sum_t [ p_t*(1 - p_t) * partial_t outer partial_t ]
      # we do that by weighting each row of Xtab by sqrt(w_t), w_t = p_t*(1 - p_t)
      w_t  <- p_vec * (1 - p_vec)
      sqrtwX <- sweep(Xtab, 1, sqrt(w_t), "*")  # each row scaled
      H0 <- crossprod(sqrtwX)
      Hll <- Hll + H0
    }
  }

  # identifiability reparam => alpha_1=1, drop that dimension
  # dimension was q+2p => free param dimension is q+2p-1
  J <- diag(q + 2*p)[-(q+1), ]  # remove row= q+1

  list(
    grad  = - J %*% gll,
    hess  = - J %*% Hll %*% t(J),
    shess =    J %*% Sll %*% t(J)
  )
}


#### ---- bilinear predictor for logistic
eta_tab_binom <- function(tab, W, X, Z)
{
  # Return the linear predictor eta_{i,j,t} = Z_{i,j,t}^T theta + alpha^T X_{i,j,t} beta
  # for logistic. We'll then do p_{i,j,t} = logistic(eta_{i,j,t}) if needed.

  p <- dim(W)[3]
  q <- dim(Z)[3]

  theta <- tab[1:q]
  alpha <- c(1, tab[(q+1):(q+p-1)])
  beta  <- tab[-(1:(q+p-1))]

  # build A, B from alpha, beta
  A <- amprod(W, matrix(alpha,1,p), 3)[,,1]
  B <- amprod(W, matrix(beta,1,p) , 3)[,,1]

  # bilinear part
  AXB <- tprod(X, list(A, B))

  # exogenous part
  ZT  <- amprod(Z, matrix(theta,1,q), 3)[,,1,]

  ZT + AXB
}


#### ---- minus log likelihood
mll_binblr <- function(tab, Y, W, X, Z)
{
  # negative log-likelihood for logistic
  #  NLL = sum_{i,j,t} [ - y_{i,j,t} * eta_{i,j,t} + log(1 + exp(eta_{i,j,t})) ]
  # We'll just build eta and then compute the sum.

  eta_arr <- eta_tab_binom(tab, W, X, Z)
  # logistic NLL
  # y * (-eta) => -y*eta
  # plus log(1 + exp(eta))

  # safe computation for log(1+exp(eta)):
  # but we'll keep it simple. For large eta, might want log1p(exp(...)) but R handles some.
  nll <- 0
  nll <- nll + sum( - Y * eta_arr, na.rm=TRUE )
  nll <- nll + sum( log(1 + exp(eta_arr)), na.rm=TRUE )

  nll
}


#### ---- fit via optim
binblr_optfit <- function(Y, W, X, Z, trace=0, start=NULL)
{
  # Minimizes the logistic NLL via BFGS using the gradient. Similar to Poisson/Gaussian approach
  p <- dim(W)[3]
  q <- dim(Z)[3]

  if(is.null(start))
  {
    # a naive start: do a logistic GLM ignoring bilinear portion
    glmData <- data.frame( Y = c(Y), apply(Z, 3, c) )
    form    <- formula(paste0("Y ~ -1 + ", paste(colnames(glmData)[-1], collapse=" + ")))
    fit0    <- glm(form, data=glmData, family=binomial())
    # get initial theta
    theta   <- coef(fit0)
    # small random for alpha[-1], beta
    start   <- c(theta, rnorm(p-1, sd=1e-3), rnorm(p, sd=1e-3))
  }

  objfun <- function(par){
    mll_binblr(par, Y, W, X, Z)
  }

  gradfun <- function(par){
    gH <- mll_gH_binom(par, Y, W, X, Z)
    as.numeric(gH$grad)
  }

  fit <- optim(par=start, fn=objfun, gr=gradfun, method="BFGS", control=list(trace=trace))
  tab <- fit$par

  list(
    tab=tab,
    value=fit$value,
    convergence=fit$convergence
  )
}


#### ---- fit via block coordinate approach (IWLS-like in sub-GLMs)
binblr_alsfit <- function(Y, W, X, Z, trace=FALSE)
{
  # We replicate the structure from the Poisson code:
  #   1) fix beta, update (theta, alpha) via a logistic GLM
  #   2) fix alpha, update (theta, beta) via logistic GLM
  #   iterate until convergence
  #
  # This approach is valid because for each subproblem, the model is linear in (theta, alpha)
  # or in (theta, beta). We can use binomial() in R's glm() for each sub-step.

  p <- dim(W)[3]
  q <- dim(Z)[3]
  m <- dim(Y)[1]

  # initial step: ignore bilinear portion => do logistic regression for theta
  glmData <- data.frame( Y=c(Y), apply(Z,3,c) )
  colnames(glmData)[2:(1+q)] <- paste0("Z",1:q)
  form <- as.formula(
    paste0("Y ~ -1 + ", paste(colnames(glmData)[-1], collapse=" + "))
  )
  fit0 <- glm(form, data=glmData, family=binomial())
  theta <- coef(fit0)
  # random small alpha, beta
  alpha <- rnorm(p, sd=0.01)
  beta  <- rnorm(p, sd=0.01)

  # track deviance
  dev_old <- Inf
  dev_new <- fit0$deviance

  # iteration
  ALPHA <- alpha
  BETA  <- beta
  THETA <- theta
  DEV   <- matrix(c(dev_old, dev_new), nrow=1)

  repeat
  {
    # 1) fix beta, update (theta, alpha)
    # build design: X * beta is the term. We flatten and do a logistic GLM
    WSbeta <- amprod(W, t(beta), 3)[,,1]
    Wbeta  <- array(dim=c(m,m,p,dim(Y)[3]))
    for(k in 1:p){
      Wbeta[,,k,] <- tprod(X, list(W[,,k], WSbeta))
    }

    glmData <- data.frame( Y=c(Y), apply(Z,3,c), apply(Wbeta, 3, c) )
    colnames(glmData)[2:(1+q)] <- paste0("Z",1:q)
    colnames(glmData)[(2+q):ncol(glmData)] <- paste0("WB",1:(p))
    form1 <- as.formula(
      paste0("Y ~ -1 + ", paste(colnames(glmData)[-1], collapse=" + "))
    )
    fitA <- glm(form1, data=glmData, family=binomial())
    coA  <- coef(fitA)
    theta <- coA[1:q]
    alpha <- coA[-(1:q)]
    names(alpha) <- paste0("alpha",1:p)

    # 2) fix alpha, update (theta, beta)
    WSalpha <- amprod(W, t(alpha), 3)[,,1]
    Walpha  <- array(dim=c(m,m,p,dim(Y)[3]))
    for(k in 1:p){
      Walpha[,,k,] <- tprod(X, list(WSalpha, W[,,k]))
    }

    glmData <- data.frame( Y=c(Y), apply(Z,3,c), apply(Walpha,3,c) )
    colnames(glmData)[2:(1+q)] <- paste0("Z",1:q)
    colnames(glmData)[(2+q):ncol(glmData)] <- paste0("WA",1:p)
    form2 <- as.formula(
      paste0("Y ~ -1 + ", paste(colnames(glmData)[-1], collapse=" + "))
    )
    fitB <- glm(form2, data=glmData, family=binomial())
    coB  <- coef(fitB)
    theta <- coB[1:q]
    beta  <- coB[-(1:q)]
    names(beta) <- paste0("beta",1:p)

    # check deviance
    dev_old <- dev_new
    dev_new <- fitB$deviance

    THETA <- rbind(THETA, theta)
    ALPHA <- rbind(ALPHA, alpha)
    BETA  <- rbind(BETA, beta)
    DEV   <- rbind(DEV, c(dev_old, dev_new))

    if(trace){
      cat("Iteration:", nrow(DEV)-1, "Dev:", dev_new, "\n")
    }
    # stopping
    if( abs(dev_old - dev_new)/abs(dev_old) < 1e-8 ) break
    if(nrow(DEV) > 50) break  # failsafe
  }

  # impose alpha_1=1
  a <- alpha[-1]/alpha[1]
  b <- beta * alpha[1]

  list(
    theta=theta,
    a=a,
    b=b,
    tab=c(theta, a, b),
    ALPHA=ALPHA,
    BETA=BETA,
    THETA=THETA,
    DEV=DEV
  )
}


#### ---- standard errors based on Hessian
se_binblr <- function(tab, Y, W, X, Z, calcSE=TRUE)
{
  # invert Hessian from identified param
  gH <- mll_gH_binom(tab, Y, W, X, Z)
  H  <- gH$hess

  se <- sqrt(diag(solve(H)))
  if(!calcSE) return(se)

  # robust sandwich
  sh <- gH$shess
  rse <- sqrt(diag(solve(H) %*% sh %*% solve(H)))
  list(se=se, rse=rse)
}


binblr <- function(Y, W, X, Z, calcSE=TRUE)
{
  # convenience wrapper:
  #  1) fit via block coord approach
  #  2) build final alpha,beta
  #  3) build A,B
  #  4) optionally compute SE

  mod <- binblr_alsfit(Y, W, X, Z, trace=FALSE)
  tab <- mod$tab

  p <- dim(W)[3]
  q <- dim(Z)[3]

  alpha <- c(1, tab[(q+1):(q+p-1)])
  beta  <- tab[-(1:(q+p-1))]

  # influence matrices
  A <- amprod(W, matrix(alpha,1,p), 3)[,,1]
  A <- A * sign(mean(diag(A)))
  diag(A) <- 0
  rownames(A) <- colnames(A) <- rownames(Y)

  B <- amprod(W, matrix(beta,1,p), 3)[,,1]
  B <- B * sign(mean(diag(B)))
  diag(B) <- 0
  rownames(B) <- colnames(B) <- rownames(Y)

  # logistic log-lik
  #   log-likelihood = sum( y*eta - log(1+exp(eta)) )
  #   ignoring sign => we do negative above
  eta_arr <- eta_tab_binom(tab, W, X, Z)
  ll <- sum( Y * eta_arr - log(1 + exp(eta_arr)), na.rm=TRUE )

  summ <- cbind(coef=tab)
  if(calcSE){
    tmp <- se_binblr(tab, Y, W, X, Z, calcSE=TRUE)
    se_  <- tmp$se
    rse_ <- tmp$rse
    t_se_  <- tab / se_
    t_rse_ <- tab / rse_
    summ <- cbind(coef=tab, se=se_, rse=rse_, t_se=t_se_, t_rse=t_rse_)
  }

  list(
    summ = summ,
    A    = A,
    B    = B,
    ll   = ll,
    ALPHA=mod$ALPHA,
    BETA=mod$BETA,
    THETA=mod$THETA,
    DEV=mod$DEV
  )
}


#### ---- explanation
# 1) This code closely parallels the Poisson or Gaussian bilinear models, but the outcome is binary.
#    We assume logistic link: Y_{ij,t} ~ Bernoulli( p_{ij,t} ), p_{ij,t} = 1/(1+exp(-eta_{ij,t})).
# 2) We maintain alpha_1=1 for identifiability of alpha,beta in the rank-1 "bilinear" part.
# 3) 'mll_gH_binom' computes gradient and Hessian for the negative log-likelihood:
#    NLL = sum( - y*eta + log(1+exp(eta)) ), plus the chain-rule trick to drop alpha_1.
# 4) 'binblr_alsfit' does block coordinate updates:
#    - fix beta, estimate (theta, alpha) via logistic GLM
#    - fix alpha, estimate (theta, beta) via logistic GLM
#    repeating until convergence.
# 5) 'binblr_optfit' uses a standard BFGS approach via 'optim'.
# 6) 'binblr' is a wrapper that defaults to the ALS routine, recovers alpha,beta, returns influence matrices, 
#    the final log-likelihood, and (optionally) standard errors (classic and robust).
# 7) The result is a bilinear logistic model for binary network edges with interpretability akin 
#    to the Poisson or Gaussian analogues.
