#### ---- functions for fitting bilinear gaussian regression models
#### ---- thoughts:
#### ---- 1. maintain the same rank-1 structure used in the Poisson code (alpha_1=1 for identifiability)
#### ---- 2. use an iterative block coordinate descent or a direct optim-based approach
#### ---- 3. for a Gaussian outcome, the negative log-likelihood is 0.5 * sum( (Y - mu)^2 ), ignoring additive constants


#### ---- gradient and Hessian of minus log likelihood
mll_gH_gauss <- function(tab, Y, W, X, Z)
{
  # This function calculates:
  #   1. The gradient of the negative log-likelihood (NLL) under a Gaussian bilinear model
  #   2. The Hessian of the NLL
  #   3. The cross-product of score vectors (S matrix) for possible robust SEs (optional)
  #
  # Model:
  #   Y_{i,j,t} ~ Normal( mu_{i,j,t}, sigma^2=1 ),
  #   mu_{i,j,t} = theta^T Z_{i,j,t} + alpha^T (some transformation of X_{i,j,t}) beta
  # Identifiability:
  #   We fix alpha_1 = 1, then store the rest in 'tab'. 
  #
  # Args:
  #   tab : parameter vector in the order (theta, alpha[-1], beta).
  #   Y   : m x m x T array of observed Gaussian responses
  #   W   : m x m x p array of nodal/dyadic influence covariates (same as Poisson code)
  #   X   : an object that encodes how Y_{i,j,t-1} or other covariates feed into the bilinear term
  #         typically an m x m x T array, then combined with W to form X_{ij,t} for alpha,beta
  #   Z   : m x m x q x T array of exogenous dyadic covariates (for theta)
  #
  # Returns a list with:
  #   grad : the gradient vector of dimension (q + 2p - 1)
  #   hess : the Hessian matrix of dimension (q + 2p - 1) x (q + 2p - 1)
  #   shess: the cross-product of individual score contributions (for robust/sandwich variance)
  #
  # The code parallels the Poisson version but uses the Gaussian negative log-likelihood.

  m <- dim(Y)[1]
  p <- dim(W)[3]   # number of influence-related covariates
  q <- dim(Z)[3]   # number of exogenous dyad-level covariates

  # parse parameter vector
  # tab is of length: q + (p-1) + p = q + 2p - 1
  theta <- tab[1:q]
  alpha <- c(1, tab[(q+1):(q + p - 1)])
  beta  <- tab[-(1:(q + p - 1))]  # last p entries

  # initialize gradient, Hessian, cross-product (for robust sandwich)
  gll   <- rep(0, q + 2*p)  # before imposing alpha_1=1 identifiability
  Hll   <- matrix(0, q + 2*p, q + 2*p)
  Sll   <- matrix(0, q + 2*p, q + 2*p)

  # loop over all i, j, ignoring diagonal (i != j)
  for(i in 1:m){
    for(j in setdiff(1:m, i))
    {
      # Y_{i,j,] : length T
      y_ij <- Y[i,j, ]

      # Z_{i,j,,]: T x q (transposed for matrix ops)
      Zij  <- t(Z[i,j,,])  # dimension: T x q

      # build the bilinear contribution for each time t:
      # analogous to Poisson code, we get X_{ij,t} from (i', j') sums
      # or from pre-multiplied arrays. We'll re-use the same "tprod" approach.
      Xij  <- tprod(X, list(t(W[i,,]), t(W[j,,])))

      # linear predictor mu_{i,j,t}:
      #   mu = Zij %*% theta + alpha^T Xij beta  (elementwise per time)
      #   dimension of Xij is typically (p, p, T)
      # same as Poisson code, but no exp(). Just identity link.
      linpred <- Zij %*% theta + c( tprod(Xij, list(matrix(alpha,1,p), matrix(beta,1,p))) )
      mu      <- linpred  # length T

      # residual r = mu - y
      r_ij    <- mu - y_ij

      # partial derivatives of mu wrt parameters
      #   * w.r.t. theta -> Z_ij
      #   * w.r.t. alpha -> X_{ij} beta
      #   * w.r.t. beta  -> X_{ij}^T alpha
      # We'll form a T x (q + p + p) matrix, but recall alpha_1=1 is "internal"
      # so effectively we have q + (p-1) + p free parameters in tab.

      # Xb: T x p, each row is X_{ij,t} * beta
      Xb <- t( amprod(Xij, matrix(beta,1,p),2 )[ ,1, ] )  # dimension T x p
      # Xa: T x p, each row is X_{ij,t}^T alpha
      Xa <- t( amprod(Xij, matrix(alpha,1,p),1 )[1,, ] )  # dimension T x p

      # design block: [Zij, Xb[,-1], Xa]
      # but we keep dimension q + p + p for the gradient. We'll fill carefully:
      #   index 1..q -> Z_ij
      #   index (q+1)..(q+p-1) -> Xb[,-1] for alpha[-1]
      #   index (q+p)..(q+2p-1) -> beta part
      # The code below just builds a T x (q+2p) then we'll handle alpha[1]=1 with J.
      # For convenience, we build a full T x (q + 2p) matrix, ignoring the eventual drop of alpha_1.
      # We'll store alpha_1 in the alpha slot but it won't appear in tab. We'll remove it with J.

      # For now, define a big matrix Xtab:
      Xtab <- cbind(Zij, Xb, Xa)  # T x (q + 2p)

      # contribution to gradient
      # derivative w.r.t. param = sum_t( r_ij[t] * partial(t) )
      eX  <- sweep(Xtab, 1, r_ij, "*")  # T x (q + 2p)
      gll <- gll + colSums(eX, na.rm=TRUE)

      # cross-product for robust or sandwich SE
      Sll <- Sll + crossprod(eX)

      # Hessian: sum over t of partial(t) partial(t)^T
      # but partial(t) is exactly the row in Xtab. We multiply each row's outer product,
      # since d^2 mu/d param^2 = 0 for a linear predictor.
      # So for each t, H += partial_t outer partial_t
      # We'll do it in a vectorized manner:
      H0 <- crossprod(Xtab)  # sum_{t} partial_t outer partial_t, ignoring r_ij
      # but for Gaussian, the negative log-lik is 0.5 * (r^2). The factor of (mu - y) in derivative
      # means effectively we get partial^2 = sum_{t} partial_t^2. 
      # Actually we need to be cautious: The derivative is (mu - y). The second derivative is
      # sum_t partial_t * partial_t (like OLS). Indeed we do:
      #   grad param_k = sum_t r_t * partial_tk
      #   Hess_{k,l} = sum_t partial_tk * partial_tl
      # so that crossprod(Xtab) is exactly what's needed. We'll add it to the total H.
      Hll <- Hll + H0
    }
  }

  # Now we incorporate the identifiability reparam:
  #   alpha_1 is fixed to 1. We do the same approach as the Poisson code, using J to drop that dimension.
  # Dimension was q + 2p. The free parameters are q + (p-1) + p = q + 2p - 1.
  # J is the same identity approach used in Poisson code:
  J <- diag(q + 2*p)[-(q+1), ]  # remove row q+1 (the alpha_1 slot)

  # return identified gradient / Hessian / sandwich
  list(
    grad  = - J %*% gll,
    hess  = - J %*% Hll %*% t(J),
    shess =    J %*% Sll %*% t(J)  # for robust sandwich if desired
  )
}


#### ---- bilinear predictor
eta_tab_gauss <- function(tab, W, X, Z)
{
  # Construct the linear predictor mu_{i,j,t} = Z_{i,j,t}^T theta + alpha^T X_{i,j,t} beta
  #
  # Args:
  #   tab:  (theta, alpha[-1], beta) with alpha_1=1 implicitly
  #   W:    m x m x p array of influence covariates
  #   X:    typically an m x m x T array capturing how lagged Y or other structure
  #         interacts with W to produce the bilinear part
  #   Z:    m x m x q x T array for exogenous dyadic covariates
  #
  # Return:
  #   A 3D array of dimension (m, m, T) containing mu_{i,j,t}.

  p <- dim(W)[3]
  q <- dim(Z)[3]

  theta <- tab[1:q]
  alpha <- c(1, tab[(q+1):(q+p-1)])
  beta  <- tab[-(1:(q+p-1))]

  # A = amprod(W, alpha, 3): m x m from alpha
  A <- amprod(W, matrix(alpha, 1, p), 3)[,,1]
  # B = amprod(W, beta, 3): m x m from beta
  B <- amprod(W, matrix(beta, 1, p), 3)[,,1]

  # Bilinear part: tprod(X, list(A,B)) => dimension (m, m, T)
  AXB <- tprod(X, list(A, B))

  # linear part from Z: ZT = amprod(Z, theta, 3) => dimension (m, m, 1, T)
  ZT  <- amprod(Z, matrix(theta, 1, q), 3)[,,1,]  # shape (m, m, T)

  # sum them
  ZT + AXB
}


#### ---- minus log likelihood
mll_gaussblr <- function(tab, Y, W, X, Z)
{
  # For Gaussian outcomes with sigma^2=1, ignoring additive constants:
  # NLL = 0.5 * sum( (Y - mu)^2 )
  #
  mu <- eta_tab_gauss(tab, W, X, Z)
  # negative log-likelihood
  0.5 * sum( (Y - mu)^2, na.rm=TRUE )
}


#### ---- fit via optim
gaussblr_optfit <- function(Y, W, X, Z, trace=0, tab=NULL)
{
  # Minimizes the Gaussian NLL = 0.5 * sum( (Y-mu)^2 ) w.r.t. tab
  # using BFGS. Similar to poisblr_optfit but no link function needed.

  p <- dim(W)[3]
  q <- dim(Z)[3]
  m <- nrow(Y)

  # if no starting tab provided, we do a simple OLS-like guess ignoring bilinear structure
  if(is.null(tab)) {
    # Flatten Y and Z for a quick OLS ignoring alpha/beta
    glmData <- data.frame( Y = c(Y), apply(Z, 3, c) )
    form    <- formula(paste0("Y ~ -1 + ", paste(colnames(glmData)[-1], collapse=" + ")))
    fit0    <- lm(form, data=glmData)
    # initial guess for theta:
    theta   <- coef(fit0)
    # small random for alpha[-1], beta
    tab     <- c( theta, rnorm(p-1, sd=1e-3), rnorm(p, sd=1e-3) )
  }

  # define objective
  objfun <- function(par) {
    mll_gaussblr(par, Y=Y, W=W, X=X, Z=Z)
  }

  # define gradient
  gradfun <- function(par) {
    gH <- mll_gH_gauss(par, Y=Y, W=W, X=X, Z=Z)
    as.numeric( gH$grad )
  }

  # run optim
  fit <- optim(tab, objfun, gradfun, method="BFGS", control=list(trace=trace))
  tab <- fit$par

  # parse final
  theta <- tab[1:q]
  a     <- tab[(q+1):(q+p-1)]
  b     <- tab[-(1:(q+p-1))]

  list(theta=theta, a=a, b=b, tab=tab, value=fit$value, conv=fit$convergence)
}


#### ---- fit via alternating reweighted least squares (coordinate descent)
gaussblr_alsfit <- function(Y, W, X, Z, trace=FALSE)
{
  # This parallels poisblr_alsfit but for Gaussian data.
  # We treat the updates for (theta, alpha) and (theta, beta)
  # each as standard linear regressions.
  # Return final parameter estimates plus iteration logs.

  p <- dim(W)[3]
  q <- dim(Z)[3]
  m <- nrow(Y)

  # Start by ignoring the bilinear portion to get a quick OLS for theta
  glmData <- data.frame( Y=c(Y), apply(Z, 3, c) )
  names(glmData)[2:(q+1)] <- paste0("Z", 1:q)
  form <- formula(paste0("Y ~ -1 + ", paste(names(glmData)[-1], collapse=" + ")))
  fit0 <- lm(form, data=glmData)
  theta <- coef(fit0)
  # random small starts for alpha, beta
  set.seed(1)
  alpha <- rnorm(p)/nrow(Y)
  beta  <- rnorm(p)/nrow(Y)

  # track iteration
  DEV <- matrix( c(Inf, sum((fit0$residuals)^2)), 1, 2,
                 dimnames=list(NULL, c("oldDev","newDev")) )
  THETA <- theta
  ALPHA <- alpha
  BETA  <- beta

  # iterative updates
  while( abs(DEV[nrow(DEV),1] - DEV[nrow(DEV),2]) / abs(DEV[nrow(DEV),2]) > 1e-9 ) {

    ## -- update theta, alpha
    # Build design: we fix beta, so X_{ij,t} beta -> say Wbeta. Then do a linear model
    WSbeta <- amprod(W, t(beta), 3)[,,1]  # an m x m matrix
    Wbeta  <- array(dim=c(m,m,p,dim(Y)[3]))
    for(k in 1:p){
      Wbeta[,,k,] <- tprod(X, list(W[,,k], WSbeta))
    }

    # Flatten
    glmData <- data.frame( Y=c(Y), apply(Z, 3, c), apply(Wbeta, 3, c) )
    names(glmData)[2:(1+q)] <- paste0("Z",1:q)
    wbn <- names(glmData)[(2+q):ncol(glmData)]
    # linear model: Y ~ -1 + Z + (X * beta)
    form1 <- formula(paste0("Y ~ -1 + ", paste(c(paste0("Z",1:q), wbn), collapse=" + ")))
    fitA  <- lm(form1, data=glmData)

    # from coefficients, extract theta, alpha
    coA   <- coef(fitA)
    theta <- coA[1:q]
    alpha <- coA[-(1:q)]
    names(alpha) <- dimnames(W)[[3]]

    ## -- update theta, beta
    WSalpha <- amprod(W, t(alpha), 3)[,,1]
    Walpha  <- array(dim=c(m,m,p,dim(Y)[3]))
    for(k in 1:p){
      Walpha[,,k,] <- tprod(X, list(WSalpha, W[,,k]))
    }

    glmData <- data.frame( Y=c(Y), apply(Z, 3, c), apply(Walpha,3,c) )
    names(glmData)[2:(1+q)] <- paste0("Z",1:q)
    wan <- names(glmData)[(2+q):ncol(glmData)]
    form2 <- formula(paste0("Y ~ -1 + ", paste(c(paste0("Z",1:q), wan), collapse=" + ")))
    fitB  <- lm(form2, data=glmData)

    coB   <- coef(fitB)
    theta <- coB[1:q]
    beta  <- coB[-(1:q)]
    names(beta) <- dimnames(W)[[3]]

    # track deviance
    DEV <- rbind(DEV, c(DEV[nrow(DEV),2], sum(fitB$residuals^2)))
    THETA <- rbind(THETA, theta)
    ALPHA <- rbind(ALPHA, alpha)
    BETA  <- rbind(BETA, beta)
    if(trace) cat("Iteration:", nrow(DEV)-1, "Dev:", DEV[nrow(DEV),2], "\n")
  }

  # impose alpha_1=1, rescale
  a <- alpha[-1]/alpha[1]
  b <- beta * alpha[1]

  # return
  list(theta=theta, a=a, b=b,
       tab=c(theta, a, b),
       ALPHA=ALPHA, BETA=BETA, THETA=THETA,
       DEV=DEV)
}


#### ---- standard errors based on Hessian
se_gaussblr <- function(tab, Y, W, X, Z, calcSE=TRUE)
{
  # Invert the Hessian from the identified parameterization,
  # returning sqrt of the diagonal for classical SE.
  # For robust/sandwich, we also have an option as in Poisson code.

  gH <- mll_gH_gauss(tab, Y, W, X, Z)
  H  <- gH$hess
  se <- sqrt(diag(solve(H)))

  if(!calcSE) return(se)

  # robust/sandwich
  sh <- gH$shess
  rse <- sqrt(diag( solve(H) %*% sh %*% solve(H) ))

  list(se=se, rse=rse)
}


gaussblr <- function(Y, W, X, Z, calcSE=TRUE)
{
  # Convenient wrapper: 
  # 1) fit with ALS
  # 2) compute predicted mu
  # 3) get log-likelihood
  # 4) optionally get SEs
  #
  # For a direct 'optim' approach, see gaussblr_optfit() similarly.

  mod <- gaussblr_alsfit(Y, W, X, Z, trace=FALSE)
  tab <- mod$tab

  p <- dim(W)[3]
  q <- dim(Z)[3]

  # Construct final alpha, beta
  alpha <- c(1, tab[(q+1):(q+p-1)])
  beta  <- tab[-(1:(q+p-1))]

  # Influence matrices A, B
  A <- amprod(W, matrix(alpha,1,p), 3)[,,1]
  A <- A * sign(mean(diag(A)))
  diag(A) <- 0
  rownames(A) <- colnames(A) <- rownames(Y)

  B <- amprod(W, matrix(beta,1,p), 3)[,,1]
  B <- B * sign(mean(diag(B)))
  diag(B) <- 0
  rownames(B) <- colnames(B) <- rownames(Y)

  # predicted mu
  mu  <- eta_tab_gauss(tab, W, X, Z)
  ll  <- -0.5 * sum( (Y - mu)^2, na.rm=TRUE ) # ignoring constant log(2*pi)

  # optional SE
  if(calcSE){
    gh <- mll_gH_gauss(tab, Y, W, X, Z)
    H  <- gh$hess
    se <- sqrt(diag(solve(H)))
    rse <- sqrt(diag( solve(H) %*% gh$shess %*% solve(H) ))

    summ <- cbind(coef=tab, se=se, rse=rse, t_se=tab/se, t_rse=tab/rse)
  } else {
    summ <- cbind(coef=tab)
  }

  list(
    summ = summ,
    A    = A,
    B    = B,
    pred = mu,
    ll   = ll,
    ALPHA= mod$ALPHA,
    BETA = mod$BETA,
    THETA= mod$THETA,
    DEV  = mod$DEV
  )
}


#### ---- explanation
# 1) This code is the Gaussian analogue of the bilinear Poisson model shown earlier.
# 2) We keep alpha_1=1 for identifiability and represent alpha in "tab" as alpha[-1].
# 3) The negative log-likelihood (NLL) for Gaussian (sigma^2=1) is 0.5 * sum( (Y - mu)^2 ).
# 4) The gradient is sum( (mu - Y) * dmu/dparam ), Hessian is sum( dmu/dparam outer dmu/dparam ).
# 5) The function 'mll_gH_gauss' parallels the Poisson 'mll_gH' but uses a Gaussian NLL.
# 6) We provide both an ALS approach ('gaussblr_alsfit') and an optim-based approach
#    ('gaussblr_optfit'). The final wrapper 'gaussblr' defaults to ALS, returns standard
#    or robust SEs, the fitted influences (A,B), and fitted mu.
# 7) This allows direct interpretability of sender/receiver influences alpha and beta, 
#    as in the Poisson case, but for continuous outcomes.
