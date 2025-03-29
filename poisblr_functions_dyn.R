#### ---- functions for fitting bilinear poisson regression models
#### ---- thoughts:
#### ---- 1. use rank 1 approx to unconstrained fit for starting value


#### ---- 



#### ---- gradient and Hessian of minus log likelihood
mll_gH_dyn<-function(tab,Y=Y,W=W,X=X,Z=Z)
{
  m<-dim(Y)[1]
  p<-dim(W)[3]
  q<-dim(Z)[3]

  theta<-tab[1:q]
  alpha<-c(1,tab[(q+1):(q+p-1)] )
  beta<-tab[-(1:(q+p-1))]

  gll<-rep(0,q+2*p)
  Hll<-Sll<-matrix(0,q+2*p,q+2*p)


  for(i in 1:m){ for(j in (1:m)[-i])
  {
    yij<-Y[i,j,]
    Zij<-t(Z[i,j,,])

    w_i_L = lapply(1:dim(W)[4], function(t){ t(W[i,,,t]) })
    w_j_L = lapply(1:dim(W)[4], function(t){ t(W[j,,,t]) })
    Xij<-tprod_dyn(X, list(w_i_L, w_j_L ))

    eta<-Zij%*%theta +  c( tprod(Xij,list(matrix(alpha,1,p),matrix(beta,1,p))))
    mu<-exp(eta)

    Xb<-t( amprod(Xij, matrix(beta,1,p),2)[,1,]  )
    Xa<-t( amprod(Xij, matrix(alpha,1,p),1)[1,,]  )
    Xtab<-cbind(Zij,Xb,Xa) 
   
    eX<-sweep( Xtab, 1, yij-mu,"*") 
    gll<- gll+ apply(eX,2,sum)

    Sll<-Sll+crossprod(eX) 

    rXij<-apply( sweep(Xij,3,(yij-mu),"*")  ,c(1,2),sum)
    H<- -t(Xtab) %*% sweep(Xtab,1,mu,"*")
    H[q+1:p,q+p+1:p]<- H[q+1:p,q+p+1:p] +rXij
    H[q+p+1:p,q+1:p]<- H[q+p+1:p,q+1:p] + t(rXij)
    Hll<-Hll+H
  }}

  J<-diag(q+2*p)[-(q+1),]
  list( grad= -J%*%gll , hess= -J%*%Hll%*%t(J), shess= J%*%Sll%*%t(J) )
}



#### ---- bilinear predictor
eta_tab_dyn<-function(tab,W,X,Z) 
{
  p<-dim(W)[3] 
  q<-dim(Z)[3] 

  theta<-tab[1:q] 
  alpha<-c(1,tab[(q+1):(q+p-1)] ) 
  beta<-tab[-(1:(q+p-1))]  

  A<-amprod(W, matrix(alpha,1,p) ,3)[,,1,] 
  A<-lapply(1:dim(A)[3], function(t){ A[,,t] })
  B<-amprod(W, matrix(beta,1,p) ,3)[,,1,] 
  B<-lapply(1:dim(B)[3], function(t){ B[,,t] })
  AXB<-tprod_dyn(X,list(A,B))
  ZT<-amprod(Z, matrix(theta,1,q) ,3)[,,1,]

  ZT+AXB
}


#### ---- minus log likelihood
mll_poisblr_dyn<-function(tab,Y,W,X,Z)
{
  ETA<-eta_tab_dyn(tab,W,X,Z)  
 -sum(Y*ETA - exp(ETA),na.rm=TRUE ) 
}


#### ---- fit via alternating reweighted least squares
poisblr_alsfit_dyn<-function(Y,W,X,Z,trace=FALSE)
{
  p<-dim(W)[3] 
  q<-dim(Z)[3] 
  n<-dim(Y)[3]
  m<-nrow(Y)

  fit<-glm( c(Y) ~ -1 + apply(Z,3,c) ,family="poisson" ) 
  theta<-fit$coef

  set.seed(1)

  THETA<-theta
  ALPHA<-alpha<-rnorm(p)/n
  BETA<-beta<-rnorm(p)/n
  DEV<-matrix( c(Inf,deviance(fit)) ,1,2)


  #### ---- block coordinate descent 
  while( abs(DEV[nrow(DEV),1]-DEV[nrow(DEV),2])/abs(DEV[nrow(DEV),2]) > 1e-9)
  {

  ## -- update theta, alpha  

  # - construct design matrix
  WSbeta<-amprod(W,t(beta),3)[,,1,]
  Wbeta<-array(dim=c(m,m,p,n))
  for(t in 1:n){
    for(k in 1:p){
      Wbeta[,,k,t] <- tprod(X[,,t],list(W[,,k,t],WSbeta[,,t]))
    } 
  }
  dimnames(Wbeta)[[3]]<-dimnames(W)[[3]]

  # - conditional fit with glm
  tastart<-c(theta,alpha)
  glmData = data.frame( cbind(Y=c(Y), apply(Z,3,c), apply(Wbeta,3,c) ) )
  names(glmData)[2:(q+1)] = paste0(names(glmData)[2:(q+1)], '_Z')
  names(glmData)[(q+2):ncol(glmData)] = paste0(names(glmData)[(q+2):ncol(glmData)], '_Wbeta')
  names(glmData) = gsub('.1', '', names(glmData))
  glmForm = formula(paste0('Y ~ -1 + ', paste(names(glmData)[-1], collapse=' + ') ))
  
  fit_alpha<-speedglm(formula=glmForm, data=glmData, family=poisson(),
    start=tastart,control=list(trace=trace) )
  
  theta<-fit_alpha$coef[  (1:dim(Z)[3]) ]
  names(theta) = paste0('(Z) ', dimnames(Z)[[3]])
  alpha<-fit_alpha$coef[ -(1:dim(Z)[3]) ]  
  names(alpha) = paste0('(alphaW) ', dimnames(W)[[3]])


  ## -- update theta, beta

  # - construct design matrix
  WSalpha<-amprod(W,t(alpha),3)[,,1,]
  alphaW<-array(dim=c(m,m,p,n))
  for(t in 1:n){
    for(k in 1:p){
       alphaW[,,k,t] <- tprod(X[,,t],list(WSalpha[,,t],W[,,k,t])) 
    } 
  }
  dimnames(alphaW)[[3]]<-dimnames(W)[[3]]

  # - conditional fit with glm
  tbstart<-c(theta,beta)
  glmData = data.frame( cbind(Y=c(Y), apply(Z,3,c), apply(alphaW,3,c) ) )
  names(glmData)[2:(q+1)] = paste0(names(glmData)[2:(q+1)], '_Z')
  names(glmData)[(q+2):ncol(glmData)] = paste0(names(glmData)[(q+2):ncol(glmData)], '_alphaW')
  names(glmData) = gsub('.1', '', names(glmData))
  glmForm = formula(paste0('Y ~ -1 + ', paste(names(glmData)[-1], collapse=' + ') ))
  
  fit_beta<-speedglm(formula=glmForm, data=glmData, family=poisson(),
    start=tbstart,control=list(trace=trace) )
  
  theta<-fit_beta$coef[  (1:dim(Z)[3]) ]
  names(theta) = paste0('(Z) ', dimnames(Z)[[3]])
  beta<-fit_beta$coef[ -(1:dim(Z)[3]) ]  
  names(beta) = paste0('(betaW) ', dimnames(W)[[3]])

  ## -- save results
  ALPHA<-rbind(ALPHA,alpha)
  BETA<-rbind(BETA,beta)
  THETA<-rbind(THETA,theta)
  DEV<-rbind(DEV,c(deviance(fit_alpha),deviance(fit_beta)))

  #cat(date(),DEV[nrow(DEV),],"\n") 
  a<-alpha[-1]/alpha[1] 
  b<-beta*alpha[1] 
  }
  a<-alpha[-1]/alpha[1] 
  b<-beta*alpha[1]  
  list(theta=theta,a=a,b=b,tab=c(theta,a,b),
    ALPHA=ALPHA, BETA=BETA, THETA=THETA, DEV=DEV) 
}

# combine fns
poisblr_dyn <- function(Y, W, X, Z, calcSE=TRUE){
  
  # get param estimates
  mod = poisblr_alsfit_dyn(Y,W,X,Z,trace=FALSE)

  # org mult effs
  tab = mod$tab ; p<-dim(W)[3] ; q<-dim(Z)[3]
  
  # sender influence
  alpha<-c(1,tab[(q+1):(q+p-1)] )
  A<-amprod(W,t(alpha),3 )[,,1,]
  for(t in 1:dim(A)[3]){
    A[,,t]<-A[,,t]*sign(mean(diag(A[,,t])))
    diag(A[,,t])<-0 }
  dimnames(A) <- dimnames(Y)

  # receiver influence
  beta<-tab[-(1:(q+p-1))]
  B<-amprod(W,t(beta),3 )[,,1,]
  for(t in 1:dim(B)[3]){
    B[,,t]<-B[,,t]*sign(mean(diag(B[,,t])))
    diag(B[,,t])<-0 }
  dimnames(B) <- dimnames(Y)

  # pred effs
  pred = exp(eta_tab_dyn(tab,W,X,Z) )
  ll = mll_poisblr_dyn(tab, Y, W, X, Z)  

  if(calcSE){
    gh = mll_gH_dyn(tab,Y,W,X,Z) 
    se = sqrt(diag(solve(gh$hess))) 
    ses = sqrt( diag( solve(gh$hess) %*% gh$shess %*% solve(gh$hess) ))
    # org output
    summ = cbind(coef=tab, se=se, ses=ses, 
      t_se=tab/se, t_ses=tab/ses) }

  if(!calcSE){ summ = cbind(coef=tab) }  

  # return
  return(list(
    summ=summ, A=A, B=B, pred=pred, ll=ll,
    ALPHA=mod$ALPHA, BETA=mod$BETA, THETA=mod$THETA, 
    DEV=mod$DEV
    ))
}