### Title:    CCPCA functions code
### Author:   Rosember Guerra-Urzola
### Created:  09-07-2021
### Modified: 

CCPCA = function(X,K,rho,alpha,Corr=NULL,W=NULL,P=NULL, comp=TRUE,maxiter= 1000, tolerance = 1e-3,
                 normalized=FALSE,progress=FALSE){
  # Function to estimate the components weigths using CCPCA
  ## Inputs:
  # X: is a data set of dimension IxJ.
  # k: is number of components. k<=min(I,J)
  # rho: is de cardinality. It must be rho=<JxK
  # W: is the initial value for W.
  # P: is the initial value for P.
  # comp: a logical value. By default TRUE,means that the cardinality is per
  #       component. FALSE means the cardinality is regarding the full weights matrix.
  # maxiter: is the maximum number of iterations
  # tolerance: is the stop criteria
  # progress: If TRUE, print the current progress.
  ## Outputs:
  # W: Sparse component weights
  
  # Define initial values #
  rho_length=length(rho)
  if(rho_length !=1 & rho_length !=K){
    stop('rho is a vector with dimension either 1 or K')
  }
  if(is.null(Corr)){
  XtX = t(X)%*%X} else{
    XtX = Corr
  }
  # if(is.null(alpha)){
  # # eigenvalues = eigen(XtX,only.values = TRUE)$values
  # alpha = max(eigen(XtX,only.values = TRUE)$values)}
  # 
  if(is.null(W) & is.null(P)){
    W = svd(X, nv=K)$v
    P = W
  } else if(is.null(W) & !is.null(P)){
    W = svd(X, nv=K)$v  
  } else if(!is.null(W) & is.null(P)){
      P = updateP(W,XTX = XtX)
    }
  lossf = rep(0,maxiter+1)
  iter = 1
  while (1) {
    # message(iter)
    iter = iter+1
    W = updateW(P=P,W=W,rho,alpha,XTX=XtX)
    P = updateP(W,XTX=XtX)
    lossf[iter] = norm(X -X%*%W%*%t(P),type = 'F')
    # message(lossf[iter])
    if(progress){
      message(paste0('Iter:',iter,', Loss:',lossf[iter]))
    }
   
    if (iter>1 & abs(lossf[iter]-lossf[(iter-1)])/abs(lossf[(iter-1)])<tolerance|| iter>maxiter){
      break
    }
    
  }
  
  vn <- dimnames(X)[[2]]
  dimnames(P) <- list(vn, paste("P", 1:dim(W)[2], sep = ""))
  dimnames(W) <- list(vn, paste("W", 1:dim(W)[2], sep = ""))
  output <- list( W = W,P = P)
  return(output)
}

updateW = function(P,W,rho,alpha,X=NULL,XTX=NULL,comp=TRUE,normalized=FALSE){
  # function to update W given P
  # P are the component loadings
  # W is the current value for W
  # rho the dacardinality of the component weights
  # X the data set
  # XTX is the correlation matrix of X.
  
  if (is.null(XTX)){
    if(is.null(X)){
      stop('Either X or XTX is needed')
    }
    XTX = t(X)%*%X
  }
  # 
  # K= dim(W)[2]
  # J= dim(W)[1]
  # B= matrix(rep(J*K),nrow = J,ncol = K)
  # for (k in 1:K) {
  #   non_zero_index_wk = W[,k]!=0
  #   B[,k]=W[,k]-(XTX[,non_zero_index_wk]%*%W[non_zero_index_wk,k]-XTX%*%P[,k])/alpha
  # }
  B=W-(XTX%*%W-XTX%*%P)/alpha # divided by the sample size?
  # W = matrix(truncatedW(b,rho),byrow = FALSE,nrow = dim(W)[1], ncol = dim(W)[2])
  W =truncated(B,rho,comp,normalized)
  return(W)
}
updateP = function(W,X=NULL,XTX=NULL){
  # function to update P given W
  # X is the data set
  if (is.null(XTX)){
    if(is.null(X)){
      stop('Either X or XTX is needed')
    }
    XTX = t(X)%*%X
  } 
  A<-t(W)%*%XTX
  svdA<-svd(A)
  P<-t(svdA$u%*%t(svdA$v))
  return(P)
}
truncated = function(W,rho, comp=TRUE,normlized=FALSE){# Check if there is a more efficent algorithm !!!
  # W is a vector or matrix
  # rho is the cardinality per component or the full matrix
  K= dim(W)[2]
  rho_length =length(rho)
  if(rho_length>K){stop('The cardinality vector length can not be greater than
                        the number of components')}
  if(comp){
    if(rho_length<K){
      rho = rep(rho,K)
    } 
    for (k in 1:K) {
      ind_order=sort(abs(W[,k]),decreasing=TRUE,index.return=TRUE)
      W[-ind_order$ix[1:rho[k]],k] =0
    }
  }else{
  ind_order=sort(abs(as.vector(W)),decreasing=TRUE,index.return=TRUE)
  W[-ind_order$ix[1:rho]] = 0 
  }
  if(normlized){
    normW=sqrt(apply(W^2, 2,sum))
    W = W%*%diag(1/normW)
  }
  return(W)
}

pev=function(W,X,totalvariance=NULL){
  # function to estimate the pev of the sparse component weights
  ## inputs
  # W: sparse component weights
  # X: data set
  # totalvariance: total variance in the data set.
  ## outputs
  # pev: vector with the pev per component
  
  
  if(is.null(totalvariance)){
    # eigenvalues = eigen(t(X)%*%X,only.values = TRUE)$values
    totalvariance = sum(eigen(t(X)%*%X,only.values = TRUE)$values)
    # totalvariance <- sum(eigenvalues)
  }
  
  
  Z <- X %*% W
  R <- qr.R(qr(Z))
  pev <- diag(R^2)/totalvariance
  return(pev)
}

TotalRecovery = function(W,supp){
  # W: the estimated component weigths
  # supp: support of the population component weights
  J=dim(W)[2]
  I=dim(W)[1]
  PerM = gtools::permutations(J,J,1:J)
  REC = rep(0,dim(PerM)[1])
  for (i in 1:dim(PerM)[1]) {
    supp_w = W[,PerM[i,]] !=0
    REC[i] = sum(supp_w == supp)/prod(I,J)
  }
  return(max(REC))
}

