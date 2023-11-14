normFact <- function(fact,X,ref,refType,k=20,t=0.5,ref2=NULL,refType2=NULL,t2=0.5,alpha,...) {
  
  if (fact=='stICA'){
    obj = unbiased_stICA(X,k,alpha=alpha)
    B=obj$B
    A=obj$A
  } else if (fact == 'SVD'){
    obj = svd(X,nu=k,nv=k)
    A = obj$u%*%diag(obj$d[1:k],k)
    B = obj$v
  } else {
    stop("Factorization method should be SVD or stICA")
  }
  
  factR2 = R2(ref,B,refType,pval=T)
  
  idx = which(factR2$allpv<t)
  
  if (t <0 | t>1){stop("t not in [0 1]")}
  
  if (!is.null(ref2)){
    if (sum(t2 <0 | t2>1)){stop("t2 not in [0 1]")}
    factR2_2 = R2(ref2,B,refType2,pval=T)
    idx_2 = c()
    if (length(t2)!=length(refType2)){
      if(length(t2)==1){
        t2=rep(t2,length(refType2))
      }
      else {
        stop("length(t2) sould be equal to 1 or length(refType2)")
      }
    }
    for (i in 1:length(refType2)){
      
      idx_2 = c(idx_2,which(factR2_2$allpv[,i]<t2[i]))
    }
    
    idx2keep =intersect(idx,idx_2)
    print(paste("Keeping",length(idx2keep), "cmpts with P value less than t2"))
    idx = setdiff(idx, idx2keep)
    
  }
  
  bestcmptA = A[,idx]
  bestcmptB = B[,idx]
  
  print(paste("Removing",length(idx),"components with P value less than",t))
  
  Xn = X - bestcmptA%*% t(bestcmptB)
  
  R2=factR2$allR2
  if (!is.null(ref2)){
    R2 = cbind(R2,factR2_2$allR2)
  }
  return(list(Xn=Xn,R2=R2,bestSV =bestcmptB,A=A,B=B))
}

unbiased_stICA <- function(X,k=10,alpha) {
  library(JADE)
  library(corpcor)
  
  jadeCummulantMatrices <- function(X) {
    n <- nrow(X)
    t <- ncol(X)
    
    M <- array(0,c(n,n,n*(n+1)/2))
    scale <- matrix(1,n,1)/t  # for convenience
    
    R <- cov(t(X)) # covariance
    k <- 1
    for (p in 1:n){
      #case q=p
      C <- ((scale %*% (X[p,]*X[p,]))*X) %*% t(X)
      E <- matrix(0,n,n)
      E[p,p] <- 1
      M[,,k] <- C - R %*% E %*% R - sum(diag(E %*% R)) * R - R %*% t(E) %*% R
      k <- k+1
      #case q<p
      if (p > 1) {
        for (q in 1:(p-1)){
          C <- ((scale %*% (X[p,]*X[q,]))*X) %*% t(X) * sqrt(2)
          E <- matrix(0,n,n)
          E[p,q] <- 1/sqrt(2)
          E[q,p] <- E[p,q]
          M[,,k] <- C - R %*% E %*% R - sum(diag(E %*% R)) * R - R %*% t(E) %*% R
          k <- k+1
        }
      }
    }
    return(M)
  }
  p <- nrow(X)
  n <- ncol(X)
  
  dimmin <- min(n,p)
  
  if (dimmin < k) {
    k <- dimmin
  }
  if (alpha <0 | alpha >1){
    stop("alpha not in [0 1]")
  }
  
  # Remove the spatiotemporal mean
  Xc <- X - matrix(rep(colMeans(X,dims=1),p),nrow = p,byrow=T);
  Xc <- Xc - matrix(rep(rowMeans(Xc,dims=1),n),nrow = p);
  
  # SVD of Xc and dimension reduction: keeping only the k first
  # components
  udv <- svd(Xc,k,k)
  D <- diag(udv$d[1:k]); if (k==1) {D <- udv$d[1]}
  U <- udv$u;
  V <- udv$v;
  
  # Estimation of the cumulant matrices
  nummat <- k*(k+1)/2;
  M <- array(0,c(k,k,2*nummat));
  Bt <- D^(1-alpha) %*% t(V)
  if (alpha == 1) { Bt <- t(V)}
  At <- D^(alpha) %*% t(U)
  if (alpha == 0) { At <- t(U)}
  M[,,1:nummat] <- jadeCummulantMatrices(Bt);
  M[,,(nummat+1):(2*nummat)] <- jadeCummulantMatrices(At)
  
  # normalization within the groups in order to allow for comparisons using
  # alpha
  M[,,1:nummat] <- alpha*M[,,1:nummat]/mean(sqrt(apply(M[,,1:nummat]*M[,,1:nummat],3,sum)));
  M[,,(nummat+1):(2*nummat)] <- (1-alpha)*M[,,(nummat+1):(2*nummat)]/mean(sqrt(apply(M[,,(nummat+1):(2*nummat)]*M[,,(nummat+1):(2*nummat)],3,sum)));
  
  # Joint diagonalization
  Worth <- rjd(M,eps = 1e-06, maxiter = 1000);
  Wo <-t (Worth$V);
  
  #     Computation of A and B
  A0 <- U %*% D^(alpha) %*% solve(Wo);
  B0 <- V%*% D^(1-alpha) %*% t(Wo);
  if (alpha == 1) { B0 <- V %*% t(Wo)}
  if (alpha == 0) { A0 <- U %*% solve(Wo)}
  
  # Add transformed means
  meanCol <- matrix(colMeans(X,dims=1),ncol =1); # spatial means
  meanRows <- matrix(rowMeans(X,dims=1),ncol = 1); # temporal means
  
  meanB <- pseudoinverse(A0) %*% (meanRows);
  meanA <- pseudoinverse(B0) %*% (meanCol);
  
  Bfin <- B0 + matrix(rep(meanB,n),nrow = n,byrow=T)
  Afin <- A0 + matrix(rep(meanA,p),nrow = p,byrow=T)
  
  return(list(A=Afin,B=Bfin,W=Wo))
}

R2 <- function(poi,V,poiType,pval = T) {
  if (is.vector(V) ){ V = matrix(V,ncol=1)}
  if (is.vector(poi)){poi = matrix(poi,nrow =length(poi))}
  
  p = nrow(V)   # number of samples
  k = ncol(V)    # number of components
  l = length(poiType)  # number of cf/poi to test
  if (is.null(l)){stop("POI type(s) neeeded")}
  p2 = nrow(poi)
  l2 = ncol(poi)
  
  if( l2 != l){ # checking poi and poiType dimensions compatiblity
    if (p2 == l){  # if poi is transposed (l*p)
      poi = t(poi)
      warning("Transposing poi to match poiType dimension")
      p2 = nrow(poi)
    } else {
      print(poi)
      print(poiType)
      stop("poi dimensions doesn't match poiType dimension")
    }
  }
  
  if (p != p2){ # checking poi and V dimensions compatiblity
    if (p2 == k){
      warnings("Transposing V to match poi dimension")
      V =t(V)
      k = p
      p = p2
    } else {
      stop("poi and V dimensions incompatible")
    }
  }
  
  R2 = rep(-1,l)
  names(R2) = colnames(poi)
  idxcorr = R2
  R2_tmp <- matrix(rep(-1,k*l),k,l,dimnames=list(colnames(V),colnames(poi)))    # r2_tmp(k,l) hold the R2 value for column k of V with poi l
  
  if (pval){
    pv = R2
    idxcorr2 = R2
    pv_tmp <- R2_tmp   # r2_tmp(k,l) hold the R2 value for column k of V with poi l
  }
  
  for (cmpt in 1:k){    # for each column of V
    cmpt2an <- V[,cmpt]
    for (ipoi in 1:l){
      idx_finite = is.finite(as.factor(poi[,ipoi]))
      poi2an = poi[idx_finite,ipoi]
      cmpt2an_finite=cmpt2an[idx_finite]
      if (poiType[ipoi] == "continuous") {  # estimation by linear regression
        coefs <- coef(lm(cmpt2an_finite~as.numeric(poi2an)))
        cmpt2an_est <- coefs[2]*as.numeric(poi2an)+coefs[1]
        nc <- 2;
      } else if (poiType[ipoi]=="categorical"){  # estimation by classe mean
        classes <- unique(poi2an)
        nc <- length(classes)
        cmpt2an_est <- rep(NA,length(cmpt2an_finite))
        for (icl in 1:length(classes) ){
          idxClasse <- which(poi2an==classes[icl])
          cmpt2an_est[idxClasse] <- mean(cmpt2an_finite[idxClasse])
        }
      } else {
        stop("Incorrect poiType. Select 'continuous' or 'categorical'. ")
      }
      sse <- sum((cmpt2an_finite-cmpt2an_est)^2)
      sst <- sum((cmpt2an_finite-mean(cmpt2an_finite))^2)
      R2_tmp[cmpt,ipoi] <-  1 - sse/sst
      if (pval){
        F <- ((sst-sse)/(nc-1))/(sse/(p-nc))
        pv_tmp[cmpt,ipoi] = 1-pf(F,nc-1,p-nc);
        if (!is.finite(pv_tmp[cmpt,ipoi])) {
          warning(paste("Non finite p-value for component ",cmpt," (pv=",pv_tmp[cmpt,ipoi],", F=",F,"), assigning NA", sep=""))
          pv_tmp[cmpt,ipoi] <- NA
        }
      }
    }
  }
  
  for (ipoi in 1:l){
    if (pval){
      pv[ipoi] <- min(pv_tmp[,ipoi])
      idxcorr2[ipoi] <- which(pv_tmp[,ipoi] == pv[ipoi])[1]   # if more than one component gives the best R2, takes the first one
    }
    R2[ipoi] <- max(R2_tmp[,ipoi])
    idxcorr[ipoi] <- which(R2_tmp[,ipoi] == R2[ipoi])[1]   # if more than one component gives the best R2, takes the first one
  }
  
  if (pval){
    return(list(R2=R2,idxcorr=idxcorr,allR2 = R2_tmp, pv=pv,idxcorr2=idxcorr2,allpv = pv_tmp))
  } else {
    return(list(R2=R2,idxcorr=idxcorr,allR2 = R2_tmp))
  }
}

WaveICA<-function(data,wf="haar",batch,group=NULL,K=20,t=0.05,t2=0.05,alpha=0){
  ### Wavelet Decomposition
  library(waveslim)
  level<-floor(log(nrow(data),2))
  if (is.null(colnames(data))){
    stop("data must have colnames")
  }
  coef<-list()
  for (k in 1:(level+1)){
    coef[[k]] <-matrix(NA,nrow(data),ncol(data))
  }
  for (j in 1:ncol(data)){
    cat(paste("######Decomposition",j,"########\n"))
    data_temp<-data[,j]
    x_modwt<-modwt(data_temp,wf=wf,n.levels =level)
    for (k in 1:(level+1)){
      coef[[k]][,j]<-x_modwt[[k]]
    }
  }
  ##### ICA
  index<-level+1
  data_wave_ICA<-list()
  for (i in (1:index)){
    cat(paste("######### ICA",i,"#############\n"))
    data_coef<-coef[[i]]
    data_coef_ICA<-normFact(fact="stICA",X=t(data_coef),ref=batch,refType ="categorical",k=K,t=t,ref2=group,refType2="categorical",t2=t2,alpha)
    data_wave_ICA[[i]]<-t(data_coef_ICA$Xn)
  }
  ### Wavelet Reconstruction
  index<-ncol(data)
  index1<-length(data_wave_ICA)
  data_coef<-matrix(NA,nrow(data_wave_ICA[[1]]),index1)
  data_wave<-matrix(NA,nrow(data_wave_ICA[[1]]),ncol(data_wave_ICA[[1]]))
  for (i in 1:index){
    cat(paste("######Reconstruction",i,"########\n"))
    for (j in 1:index1){
      data_coef[,j]<-data_wave_ICA[[j]][,i]
    }
    data_temp<-data[,i]
    data_coef<-as.data.frame(data_coef)
    colnames(data_coef)<-c(paste("d",1:(index1-1),sep=""),paste("s",(index1-1),sep=""))
    y<-as.list(data_coef)
    attributes(y)$class<-"modwt"
    attributes(y)$wavelet<-wf
    attributes(y)$boundary<-"periodic"
    data_wave[,i]<-imodwt(y)+mean(data_temp)
  }
  rownames(data_wave)<-rownames(data)
  colnames(data_wave)<-colnames(data)
  return(list(data_wave=data_wave))
}