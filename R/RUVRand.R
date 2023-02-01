# Y: A metabolomics data matrix with samples in rows and metabolites in columns
# ctl: A logical vector indicating the controls
# lambda: regularization parameter for the unwanted variation.
# k: rank of the unwanted variation term.

RUVRand <- function(Y, ctl,lambda=NULL, k=NULL,plotk=FALSE,...){
  
  Yc<-Y[, ctl]
  svdYc <- svd(Yc)
  fullW <- svdYc$u %*% diag(svdYc$d)  
  if(plotk){
    barplot(prcomp(Yc,scale. =T)$sdev^2/
              sum(prcomp(Yc,scale. =T)$sdev^2),
#            xlim=c(0,10),
#            names.arg =c(1:9),
            ylim=c(0,1),
            xlab="k",
            ylab="proportion of the variance",
            cex.lab=1.2,cex.axis=1.2,las=2)
    
  }
  if(is.null(k))
    stop('k must be entered')    
  if (!is.null(k) & is.null(lambda)){ 
    optklambda <- opt(ktry=k, W=fullW,Yc=Yc)
    lambda <- optklambda$optmat[,3]     
  } else optklambda <- NULL
  
  W<-fullW[,1:k,drop=FALSE]
  alpha<-solve(t(W)%*%W + lambda*diag(k), t(W) %*% Y)
  uvcomp<-W %*% alpha
  newY <- Y - uvcomp 
return(list(unadjY=Y,newY=newY,UVcomp=uvcomp,W=W,alpha= alpha,opt=optklambda,
            k=k,lambda=lambda,ctl=ctl))  
}

#kvec : A numerical vector with values of k for which lambda needs to be estimated
# W: A data matrix of unwanted variation with samples in rows and factors in columns
# Yc: A data matrix with samples in rows and quality control metabolites in columns
opt <- function(ktry,W,Yc){
  opt<-list()
  optmat<-matrix(NA,nrow=1,ncol=8)
  colnames(optmat)<-c("sigma2.a","sigma2.e","nu",
                      "lower_sigma2.a","upper_sigma2.a",
                      "lower_sigma2.e","upper_sigma2.e",
                      "convergence")
  
  opt<-optim(c(0.1,0.1),
             loglik,
             Y=t(Yc),
             W=W[,1:ktry,drop=FALSE],
             hessian=T)    
  fisher_info<-solve(opt$hessian/2)
  se<-sqrt(diag(fisher_info))
  upper_par1<-opt$par[1]+1.96*se[1]
  lower_par1<-opt$par[1]-1.96*se[1]
  upper_par2<-opt$par[2]+1.96*se[2]
  lower_par2<-opt$par[2]-1.96*se[2]
  
  optmat[1,]<-c(opt$par[1],
                opt$par[2],
                opt$par[2]/opt$par[1],
                lower_par1, upper_par1,
                lower_par2, upper_par2,                  
                opt$convergence)
  
  rownames(optmat)<-ktry
  return(list(optmat=optmat, opt=opt))
}

# Y: A metabolomics data matrix with samples in rows and metabolites in columns
# W: A data matrix of unwanted variation with samples in rows and factors in columns
loglik<- function (par, Y,W) 
{
  m <- ncol(Y)
  n<-nrow(Y)
  if (!is.null(W)){
    sigma2.a<-par[1]
    sigma2.e<-par[2]  
    if ((sigma2.a<0)|(sigma2.e<0))
      return(1e6)
    Sigma<-sigma2.a*(W%*%t(W))+sigma2.e*diag(m)
  } else{
    Sigma <- diag(m) 
    Sigma[upper.tri(Sigma, diag=TRUE)] <- par 
    Sigma <- Sigma + t(Sigma) - diag(diag(Sigma)) 
  }
  ed = eigen(Sigma, symmetric = TRUE)
  ev = ed$values
  if (!all(ev >= -1e-06 * abs(ev[1])))
    return(1e6)
  mu<-rep(0,m)
  centeredx<-sweep(Y,2,mu,"-")
  ssnew<-  t(centeredx)%*%(centeredx)
  if (is.null(tryCatch(solve(Sigma), error=function(e) NULL)))
    return(1e6)
  else
    inv.Sigma<-solve(Sigma) 
  Sigmainvss<-inv.Sigma%*%ssnew
  return(n*determinant(Sigma,logarithm=T)$mod+sum(diag(Sigmainvss)))
  
}
