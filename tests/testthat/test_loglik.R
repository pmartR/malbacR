context('Run ruvrand to run ruv random on edata')

test_that('loglik function in RUVrand returns correct value',{
  
  # Load the reduced peptide data frames ---------------------------------------
  load(system.file('testdata',
                   'mini_udn.Rdata',
                   package = 'malbacR'))
  
  # Run as.metabData with agreeable data frames --------------------------------
  
  # Construct a metabData object with the edata, fdata, and emeta data frames.
  mdata <- pmartR::as.metabData(e_data = edata,
                                f_data = fdata,
                                e_meta = emeta,
                                edata_cname = 'Metabolite',
                                fdata_cname = 'SampleID',
                                emeta_cname = 'Metabolite',
                                data_scale = 'log2')
  molfilt <- pmartR::molecule_filter(mdata)
  mdata <- pmartR::applyFilt(molfilt,mdata)
  impObj <- imputation(mdata)
  mdata <- apply_imputation(impObj,mdata)
  # data manipulation for setting up ruv-random
  edat <- as.matrix(mdata$e_data[,-1]) %>%
    t()
  molecules <- mdata$e_data[,1]
  
  # find the parameter ctl (the negative controls)
  edat_cname = pmartR::get_edata_cname(mdata)
  proper_order_ctl <- mdata$e_data %>%
    dplyr::select(dplyr::all_of(edat_cname)) %>%
    dplyr::left_join(mdata$e_meta)
  nc_cname = "IS"; nc_val = "IS"; k = 3
  nc_cnameNum <- which(colnames(proper_order_ctl) == nc_cname)
  ctlRUV <- proper_order_ctl[,nc_cnameNum] == nc_val
  # set up parameters for the ruv-random and opt functions
  Y = edat; ctl = ctlRUV; k = k; plotk = FALSE; lambda = NULL
  Yc<-Y[, ctl]
  svdYc <- svd(Yc)
  fullW <- svdYc$u %*% diag(svdYc$d)  
  W = fullW; ktry = k; Yc = Yc

  # function version  
  loglik_function <- malbacR:::loglik(par = c(0.1,0.1),Y = t(Yc),W = W[,1:ktry,drop = FALSE])
  
  # manual version
  Y = t(Yc); W = W[,1:ktry, drop = FALSE]; par = c(0.1,0.1)
  m <- ncol(Y)
  n<-nrow(Y)
  sigma2.a<-par[1]
  sigma2.e<-par[2]  
  Sigma<-sigma2.a*(W%*%t(W))+sigma2.e*diag(m)
  
  ed = eigen(Sigma, symmetric = TRUE)
  ev = ed$values
  mu<-rep(0,m)
  centeredx<-sweep(Y,2,mu,"-")
  ssnew<-  t(centeredx)%*%(centeredx)
  inv.Sigma<-solve(Sigma) 
  Sigmainvss<-inv.Sigma%*%ssnew
  loglik_manual <- (n*determinant(Sigma,logarithm=T)$mod+sum(diag(Sigmainvss)))
  
  # compare function vs manual
  expect_equal(loglik_function,loglik_manual)
})
