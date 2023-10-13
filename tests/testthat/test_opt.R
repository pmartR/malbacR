context('Run ruvrand to run ruv random on edata')

test_that('opt function in RUVrand returns list of',{
  
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
  molfilt <- molecule_filter(mdata)
  mdata <- applyFilt(molfilt,mdata)
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
  nc_cnameNum <- which(colnames(proper_order_ctl) == "IS")
  ctlRUV <- proper_order_ctl[,nc_cnameNum] == "IS"
  
  # set up parameters for the opt function
  Y = edat; ctl = ctlRUV; k = 3
  
  # now set up opt function for use
  Yc<-Y[, ctl]
  svdYc <- svd(Yc)
  fullW <- svdYc$u %*% diag(svdYc$d)  
  ktry=k; W=fullW; Yc=Yc
  
  # function version
  optklambda <- malbacR:::opt(ktry=k, W=fullW,Yc=Yc)
  
  # manual version
  opt<-list()
  optmat<-matrix(NA,nrow=1,ncol=8)
  colnames(optmat)<-c("sigma2.a","sigma2.e","nu",
                      "lower_sigma2.a","upper_sigma2.a",
                      "lower_sigma2.e","upper_sigma2.e",
                      "convergence")
  
  opt<-optim(c(0.1,0.1),
             malbacR:::loglik,
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
  opt_manual <- list(optmat=optmat, opt=opt)
  
  # compare function vs manual
  expect_equal(optklambda,opt_manual)
})

