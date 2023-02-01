context('Run ruvrand to run ruv random on edata')

test_that('RUVrand returns the same abundance values as bc_ruvRandom',{
  
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
                                emeta_cname = 'Metabolite')
  # already on log scale
  attributes(mdata)$data_info$data_scale_orig <- "log2"
  attributes(mdata)$data_info$data_scale <- "log2"
  
  # remove molecules with too much missingness
  molfilt <- pmartR::molecule_filter(mdata)
  mdata2 <- pmartR::applyFilt(molfilt,mdata,min_num = 2)
  impObj <- imputation(mdata2)
  mdataImp <- apply_imputation(impObj,mdata2)

  # find the values needed for the normalize function
  edat <- as.matrix(mdataImp$e_data[,-1]) %>%
    t()
  molecules <- mdataImp$e_data[,1]
  
  # find the parameter ctl (the negative controls)
  ctlRUV <- mdataImp$e_meta$IS == "IS"
  
  # now create the ruv-random abundance values!
  # we set lambda to be null, plotk to be FALSE
  ruvObj <- malbacR:::RUVRand(Y = edat,
                             ctl = ctlRUV,
                             k = 3,
                             plotk = FALSE,
                             lambda = NULL)
  # to match bc_ruvRandom setup, transpose the data
  ruvObjT <- ruvObj$newY %>%
    t()
  # remove columns with internal standards
  ruvNoIS <- ruvObjT[!ctlRUV,]
  
  # now run using bc_ruvRandom for comparison
  mdata_ruv <- bc_ruvRandom(mdataImp,"IS","IS")
  # create it as a matrix to compare with manual version
  bc_ruv <- as.matrix(mdata_ruv$e_data[,-1])
  
  # run the test
  expect_equal(ruvNoIS,bc_ruv)
  # make sure the data is of class list
  expect_equal(class(ruvObj),'list')
})
