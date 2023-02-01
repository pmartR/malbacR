context('Run get_params which is a subsection of bc_qcrlsc')

test_that('get_params obtains the parameters for QC-RSLC signal correction', {
  
  # Load the reduced metabolite data frames ------------------------------------
  load(system.file('testdata',
                   'mini_udn.Rdata',
                   package = 'malbacR'))
  
  # Run as.metabData with agreeable data frames ----------------------------------
  
  # Construct a metabData object with the edata, fdata, and emeta data frames.
  mdata <- pmartR::as.metabData(e_data = edata,
                                f_data = fdata,
                                e_meta = emeta,
                                edata_cname = 'Metabolite',
                                fdata_cname = 'SampleID',
                                emeta_cname = "Metabolite")
  attributes(mdata)$data_info$data_scale_orig <- "log2"
  attributes(mdata)$data_info$data_scale <- "log2"
  
  # Run through the different warnings! ----------------------------------------
  # block_cname needs to be a column in the fdata
  expect_error(malbacR:::get_params(omicsData = mdata,
                         block_cname = "BatchData",
                         qc_cname = "QC",
                         qc_ind = "QC.NIST",
                         order_cname = "RunOrderOverall",
                         missing_thresh = 0.5,
                         rsd_thresh = 0.3),
               "The f_data component of omicsData must contain a column for 'Batch'")
  
  # block_cname needs to be either numeric or factor
  expect_error(malbacR:::get_params(omicsData = mdata,
                         block_cname = "BatchName",
                         qc_cname = "QC",
                         qc_ind = "QC.NIST",
                         order_cname = "RunOrderOverall",
                         missing_thresh = 0.5,
                         rsd_thresh = 0.3),
               "Values in block_cname column")
  
  # qc_cname needs to be a column in the fdata
  expect_error(malbacR:::get_params(omicsData = mdata,
                          block_cname = "BatchName",
                          qc_cname = "QC_data",
                          qc_ind = "QC.NIST",
                          order_cname = "RunOrderOverall",
                          missing_thresh = 0.5,
                          rsd_thresh = 0.3),
               "The f_data component of omicsData must contain a column for 'QCtype'")
  
  # qc_val needs to be a value in qc_name
  expect_error(malbacR:::get_params(omicsData = mdata,
                          block_cname = "BatchName",
                          qc_cname = "QC",
                          qc_ind = "QC.Pooled",
                          order_cname = "RunOrderOverall",
                          missing_thresh = 0.5,
                          rsd_thresh = 0.3),
               "Input parameter qc_ind must be a value found in the qc_cname")
  
  # order_cname must be a column in the fdata
  expect_error(malbacR:::get_params(omicsData = mdata,
                          block_cname = "BatchName",
                          qc_cname = "QC",
                          qc_ind = "QC.NIST",
                          order_cname = "RunOrder",
                          missing_thresh = 0.5,
                          rsd_thresh = 0.3),
               "The f_data component of omicsData must contain a column for 'RunOrder'")
  
  # order_cname must be numeric
  expect_error(malbacR:::get_params(omicsData = mdata,
                          block_cname = "BatchNum",
                          qc_cname = "QC",
                          qc_ind = "QC.NIST",
                          order_cname = "BatchName",
                          missing_thresh = 0.5,
                          rsd_thresh = 0.3),
               "Values in order_cname column")
  
  # we need to load the bigger data set 
  # Load the reduced metabolite data frames ------------------------------------
  load(system.file('testdata',
                   'mini_udn_four.Rdata',
                   package = 'malbacR'))
  
  # Run as.metabData with agreeable data frames ----------------------------------
  
  # Construct a metabData object with the edata, fdata, and emeta data frames.
  mdata4 <- pmartR::as.metabData(e_data = edata,
                                f_data = fdata,
                                e_meta = emeta,
                                edata_cname = 'Metabolite',
                                fdata_cname = 'SampleID',
                                emeta_cname = 'Metabolite')
  attributes(mdata4)$data_info$data_scale_orig <- "log2"
  attributes(mdata4)$data_info$data_scale <- "log2"
  
  # run get_params function ----------------------------------------------------
  params_obj = malbacR:::get_params(omicsData = mdata4,
                      block_cname = "BatchNum",
                      qc_cname = "QC",
                      qc_ind = "QC.NIST",
                      order_cname = "RunOrderOverall",
                      missing_thresh = 0.5,
                      rsd_thresh = 0.3)
  
  # Dimension check time -------------------------------------------------------
  # check that we have a list output of length 2 for params_obj
  expect_equal(length(params_obj),2)
  expect_equal(is.list(params_obj),TRUE)
  
  # now let's look at dimension of params_obj$bad_feats
  # find the bad features by hand
  qc_cols = as.character(mdata4$f_data[which(mdata4$f_data[, "QC"] == "QC.NIST"), "SampleID"])
  qc_dat <- mdata4$e_data[, which(names(mdata4$e_data) %in% c("Metabolite", qc_cols))]
  names(qc_dat)[1] = "Metabolite"
  frac_missing = apply(is.na(qc_dat[, -1]), 1, function(x) sum(x) / length(x))
  bad_feats = qc_dat[which(frac_missing > 0.5), which(names(qc_dat) =="Metabolite")]
  
  # compare code with by hand
  expect_equal(bad_feats,params_obj$bad_feats)
  
  # now we look at dimensions of params_obj$final_ests
  expect_equal(colnames(params_obj$final_ests),c("Block","Metabolite","Span","Poly_Degree","MSE"))
  numBatch <- length(unique(mdata4$f_data$BatchNum))
  expect_equal(nrow(params_obj$final_ests),(nrow(mdata4$e_data)-length(bad_feats))*numBatch)
})

