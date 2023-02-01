context('Run normalize_qcrlsc which is a subsection of bc_qcrlsc')

test_that('normalize_qcrlsc runs QC-RLSC on the data', {
  
  # Load the reduced metabolite data frames ------------------------------------
  load(system.file('testdata',
                   'mini_udn_four.Rdata',
                   package = 'malbacR'))
  
  # Run as.metabData with agreeable data frames --------------------------------
  
  # Construct a metabData object with the edata, fdata, and emeta data frames.
  mdata4 <- pmartR::as.metabData(e_data = edata,
                                 f_data = fdata,
                                 e_meta = emeta,
                                 edata_cname = 'Metabolite',
                                 fdata_cname = 'SampleID',
                                 emeta_cname = 'Metabolite')
  attributes(mdata4)$data_info$data_scale_orig <- "log2"
  attributes(mdata4)$data_info$data_scale <- "log2"
  
  # run the get_params function ------------------------------------------------
  params_obj = malbacR:::get_params(omicsData = mdata4,
                          block_cname = "BatchNum",
                          qc_cname = "QC",
                          qc_ind = "QC.NIST",
                          order_cname = "RunOrderOverall",
                          missing_thresh = 0.5,
                          rsd_thresh = 0.3)
  
  # need to remove the bad features
  cfilt <- pmartR::custom_filter(mdata4,e_data_remove = params_obj$bad_feats)
  mdataFilt <- pmartR::applyFilt(cfilt,mdata4)
  
  # check assumptions ----------------------------------------------------------
  # check that optimal params is data.frame with expected elements
  expect_error(malbacR:::normalize_qcrlsc(omicsData = mdataFilt,
                                optimal_params = params_obj,
                                block_cname = "BatchNum",
                                qc_cname = "QC",
                                qc_ind = "QC.NIST",
                                backtransform = FALSE),
               "Input parameter optimal_params must be a data.frame")
  
  # check that the data.frame is the right object based on column names
  # remove span column just for demonstration of error message
  no_span <- params_obj$final_ests
  no_span <- no_span[,-3]
  expect_error(malbacR:::normalize_qcrlsc(omicsData = mdataFilt,
                                optimal_params = no_span,
                                block_cname = "BatchNum",
                                qc_cname = "QC",
                                qc_ind = "QC.NIST",
                                backtransform = FALSE),
               "Input parameter optimal_params must be the final_ests element from the output from 'get_params'")
  
  # block_cname needs to be a column in the fdata
  expect_error(malbacR:::normalize_qcrlsc(omicsData = mdataFilt,
                                optimal_params = params_obj$final_ests,
                                block_cname = "BatchData",
                                qc_cname = "QC",
                                qc_ind = "QC.NIST",
                                backtransform = FALSE),
               "The f_data component of omicsData_filt must contain a column for 'Batch'")
  
  # block_cname needs to be either numeric or factor
  expect_error(malbacR:::normalize_qcrlsc(omicsData = mdataFilt,
                                optimal_params = params_obj$final_ests,
                                block_cname = "BatchName",
                                qc_cname = "QC",
                                qc_ind = "QC.NIST",
                                backtransform = FALSE),
               "Values in block_cname column")

  # qc_cname needs to be a column in the fdata
  expect_error(malbacR:::normalize_qcrlsc(omicsData = mdataFilt,
                                optimal_params = params_obj$final_ests,
                                block_cname = "BatchNum",
                                qc_cname = "QC_data",
                                qc_ind = "QC.NIST",
                                backtransform = FALSE),
               "The f_data component of omicsData_filt must contain a column for 'QCtype'")
  
  # qc_val needs to be a value in qc_name
  expect_error(malbacR:::normalize_qcrlsc(omicsData = mdataFilt,
                                optimal_params = params_obj$final_ests,
                                block_cname = "BatchNum",
                                qc_cname = "QC",
                                qc_ind = "QC.Pooled",
                                backtransform = FALSE),
               "The value for qc_ind is not present in the column specified by qc_cname")
  
  # backtransform needs to be logical
  expect_error(malbacR:::normalize_qcrlsc(omicsData = mdataFilt,
                                optimal_params = params_obj$final_ests,
                                block_cname = "BatchNum",
                                qc_cname = "QC",
                                qc_ind = "QC.NIST",
                                backtransform = 2),
               "Input parameter backtransform must be of class 'logical'")
  
  # run normalize_qcrlsc function ----------------------------------------------
  mdata_qcrlsc <- malbacR:::normalize_qcrlsc(omicsData = mdataFilt,
                                   optimal_params = params_obj$final_ests,
                                   block_cname = "BatchNum",
                                   qc_cname = "QC",
                                   qc_ind = "QC.NIST",
                                   backtransform = FALSE)
  
  # dimension check time -------------------------------------------------------
  # check that we have a list of length 3 like we do with pmart object
  expect_equal(length(mdataFilt), length(mdata_qcrlsc))
  expect_equal(is.list(mdata_qcrlsc),TRUE)
  
  # now let's look at dimension of mdata_qcrlsc
  # we have removed the QC samples so find how many QC samples should have been removed
  numQC <- sum(mdataFilt$f_data$QC == "QC.NIST")
  # edata dimensions
  expect_equal(nrow(mdataFilt$e_data),nrow(mdata_qcrlsc$e_data))
  expect_equal(ncol(mdataFilt$e_data),ncol(mdata_qcrlsc$e_data)+numQC)
  # fdata dimensions
  expect_equal(nrow(mdataFilt$f_data),nrow(mdata_qcrlsc$f_data)+numQC)
  expect_equal(ncol(mdataFilt$f_data),ncol(mdata_qcrlsc$f_data))
  # emeta dimensions
  expect_equal(nrow(mdataFilt$e_meta),nrow(mdata_qcrlsc$e_meta))
  expect_equal(ncol(mdataFilt$e_meta),ncol(mdata_qcrlsc$e_meta))
  
  # check check.names
  expect_equal(attr(mdataFilt,"check.names"),attr(mdata_qcrlsc,"check.names"))
  # check meta_info
  expect_equal(attr(mdataFilt,"meta_info"),attr(mdata_qcrlsc,"meta_info"))
  # check filters
  # we should have one more filter now in the mdata_qcrlsc pmart object than in mdataFilt
  expect_equal(length(attr(mdataFilt,"filters"))+1,length(attr(mdata_qcrlsc,"filters")))
  # check group_DF
  expect_equal(attr(mdataFilt,"group_DF"),attr(mdata_qcrlsc,"group_DF"))
})

