context('Run QC-RLSC batch correction')

test_that('bc_qcrlsc returns the correct data frame and attributes', {
  
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
  expect_error(bc_qcrlsc(omicsData = mdata,
                         block_cname = "BatchData",
                         qc_cname = "QC",
                         qc_val = "QC.NIST",
                         order_cname = "RunOrderOverall"),
               "The f_data component of omicsData must contain a column for 'Batch'")
  
  # block_cname needs to be either numeric or factor
  expect_error(bc_qcrlsc(omicsData = mdata,
                         block_cname = "BatchName",
                         qc_cname = "QC",
                         qc_val = "QC.NIST",
                         order_cname = "RunOrderOverall"),
               "Values in block_cname column")
  
  # qc_cname needs to be a column in the fdata
  expect_error(bc_qcrlsc(omicsData = mdata,
                         block_cname = "BatchNum",
                         qc_cname = "QC_data",
                         qc_val = "QC.NIST",
                         order_cname = "RunOrderOverall"),
               "The f_data component of omicsData must contain a column for 'QCtype'")
  
  # qc_val needs to be a value in qc_name
  expect_error(bc_qcrlsc(omicsData = mdata,
                         block_cname = "BatchNum",
                         qc_cname = "QC",
                         qc_val = "QC.Pooled",
                         order_cname = "RunOrderOverall"),
               "Input parameter qc_val must be a value found in the qc_cname")
  
  # order_cname must be a column in the fdata
  expect_error(bc_qcrlsc(omicsData = mdata,
                         block_cname = "BatchNum",
                         qc_cname = "QC",
                         qc_val = "QC.NIST",
                         order_cname = "RunOrder"),
               "The f_data component of omicsData must contain a column for 'RunOrder'")
  
  # order_cname must be numeric
  expect_error(bc_qcrlsc(omicsData = mdata,
                         block_cname = "BatchNum",
                         qc_cname = "QC",
                         qc_val = "QC.NIST",
                         order_cname = "BatchName"),
               "Values in order_cname column")
  
  # keep_qc must be of length 1
  expect_error(bc_qcrlsc(omicsData = mdata,
                         block_cname = "BatchNum",
                         qc_cname = "QC",
                         qc_val = "QC.NIST",
                         order_cname = "RunOrderOverall",
                         keep_qc = c(TRUE,FALSE)),
               "Input parameter qc_val must be of length 1")
  # keep_qc must be logical
  expect_error(bc_qcrlsc(omicsData = mdata,
                         block_cname = "BatchNum",
                         qc_cname = "QC",
                         qc_val = "QC.NIST",
                         order_cname = "RunOrderOverall",
                         keep_qc = "True"),
               "Input parameter keep_qc must be logical")
  
  # error out if we do not have enough QC samples
  expect_error(bc_qcrlsc(omicsData = mdata,
                         block_cname = "BatchNum",
                         qc_cname = "QC",
                         qc_val = "QC.NIST",
                         order_cname = "RunOrderOverall"),
               "The following molecules have too few non-missing QC data points in at least one")
  
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
  
  # error out if we do not have a QC at the beginning and end of each batch
  # update some QC values so that we error out
  mdata4$f_data$QC[1] = "NotQC"
  expect_error(bc_qcrlsc(omicsData = mdata4,
                         block_cname = "BatchNum",
                         qc_cname = "QC",
                         qc_val = "QC.NIST",
                         order_cname = "RunOrderOverall"),
               "The first and last sample run for each batch must be a QC sample")
  # put the QC value back to normal
  mdata4$f_data$QC[1] = "QC.NIST"
  
  # check the dimensions of the datasets ---------------------------------------
  # run the batch correction!
  udn_QC <- bc_qcrlsc(omicsData = mdata4,
                      block_cname = "BatchNum",
                      qc_cname = "QC",
                      qc_val = "QC.NIST",
                      order_cname = "RunOrderOverall")
  
  # Dimension check time -------------------------------------------------------
  # Investigate the e_data, f_data, and e_meta data frames.
  # the data will not have the NIST samples anymore so we need to account for that
  numNIST <- sum(grepl("NIST",(colnames(mdata4$e_data)[-1])))
  # the function also removes bad molecules so we need to account for that too
  numRemoved <- length(attr(udn_QC,"filters")[[1]]$filtered$e_data_remove)

  # check the dimensions make sense
  expect_equal(nrow(mdata4$e_data),
               nrow(udn_QC$e_data) + numRemoved)
  expect_equal(ncol(mdata4$e_data),
               ncol(udn_QC$e_data)+numNIST)
  expect_equal(nrow(mdata4$f_data),
               nrow(udn_QC$f_data)+numNIST)
  expect_equal(ncol(mdata4$f_data),
               ncol(udn_QC$f_data))
  
  # Attribute check time -------------------------------------------------------
  og_group_without_QC = attr(mdata4,"group_DF")
  expect_identical(data.frame(og_group_without_QC),
                   data.frame(attr(udn_QC,"group_DF")))
  expect_identical(attr(mdata4, 'cnames'),
                   attr(udn_QC, 'cnames'))
  # all information other than batch info should stay the same
  expect_identical(attributes(mdata4)$data_info[1:2],
                   attributes(udn_QC)$data_info[1:2])
  expect_equal(attributes(udn_QC)$data_info$norm_info$is_norm,TRUE)
  expect_equal(attr(udn_QC,"data_info")$num_edata, nrow(udn_QC$e_data))
  expect_equal(attr(udn_QC,"data_info")$num_miss_obs, sum(is.na(udn_QC$e_data)))
  expect_equal(attr(udn_QC,"data_info")$prop_missing, sum(is.na(udn_QC$e_data))/(nrow(udn_QC$e_data)*(ncol(udn_QC$e_data)-1)))
  expect_equal(attr(udn_QC,"data_info")$num_samps, nrow(udn_QC$f_data))
  
  # batch info shold be updated
  expect_identical(attributes(udn_QC)$data_info$batch_info,
                   list(is_bc = TRUE,bc_method = "bc_qcrlsc", params = list(block_cname = "BatchNum",
                                                                                 qc_cname = "QC",
                                                                                 qc_val = "QC.NIST",
                                                                                 order_cname = "RunOrderOverall",
                                                                                 missing_thresh = 0.5,
                                                                                 rsd_thresh = 0.3,
                                                                                 backtransform = FALSE,
                                                                                 keep_qc = FALSE)))
  expect_identical(attr(mdata4, 'meta_info')$meta_data,
                   attr(udn_QC, 'meta_info')$meta_data) 
  expect_identical(attr(mdata4, 'meta_info')$num_emeta,
                   attr(udn_QC, 'meta_info')$num_emeta + numRemoved)
  # we should have filters now unlike at the beginning
  expect_false(is.null(attr(udn_QC,"filters")))
  
  # do the same thing but with emeta data too ----------------------------------
  # create an emeta column to ensure that emeta also runs as it should
  emet <- data.frame(Metabolite = mdata4$e_data$Metabolite,
                     Ref_ID = seq(1:nrow(mdata4$e_data)))
  
  udn_with_emet <- pmartR::as.metabData(e_data = mdata4$e_data,
                                        f_data = mdata4$f_data,
                                        e_meta = emet,
                                        edata_cname = "Metabolite",
                                        fdata_cname = "SampleID",
                                        emeta_cname = "Metabolite")
  attr(udn_with_emet,"data_info") <- attr(mdata4,"data_info")
  
  # run the batch correction!
  udn_QC2 <- bc_qcrlsc(omicsData = udn_with_emet,
                      block_cname = "BatchNum",
                      qc_cname = "QC",
                      qc_val = "QC.NIST",
                      order_cname = "RunOrderOverall")
  
  # Dimension check time (with emeta) ------------------------------------------
  # Investigate the e_data, f_data, and e_meta data frames.
  # the data will not have the NIST samples anymore so we need to account for that
  numNIST <- sum(grepl("NIST",(colnames(udn_with_emet$e_data)[-1])))
  # the function also removes bad molecules so we need to account for that too
  numRemoved <- length(attr(udn_QC2,"filters")[[1]]$filtered$e_data_remove)
  
  # check the dimensions make sense
  expect_equal(nrow(udn_with_emet$e_data),
               nrow(udn_QC2$e_data) + numRemoved)
  expect_equal(ncol(udn_with_emet$e_data),
               ncol(udn_QC2$e_data)+numNIST)
  expect_equal(nrow(udn_with_emet$f_data),
               nrow(udn_QC2$f_data)+numNIST)
  expect_equal(ncol(udn_with_emet$f_data),
               ncol(udn_QC2$f_data))
  
  # Attribute check time (with emeta) ------------------------------------------
  # inspect attributes of bc_qcrlsc data frame
  og_group_without_QC = attr(udn_with_emet,"group_DF")
  expect_identical(data.frame(og_group_without_QC),
                   data.frame(attr(udn_QC2,"group_DF")))
  expect_identical(attr(udn_with_emet, 'cnames'),
                   attr(udn_QC2, 'cnames'))
  # look at data_info
  expect_identical(attributes(mdata4)$data_info[1:2],
                   attributes(udn_QC2)$data_info[1:2])
  expect_equal(attributes(udn_QC2)$data_info$norm_info$is_norm,TRUE)
  expect_equal(attr(udn_QC2,"data_info")$num_edata, nrow(udn_QC2$e_data))
  expect_equal(attr(udn_QC2,"data_info")$num_miss_obs, sum(is.na(udn_QC2$e_data)))
  expect_equal(attr(udn_QC2,"data_info")$prop_missing, sum(is.na(udn_QC2$e_data))/(nrow(udn_QC2$e_data)*(ncol(udn_QC2$e_data)-1)))
  expect_equal(attr(udn_QC2,"data_info")$num_samps, nrow(udn_QC2$f_data))
  # batch info shold be updated
  expect_identical(attributes(udn_QC2)$data_info$batch_info,
                   list(is_bc = TRUE,bc_method = "bc_qcrlsc", params = list(block_cname = "BatchNum",
                                                                            qc_cname = "QC",
                                                                            qc_val = "QC.NIST",
                                                                            order_cname = "RunOrderOverall",
                                                                            missing_thresh = 0.5,
                                                                            rsd_thresh = 0.3,
                                                                            backtransform = FALSE,
                                                                            keep_qc = FALSE)))
  expect_identical(attr(mdata4, 'meta_info')$meta_data,
                   attr(udn_QC2, 'meta_info')$meta_data) 
  expect_identical(attr(mdata4, 'meta_info')$num_emeta,
                   attr(udn_QC2, 'meta_info')$num_emeta + numRemoved)
  # we should have filters now unlike at the beginning
  expect_false(is.null(attr(udn_QC2,"filters")))
})

