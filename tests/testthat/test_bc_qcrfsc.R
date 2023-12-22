context('Run QC-RFSC batch correction')

test_that('bc_qcrfsc returns the correct data frame and attributes', {
  
  # Load the reduced metabolite data frames ------------------------------------
  load(system.file('testdata',
                   'mini_udn_four.Rdata',
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
  # must be run on abundance scale
  # must be a character value
  expect_error(bc_qcrfsc(omicsData = mdata,
                         qc_cname = "QC",
                         qc_val = "QC.NIST",
                         group_cname = "Age",
                         order_cname = "RunOrderOverall"),
               "For QCRFSC, omicsData must be ran with the scale 'abundance'")
  # so we fix this for future process
  mdata <- pmartR::edata_transform(mdata,"abundance")
  
  # must undergo group designation
  expect_error(bc_qcrfsc(omicsData = mdata,
                         qc_cname = "QC",
                         qc_val = "QC.NIST",
                         group_cname = "Age",
                         order_cname = "RunOrderOverall"),
               "omicsData must have undergone group designation")
  # so add in group designation
  mdata <- pmartR::group_designation(mdata,main_effects = "Age")
  
  # QCCNAME
  # must be a character value
  expect_error(bc_qcrfsc(omicsData = mdata,
                         qc_cname = mdata$e_data,
                         qc_val = "QC.NIST",
                         group_cname = "Age",
                         order_cname = "RunOrderOverall"),
               "Input parameter qc_cname must be of class 'character'.")
  # must be of length 1
  expect_error(bc_qcrfsc(omicsData = mdata,
                         qc_cname = c("QC_data","QC"),
                         qc_val = "QC.NIST",
                         group_cname = "Age",
                         order_cname = "RunOrderOverall"),
               "Input parameter qc_cname must be of length 1")
  # qc_cname needs to be a column in the fdata
  expect_error(bc_qcrfsc(omicsData = mdata,
                         qc_cname = "QC_data",
                         qc_val = "QC.NIST",
                         group_cname = "Age",
                         order_cname = "RunOrderOverall"),
               "Input parameter qc_cname must be a column found in f_data of omicsData.")
  # QCVAL
  # must be a character
  expect_error(bc_qcrfsc(omicsData = mdata,
                         qc_cname = "QC",
                         qc_val = mdata$e_data,
                         group_cname = "Age",
                         order_cname = "RunOrderOverall"),
               "Input parameter qc_val must be of class 'character'")
  # must be a value in qc_cname
  expect_error(bc_qcrfsc(omicsData = mdata,
                         qc_cname = "QC",
                         qc_val = "QC",
                         group_cname = "Age",
                         order_cname = "RunOrderOverall"),
               "Input parameter qc_val must be a value in qc_cname column")
  # keep_qc must be of length 1
  expect_error(bc_qcrfsc(omicsData = mdata,
                         qc_cname = "QC",
                         qc_val = "QC.NIST",
                         group_cname = "Age",
                         order_cname = "RunOrderOverall",
                         keep_qc = c(TRUE,FALSE)),
               "Input parameter qc_val must be of length 1")
  # keep_qc must be logical
  expect_error(bc_qcrfsc(omicsData = mdata,
                         qc_cname = "QC",
                         qc_val = "QC.NIST",
                         group_cname = "Age",
                         order_cname = "RunOrderOverall",
                         keep_qc = "True"),
               "Input parameter keep_qc must be logical")
  # ORDERCNAME
  # must be a character
  expect_error(bc_qcrfsc(omicsData = mdata,
                         qc_cname = "QC",
                         qc_val = "QC.NIST",
                         group_cname = "Age",
                         order_cname = mdata$e_data),
               "Input parameter order_cname must be of class 'character'")
  # must be length 1
  expect_error(bc_qcrfsc(omicsData = mdata,
                         qc_cname = "QC",
                         qc_val = "QC.NIST",
                         group_cname = "Age",
                         order_cname = c("RunOrder","RunTime")),
               "Input parameter order_cname must be of length 1")
  # order_cname must be a column in the fdata
  expect_error(bc_qcrfsc(omicsData = mdata,
                         qc_cname = "QC",
                         qc_val = "QC.NIST",
                         group_cname = "Age",
                         order_cname = "RunOrder"),
               "Input parameter order_cname must be a column found in f_data of omicsData.")
  # order_cname must be numeric
  expect_error(bc_qcrfsc(omicsData = mdata,
                         qc_cname = "QC",
                         qc_val = "QC.NIST",
                         order_cname = "BatchName"),
               "Input parameter order_cname must contain integer values for run order")
  
  # GROUP_CNAME
  # must be a character
  expect_error(bc_qcrfsc(omicsData = mdata,
                         qc_cname = "QC",
                         qc_val = "QC.NIST",
                         group_cname = mdata$f_data,
                         order_cname = "RunOrderOverall"),
               "Input parameter group_cname must be of class 'character'")
  # must be length 1
  expect_error(bc_qcrfsc(omicsData = mdata,
                         qc_cname = "QC",
                         qc_val = "QC.NIST",
                         group_cname = c("Age","Sex","Health"),
                         order_cname = "RunOrderOverall"),
               "Input parameter group_cname must be of length 1 or 2")
  # order_cname must be a column in the fdata
  expect_error(bc_qcrfsc(omicsData = mdata,
                         qc_cname = "QC",
                         qc_val = "QC.NIST",
                         group_cname = "Health",
                         order_cname = "RunOrderOverall"),
               "Input parameter group_cname must be a column found in f_data of omicsData.")
  # we cannot have duplicate run orders
  mdata2 <- mdata
  mdata2$f_data$RunOrderOverall[3] <- as.integer(6)
  expect_error(bc_qcrfsc(omicsData = mdata2,
                         qc_cname = "QC",
                         qc_val = "QC.NIST",
                         group_cname = "Age",
                         order_cname = "RunOrderOverall"),
               "Input parameter order_cname cannot contain duplicate values for the overall run order")
  # error out if we have missing values
  expect_error(bc_qcrfsc(omicsData = mdata,
                         qc_cname = "QC",
                         qc_val = "QC.NIST",
                         group_cname = "Age",
                         order_cname = "RunOrderOverall"),
               "QCRFSC requires no missing observations. Remove molecules with missing samples.")
  
  # apply imputation to remove missing values
  molfilt <- pmartR::molecule_filter(mdata,use_groups = TRUE)
  mdata <- pmartR::applyFilt(molfilt,mdata)
  mdata_log <- pmartR::edata_transform(mdata,"log2")
  impObj <- imputation(mdata_log)
  mdata_imp <- apply_imputation(impObj,mdata_log)
  mdata_imp <- pmartR::edata_transform(mdata_imp,"abundance")
  
  # error out if we do not have a QC at the beginning and end of each batch
  # update some QC values so that we error out
  mdata_imp$f_data$QC[1] = "NotQC"
  expect_error(bc_qcrfsc(omicsData = mdata_imp,
                         qc_cname = "QC",
                         qc_val = "QC.NIST",
                         group_cname = "Age",
                         order_cname = "RunOrderOverall"),
               "QCRFSC requires that the first and last sample ran are QC samples")
  # put the QC value back to normal
  mdata_imp$f_data$QC[1] = "QC.NIST"
  
  # check the dimensions of the datasets ---------------------------------------
  # run the batch correction!
  udn_QC <- bc_qcrfsc(omicsData = mdata_imp,
                      qc_cname = "QC",
                      qc_val = "QC.NIST",
                      group_cname = "Age",
                      order_cname = "RunOrderOverall")
  
  # Dimension check time -------------------------------------------------------
  # Investigate the e_data, f_data, and e_meta data frames.
  # the data will not have the NIST samples anymore so we need to account for that
  numNIST <- sum(grepl("NIST",(colnames(mdata$e_data)[-1])))

  # check the dimensions make sense
  expect_equal(nrow(mdata_imp$e_data),
               nrow(udn_QC$e_data))
  expect_equal(ncol(mdata_imp$e_data),
               ncol(udn_QC$e_data)+numNIST)
  expect_equal(nrow(mdata_imp$f_data),
               nrow(udn_QC$f_data)+numNIST)
  expect_equal(ncol(mdata_imp$f_data),
               ncol(udn_QC$f_data))
  
  # Attribute check time -------------------------------------------------------
  # inspect attributes of bc_qcrfsc data frame
  og_group_without_QC = attr(mdata,"group_DF") %>% dplyr::filter(Group != "QC.NIST")
  expect_identical(data.frame(og_group_without_QC),
               data.frame(attr(udn_QC,"group_DF")))
  expect_identical(attr(mdata, 'cnames'),
                   attr(udn_QC, 'cnames'))
  # all information other than batch info should stay the same
  expect_identical(attributes(mdata)$data_info[1:2],
                   attributes(udn_QC)$data_info[1:2])
  expect_equal(attributes(udn_QC)$data_info$norm_info$is_norm,TRUE)
  expect_equal(attr(udn_QC,"data_info")$num_edata, nrow(udn_QC$e_data))
  expect_equal(attr(udn_QC,"data_info")$num_miss_obs, sum(is.na(udn_QC$e_data)))
  expect_equal(attr(udn_QC,"data_info")$prop_missing, sum(is.na(udn_QC$e_data))/(nrow(udn_QC$e_data)*(ncol(udn_QC$e_data)-1)))
  expect_equal(attr(udn_QC,"data_info")$num_samps, nrow(udn_QC$f_data))
  
  # batch info shold be updated
  expect_identical(attributes(udn_QC)$data_info$batch_info,
                   list(is_bc = TRUE,bc_method = "bc_qcrfsc", params = list(qc_cname = "QC",
                                                                            qc_val = "QC.NIST",
                                                                            order_cname = "RunOrderOverall",
                                                                            group_cname = "Age",
                                                                            ntree = 500,
                                                                            keep_qc = FALSE)))
  expect_identical(attr(mdata, 'meta_info')$meta_data,
                   attr(udn_QC, 'meta_info')$meta_data) 
  expect_identical(attr(mdata, 'meta_info')$num_emeta,
                   attr(udn_QC, 'meta_info')$num_emeta)
  # we should have filters now unlike at the beginning
  expect_false(is.null(attr(udn_QC,"filters")))
  
  # do the same thing but with emeta data too ----------------------------------
  # create an emeta column to ensure that emeta also runs as it should
  emet <- data.frame(Metabolite = mdata$e_data$Metabolite,
                     Ref_ID = seq(1:nrow(mdata$e_data)))
  
  udn_with_emet <- pmartR::as.metabData(e_data = mdata$e_data,
                                        f_data = mdata$f_data,
                                        e_meta = emet,
                                        edata_cname = "Metabolite",
                                        fdata_cname = "SampleID",
                                        emeta_cname = "Metabolite")
  attr(udn_with_emet,"data_info") <- attr(mdata,"data_info")
  
  # apply imputation to remove missing values
  udn_with_emet <- group_designation(udn_with_emet,main_effects = "Age")
  molfilt <- molecule_filter(udn_with_emet,use_groups = TRUE)
  udn_with_emet <- applyFilt(molfilt,udn_with_emet)
  udn_with_emet_log <- pmartR::edata_transform(udn_with_emet,"log2")
  impObj <- imputation(udn_with_emet_log)
  udn_with_emet <- apply_imputation(impObj,udn_with_emet_log)
  udn_with_emet <- pmartR::edata_transform(udn_with_emet,"abundance")
  
  # run the batch correction!
  udn_QC2 <- bc_qcrfsc(omicsData = udn_with_emet,
                      qc_cname = "QC",
                      qc_val = "QC.NIST",
                      group_cname = "Age",
                      order_cname = "RunOrderOverall")
  
  # Dimension check time (with emeta) ------------------------------------------
  # Investigate the e_data, f_data, and e_meta data frames.
  # the data will not have the NIST samples anymore so we need to account for that
  numNIST <- sum(grepl("NIST",(colnames(udn_with_emet$e_data)[-1])))

  # check the dimensions make sense
  expect_equal(nrow(udn_with_emet$e_data),
               nrow(udn_QC2$e_data))
  expect_equal(ncol(udn_with_emet$e_data),
               ncol(udn_QC2$e_data)+numNIST)
  expect_equal(nrow(udn_with_emet$f_data),
               nrow(udn_QC2$f_data)+numNIST)
  expect_equal(ncol(udn_with_emet$f_data),
               ncol(udn_QC2$f_data))
  
  # Attribute check time (with emeta) ------------------------------------------
  # inspect attributes of bc_qcrfsc data frame
  og_group_without_QC = attr(udn_with_emet,"group_DF") %>% dplyr::filter(Group != "QC.NIST")
  expect_identical(data.frame(og_group_without_QC),
                   data.frame(attr(udn_QC2,"group_DF")))
  expect_identical(attr(udn_with_emet, 'cnames'),
                   attr(udn_QC2, 'cnames'))
  # look at data_info
  expect_identical(attributes(mdata)$data_info[1:2],
                   attributes(udn_QC2)$data_info[1:2])
  expect_equal(attributes(udn_QC2)$data_info$norm_info$is_norm,TRUE)
  expect_equal(attr(udn_QC2,"data_info")$num_edata, nrow(udn_QC2$e_data))
  expect_equal(attr(udn_QC2,"data_info")$num_miss_obs, sum(is.na(udn_QC2$e_data)))
  expect_equal(attr(udn_QC2,"data_info")$prop_missing, sum(is.na(udn_QC2$e_data))/(nrow(udn_QC2$e_data)*(ncol(udn_QC2$e_data)-1)))
  expect_equal(attr(udn_QC2,"data_info")$num_samps, nrow(udn_QC2$f_data))
  # batch info shold be updated
  expect_identical(attributes(udn_QC2)$data_info$batch_info,
                   list(is_bc = TRUE,bc_method = "bc_qcrfsc", params = list(qc_cname = "QC",
                                                                            qc_val = "QC.NIST",
                                                                            order_cname = "RunOrderOverall",
                                                                            group_cname = "Age",
                                                                            ntree = 500,
                                                                            keep_qc = FALSE)))
  expect_identical(attr(mdata, 'meta_info')$meta_data,
                   attr(udn_QC2, 'meta_info')$meta_data) 
  expect_identical(attr(mdata, 'meta_info')$num_emeta,
                   attr(udn_QC2, 'meta_info')$num_emeta)
  # we should have filters now unlike at the beginning
  expect_false(is.null(attr(udn_QC2,"filters")))
})

