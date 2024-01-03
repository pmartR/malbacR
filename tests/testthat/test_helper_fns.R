context('test helper functions')

test_that('helper functions pull from correct attributes',{
  
  load(system.file('testdata',
                   'mini_udn.Rdata',
                   package = 'malbacR'))
  
  # Construct a pepData object with the edata, fdata, and emeta data frames.
  mdata <- pmartR::as.metabData(e_data = edata,
                                f_data = fdata,
                                e_meta = emeta,
                                edata_cname = 'Metabolite',
                                fdata_cname = 'SampleID',
                                emeta_cname = 'Metabolite',
                                data_scale = "log2")
  
  expect_identical(get_batch_method(mdata), attributes(mdata)$data_info$batch_info$bc_method)
  expect_identical(get_batch_parameters(mdata), attributes(mdata)$data_info$batch_info$params)
  expect_identical(get_data_batch(mdata),attributes(mdata)$data_info$batch_info$is_bc)
  
  # attempt some batch correction methods
  mdata <- pmartR::group_designation(mdata,main_effects = "Age",batch_id = "BatchNum")
  molfilt <- pmartR::molecule_filter(mdata,use_groups=TRUE)
  mdata <- pmartR::applyFilt(molfilt,mdata)
  impObj <- imputation(mdata)
  mdata <- apply_imputation(impObj,mdata)
  mdata_abundance <- pmartR::edata_transform(mdata,"abundance")
  
  # run batch correction
  mdata_nomis <- bc_nomis(mdata_abundance,"IS","IS")
  mdata_eigenMS <- bc_eigenMS(mdata)
  mdata_pareto <- bc_pareto(mdata)
  mdata_qcrfsc <- bc_qcrfsc(mdata_abundance,"QC","QC.NIST","RunOrderOverall","Age")
  # combat requires data to be normalized
  mdata_norm <- pmartR::normalize_global(mdata,"all","median",apply_norm = TRUE,backtransform = TRUE)
  mdata_qcrfsc <- bc_combat(mdata_norm)

  # check bc method of batch corrected data
  expect_identical(get_batch_method(mdata_nomis), attributes(mdata_nomis)$data_info$batch_info$bc_method)
  expect_identical(get_batch_method(mdata_qcrfsc), attributes(mdata_qcrfsc)$data_info$batch_info$bc_method)
  expect_identical(get_batch_method(mdata_qcrfsc), attributes(mdata_qcrfsc)$data_info$batch_info$bc_method)
  expect_identical(get_batch_method(mdata_eigenMS), attributes(mdata_eigenMS)$data_info$batch_info$bc_method)
  expect_identical(get_batch_method(mdata_pareto), attributes(mdata_pareto)$data_info$batch_info$bc_method)
  
  # check parameters of batch corrected data
  expect_identical(get_batch_parameters(mdata_nomis), attributes(mdata_nomis)$data_info$batch_info$params)
  expect_identical(get_batch_parameters(mdata_qcrfsc), attributes(mdata_qcrfsc)$data_info$batch_info$params)
  expect_identical(get_batch_parameters(mdata_qcrfsc), attributes(mdata_qcrfsc)$data_info$batch_info$params)
  expect_identical(get_batch_parameters(mdata_eigenMS), attributes(mdata_eigenMS)$data_info$batch_info$params)
  expect_identical(get_batch_parameters(mdata_pareto), attributes(mdata_pareto)$data_info$batch_info$params)
  
  # check true/false if batch correction has been done
  # check parameters of batch corrected data
  expect_identical(get_data_batch(mdata_nomis), attributes(mdata_nomis)$data_info$batch_info$is_bc)
  expect_identical(get_data_batch(mdata_qcrfsc), attributes(mdata_qcrfsc)$data_info$batch_info$is_bc)
  expect_identical(get_data_batch(mdata_qcrfsc), attributes(mdata_qcrfsc)$data_info$batch_info$is_bc)
  expect_identical(get_data_batch(mdata_eigenMS), attributes(mdata_eigenMS)$data_info$batch_info$is_bc)
  expect_identical(get_data_batch(mdata_pareto), attributes(mdata_pareto)$data_info$batch_info$is_bc)
})
  