context('Run scale_method to obtain scaled version of the data matrix')

test_that('scale_method scales the data based on input supplied',{
  
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
  attributes(mdata)$data_info$data_scale_orig <- "log2"
  attributes(mdata)$data_info$data_scale <- "log2"
  
  # remove molecules with too much missingness
  molfilt <- pmartR::molecule_filter(mdata)
  mdata2 <- pmartR::applyFilt(molfilt,mdata,min_num = 2)

  # Run through the different scaling options ----------------------------------
  
  # simplify the name of the data frame
  dat <- mdata2$e_data[,-1]

  #### AUTO --------------------------------------------------------------------
  # dance R method
  auto_dance <- malbacR:::scale_method(dat,"auto")
  # manual method
  # baseline row
  mean_info <- mean(as.numeric(dat[,1]),na.rm = TRUE)
  sd_info <- sd(as.numeric(dat[,1]),na.rm = TRUE)
  auto_manual <- (dat[,1] - mean_info)/sd_info
  # for loop for the rest of the loops
  for(i in 2:ncol(dat)){
    mean_info <- mean(as.numeric(dat[,i]),na.rm = TRUE)
    sd_info <- sd(as.numeric(dat[,i]),na.rm = TRUE)
    new_row <- (dat[,i] - mean_info)/sd_info
    auto_manual <- cbind(auto_manual,new_row)
  }
  rownames(auto_dance) <- NULL
  colnames(auto_dance) <- NULL
  rownames(auto_manual) <- NULL
  colnames(auto_manual) <- NULL
  
  # run the test
  expect_equal(as.matrix(auto_dance),as.matrix(auto_manual))
  
  #### RANGE -------------------------------------------------------------------
  # dance R method
  range_dance <- malbacR:::scale_method(dat,"range")

  # manual method
  # baseline row
  mean_info <- mean(as.numeric(dat[,1]),na.rm = TRUE)
  range_info <- range(as.numeric(dat[,1]),na.rm = TRUE)
  range_manual <- (dat[,1] - mean_info)/(range_info[2] - range_info[1])
  # for loop for the rest of the loops
  for(i in 2:ncol(dat)){
    mean_info <- mean(as.numeric(dat[,i]),na.rm = TRUE)
    range_info <- range(as.numeric(dat[,i]),na.rm = TRUE)
    new_row <- (dat[,i] - mean_info)/(range_info[2] - range_info[1])
    range_manual <- cbind(range_manual,new_row)
  }
  rownames(range_dance) <- NULL
  colnames(range_dance) <- NULL
  rownames(range_manual) <- NULL
  colnames(range_manual) <- NULL
  # run the test
  expect_equal(as.matrix(range_dance),as.matrix(range_manual))
  
  #### PARETO ------------------------------------------------------------------
  # dance R method
  pareto_dance <- malbacR:::scale_method(dat,"pareto")
                                                         
  # manual method
  # baseline row
  mean_info <- mean(as.numeric(dat[,1]),na.rm = TRUE)
  sd_info <- sd(as.numeric(dat[,1]),na.rm = TRUE)
  pareto_manual <- (dat[,1] - mean_info)/sqrt(sd_info)
  # for loop for the rest of the loops
  for(i in 2:ncol(dat)){
    mean_info <- mean(as.numeric(dat[,i]),na.rm = TRUE)
    sd_info <- sd(as.numeric(dat[,i]),na.rm = TRUE)
    new_row <- (dat[,i] - mean_info)/sqrt(sd_info)
    pareto_manual <- cbind(pareto_manual,new_row)
  }
  rownames(pareto_dance) <- NULL
  colnames(pareto_dance) <- NULL
  rownames(pareto_manual) <- NULL
  colnames(pareto_manual) <- NULL
  # run the test
  expect_equal(as.matrix(pareto_dance),as.matrix(pareto_manual))
  
  #### VAST --------------------------------------------------------------------
  # dance R method
  vast_dance <- malbacR:::scale_method(dat,"vast")
  
  # manual method
  # baseline row
  mean_info <- mean(as.numeric(dat[,1]),na.rm = TRUE)
  sd_info <- sd(as.numeric(dat[,1]),na.rm = TRUE)
  vast_manual <- mean_info * (dat[,1] - mean_info)/(sd_info^2)
  # for loop for the rest of the loops
  for(i in 2:ncol(dat)){
    mean_info <- mean(as.numeric(dat[,i]),na.rm = TRUE)
    sd_info <- sd(as.numeric(dat[,i]),na.rm = TRUE)
    new_row <- mean_info * (dat[,i] - mean_info)/(sd_info^2)
    vast_manual <- cbind(vast_manual,new_row)
  }
  rownames(vast_dance) <- NULL
  colnames(vast_dance) <- NULL
  rownames(vast_manual) <- NULL
  colnames(vast_manual) <- NULL
  # run the test
  expect_equal(as.matrix(vast_dance),as.matrix(vast_manual))
  
  #### LEVEL -------------------------------------------------------------------
  # dance R method
  level_dance <- malbacR:::scale_method(dat,"level")

  # manual method
  # baseline row
  mean_info <- mean(as.numeric(dat[,1]),na.rm = TRUE)
  level_manual <- (dat[,1] - mean_info)/mean_info
  # for loop for the rest of the loops
  for(i in 2:ncol(dat)){
    mean_info <- mean(as.numeric(dat[,i]),na.rm = TRUE)
    new_row <- (dat[,i] - mean_info)/mean_info
    level_manual <- cbind(level_manual,new_row)
  }
  rownames(level_dance) <- NULL
  colnames(level_dance) <- NULL
  rownames(level_manual) <- NULL
  colnames(level_manual) <- NULL
  # run the test
  expect_equal(as.matrix(level_dance),as.matrix(level_manual))
  
  #### POWER -------------------------------------------------------------------
  # dance R method
  power_dance <- malbacR:::scale_method(dat,"power")
  
  # manual method
  # baseline row
  vals = sqrt(dat[,1])
  power_manual <- sqrt(dat[,1]) - mean(as.numeric(vals),na.rm=T)
  # for loop for the rest of the loops
  for(i in 2:ncol(dat)){
    vals = sqrt(dat[,i])
    new_row <- sqrt(dat[,i]) - mean(as.numeric(vals),na.rm=T)
    power_manual <- cbind(power_manual,new_row)
  }
  rownames(power_dance) <- NULL
  colnames(power_dance) <- NULL
  rownames(power_manual) <- NULL
  colnames(power_manual) <- NULL
  # run the test
  expect_equal(as.matrix(power_dance),as.matrix(power_manual))
})
