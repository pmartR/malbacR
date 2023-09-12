#' SERRF batch correction
#' 
#' This function returns a pmart object that has been undergone SERRF
#'  batch effect correction
#' 
#' @param omicsData an object of the class 'pepData', 'proData', 'metabData',
#'   'lipidData', or 'nmrData', usually created by \code{\link{as.pepData}},
#'   \code{\link{as.proData}}, \code{\link{as.metabData}},
#'   \code{\link{as.lipidData}}, or \code{\link{as.nmrData}}, respectively.
#' @param sampletype_cname character string giving the name of the column in omicsData$f_data
#'  that contains the sample type information (such as quality control samples)
#' @param test_val character string giving the name of the value within the column sampletype_cname to be used
#'  as the testing value for SERRF
#'  
#' @return Object of same class as omicsData that has been undergone
#'   SERRF normalization
#' 
#' @examples
#' library(malbacR)
#' library(pmartR)
#' data("pmart_amide")
#' pmart_amide <- edata_transform(pmart_amide,"log2")
#' pmart_amide <- group_designation(pmart_amide,main_effects = "group",batch_id = "batch")
#' pmart_amide <- normalize_global(pmart_amide,subset_fn = "all",norm_fn = "median",
#'                                apply_norm = TRUE,backtransform = TRUE)
#' impObj <- imputation(omicsData = pmart_amide)
#' amide_imp <- apply_imputation(imputeData = impObj, omicsData = pmart_amide)
#' amide_serrf <- bc_serrf(omicsData = amide_imp,sampletype_cname = "group",test_val = "QC")
#' 
#' @author Damon Leach
#' 
#' @export
#' 
bc_serrf <- function(omicsData, sampletype_cname, test_val){
  # run through checks ---------------------------------------------------------
  
  # check that omicsData is of appropriate class #
  if (!inherits(omicsData, c("pepData", "proData", "metabData", "lipidData",
                             "nmrData"))) {
    
    stop (paste("omicsData must be of class 'pepData', 'proData', 'metabData',",
                "'lipidData', or 'nmrData'",
                sep = ' '))
  }
  
  # check that omicsData has batch id information
  if (is.null(attributes(attr(omicsData,"group_DF"))$batch_id)){
    stop (paste("omicsData must have batch_id attribute for batch correction",
                sep = ' '))
  }
  
  # sampletype_cname - type of each sample
  if (class(sampletype_cname) != "character") {
    stop("Input parameter sampletype_cname must be of class 'character'.")
  }
  
  if (length(sampletype_cname) != 1) {
    stop("Input parameter sampletype_cname must be of length 1 (e.g. vector containing a single element")
  }
  
  if (!any(names(omicsData$f_data) == sampletype_cname)) {
    stop("Input parameter sampletype_cname must be a column found in f_data of omicsData.")
  }
  
  # test_val - the value that will be used to determine what is the testing and what is the training
  if (class(test_val) != "character") {
    stop("Input parameter test_val must be of class 'character'.")
  }
  if (length(test_val) != 1) {
    stop("Input parameter test_val must be of length 1 (e.g. vector containing a single element)")
  }
  if (!test_val %in% omicsData$f_data[,sampletype_cname]){
    stop("Input parameter test_val must be a value in sampletype_cname column in omicsData$f_data")
  }
  
  # check that omicsData has batch id information
  if (is.null(attributes(attr(omicsData,"group_DF"))$batch_id)){
    stop (paste("omicsData must have batch_id attribute for batch correction",
                sep = ' '))
  }
  
  # check that omicsData has batch id information
  if (is.null(attributes(attr(omicsData,"group_DF"))$batch_id)){
    stop (paste("omicsData must have batch_id attribute for batch correction",
                sep = ' '))
  }
  
  # useful information
  edata_cname <- pmartR::get_edata_cname(omicsData)
  fdata_cname <- pmartR::get_fdata_cname(omicsData)
  edata_cnameCol <- which(colnames(omicsData$e_data) == edata_cname)
  fdata_cnameCol <- which(colnames(omicsData$f_data) == fdata_cname)
  
  # SERRF does not like missing values
  if (sum(is.na(omicsData$e_data[,-edata_cnameCol])) != 0) {
    stop ("SERRF requires no missing observations. Remove molecules with missing samples.")
  }
  
  # we cannot have negative values
  if(sum(omicsData$e_data < 0,na.rm=TRUE) > 0){
    stop("SERRF cannot run with expression data that has negative values (likely
         due to normalization without a backtransform")
  }
  
  # run the SERRF calculations -------------------------------------------------
  # set up row names to be the edata
  rownames(omicsData$e_data) <- omicsData$e_data[,edata_cnameCol]
  edata <- omicsData$e_data[,-edata_cnameCol]
  fdata <- omicsData$f_data
  
  
  
  # separate the edata based on QC or not QC
  qc_sampNames = fdata[fdata[,sampletype_cname] == test_val,][[fdata_cnameCol]]
  nonqc_sampNames = fdata[fdata[,sampletype_cname] != test_val,][[fdata_cnameCol]]
  edataQC <- edata %>% dplyr::select(dplyr::all_of(qc_sampNames))
  edata_noQC <- edata %>% dplyr::select(dplyr::all_of(nonqc_sampNames))
  
  # now add in fdata information to the data for QC and non QC samples
  edataQC_t <- edataQC %>%
    t() %>%
    data.frame() %>%
    tibble::rownames_to_column(var = fdata_cname)
  colnames(edataQC_t) <- c(fdata_cname,rownames(edata))
  edata_noQC_t <- edata_noQC %>%
    t() %>%
    data.frame() %>%
    tibble::rownames_to_column(var = fdata_cname)
  colnames(edata_noQC_t) <- c(fdata_cname,rownames(edata))
  training_dat <- fdata %>%
    dplyr::right_join(edataQC_t, by = fdata_cname)
  testing_dat <- fdata %>%
    dplyr::right_join(edata_noQC_t, by = fdata_cname)
  
  # obtain the batch information
  batch_info <- attributes(attr(omicsData,"group_DF"))$batch_id
  batch_sampleID <- which(colnames(batch_info) == pmartR::get_fdata_cname(omicsData))
  batchName <- colnames(batch_info)[-batch_sampleID]
  
  # list the columns of fdata to remove
  badfdata <- colnames(fdata)
  # we want to keep batch number and sample ID in the datasets
  batchNum <- which(badfdata == batchName)
  sampNum <- which(badfdata == fdata_cname)
  badfdata <- badfdata[-c(batchNum,sampNum)]
  
  trainBatchNum <- which(colnames(training_dat) == batchName)
  colnames(training_dat)[trainBatchNum] <- "batch"
  testBatchNum <- which(colnames(testing_dat) == batchName)
  colnames(testing_dat)[testBatchNum] <- "batch"
  
  # set up correlation matrix lists
  corrs_train = list()
  corrs_target = list()
  
  # find the correlation matrices for training and testing data
  for(b in 1:length(unique(batch_info[[batchName]]))){
    #b = 2
    # find current batch
    current_batch = unique(batch_info[[batchName]])[b]
    
    train_one_batch = subset(training_dat,training_dat$batch == current_batch)
    test_one_batch = subset(testing_dat,testing_dat$batch == current_batch)
    
    # scale the training/QC data
    # subset to only training data and transpose the data so we have molecule x sample
    rownames(train_one_batch) <- NULL
    train_one_batch <- tibble::column_to_rownames(train_one_batch,var = fdata_cname)
    train_one_batch <- train_one_batch %>%
      t() %>%
      data.frame()
    train_one_batch <- subset(train_one_batch,rownames(train_one_batch) %in% omicsData$e_data[[edata_cname]])
    train_one_batch <- apply(train_one_batch,2,as.numeric)
    
    # do this for test data too
    rownames(test_one_batch) <- NULL
    test_one_batch <- tibble::column_to_rownames(test_one_batch,var = fdata_cname)
    test_one_batch <- test_one_batch %>%
      t() %>%
      data.frame()
    test_one_batch <- subset(test_one_batch,rownames(test_one_batch) %in% omicsData$e_data[[edata_cname]])
    test_one_batch <- apply(test_one_batch,2,as.numeric)
    
    # scale the data (different for target and training it appears)
    train_scale = t(apply(train_one_batch,1,scale))
    #train_scale = t(scale(train_one_batch))
    test_scale = scale(test_one_batch)
    
    # find the correlation values for both target and training data (scaled)
    # e_current_batch = all_scale
    corrs_train[[b]] = cor(t(train_scale), method = "spearman")
    corrs_target[[b]] = cor(t(test_scale), method = "spearman")
  }
  
  # set up nested data so we can do analyses by batch
  # first we set up training data
  training_nest <- training_dat %>%
    # remove non essential fdata
    dplyr::select(-dplyr::any_of(badfdata)) %>%
    # nest by batch
    dplyr::nest_by(batch) %>%
    # rename to be training data
    dplyr::rename(train_data = data) %>%
    # ungroup to avoid complications
    dplyr::ungroup()
  
  testing_nest <- testing_dat %>%
    # remove non essential data
    dplyr::select(-dplyr::any_of(badfdata)) %>%
    # nest by batch
    dplyr::nest_by(batch) %>%
    # rename to be testing data
    dplyr::rename(test_data = data) %>%
    # ungroup to avoid complications
    dplyr::ungroup()
  
  # left join the training and testing data
  dat_nest <- training_nest %>%
    dplyr::left_join(testing_nest)
  
  # get information regarding all the nonQC samples intensity before scaling
  fdataCols <- which(colnames(testing_dat) %in% colnames(fdata))
  fdataCols = c(fdataCols,which(colnames(testing_dat) == "batch"))
  all_nonQC <- unlist(matrix(testing_dat[,-fdataCols]))
  
  
  dat_nest <- dat_nest %>%
    dplyr::mutate(norm_data = purrr::map2(train_data,test_data,function(trn,tst){
      #trn <- dat_nest$train_data[[1]]
      #tst <- dat_nest$test_data[[1]]
      # keep sample ID information as rownames
      trn <- trn %>%
        tibble::column_to_rownames(var = fdata_cname)
      tst <- tst %>%
        tibble::column_to_rownames(var = fdata_cname)
      
      
      
      # create empty data matrix for new normalization methods
      tbl_colnames <- colnames(tst)
      norm_dat <- matrix(nrow=nrow(tst),ncol=ncol(tst))
      colnames(norm_dat) <- tbl_colnames
      
      doParallel::registerDoParallel(parallel::detectCores()-1)
      # for each molecule normalize the data
      nd <- foreach::foreach(i = 1:ncol(norm_dat),.combine = cbind) %dopar% {
        # find what variables we are working with num_cor
        num_cor = 10
        for(b in 1:length(unique(batch_info[[batchName]]))){
          #b = 1; i = 1
          # find current batch
          current_batch = unique(batch_info[[batchName]])[b]
          
          # select only the batch we are working with
          corr_train_one_batch = corrs_train[[b]]
          corr_test_one_batch = corrs_target[[b]]
          
          # order the correlation as it pertains to the molecule of interest j
          corr_train_order = order(abs(corr_train_one_batch[,i]),decreasing = TRUE)
          corr_test_order = order(abs(corr_test_one_batch[,i]),decreasing = TRUE)
          
          # initialize sel_var vector and set l to be default number to start with
          sel_var = c()
          l = num_cor
          
          # want to find the num (10) most correlated molecules to our current molecule
          # of interest across both training and testing data
          while(length(sel_var)<(num_cor)){
            # find the intersection of the training and target correlation values
            # until we find 10 that match (or other number if not 10 as default)
            sel_var = intersect(corr_train_order[1:l], corr_test_order[1:l])
            # we cannot use our own molecule
            sel_var = sel_var[!sel_var == i]
            # keep adding one to l until we get to 10 molecules that intersect
            l = l+1
          }
        }
        
        #i = 1
        # identify the current variable
        current_var = colnames(norm_dat)[i]
        # find which columns correspond to that variable
        trn_current_var_col = which(colnames(trn) == current_var)
        tst_current_var_col = which(colnames(tst) == current_var)
        # find the y response that molecule
        train_data_y = scale(trn[,trn_current_var_col],scale = F)
        # scale the QC training data
        train_data_x = apply(trn[,sel_var],2,scale)
        train_data = data.frame(train_data_y,train_data_x )
        og_train_data_colnames = c("y",colnames(norm_dat)[sel_var])
        # scale the testing data
        test_data = apply(tst[,sel_var],2,scale)
        test_data = data.frame(test_data)
        # ensure that column names are the molecules
        og_test_data_colnames = colnames(norm_dat)[sel_var]
        #colnames(test_data)[tst_current_var_col] <- "y"
        
        # run random forest on this dataset
        # r sometimes gets angry at the molecule names so 
        # we simplify them here for the model
        predictor_vars = paste0("var",seq(from = 1, to = length(colnames(norm_dat)[sel_var]),by = 1))
        colnames(train_data) = c("y",predictor_vars)
        colnames(test_data) = predictor_vars
        model = ranger::ranger(y ~ ., data = train_data)
        
        # begin adjusting the normalized values
        # find the important values
        # for QC samples
        pred_bqc = predict(model,data = train_data)$predictions
        abund_val_bqc = trn[,i]
        abund_val_bqc_mean = mean(trn[,i],na.rm=T)
        all_data_var_num_qc = which(colnames(training_dat) == colnames(trn)[i])
        abund_val_qc_median = median(training_dat[,all_data_var_num_qc],na.rm=T)
        # for nonQC samples
        pred_bnq = predict(model, data = test_data)$predictions
        abund_val_bnq = tst[,i]
        abund_val_bnq_mean = mean(tst[,i],na.rm=T)
        all_data_var_num_nq = which(colnames(training_dat) == colnames(trn)[i])
        abund_val_nq_median = median(testing_dat[,all_data_var_num_nq],na.rm=T)
        
        # now calculate serff adjsuted values
        norm_vals_bqc = abund_val_bqc/((pred_bqc + abund_val_bqc_mean)/abund_val_qc_median)
        norm_vals_bnq = abund_val_bnq/((pred_bnq + abund_val_bnq_mean)/abund_val_nq_median)
        serrf_bqc = norm_vals_bqc/(median(norm_vals_bqc,na.rm=T)/abund_val_qc_median)
        serrf_bnq = norm_vals_bnq/(median(norm_vals_bnq,na.rm=T)/abund_val_nq_median)
        
        norm_dat[,trn_current_var_col] <- serrf_bnq
      }
      
      # now make it a datframe and add back in rownname information
      norm_dat <- data.frame(nd)
      colnames(norm_dat) = colnames(tst)
      rownames(norm_dat) <- rownames(tst)
      return(norm_dat)
    }))
  
  # combine all the different batch normalizations together
  all_dat <- do.call("rbind", dat_nest$norm_data)
  # put back into molecule x sample order
  edata_serrf <- all_dat %>%
    t() %>%
    data.frame() %>%
    tibble::rownames_to_column(var = edata_cname)
  
  # put back into proper ordering (after removing QC samples)
  og_edata <- omicsData$e_data %>%
    dplyr::select(-dplyr::all_of(qc_sampNames))
  og_edata_ordering <- colnames(og_edata)
  edata_serrf <- edata_serrf %>% dplyr::select(dplyr::all_of(og_edata_ordering))
  edata_serrf[,edata_cname] <- edata_serrf[order(edata_serrf[,edata_cname],omicsData$e_data[,edata_cname]),][[edata_cname]]
  
  # filter out the old QC samples
  fdata_serrf <- fdata[fdata[,sampletype_cname] != test_val,]
  
  # move back into pmart object ------------------------------------------------
  
  # assume emeta is NULL unless otherwise stated
  emeta_cname = NULL
  emet = NULL
  if(!is.null(omicsData$e_meta)){
    emeta_cname = pmartR::get_emeta_cname(omicsData)
    emet = omicsData$e_meta
  }
  
  # create the pmart object #
  if(inherits(omicsData,"isobaricpepData") & inherits(omicsData,"pepData")){
    pmartObj = pmartR::as.isobaricpepData(e_data = edata_serrf,
                                          edata_cname = edata_cname,
                                          f_data = fdata_serrf,
                                          fdata_cname = fdata_cname,
                                          e_meta = emet,
                                          emeta_cname = emeta_cname)
  }
  else if(inherits(omicsData,"pepData")){
    pmartObj = pmartR::as.pepData(e_data = edata_serrf,
                                  edata_cname = edata_cname,
                                  f_data = fdata_serrf,
                                  fdata_cname = fdata_cname,
                                  e_meta = emet,
                                  emeta_cname = emeta_cname)
  }
  else if(inherits(omicsData,"proData")){
    pmartObj = pmartR::as.proData(e_data = edata_serrf,
                                  edata_cname = edata_cname,
                                  f_data = fdata_serrf,
                                  fdata_cname = fdata_cname,
                                  e_meta = emet,
                                  emeta_cname = emeta_cname)
  }
  else if(inherits(omicsData,"metabData")){
    pmartObj = pmartR::as.metabData(e_data = edata_serrf,
                                    edata_cname = edata_cname,
                                    f_data = fdata_serrf,
                                    fdata_cname = fdata_cname,
                                    e_meta = emet,
                                    emeta_cname = emeta_cname)
  }
  else if(inherits(omicsData,"lipidData")){
    pmartObj = pmartR::as.lipidData(e_data = edata_serrf,
                                    edata_cname = edata_cname,
                                    f_data = fdata_serrf,
                                    fdata_cname = fdata_cname,
                                    e_meta = emet,
                                    emeta_cname = emeta_cname)
  }
  else if(inherits(omicsData,"nmrData")){
    pmartObj = pmartR::as.nmrData(e_data = edata_serrf,
                                  edata_cname = edata_cname,
                                  f_data = fdata_serrf,
                                  fdata_cname = fdata_cname,
                                  e_meta = emet,
                                  emeta_cname = emeta_cname)
  }
  
  # Update the data_info attribute.
  attr(pmartObj, 'data_info') <- pmartR:::set_data_info(
    e_data = pmartObj$e_data,
    edata_cname = pmartR::get_edata_cname(omicsData),
    data_scale_orig = pmartR::get_data_scale_orig(omicsData),
    data_scale = pmartR::get_data_scale(omicsData),
    data_types = pmartR::get_data_info(omicsData)$data_types,
    norm_info = pmartR::get_data_info(omicsData)$norm_info,
    is_normalized = pmartR::get_data_info(omicsData)$norm_info$is_normalized,
    batch_info = pmartR::get_data_info(omicsData)$batch_info,
    is_bc = pmartR::get_data_info(omicsData)$batch_info$is_bc
  )
  
  # Add the group information to the group_DF attribute in the omicsData object.
  attr(pmartObj, "group_DF") = attr(omicsData,"group_DF")
  
  # Update the data_info attribute.
  attributes(pmartObj)$data_info$batch_info <- list(
    is_bc = TRUE,
    bc_method = "serrf",
    params = list()
  )
  
  # Update the meta_info attribute.
  attr(pmartObj, 'meta_info') <- pmartR:::set_meta_info(
    e_meta = omicsData$e_meta,
    emeta_cname = pmartR::get_emeta_cname(omicsData)
  )
  
  # Update the filters attribute (but keep all the filtering)
  attr(pmartObj, 'filters') <- attr(omicsData,'filters')
  
  # return pmart object
  return(pmartObj)
  
}
