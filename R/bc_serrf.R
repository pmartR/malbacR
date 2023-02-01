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
#' \dontrun{
#' library(malbacR)
#' data("pmart_amide")
#' pmart_amide <- edata_transform(pmart_amide,"log2")
#' pmart_amide <- group_designation(pmart_amide,main_effects = "group",batch_id = "batch")
#' pmart_amide <- normalize_global(pmart_amide,subset_fn = "all",norm_fn = "median",
#'                                apply_norm = TRUE,backtransform = TRUE)
#' amide_combat <- bc_combat(omicsData = pmart_amide, use_groups = FALSE)
#' impObj <- imputation(omicsData = pmart_amide)
#' amide_imp <- apply_imputation(imputeData = impObj, omicsData = pmart_amide)
#' amide_serrf <- bc_serrf(omicsData = amide_imp,sampletype_cname = "group",test_val = "QC")
#' 
#' }
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
  
  # run the SERRF calculations -------------------------------------------------
  # set up row names to be the edata
  rownames(omicsData$e_data) <- omicsData$e_data[,edata_cnameCol]
  edata <- omicsData$e_data[,-edata_cnameCol]
  fdata <- omicsData$f_data
  
  # separate the edata based on QC or not QC
  edataQC <- edata[,fdata[,sampletype_cname] == test_val]
  edata_noQC <- edata[,fdata[,sampletype_cname] != test_val]
  
  # now add in fdata information to the data for QC and non QC samples
  edataQC_t <- edataQC %>%
    t() %>%
    data.frame() %>%
    tibble::rownames_to_column(var = fdata_cname)
  edata_noQC_t <- edata_noQC %>%
    t() %>%
    data.frame() %>%
    tibble::rownames_to_column(var = fdata_cname)
  training_dat <- fdata %>%
    dplyr::right_join(edataQC_t, by = fdata_cname)
  testing_dat <- fdata %>%
    dplyr::right_join(edata_noQC_t, by = fdata_cname)
  
  # obtain the batch information
  batch_info <- attributes(attr(omicsData,"group_DF"))$batch_id
  batchName <- colnames(batch_info)[2]
  
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
        #i = 1
        # identify the current variable
        current_var = colnames(norm_dat)[i]
        # find which columns correspond to that variable
        trn_current_var_col = which(colnames(trn) == current_var)
        tst_current_var_col = which(colnames(tst) == current_var)
        # find the y response that molecule
        train_data_y = scale(trn[,trn_current_var_col],scale = F)
        # scale the QC training data
        train_data_x = t(apply(trn[,-trn_current_var_col],1,scale))
        train_data = data.frame(train_data_y,train_data_x )
        colnames(train_data) = c("y",colnames(norm_dat)[-trn_current_var_col])
        # scale the testing data
        test_data = t(apply(tst,1,scale))
        test_data = data.frame(test_data)
        # ensure that column names are the molecules
        colnames(test_data) = colnames(norm_dat)
        colnames(test_data)[tst_current_var_col] <- "y"
        
        # run random forest
        model = ranger::ranger(y~.,data= train_data)
        
        # current values for this batch and molecule
        ii = as.numeric(unlist(tst[,trn_current_var_col]))
        
        # systematic error
        si = (predict(model,data = test_data)$predictions  + mean(ii,na.rm=TRUE))/(median(all_nonQC,na.rm = TRUE))
        norm_val = ii/si
        
        # multiply by iibar
        iibar = median(all_nonQC,na.rm=TRUE)/median(norm_val,na.rm = TRUE)
        norm_val = norm_val*iibar
        
        # add to the data frame of normalized values
        norm_dat[,trn_current_var_col] <- norm_val
      }
      
      # now make it a datframe and add back in rownname information
      norm_dat <- data.frame(nd)
      colnames(norm_dat) = colnames(norm_dat)
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
  adjustRank <- rank(match(colnames(edata_serrf), colnames(omicsData$e_data)))
  edata_serrf <- edata_serrf[,adjustRank]
  edata_serrf[,edata_cname] <- omicsData$e_data[,edata_cname]

  
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
