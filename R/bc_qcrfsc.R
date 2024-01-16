#' QCRFSC batch correction
#' 
#' This function returns a pmart object that has been undergone QCRFSC
#'  batch effect correction
#' 
#' @param omicsData an object of the class 'pepData', 'proData', 'metabData',
#'   'lipidData', or 'nmrData', usually created by \code{\link{as.pepData}},
#'   \code{\link{as.proData}}, \code{\link{as.metabData}},
#'   \code{\link{as.lipidData}}, or \code{\link{as.nmrData}}, respectively.
#' @param qc_cname character string giving name of column in omicsData$f_data
#'  that contains the factor variable indicating whether sample is QC or not
#' @param qc_val character string giving the value from the qc_cname column that
#'  indicates a QC sample
#' @param order_cname character string giving name of column in omicsData$f_data
#'  that contains the run order
#' @param group_cname character string giving the name of the column in omicsData$f_data
#'  that contians the group information
#' @param ntree number of trees to grow in random forest model (default is set to 500)
#' @param keep_qc logical value to determine whether or not to include QC samples in the final output
#'   of the data (default is set to FALSE)
#'  
#' @return Object of same class as omicsData that has been undergone
#'   QCRFSC normalization
#'   
#' @details QCRFSC is ran on the raw abundance values. However, it is recommended to 
#' run imputation on the log2 scale. Therefore, when using QCRFSC, it is encouraged to
#' transform the data to log2 scale for imputation, and then transform the data back to
#' a raw abundance scale for bc_qcrfsc
#' 
#' @examples
#' library(malbacR)
#' library(pmartR)
#' data("pmart_amide")
#' pmart_amide <- edata_transform(pmart_amide,"log2")
#' impObj <- imputation(omicsData = pmart_amide)
#' amide_imp <- apply_imputation(imputeData = impObj, omicsData = pmart_amide)
#' amide_imp_abund <- edata_transform(amide_imp,"abundance")
#' amide_imp_abund <- group_designation(amide_imp_abund,main_effects = "group")
#' amide_qcrfsc <- bc_qcrfsc(omicsData = amide_imp_abund,qc_cname = "group",qc_val = "QC",order_cname = "Injection_order",group_cname = "group", ntree = 500,keep_qc = FALSE)
#' 
#' @author Damon Leach
#' 
#' @export
#' 
bc_qcrfsc <- function(omicsData,qc_cname,qc_val,order_cname,group_cname,ntree = 500,keep_qc = FALSE){
  # run through checks ---------------------------------------------------------
  
  # check that omicsData is of appropriate class #
  if (!inherits(omicsData, c("pepData", "proData", "metabData", "lipidData",
                             "nmrData"))) {
    
    stop (paste("omicsData must be of class 'pepData', 'proData', 'metabData',",
                "'lipidData', or 'nmrData'",
                sep = ' '))
  }
  
  # check that data is on abundance scale
  if(attributes(omicsData)$data_info$data_scale != "abundance"){
    stop ("QC-RFSC must be ran with raw abundance values. Please transform your data to 'abundance'.")
  }
  
  if(is.null(attributes(omicsData)$group_DF)){
    stop (paste("omicsData must have undergone group designation"))
  }
  
  # qc_cname - type of each sample
  if (class(qc_cname) != "character") {
    stop("Input parameter qc_cname must be of class 'character'.")
  }
  
  if (length(qc_cname) != 1) {
    stop("Input parameter qc_cname must be of length 1 (e.g. vector containing a single element")
  }
  
  if (!any(names(omicsData$f_data) == qc_cname)) {
    stop("Input parameter qc_cname must be a column found in f_data of omicsData.")
  }
  
  # qc_val - value for QC that must be in the column qc_cname
  # test_val - the value that will be used to determine what is the testing and what is the training
  if (class(qc_val) != "character") {
    stop("Input parameter qc_val must be of class 'character'.")
  }
  if (length(qc_val) != 1) {
    stop("Input parameter qc_val must be of length 1 (e.g. vector containing a single element)")
  }
  if (!qc_val %in% omicsData$f_data[,qc_cname]){
    stop("Input parameter qc_val must be a value in qc_cname column in omicsData$f_data")
  }
  
  # order_cname - the OVERALL run order
  if (class(order_cname) != "character") {
    stop("Input parameter order_cname must be of class 'character'.")
  }
  if (length(order_cname) != 1) {
    stop("Input parameter order_cname must be of length 1 (e.g. vector containing a single element)")
  }
  
  if (!any(names(omicsData$f_data) == order_cname)) {
    stop("Input parameter order_cname must be a column found in f_data of omicsData.")
  }
  
  if (!is.integer(omicsData$f_data[,order_cname])){
    stop ("Input parameter order_cname must contain integer values for run order")
  }
  
  if(length(omicsData$f_data[,order_cname]) != length(unique(omicsData$f_data[,order_cname]))){
    stop ("Input parameter order_cname cannot contain duplicate values for the overall run order")
  }
  
  # group_cname - type of each sample
  if (class(group_cname) != "character") {
    stop("Input parameter group_cname must be of class 'character'.")
  }
  
  if (length(group_cname) > 2) {
    stop("Input parameter group_cname must be of length 1 or 2 (e.g. vector containing a one or two elements")
  }
  
  if (!any(names(omicsData$f_data) == group_cname)) {
    stop("Input parameter group_cname must be a column found in f_data of omicsData.")
  }
  
  # check that value in keep_qc is logical
  if (!is.logical(keep_qc)) {
    stop("Input parameter keep_qc must be logical (either TRUE or FALSE)")
  }
  
  # check that value in keep_qc is length 1
  if (length(keep_qc) != 1) {
    stop("Input parameter qc_val must be of length 1 (e.g. vector containing a single element).")
  }
  
  # check for missing values
  # qcrfsc does not like missing values
  if (attributes(omicsData)$data_info$num_miss_obs != 0) {
    stop ("QCRFSC requires no missing observations. Remove molecules with missing samples.")
  }
  
  # QCRFSC requires that first and last sample be a QC sample
  minOrder <- min(omicsData$f_data[,order_cname],na.rm = T)
  maxOrder<- max(omicsData$f_data[,order_cname],na.rm = T)
  if(!(omicsData$f_data[,qc_cname][which(omicsData$f_data[[order_cname]] == minOrder)] == qc_val & 
       omicsData$f_data[,qc_cname][which(omicsData$f_data[[order_cname]] == maxOrder)] == qc_val)){
    stop ("QCRFSC requires that the first and last sample ran are QC samples")
  }
  
  # make sure that we have enough observed values prior to imputation
  # a bit challenging to check but alas we make do
  if(!is.null(attributes(omicsData)$imputation_info)){
    orig_edata = attributes(omicsData)$original_edata
    if(sum(dim(orig_edata) == dim(omicsData$e_data)) != 2){
      stop (paste("Imputation was applied prior to filtering. Please filter the data and then
                  apply imputation on the data."))
    }
  }
  
  if(!is.null(attributes(omicsData)$imputation_info)){
    orig_edata = attributes(omicsData)$original_edata
    groupDat <- attr(omicsData,"group_DF")
    id_col <- which(names(orig_edata) == pmartR::get_edata_cname(omicsData))
    ordering = omicsData$e_data[,id_col]
    # Create a data frame with the ID columns and the minimum number of non-missing values per grouping
    output <- orig_edata %>%
      tidyr::pivot_longer(cols = -tidyselect::all_of(id_col), names_to = names(groupDat)[1], values_to = "value") %>%
      dplyr::left_join(groupDat, by = pmartR::get_fdata_cname(omicsData)) %>%
      dplyr::group_by(dplyr::across(tidyselect::all_of(id_col)), Group) %>%
      {
        if(inherits(omicsData, "seqData")) {
          dplyr::summarise(., num_obs = sum(value != 0), .groups = "keep")
        } else {
          dplyr::summarise(., num_obs = sum(!is.na(value)),.groups = "keep")
        }
      } %>% 
      dplyr::group_by(dplyr::across(tidyselect::all_of(id_col))) %>%
      dplyr::summarise(min_num_obs = as.numeric(min(num_obs)),.groups = "keep") %>%
      dplyr::ungroup() %>%
      dplyr::rename(molecule = tidyselect::all_of(id_col)) %>%
      dplyr::arrange(match(molecule,ordering)) %>%
      data.frame()
    if(min(output$min_num_obs,na.rm=T) < 2){
      stop (paste("There is at least one molecule that did not have sufficient observed values per group
                  prior to imputation. Please run molecule_filter using groups = TRUE prior to running
                  imputation."))
    }
  }
  
  # run QCRSFC -----------------------------------------------------------------
  # retain seed after  running code
  if (!exists(".Random.seed")) runif(1)
  old_seed <- .Random.seed
  on.exit(.Random.seed <- old_seed)
  
  # obtain important information
  edata = omicsData$e_data
  fdata = omicsData$f_data
  edata_cname = pmartR::get_edata_cname(omicsData)
  fdata_cname = pmartR::get_fdata_cname(omicsData)
  edata_cname_colNum = which(colnames(edata) == edata_cname)
  
  # arrange the data by order number
  original_ordering = fdata[[order_cname]]
  vector_of_vars <- rlang::quos(order_cname)
  
  # order the data by run order
  fdata_ordered <- fdata %>%
    dplyr::arrange_at(order_cname)
  edata_ordered <- edata %>%
    dplyr::select()
  
  
  # set up random forest parameters
  numX <- 1:nrow(fdata_ordered)
  edata_df <- edata[,-edata_cname_colNum]
  edata_df <- edata_df %>%
    dplyr::select(fdata_ordered[[fdata_cname]])
  cn = colnames(edata_df)
  edata_matrix <- as.matrix(edata_df)
  
  # find which samples are the QC samples
  qcid <- which(fdata_ordered[,qc_cname] == qc_val)
  
  # run random forest for each molecule
  for(i in 1:nrow(edata)){
    set.seed(2023)
    temp <- randomForest::randomForest(data.frame(qcid), as.numeric(edata_matrix[i, qcid]), ntree = ntree)
    y <- data.frame(numX)
    colnames(y) <- "qcid"
    # predict based on RF results
    rfP <- predict(temp,y)
    # add to matrix information
    edata_matrix[i,] <- as.numeric(edata_matrix[i,])/rfP
  }
  # backtransform to 1000
  edata_matrix <- edata_matrix * 1000
  rownames(edata_matrix) <- edata[[edata_cname]]

  # remove QC samples from edata
  edata_loess = edata_matrix
  fdata_qcrfsc <- fdata_ordered
  
  if(keep_qc == FALSE){
    # filter out the old QC samples
    edata_loess = edata_matrix[,-qcid]
    fdata_qcrfsc <- fdata_ordered[fdata_ordered[,qc_cname] != qc_val,]
    
    # put fdata back into proper ordering
    fdata_subset = fdata[fdata[[order_cname]] %in% fdata_qcrfsc[[order_cname]],]
    fdata_qcrfsc <- fdata_qcrfsc[match(fdata_subset[[order_cname]],fdata_qcrfsc[[order_cname]]),]
  }
 
  edata_qcrfsc <- edata_loess %>% data.frame()
  # put edata back into proper ordering from initial running
  edata_qcrfsc <- tibble::rownames_to_column(edata_qcrfsc, var = edata_cname)
  edata_ordering <- rank(match(colnames(edata_qcrfsc),colnames(omicsData$e_data)))
  edata_qcrfsc <- edata_qcrfsc %>%
    dplyr::select(colnames(edata_qcrfsc)[edata_ordering])
  
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
    pmartObj = pmartR::as.isobaricpepData(e_data = edata_qcrfsc,
                                          edata_cname = edata_cname,
                                          f_data = fdata_qcrfsc,
                                          fdata_cname = fdata_cname,
                                          e_meta = emet,
                                          emeta_cname = emeta_cname)
  }
  else if(inherits(omicsData,"pepData")){
    pmartObj = pmartR::as.pepData(e_data = edata_qcrfsc,
                                  edata_cname = edata_cname,
                                  f_data = fdata_qcrfsc,
                                  fdata_cname = fdata_cname,
                                  e_meta = emet,
                                  emeta_cname = emeta_cname)
  }
  else if(inherits(omicsData,"proData")){
    pmartObj = pmartR::as.proData(e_data = edata_qcrfsc,
                                  edata_cname = edata_cname,
                                  f_data = fdata_qcrfsc,
                                  fdata_cname = fdata_cname,
                                  e_meta = emet,
                                  emeta_cname = emeta_cname)
  }
  else if(inherits(omicsData,"metabData")){
    pmartObj = pmartR::as.metabData(e_data = edata_qcrfsc,
                                    edata_cname = edata_cname,
                                    f_data = fdata_qcrfsc,
                                    fdata_cname = fdata_cname,
                                    e_meta = emet,
                                    emeta_cname = emeta_cname)
  }
  else if(inherits(omicsData,"lipidData")){
    pmartObj = pmartR::as.lipidData(e_data = edata_qcrfsc,
                                    edata_cname = edata_cname,
                                    f_data = fdata_qcrfsc,
                                    fdata_cname = fdata_cname,
                                    e_meta = emet,
                                    emeta_cname = emeta_cname)
  }
  else if(inherits(omicsData,"nmrData")){
    pmartObj = pmartR::as.nmrData(e_data = edata_qcrfsc,
                                  edata_cname = edata_cname,
                                  f_data = fdata_qcrfsc,
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

  # check group designation
  # since we are removing samples if keep_qc != TRUE
  if(!is.null(attributes(attr(omicsData,"group_DF"))$batch_id)){
    batch_id_col = which(colnames(attributes(attr(omicsData,"group_DF"))$batch_id) != fdata_cname)
    batch_id_name = colnames(attributes(attr(omicsData,"group_DF"))$batch_id)[batch_id_col]
    pmartObj <- pmartR::group_designation(pmartObj,main_effects = group_cname,
                                          batch_id = batch_id_name)
  } else {
    pmartObj <- pmartR::group_designation(pmartObj,main_effects = group_cname)
  }
  if(keep_qc == TRUE){
    attr(pmartObj,"group_DF") = attr(omicsData,"group_DF")
  }

  
  # Update the data_info attribute for batch
  attributes(pmartObj)$data_info$batch_info <- list(
    is_bc = TRUE,
    bc_method = "bc_qcrfsc",
    params = list(qc_cname = qc_cname,
                  qc_val = qc_val,
                  order_cname = order_cname,
                  group_cname = group_cname,
                  ntree = ntree,
                  keep_qc = keep_qc)
  )
  # update normalization as well 
  attributes(pmartObj)$data_info$norm_info <- list(
    is_normalized = TRUE,
    norm_type = "bc_qcrfsc"
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
