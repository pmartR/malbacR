#' NOMIS normalized adjusted object
#' 
#' This function returns a pmart object that has been undergone NOMIS
#' normalization regarding batch effect correction
#' 
#' @param omicsData an object of the class 'pepData', 'proData', 'metabData',
#'   'lipidData', or 'nmrData', usually created by \code{\link{as.pepData}},
#'   \code{\link{as.proData}}, \code{\link{as.metabData}},
#'   \code{\link{as.lipidData}}, or \code{\link{as.nmrData}}, respectively.
#' @param is_cname character string giving the name of the column in omicsData$e_meta
#'  that contains the factor variable indicating whether the molecule is internal
#'  standards or not
#' @param is_val value giving the name/value in the is_cname column that indicates
#'  that the sample is an Internal Standards sample
#' @param num_pc numeric variable providing the number of principal components for NOMIS
#'  to consider when running the batch correction method, defaults to 2
#'  
#' @return Object of same class as omicsData that has been undergone
#'   NOMIS normalization
#' 
#' @examples
#' library(malbacR)
#' library(pmartR)
#' data("pmart_mix")
#' pmart_mix <- edata_transform(pmart_mix,"log2")
#' pmart_mix <- group_designation(pmart_mix,main_effects = "BatchNum",batch_id = "BatchNum")
#' pmart_mix <- normalize_global(pmart_mix,subset_fn = "all",norm_fn = "median",
#'                               apply_norm = TRUE,backtransform = TRUE)
#' mix_nomis <- bc_nomis(omicsData = pmart_mix, is_cname = "tag", is_val = "IS", num_pc = 2)
#' 
#' @author Damon Leach
#' 
#' @export
#' 
bc_nomis <- function(omicsData,is_cname,is_val,num_pc = 2){
  
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
  
  # is_cname
  if (class(is_cname) != "character") {
    stop("Input parameter is_cname must be of class 'character'.")
  }
  
  if (length(is_cname) != 1) {
    stop("Input parameter is_cname must be of length 1 (e.g. vector containing a single element")
  }
  
  if (!any(names(omicsData$e_meta) == is_cname)) {
    stop("Input parameter is_cname must be a column found in e_meta of omicsData.")
  }
  
  # is_val
  if (class(is_val) != "character") {
    stop("Input parameter is_val must be of class 'character'.")
  }
  
  if (length(is_val) != 1) {
    stop("Input parameter is_val must be of length 1 (e.g. vector containing a single element")
  }
  
  if (!(is_val %in% omicsData$e_meta[, is_cname])) {
    stop(
      "Input parameter is_val must be a value found in the is_cname column from the e_meta of omicsData."
    )
  }
  
  # num_pc
  if (class(num_pc) != "numeric"){
    stop("Input parameter num_pc must be of class 'numeric'.")
  }
  
  # NOMIS requires no missing values
  # which column has the edata cname
  cnameCol <- which(colnames(omicsData$e_data) == pmartR::get_edata_cname(omicsData))
  if (sum(is.na(omicsData$e_data[,-cnameCol])) != 0) {
    stop ("NOMIS requires no missing observations. Remove molecules with missing samples.")
  }
  
  # we cannot have negative values
  if(sum(omicsData$e_data < 0,na.rm=TRUE) > 0){
    stop("NOMIS cannot run with expression data that has negative values (likely
         due to normalization without a backtransform")
  }
  
  # run the nomis calculations -------------------------------------------------
  
  # find the values needed for the normalize function
  # set up the object parameter (the e_data)
  edat <- as.matrix(omicsData$e_data[,-1])
  rownames(edat) <- omicsData$e_data[,1]
  rowNamEdat = rownames(edat)

  # set up the factors parameter which is model.matrix with batch information
  # obain batch information
  batch_attributes <- attributes(attr(omicsData,"group_DF"))$batch_id
  # find proper order for the samples in case edata and fdata are mixed around
  info = pmartR::get_fdata_cname(omicsData)
  batch_info <- data.frame(samples = colnames(edat))
  # update the attribute order to match edata
  batch_info <- batch_info%>%
    dplyr::left_join(batch_attributes, by = c("samples" = info))
  
  batch <- as.factor(batch_info[,2])
  modMatrix <- stats::model.matrix(~-1+batch)
  
  # set us the standards parameter which is the internal standards information
  # make sure the emeta and edata are ordered in the correct order
  edat_cname = pmartR::get_edata_cname(omicsData)
  proper_order_is <- omicsData$e_data %>%
    dplyr::select(dplyr::all_of(edat_cname)) %>%
    dplyr::left_join(omicsData$e_meta)
  
  is_cnameNum <- which(colnames(proper_order_is) == is_cname)
  isInternal <- proper_order_is[,is_cnameNum] == is_val
  
  # now create the nomis abundance values!
  # we default the number of principal components to be 2
  # we set method to be NOMIS
  edatNomis <- crmn::normalize(object = edat,
                         method = "nomis",
                         factor = modMatrix,
                         standards = isInternal,
                         ncomp = num_pc)
  
  # move back into pmart object ------------------------------------------------
  
  # put edata into proper format
  edatNomis <- data.frame(edatNomis)
  edatNomis <- tibble::rownames_to_column(edatNomis, var = pmartR::get_edata_cname(omicsData))
  edat_colNames <- colnames(omicsData$e_data)
  colnames(edatNomis) <- edat_colNames
  
  # find the important values for pmart creation
  edat_cname = pmartR::get_edata_cname(omicsData)
  fdata_cname = pmartR::get_fdata_cname(omicsData)
  fdat = omicsData$f_data[omicsData$f_data[[fdata_cname]] %in% colnames(edatNomis),]
  fdat_cname = pmartR::get_fdata_cname(omicsData)
  
  # assume emeta is NULL unless otherwise stated
  emet_cname = NULL
  emet = NULL
  if(!is.null(omicsData$e_meta)){
    # need to update the e_meta if it exists as NOMIS removes the IS values
    emet_cname = pmartR::get_emeta_cname(omicsData)
    emet_cnameNum = which(colnames(omicsData$e_meta) == emet_cname)
    edat_cnameNum = which(colnames(omicsData$e_data) == edat_cname)
    emet = omicsData$e_meta[omicsData$e_meta[,emet_cnameNum] %in% edatNomis[,edat_cnameNum],]
  }
  
  # create the pmart object #
  if(inherits(omicsData,"isobaricpepData") & inherits(omicsData,"pepData")){
    pmartObj = pmartR::as.isobaricpepData(e_data = edatNomis,
                                          edata_cname = edat_cname,
                                          f_data = fdat,
                                          fdata_cname = fdat_cname,
                                          e_meta = emet,
                                          emeta_cname = emet_cname)
  }
  else if(inherits(omicsData,"pepData")){
    pmartObj = pmartR::as.pepData(e_data = edatNomis,
                                  edata_cname = edat_cname,
                                  f_data = fdat,
                                  fdata_cname = fdat_cname,
                                  e_meta = emet,
                                  emeta_cname = emet_cname)
  }
  else if(inherits(omicsData,"proData")){
    pmartObj = pmartR::as.proData(e_data = edatNomis,
                                  edata_cname = edat_cname,
                                  f_data = fdat,
                                  fdata_cname = fdat_cname,
                                  e_meta = emet,
                                  emeta_cname = emet_cname)
  }
  else if(inherits(omicsData,"metabData")){
    pmartObj = pmartR::as.metabData(e_data = edatNomis,
                                    edata_cname = edat_cname,
                                    f_data = fdat,
                                    fdata_cname = fdat_cname,
                                    e_meta = emet,
                                    emeta_cname = emet_cname)
  }
  else if(inherits(omicsData,"lipidData")){
    pmartObj = pmartR::as.lipidData(e_data = edatNomis,
                                    edata_cname = edat_cname,
                                    f_data = fdat,
                                    fdata_cname = fdat_cname,
                                    e_meta = emet,
                                    emeta_cname = emet_cname)
  }
  else if(inherits(omicsData,"nmrData")){
    pmartObj = pmartR::as.nmrData(e_data = edatNomis,
                                  edata_cname = edat_cname,
                                  f_data = fdat,
                                  fdata_cname = fdat_cname,
                                  e_meta = emet,
                                  emeta_cname = emet_cname)
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
    bc_method = "nomis",
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
