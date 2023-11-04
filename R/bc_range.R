#' Range Scaling adjusted object
#' 
#' This function returns a pmart object that has been undergone Range Scaling
#' normalization
#' 
#' @param omicsData an object of the class 'pepData', 'proData', 'metabData',
#'   'lipidData', or 'nmrData', usually created by \code{\link{as.pepData}},
#'   \code{\link{as.proData}}, \code{\link{as.metabData}},
#'   \code{\link{as.lipidData}}, or \code{\link{as.nmrData}}, respectively.
#'
#' @return Object of same class as omicsData that has been undergone
#'   Range Scaling Normalization as batch effect correction
#' 
#' @examples
#' library(malbacR)
#' library(pmartR)
#' data("pmart_amide")
#' pmart_amide <- edata_transform(pmart_amide,"log2")
#' pmart_amide <- group_designation(pmart_amide,main_effects = "group",batch_id = "batch")
#' amide_range <- bc_range(omicsData = pmart_amide)
#' 
#' @author Damon Leach
#' 
#' @export
#' 
bc_range <- function(omicsData) {
  
  # run checks -----------------------------------------------------------------
  
  # check that omicsData is of appropriate class #
  if (!inherits(omicsData, c("pepData", "proData", "metabData", "lipidData",
                             "nmrData"))) {
    
    stop (paste("omicsData must be of class 'pepData', 'proData', 'metabData',",
                "'lipidData', or 'nmrData'",
                sep = ' '))
  }
  
  # ensure that we have at least two abundance values for each biomolecule
  # find the edata_cname column
  edat_id <- which(colnames(omicsData$e_data) == pmartR::get_edata_cname(omicsData))
  # find the number of non missing values for each biomolecule
  missing_val_sum <- rowSums(!is.na(omicsData$e_data[,-edat_id]))
  # if there are 0 or 1 abundance values for the biomolecule, we have an issue
  if(min(missing_val_sum) < 2) {
    stop (paste("Range Scaling requires that each biomolecule be present in at least 2 samples"))
  }
  
  # important info for compiling pmart object later ----------------------------
  # find the individual data sets
  edat <- omicsData$e_data
  fdat <- omicsData$f_data
  emet <- omicsData$e_meta
  
  # get the variables we need to create pmart object
  edat_cname = pmartR::get_edata_cname(omicsData)
  fdat_cname = pmartR::get_fdata_cname(omicsData)
  emet_cname = pmartR::get_emeta_cname(omicsData)
  
  # run scaling on data --------------------------------------------------------
  # remove the name
  cname <- pmartR::get_edata_cname(omicsData)
  
  # range scaling calculation
  edatScaled <- edat %>%
    dplyr::select(-cname) %>%
    scale_method(methods = "range") %>%
    t()
  
  # create pmart object --------------------------------------------------------
  
  # add back in edata_cname to edata
  edatScaled <- edat %>%
    dplyr::select(cname) %>%
    cbind(edatScaled)
  
  # create the pmart object #
  if(inherits(omicsData,"isobaricpepData") & inherits(omicsData,"pepData")){
    pmartObj = pmartR::as.isobaricpepData(e_data = edatScaled,
                                          edata_cname = edat_cname,
                                          f_data = fdat,
                                          fdata_cname = fdat_cname,
                                          e_meta = emet,
                                          emeta_cname = emet_cname)
  }
  else if(inherits(omicsData,"pepData")){
    pmartObj = pmartR::as.pepData(e_data = edatScaled,
                                  edata_cname = edat_cname,
                                  f_data = fdat,
                                  fdata_cname = fdat_cname,
                                  e_meta = emet,
                                  emeta_cname = emet_cname)
  }
  else if(inherits(omicsData,"proData")){
    pmartObj = pmartR::as.proData(e_data = edatScaled,
                                  edata_cname = edat_cname,
                                  f_data = fdat,
                                  fdata_cname = fdat_cname,
                                  e_meta = emet,
                                  emeta_cname = emet_cname)
  }
  else if(inherits(omicsData,"metabData")){
    pmartObj = pmartR::as.metabData(e_data = edatScaled,
                                    edata_cname = edat_cname,
                                    f_data = fdat,
                                    fdata_cname = fdat_cname,
                                    e_meta = emet,
                                    emeta_cname = emet_cname)
  }
  else if(inherits(omicsData,"lipidData")){
    pmartObj = pmartR::as.lipidData(e_data = edatScaled,
                                    edata_cname = edat_cname,
                                    f_data = fdat,
                                    fdata_cname = fdat_cname,
                                    e_meta = emet,
                                    emeta_cname = emet_cname)
  }
  else if(inherits(omicsData,"nmrData")){
    pmartObj = pmartR::as.nmrData(e_data = edatScaled,
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
  
  # Update the data_info attribute for batch
  attributes(pmartObj)$data_info$batch_info <- list(
    is_bc = TRUE,
    bc_method = "bc_range",
    params = list()
  )
  # update normalization as well 
  attributes(pmartObj)$data_info$norm_info <- list(
    is_normalized = TRUE,
    norm_type = "bc_range"
  )
  
  # Update the meta_info attribute.
  attr(pmartObj, 'meta_info') <- pmartR:::set_meta_info(
    e_meta = omicsData$e_meta,
    emeta_cname = pmartR::get_emeta_cname(omicsData)
  )
  
  # Update the filtesr attribute
  attr(pmartObj, 'filters') <- attr(omicsData,'filters')
  
  # return pmart object
  return(pmartObj)
}
