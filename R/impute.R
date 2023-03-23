#' imputation
#' 
#' This function returns a pmart filter that has imputed missing values by "em" method
#' 
#' @param omicsData an object of the class 'pepData', 'proData', 'metabData',
#'   'lipidData', or 'nmrData', usually created by \code{\link{as.pepData}},
#'   \code{\link{as.proData}}, \code{\link{as.metabData}},
#'   \code{\link{as.lipidData}}, or \code{\link{as.nmrData}}, respectively.
#'  
#' @return object of class imputeData and dataframe with the imputated e_data 
#' 
#' @examples
#' library(malbacR)
#' data("pmart_amide")
#' impFilt <- imputation(pmart_amide)
#' 
#' @author Evan Martin, Damon Leach
#' 
#' @export
imputation <- function (omicsData) {

  # Make sure the input is the correct class.
  if (!inherits(omicsData, c("pepData", "proData", "metabData", "lipidData",
                             "nmrData"))) {
    
    stop (paste("omicsData must be of class 'pepData', 'proData', 'metabData',",
                "'lipidData', or 'nmrData'",
                sep = ' '))
  }

  # # Ensure the imputation method is an accepted method.
  # if (!method %in% c("em")) {
  # 
  #   # You shall not pass!!!
  #   stop ("The input to method does not match any of the accepted values.")
  # 
  # }
  method = "em"
  
  # get edata name information
  edat_cname = pmartR::get_edata_cname(omicsData)
  id_col <- which(colnames(omicsData$e_data) == edat_cname)
  
  # make sure there are no rows with completely NA values
  if(sum(rowSums(!is.na(omicsData$e_data[,-id_col])) == 0)){
    stop ("At least one row has all NA values. Run 'molecule_filter' before running any imputation.")
  }
  
  # make sure we actual have missing values
  if(sum(is.na(omicsData$e_data[,-id_col])) == 0){
    stop ("There are no missing values in this dataset. Imputation is not necessary.")
  }
  
  # tranpose the data to be molecules as columns and rows are samples
  tdat <- t(omicsData$e_data[,-id_col]) %>%
    data.frame()
  
    
  # Impute data using the impute_EM function
  imputed <- mvdalab::imputeEM(tdat,impute.ncomps = 1)
  imputed_t <- data.frame(imputed$Imputed.DataFrames)

  # transpose back
  imputed_data <- t(imputed_t)

  # Return an imputation object with the imputed_values and missing_by_sample
  # attributes.
  return (
    structure(
      as.data.frame(imputed_data, check.names = FALSE),
      # Add in imputeData attributes.
      original_edata = omicsData$e_data,
      missing_by_sample = omicsData$e_data[, -id_col] %>%
        purrr::map_df(~sum(is.na(.)) / length(.)) %>%
        `class<-`("data.frame"),
      imputation_method = method,
      class = c("imputeData", "data.frame")
    )
  )

}

#' apply_imputation
#' 
#' This function returns a pmart object that has undergone imputation
#' 
#' @param imputeData object of class imputeData that contains the imputed e_data
#'   determined using the function 'imputation'
#' @param omicsData an object of the class 'pepData', 'proData', 'metabData',
#'   'lipidData', or 'nmrData', usually created by \code{\link{as.pepData}},
#'   \code{\link{as.proData}}, \code{\link{as.metabData}},
#'   \code{\link{as.lipidData}}, or \code{\link{as.nmrData}}, respectively.
#'  
#' @return object of pmartR object containing imputated e_data 
#' 
#' @examples
#' library(malbacR)
#' data("pmart_amide")
#' impFilt <- imputation(pmart_amide)
#' imputed_data <- apply_imputation(impFilt,pmart_amide)
#' 
#' @author Evan Martin, Damon Leach
#' 
#' @export
apply_imputation <- function (imputeData, omicsData) {

  # Preliminaries --------------------------------------------------------------

  # We have to carry out all these preliminary checks because we were asked:
  # "When using this software, what would an idiot do?"

  if (!inherits(imputeData, "imputeData")) {
    stop ("The input to the imputeData argument must be an imputeData object.")
  }

  # Make sure the input is the correct class.
  if (!inherits(omicsData, c("pepData", "proData", "metabData", "lipidData",
                             "nmrData"))) {
    
    stop (paste("omicsData must be of class 'pepData', 'proData', 'metabData',",
                "'lipidData', or 'nmrData'",
                sep = ' '))
  }

  # Apply imputation -----------------------------------------------------------
  
  edat_cname = pmartR::get_edata_cname(omicsData)
  id_col <- which(colnames(omicsData$e_data) == edat_cname)
  
  # Replace the original values in e_data with the imputed values.
  omicsData$e_data[, -id_col] <- imputeData

  # Extract the data_info attribute from the omicsData object. This will be used to
  # update the data_info attribute with the imputed data.
  the_info <- pmartR::get_data_info(omicsData)

  # Oh look! Here is the end of the function already.
  #
  # The actual imputation was carried out in the imputation function. All we
  # need to do here is add the omicsData attributes to the data frame with the
  # imputed values. The data_info attribute will be updated but all other
  # attributes will remain the same.
  return (
    structure(
      omicsData,
      data_info = pmartR:::set_data_info(
        e_data = omicsData$e_data,
        edata_cname = pmartR::get_edata_cname(omicsData),
        data_scale_orig = the_info$data_scale_orig,
        data_scale = the_info$data_scale,
        data_types = the_info$data_types,
        norm_info = the_info$norm_info,
        is_normalized = the_info$norm_info$is_normalized,
        batch_info = the_info$batch_info,
        is_bc = the_info$batch_info$is_bc
      ),
      original_edata = attr(imputeData, "original_edata"),
      imputation_info = list(
        imputation_method = attr(imputeData, "imputation_method"),
        prop_imputed = the_info$prop_missing
      )
    )
  )

}
