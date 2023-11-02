#' Fetch the Batch Effect Parameters
#' 
#' Retrieves the values from the batch effect parameters attributes
#' 
#' @param omicsData An object of class pepData, proData, metabData, lipidData,
#'                  or nmrData.
#' 
#' @return A list for each of the different parameters used to run batch correction
#' 
#' @export
#' @name get_batch_parameters
#' 
get_batch_parameters <- function (omicsData) {
  # Check class of omicsData.
  if (!inherits(omicsData, c("pepData",
                             "proData",
                             "metabData",
                             "lipidData",
                             "nmrData",
                             "seqData"))) {
    
    # Lay down an error in the console.
    stop (paste("omicsData must be of class 'pepData', 'proData',",
                "'metabData', 'lipidData', 'nmrData', or 'seqData'",
                sep = " "))
    
  }
  
  # Extract and return the group_DF attribute.
  return (attributes(omicsData)$data_info$batch_info$params)
}

#' Fetch the Batch Effect Method Name
#' 
#' Retrieves the name of the batch effect
#' 
#' @param omicsData An object of class pepData, proData, metabData, lipidData,
#'                  or nmrData.
#' 
#' @return A character vector specifying the batch correction method used (NULL if no batch correction applied)
#' 
#' @export
#' @name get_batch_method
#' 
get_batch_method <- function (omicsData) {
  # Check class of omicsData.
  if (!inherits(omicsData, c("pepData",
                             "proData",
                             "metabData",
                             "lipidData",
                             "nmrData",
                             "seqData"))) {
    
    # Lay down an error in the console.
    stop (paste("omicsData must be of class 'pepData', 'proData',",
                "'metabData', 'lipidData', 'nmrData', or 'seqData'",
                sep = " "))
    
  }
  
  # Extract and return the group_DF attribute.
  return (attributes(omicsData)$data_info$batch_info$bc_method)
}
