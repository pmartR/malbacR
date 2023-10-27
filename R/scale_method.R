#' Perform scale method on data
#' 
#' This function is a modified version of the function ScalingMethods from the 
#' package DiffCorr, but it allows for missing values for Power and Range scaling
#' 
#' @param data a data matrix or data frame object where each row is a biomolecule
#' and each column each sample or replicate
#' @param methods a string character specifying the pre-treatment method and can be 
#' "auto", "range", "pareto", "vast", "level", or "power"
#'
#' @return a data matrix or data frame object (depending on what was the input) that
#' has undergone pre-treatment method
#' 
#' @examples
#' library(malbacR)
#' data(pmart_amide)
#' malbacR:::scale_method(pmart_amide$e_data[,-1])
#' 
#' @author Damon Leach
#' 
#' 
scale_method <- function (data, methods = c("auto", "range", "pareto", "vast", 
                                       "level", "power")) 
{
  methods <- match.arg(methods)
  if (ncol(data) > 1) {
    switch(methods, auto = {
      res <- apply(data, 2, function(x) (x - mean(x, na.rm = TRUE))/sd(x, 
                                                                       na.rm = TRUE))
      return(data.frame(res, check.names = FALSE))
    }, range = {
      res <- apply(data, 2, function(x) (x - mean(x, na.rm = TRUE))/(range(x,na.rm=T)[2] - 
                                                                       range(x,na.rm=T)[1]))
      return(data.frame(res, check.names = FALSE))
    }, pareto = {
      res <- apply(data, 2, function(x) (x - mean(x, na.rm = TRUE))/sqrt(sd(x, 
                                                                            na.rm = TRUE)))
      return(data.frame(res, check.names = FALSE))
    }, vast = {
      res <- apply(data, 2, function(x) mean(x, na.rm = TRUE) * 
                     (x - mean(x, na.rm = TRUE))/(sd(x, na.rm = TRUE)^2))
      return(data.frame(res, check.names = FALSE))
    }, level = {
      res <- apply(data, 2, function(x) (x - mean(x, na.rm = TRUE))/mean(x, 
                                                                         na.rm = TRUE))
      return(data.frame(res, check.names = FALSE))
    }, power = {
      res <- apply(data, 2, function(x) sqrt(x) - mean(sqrt(x),na.rm=T))
      return(data.frame(res, check.names = FALSE))
    })
  }
}
