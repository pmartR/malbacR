#' Get optimal span value
#'
#' Get optimal span value for QC-RLSC normalization
#'
#' @param qc_data data.frame of only the QC samples from one block/batch, where
#'   the first column gives the run order and the second column gives the
#'   abundance value
#' @param lam vector of (numeric) possible polynomial degrees to test in the
#'   loess fitting; defaults to 1,
#' @return list containing elements for: Lambda, Alpha, MSE
#' @details This internal function is called by the \code{get_params} function
#' 
#' @author Lisa Bramer, Kelly Stratton
#' @rdname optimal_span
#' @name optimal_span
#' @export
optimal_span <- function(qc_data, lam = c(1, 2)) {
  ### initial checks ###

  # check that qc_data is data.frame with 2 columns #
  if (class(qc_data) != "data.frame") {
    stop("qc_data must be a data.frame")
  }
  if (ncol(qc_data) != 2) {
    stop("qc_data must have 2 columns")
  }

  # check that lam is numeric vector #
  if (!is.numeric(lam)) {
    stop("lam must be a numeric vector")
  }
  if (!is.vector(lam)) {
    stop("lam must be a numeric vector")
  }

  ### end of initial checks ###

  qc_data = qc_data[!is.na(qc_data[, 2]), ] # remove QCs that have NA values for abundance

  # determine number of QC samples #
  n_qc = nrow(qc_data) - 1

  # # initialize results list #
  # all_res = list()
  
  # create a tibble for the different data we need to map
  span_list = tibble::tibble(cur_lam = as.list(lam))
  
  span_list <- span_list %>%
    dplyr::mutate(
      # create a sequence of possible alpha values to be considered
      alp_vals = purrr::map(cur_lam,function(cl){
        alp_vals = seq(from = (cl + 1),
                       to = n_qc,
                       by = 1)/n_qc
        
        alp_vals <- alp_vals[alp_vals*nrow(qc_data) >= 5]
        return(alp_vals)
        }),
      
      # set up a res matrix
      res = purrr::map(alp_vals,function(av){
        res = matrix(NA,nrow = (n_qc+1), ncol = length(av))
        return(res)
      }),
      
      # set up current alpha value
      cur_alp = purrr::map(alp_vals, function(av){
        cur_alp = as.list(av)
        return(cur_alp)
      }),
      
      # set up k value
      k_val = purrr::map(cur_alp,function(cl){
        purrr::map(cl,function(clmini){
          return(as.list(seq(1:(n_qc + 1))))
        })
      }))
  
  # this portion was just used for debugging
  # k_val = span_list$k_val
  # cur_alp = span_list$cur_alp
  # cur_lam = span_list$cur_lam
  # qc_data
  # qc6 <- list(k_val = k_val,cur_alp = cur_alp,cur_lam = cur_lam,qc_data = qc_data)
  # saveRDS(qc6,"qc6.rds")
  # qc_data <- qc4$qc_data
  # for(i in 1:length(k_val)){
  #   i = 2
  #   for(j in 1:length(k_val[[1]])){
  #     for(k in 1:length(k_val[[1]][[1]])){
  #       k = 25
  #       data = qc_data[-(k_val[[i]][[j]][[k]]),]
  #       model = stats::loess(
  #         data[,2] ~ data[,1],
  #         span = cur_alp[[i]][[j]],
  #         # the issue with degree 2 because polynomial
  #         degree = cur_lam[[i]],
  #         control = loess.control(surface = "direct"))
  #       pred = predict(model,newdata = qc_data[k_val[[i]][[j]][[k]],1])
  #     }
  #   }
  # }
  
  span_list <- span_list %>%
    dplyr::mutate(
      # create a loess model and find the predicted value
      pred = purrr::pmap(list(k_val,cur_alp,cur_lam),function(kv,ca,cl){
        purrr::map2(kv,ca,function(kvmini,camini){
          purrr::map(kvmini,function(kvmini2){
            data = qc_data[-(kvmini2),]
            model = suppressWarnings(stats::loess(
              data[,2] ~ data[,1],
              span = camini,
              degree = cl,
              control = loess.control(surface = "direct")))
            pred = predict(model,newdata = qc_data[kvmini2,1])
            return(pred)
          })
        })
      }),
      potential_warning = purrr::pmap(list(k_val,cur_alp,cur_lam),function(kv,ca,cl){
        purrr::map2(kv,ca,function(kvmini,camini){
          purrr::map(kvmini,function(kvmini2){
            data = qc_data[-(kvmini2),]
            
            potential_warn = tryCatch(
              {model = stats::loess(
                data[,2] ~ data[,1],
                span = camini,
                # the issue with degree 2 because polynomial
                degree = cl,
                control = loess.control(surface = "direct"))},
              warning = function(w) {c(lambda = camini,alpha = cl)}
            )
            if(!is.numeric(potential_warn)){potential_warn = c(lambda = NA,alpha = NA)}
            return(potential_warn)
          })
        })
      }),
      # update the res matrix
      res = purrr::map2(res,pred,function(resobj,predobj){
        for(j in 1:ncol(resobj)) {
          resobj[,j] = unlist(predobj[[j]])
        }
        return(resobj)
      })
  )
  
  span_list$temp = lapply(span_list$res, function(x) {
    apply((x - qc_data[,2]) ^2, 2, mean)
  })
  
  span_list <- span_list %>%
    dplyr::mutate(res2 = purrr::pmap(list(temp,cur_lam),function(tobj,cl){
      alp = seq(from = cl + 1,
                to = n_qc,
                by = 1) / n_qc
      alp = alp[alp*nrow(qc_data) >= 5]

      res = data.frame(Lambda = rep(cl,length(alp)),
                       Alpha = alp,
                       MSE = tobj)
      return(res)
    }))
  
  # if we are using the one with errors let the user know
  potWarn = unlist(span_list$potential_warning)
  alphaVal = potWarn[names(potWarn) == "alpha"]
  alphaVal = alphaVal[!is.na(alphaVal)]
  lambdaVal = potWarn[names(potWarn) == "lambda"]
  lambdaVal = lambdaVal[!is.na(lambdaVal)]

  final = do.call(rbind, span_list$res2)
  final2 = final[which.min(final$MSE),]
  final3 <- list(
    Lambda = final2$Lambda,
    Alpha = final2$Alpha,
    MSE = final2$MSE,
    all_res = final
  )
  
  if(final3$Alpha %in% alphaVal & final3$Lambda %in% lambdaVal){
    warning("The optimal spans resulted in NA for standard errors. Proceed with caution")
  }
  return(final3)
}
