#' Internal.compute_cor
#' 
#' This function is an internal function to compute correlation
#'  
#' @return return correlation computed
#' 
Internal.compute_cor <- function(train_num, test_num = NULL,
                                 selectVar_corType   = c("pcor", "cor"),
                                 selectVar_corMethod = c("spearman", "pearson")) {
  
  if (selectVar_corType == "pcor") {
    train_num_noNA <- Internal.remove_NA(train_num, data_label = "training data")
    if (!is.null(test_num)) test_num_noNA <- Internal.remove_NA(test_num, data_label = "test data")
    
    if (ncol(train_num_noNA) > 500) message("      Your data have more than 500 variables. It may take some time to process large datasets.")
    
    train_cor <- data.frame(ppcor::pcor(train_num_noNA, method = selectVar_corMethod)$estimate)
    names(train_cor) <- names(train_num_noNA)
    row.names(train_cor) <- names(train_num_noNA)
    
    if (!is.null(test_num)) {
      test_cor <- data.frame(ppcor::pcor(test_num_noNA, method = selectVar_corMethod)$estimate)
      names(test_cor) <- names(test_num_noNA)
      row.names(test_cor) <- names(test_num_noNA)
    } else test_cor <- NULL
    
  } else {
    
    train_cor <- data.frame(cor(train_num, method = selectVar_corMethod, use = "complete.obs"))
    
    if (!is.null(test_num)) {
      test_cor <- data.frame(cor(test_num, method = selectVar_corMethod, use = "complete.obs"))
    } else test_cor <- NULL
  }
  
  cor_info <- list(variable_name = names(train_cor),
                   train_cor = train_cor, test_cor = test_cor)
}

Internal.select_variable <- function(cor_info, selectVar_minNum = NULL,
                                     selectVar_maxNum = NULL) {
  
  variable_name <- cor_info$variable_name
  
  train_cor <- cor_info$train_cor
  test_cor  <- cor_info$test_cor
  
  train_cor[is.na(train_cor)] <- 0
  if (!is.null(test_cor)) test_cor[is.na(test_cor)] <- 0
  
  pb <- pbapply::timerProgressBar(min = 0, max = length(variable_name),
                                  initial = 0, style = 3, width = 70,
                                  min_time = 20)
  selected_var <- lapply(1:length(variable_name), function(var_idx,
                                                           variable_name,
                                                           train_cor, test_cor,
                                                           selectVar_minNum, selectVar_maxNum) {
    pbapply::setTimerProgressBar(pb, var_idx)
    
    input_one_variable_name <- variable_name[[var_idx]]
    correlated_train <- abs(train_cor[input_one_variable_name])
    correlated_train <- correlated_train[!row.names(correlated_train) %in% input_one_variable_name,,drop = FALSE]
    correlated_train <- correlated_train[order(correlated_train[[1]], decreasing = TRUE),,drop = FALSE]
    
    correlated_train_name <- row.names(correlated_train)
    candidate_train_name  <- correlated_train_name[correlated_train[[1]] > 0.5]
    candidate_var_name    <- candidate_train_name
    
    if (!is.null(test_cor)) {
      correlated_test <- abs(test_cor[input_one_variable_name])
      correlated_test <- correlated_test[!row.names(correlated_test) %in% input_one_variable_name,,drop = FALSE]
      correlated_test <- correlated_test[order(correlated_test[[1]], decreasing = TRUE),,drop = FALSE]
      
      correlated_test_name <- row.names(correlated_test)
      candidate_test_name  <- correlated_test_name[correlated_test[[1]] > 0.5]
      candidate_var_name   <- intersect(candidate_var_name, candidate_test_name)
    }
    
    if (length(candidate_var_name) < selectVar_minNum) {
      
      if (is.null(test_cor)) {
        selected_var_name <- correlated_train_name[1:selectVar_minNum]
      } else {
        current_upper_limit <- min(length(candidate_train_name), length(candidate_test_name))
        current_upper_limit <- max(current_upper_limit, selectVar_minNum)
        
        while (1) {
          candidate_test_name_tmp  <- correlated_test_name[1:current_upper_limit]
          candidate_train_name_tmp <- correlated_train_name[1:current_upper_limit]
          candidate_var_name_tmp <- intersect(candidate_train_name_tmp, candidate_test_name_tmp)
          if (length(candidate_var_name_tmp) < selectVar_minNum) {
            current_upper_limit <- current_upper_limit + 1
          } else {
            selected_var_name <- candidate_var_name_tmp
            break
          }
        }
        
        if (length(selected_var_name) > selectVar_maxNum) {
          selected_var_name <- selected_var_name[1:selectVar_maxNum]
        }
        
      }
      
    } else if (length(candidate_var_name) > selectVar_maxNum) {
      
      selected_var_name <- candidate_var_name[1:selectVar_maxNum]
      
    } else {
      selected_var_name <- candidate_var_name
    }
    selected_var_name
  },
  variable_name = variable_name,
  train_cor = train_cor, test_cor = test_cor,
  selectVar_minNum = selectVar_minNum, selectVar_maxNum = selectVar_maxNum)
  pbapply::closepb(pb)
  names(selected_var) <- variable_name
  selected_var
}

#' Internal.compute_targetVal
#' 
#' This function is an internal function to compute target values
#'  
#' @return return computed target values
#' 
Internal.compute_targetVal <- function(input_col, targetVal_method = c("median", "mean"),
                                       targetVal_removeOutlier = TRUE) {
  
  input_col <- input_col[is.finite(input_col)]
  if (targetVal_removeOutlier) outlier_res <- grDevices::boxplot.stats(input_col, coef = 1.5)$out
  
  if (targetVal_method == "mean") {
    if (targetVal_removeOutlier) input_col[input_col %in% outlier_res] <- mean(input_col[!input_col %in% outlier_res], na.rm = TRUE)
    res <- mean(input_col, na.rm = TRUE)
  } else {
    if (targetVal_removeOutlier) input_col[input_col %in% outlier_res] <- median(input_col[!input_col %in% outlier_res], na.rm = TRUE)
    res <- median(input_col, na.rm = TRUE)
  }
  res
}

#' Internal.compute_errorRatio
#' 
#' This function is an internal function to compute error ratio
#'  
#' @return return computed error ratio
#' 
Internal.compute_errorRatio <- function(rawVal, sampleType,
                                        targetVal) {
  
  out <- sapply(1:length(rawVal), function(val_idx, rawVal, sampleType, targetVal) {
    set.seed(1)
    current_targetVal <- targetVal[row.names(targetVal) == sampleType[val_idx],] [[1]]
    current_rawVal <- rawVal[val_idx]
    errorRatio <- ifelse(current_targetVal == 0, 0, (current_rawVal - current_targetVal) / current_targetVal)
    out <- c(y = errorRatio, y_target = current_targetVal, y_raw = current_rawVal)
    out
  }, rawVal = rawVal, sampleType = sampleType, targetVal = targetVal)
  
  out_df <- data.frame(t(out))
  
  if (all(out_df$y == 0)) {
    out_df$y <- rnorm(length(out_df$y), sd = 0.001, mean = 0.0005)
  }
  out_df
}

#' compute_targetVal
#' 
#' This function computes the target values for ensemble learning architecture
#' 
#' @param QC_num a numeric data.frame including the metabolite values of quality control (QC) samples. Missing values and infinite values will not be taken into account. Row: sample. Column: metabolite variable. See Examples.
#' @param sampleType a vector corresponding to \code{QC_num} to specify the type of each QC sample. QC samples of the \strong{same type} should have the \strong{same type name}. See Examples.
#' @param batchID a vector corresponding to \code{QC_num} to specify the batch of each sample. Ignored if \code{targetVal_batchWise = FALSE}. See Examples.
#' @param targetVal_method a character string specifying how the target values are computed. Can be \code{"mean"} (default) or \code{"median"}. See Details.
#' @param targetVal_batchWise logical. If \code{TRUE}, the target values will be computed based on each batch, otherwise, based on the whole dataset. Setting \code{TRUE} might be useful if your dataset has very obvious batch effects, but this may also make the algorithm less robust. See Details. Default: \code{FALSE}.
#' @param targetVal_removeOutlier logical. If \code{TRUE}, outliers will be removed before the computation. Outliers are determined with 1.5 * IQR (interquartile range) rule. We recommend turning this off when the target values are computed based on batches. See Details. Default: \code{!targetVal_batchWise}.
#' @param coerce_numeric logical. If \code{TRUE}, values in \code{QC_num} will be coerced to numeric before the computation. The columns cannot be coerced will be removed (with warnings). See Examples. Default: \code{FALSE}.
#'  
#' @return the function returns a list of length one containing the target values computed on the whole dataset.
#' 
compute_targetVal <- function(QC_num, sampleType, batchID = NULL,
                              targetVal_method = c("mean", "median"),
                              targetVal_batchWise = FALSE,
                              targetVal_removeOutlier = !targetVal_batchWise,
                              coerce_numeric = FALSE) {
  
  message("+ Computing target values...   ", Sys.time())
  
  targetVal_method <- match.arg(targetVal_method)
  
  sampleType <- as.factor(sampleType)
  batchID    <- as.factor(batchID)
  
  if(coerce_numeric) {
    QC_num <- as.data.frame(sapply(QC_num, as.numeric))
    idx_NA <- sapply(QC_num, function(x) {
      all(is.na(x))
    })
    
    if (sum(idx_NA) > 0) {
      QC_num <- QC_num[,!idx_NA]
      warning("  ", sum(idx_NA), " column(s) in QC_num removed due to non-numeric values." )
    }
  } else {
    if (!all(sapply(QC_num, is.numeric))) stop("  The values in QC_num should be numeric!")
  }
  
  if (targetVal_batchWise) {
    target_values <- aggregate(QC_num, by = list(batch = batchID, sample = sampleType),
                               FUN = Internal.compute_targetVal, targetVal_method = targetVal_method,
                               targetVal_removeOutlier = targetVal_removeOutlier)
    batchID <- target_values$batch
    target_values$batch <- NULL
    target_values_list <- split(target_values, f = batchID)
  } else {
    target_values_list <- list(wholeDataset = aggregate(QC_num,
                                                        by = list(sample = sampleType),
                                                        FUN =  Internal.compute_targetVal,
                                                        targetVal_method = targetVal_method,
                                                        targetVal_removeOutlier = targetVal_removeOutlier))
  }
  
  target_values <- lapply(target_values_list, function(x) {
    row.names(x) <- x$sample
    x$sample <- NULL
    x
  })
  target_values
}

#' Internal.run_ensemble
#' 
#' This function is an internal function to run ensemble learning
#'  
#' @return return ensemble learning
#' 
Internal.run_ensemble <- function(trainSet, testSet,
                                  mtry_percent = seq(0.2, 0.8, 0.2),
                                  nodesize_percent = seq(0.2, 0.8, 0.2),
                                  ..., return_base_res = FALSE) {
  
  if (!is.null(mtry_percent)) mtry <- round(mtry_percent * (ncol(trainSet) - 3))
  if (!is.null(nodesize_percent)) nodesize <- round(nodesize_percent * nrow(trainSet))
  
  rf_hyperparams <- expand.grid(mtry = unique(mtry),
                                nodesize = unique(nodesize),
                                ... = ...)
  
  pred_ensemble <- lapply(1:nrow(rf_hyperparams), function(idx) {
    set.seed(1)
    current_hyperparams <- as.list(rf_hyperparams[idx,])
    
    folds_train <- caret::createFolds(1:length(trainSet$y), k = 5, returnTrain = TRUE)
    
    res_folds <- lapply(folds_train, function(train_idx, rf_params) {
      
      train_fold     <- trainSet[train_idx,]
      validate_fold  <- trainSet[-train_idx,]
      
      fold_formula <- c(formula = as.formula(y ~ .),
                        data = list(train_fold[!names(train_fold) %in% c("y_target", "y_raw")]),
                        current_hyperparams)
      RF_fold_mod  <- do.call(randomForest::randomForest, fold_formula)
      
      pred_fold <- predict(RF_fold_mod, validate_fold)
      
      pred_convert <- validate_fold$y_raw / (pred_fold + 1)
      pred_loss <- abs(pred_convert - validate_fold$y_target) / validate_fold$y_target
      pred_loss
    })
    
    mean_loss <- mean(unlist(res_folds), na.rm = TRUE)
    mod_weight <- 1/exp(mean_loss)
    
    base_formula <- c(formula = as.formula(y ~ .), data = list(trainSet[!names(trainSet) %in% c("y_target", "y_raw")]),
                      current_hyperparams)
    RF_base_mod <- do.call(randomForest::randomForest, base_formula)
    
    pred_test   <- predict(RF_base_mod, testSet)
    pred_test_convert <- testSet$y_raw / (pred_test + 1)
    
    if (return_base_res) {
      out <- list(mod_weight = mod_weight, pred_test_convert = pred_test_convert, RF_base_mod = RF_base_mod)
    } else {
      out <- list(mod_weight = mod_weight, pred_test_convert = pred_test_convert)
    }
    out
  })
  mod_weights <- sapply(pred_ensemble, function(x) x$mod_weight)
  mod_weights_norm <- mod_weights / sum(mod_weights, na.rm = TRUE)
  
  pred_test <- sapply(pred_ensemble, function(x) x$pred_test_convert)
  pred_norm <- apply(pred_test, 1, function(x) sum(x * mod_weights_norm, na.rm = TRUE))
  
  if (return_base_res) {
    output <- list(pred_norm = pred_norm, base_res = pred_ensemble, base_weights = mod_weights_norm, rf_hyperparams = rf_hyperparams)
  } else {
    output <- pred_norm
  }
  output
}

#' run_TIGER_internal
#' 
#' Use TIGER algorithm to eliminate the technical variation in metabolomics data. TIGER supports targeted and untargeted metabolomics data and is competent to perform both intra- and inter-batch technical variation removal.
#' 
#' @param test_samples (required) a data.frame containing the samples to be corrected (for example, subject samples). This data.frame should contain columns of
#' \itemize{
#' \item sample ID (required): name or label for each sample,
#' \item sample type (required): indicating the type of each sample,
#' \item batch ID (required): the batch of each sample,
#' \item order information (optional): injection order or temporal information of each sample,
#' \item position information (optional): well position of each sample,
#' \item metabolite values (required): values to be normalised. Infinite values are not allowed.
#' }
#' Row: sample. Column: variable. See Examples.
#' @param train_samples (required) a data.frame containing the quality control (QC) samples used for model training. The columns in this data.frame should correspond to the columns in \code{test_samples}. And \code{test_samples} and \code{train_samples} should have the identical column names.
#' @param col_sampleID  (required) a character string indicating the name of the column that specifies the sample ID of each sample. The values in this column will not affect the data correction process but can act as labels for different samples. See Examples.
#' @param col_sampleType (required) a character string indicating the name of the column that specifies the type (such as QC1, QC2, subject) of each sample. This column can be used to indicate different kinds of QC samples in \code{train_samples}. QC samples of the \strong{same type} should have the \strong{same type name}. See Examples.
#' @param col_batchID (required) a character string indicating the name of the column that specifies the batch ID of each sample. See Examples.
#' @param col_order (optional) \code{NULL} or a character string indicating the name of the column that contains the injection order or temporal information (numeric values). This can explicitly ask the algorithm capture the technical variation introduced by injection order, which might be useful when your data have very obvious temporal drifts. If \code{NULL} (default), \code{train_samples} and \code{test_samples} should have \strong{No} column contains injection order information.
#' @param col_position (optional) \code{NULL} or a character string indicating the name of the column that contains the well position information (numeric values). This can explicitly ask the algorithm capture the technical variation introduced by well position, which might be useful when the well position has a great impact during data acquisition. If \code{NULL} (default), \code{train_samples} and \code{test_samples} should have \strong{No} column contains well position information.
#' @param targetVal_external (optional) a list generated by function \code{\link{compute_targetVal}}. See Details.
#' @param targetVal_method a character string specifying how target values are to be computed. Can be \code{"mean"} (default) or \code{"median"}. Ignored if a list of external target values has been assigned to  \code{targetVal_external}.
#' @param targetVal_batchWise logical. If \code{TRUE}, the target values will be computed based on each batch, otherwise, based on the whole dataset. Setting \code{TRUE} might be useful if your dataset has very obvious batch effects, but this may also make the algorithm less robust. Default: \code{FALSE}. Ignored if a list of external target values has been assigned to  \code{targetVal_external}.
#' @param targetVal_removeOutlier logical. If \code{TRUE}, outliers will be removed before the computation. Outliers are determined with 1.5 * IQR (interquartile range) rule. We recommend turning this off when the target values are computed based on batches. Default: \code{!targetVal_batchWise}. Ignored if a list of external target values has been assigned to  \code{targetVal_external}.
#' @param selectVar_external (optional) a list generated by function \code{\link{select_variable}}. See Details.
#' @param selectVar_corType a character string indicating correlation (\code{"cor"}, default) or partial correlation (\code{"pcor"}) is to be used. Can be abbreviated. Ignored if a list of selected variables has been assigned to \code{selectVar_external}. \strong{Note}: computing partial correlations of a large dataset can be very time-consuming.
#' @param selectVar_corMethod a character string indicating which correlation coefficient is to be computed. One of \code{"spearman"} (default) or \code{"pearson"}. Can be abbreviated. Ignored if a list of selected variables has been assigned to \code{selectVar_external}.
#' @param selectVar_minNum an integer specifying the minimum number of the selected metabolite variables (injection order and well position are not regarded as metabolite variables). If \code{NULL}, no limited, but 1 at least. Default: \code{5}. Ignored if a list of selected variables has been assigned to \code{selectVar_external}.
#' @param selectVar_maxNum an integer specifying the maximum number of the selected metabolite variables (injection order and well position are not regarded as metabolite variables). If \code{NULL}, no limited, but no more than the number of all available metabolite variables. Default: \code{10}. Ignored if a list of selected variables has been assigned to \code{selectVar_external}.
#' @param selectVar_batchWise (advanced) logical. Specify whether the variable selection should be performed based on each batch. Default: \code{FALSE}. Ignored if a list of selected variables has been assigned to \code{selectVar_external}. \strong{Note}: the support of batch-wise variable selection is provided for data requiring special processing (for example, data with strong batch effects). But in most case, batch-wise variable selection is not necessary. Setting \code{TRUE} can make the algorithm less robust.
#' @param mtry_percent (advanced) a numeric vector indicating the percentages of selected variables randomly sampled as candidates at each split when training random forest models (base learners). \strong{Note}: providing more arguments will include more base learners into the ensemble model, which will increase the processing time. Default: \code{seq(0.2, 0.8, 0.2)}.
#' @param nodesize_percent (advanced) a numeric vector indicating the percentages of sample size used as the minimum sizes of the terminal nodes in random forest models (base learners). \strong{Note}: providing more arguments will include more base learners into the ensemble model, which will increase the processing time. Default: \code{seq(0.2, 0.8, 0.2)}.
#' @param ... (advanced) optional arguments (except \code{mtry} and \code{nodesize}) to be passed to \code{\link[randomForest]{randomForest}} for model training. Arguments \code{mtry} and \code{nodesize} are determined by \code{mtry_percent} and \code{nodesize_percent}. See \code{\link[randomForest]{randomForest}} and Examples. \strong{Note}: providing more arguments will include more base learners into the ensemble model, which will increase the processing time.
#' @param parallel.cores an integer (== -1 or >= 1) specifying the number of cores for parallel computation. Setting \code{-1} to run with all cores. Default: \code{2}.
#'  
#' @return return tiger batch corrected results
#' 
run_TIGER_internal <- function(test_samples, train_samples,
                       col_sampleID, col_sampleType, col_batchID,
                       col_order = NULL, col_position = NULL,
                       
                       targetVal_external = NULL, targetVal_method = c("mean", "median"),
                       targetVal_batchWise = FALSE, targetVal_removeOutlier = !targetVal_batchWise,
                       
                       selectVar_external = NULL, selectVar_corType = c("cor", "pcor"),
                       selectVar_corMethod = c("pearson", "spearman"),
                       selectVar_minNum = 5, selectVar_maxNum = 10,
                       selectVar_batchWise = FALSE,
                       
                       mtry_percent = seq(0.2, 0.8, 0.2),
                       nodesize_percent = seq(0.2, 0.8, 0.2),
                       ..., parallel.cores = 2) {
  require(caret)
  
  message("+ Initialising...   ", Sys.time())
  
  selectVar_corType   <- match.arg(selectVar_corType)
  selectVar_corMethod <- match.arg(selectVar_corMethod)
  targetVal_method    <- match.arg(targetVal_method)
  
  for (col_idx in c(col_sampleID, col_sampleType, col_batchID)) {
    test_samples[[col_idx]]  <- as.character(test_samples[[col_idx]])
    train_samples[[col_idx]] <- as.character(train_samples[[col_idx]])
  }
  
  if (!is.null(col_order)) {
    if (anyNA(test_samples[[col_order]])  | any(!is.finite(test_samples[[col_order]]))  ) stop("  test samples: col_order should be numeric only!")
    if (anyNA(train_samples[[col_order]]) | any(!is.finite(train_samples[[col_order]])) ) stop("  train samples: col_order should be numeric only!")
  }
  
  if (!is.null(col_position)) {
    if (anyNA(test_samples[[col_position]])  | any(!is.finite(test_samples[[col_position]]))  ) stop("  test samples: col_position should be numeric only!")
    if (anyNA(train_samples[[col_position]]) | any(!is.finite(train_samples[[col_position]])) ) stop("  train samples: col_position should be numeric only!")
  }
  
  if (!all(sapply(train_samples[!names(train_samples) %in% c(col_sampleID, col_sampleType, col_batchID)], is.numeric))) {
    stop("  The values of train samples (except sampleType and batchID) should be numeric!")
  }
  
  if (!all(sapply(test_samples[!names(test_samples) %in%c(col_sampleID, col_sampleType, col_batchID)], is.numeric))) {
    stop("  The values of test samples (except sampleType and batchID) should be numeric!")
  }
  
  # To check infinite values. - Deprecated. Infinite values are not allowed.
  # train_samples_list  <- split(train_samples, f = train_samples[[col_sampleType]])
  # train_samples_check <- lapply(train_samples_list, Internal.impute_infinite)
  # train_samples <- do.call("rbind", train_samples_check)
  # test_samples  <- Internal.impute_infinite(test_samples)
  
  test_samples_bak <- test_samples
  
  batchID_train <- unique(train_samples[[col_batchID]])
  batchID_test  <- unique(test_samples[[col_batchID]])
  
  if(!all(batchID_test %in% batchID_train)) stop("  The batchID in train and test samples cannot match!")
  batchID <- batchID_test
  
  train_num <- train_samples[!names(train_samples) %in% c(col_sampleID, col_sampleType, col_batchID, col_order, col_position)]
  test_num  <- test_samples[!names(test_samples)   %in% c(col_sampleID, col_sampleType, col_batchID, col_order, col_position)]
  
  var_names <- names(test_num)
  
  # Target value computation
  if (is.null(targetVal_external)) {
    targetVal_list <- compute_targetVal(QC_num = train_num,
                                        sampleType = train_samples[[col_sampleType]],
                                        batchID    = train_samples[[col_batchID]],
                                        targetVal_method = targetVal_method,
                                        targetVal_batchWise = targetVal_batchWise,
                                        targetVal_removeOutlier = targetVal_removeOutlier,
                                        coerce_numeric = FALSE)
  } else {
    message("+ External target values loaded.   ", Sys.time())
    targetVal_list <- targetVal_external
    if (length(targetVal_list) == 1 & names(targetVal_list[1]) == "wholeDataset") {
      targetVal_batchWise <- FALSE
    } else {
      if (all(names(targetVal_list) %in% batchID_train)) {
        targetVal_batchWise <- TRUE
      } else stop("  Batch ID of targetVal_external cannot match your training data!")
    }
  }
  
  # Variable selection
  if (is.null(selectVar_external)) {
    var_selected_list <- select_variable(train_num = train_num, test_num = test_num,
                                         train_batchID = train_samples[[col_batchID]],
                                         test_batchID  = test_samples[[col_batchID]],
                                         selectVar_batchWise = selectVar_batchWise,
                                         selectVar_corType   = selectVar_corType,
                                         selectVar_corMethod = selectVar_corMethod,
                                         selectVar_minNum = selectVar_minNum,
                                         selectVar_maxNum = selectVar_maxNum,
                                         coerce_numeric = FALSE)
  } else {
    message("+ External selected variables loaded.   ", Sys.time())
    var_selected_list <- selectVar_external
    if (length(var_selected_list) == 1 & names(var_selected_list[1]) == "wholeDataset") {
      selectVar_batchWise <- FALSE
    } else {
      if (all(names(var_selected_list) %in% batchID_test)) {
        selectVar_batchWise <- TRUE
      } else stop("  Batch ID of selectVar_external cannot match your test data!")
    }
  }
  
  idx_test_na  <- is.na(test_samples)
  idx_train_na <- is.na(train_samples)
  if (any(idx_test_na))  test_samples[idx_test_na]   <- 0
  if (any(idx_train_na)) train_samples[idx_train_na] <- 0
  
  message("+ Data correction started.   ", Sys.time())
  parallel.cores <- ifelse(parallel.cores == -1, parallel::detectCores(), parallel.cores)
  cl <- parallel::makeCluster(parallel.cores)
  parallel::clusterExport(cl = cl, varlist = c("Internal.compute_errorRatio", "Internal.run_ensemble","createFolds"), envir = environment())
  pbapply::pboptions(type = "timer", style = 3, char = "=", txt.width = 70)
  
  # Original sample order backup
  test_samples <- cbind(original_idx = 1:nrow(test_samples), test_samples)
  
  res_var <- pbapply::pblapply(var_names, function(current_var, var_selected_list, targetVal_list,
                                                   targetVal_batchWise, selectVar_batchWise,
                                                   train_samples, test_samples, col_sampleID, col_sampleType,
                                                   col_batchID, col_order, col_position, batchID, mtry_percent,
                                                   targetVal_method, nodesize_percent, ...) {
    if (!targetVal_batchWise) {
      train_y_all <- Internal.compute_errorRatio(rawVal     = train_samples[[current_var]],
                                                 sampleType = train_samples[[col_sampleType]],
                                                 targetVal  = targetVal_list$wholeDataset[current_var])
      
    }
    
    test_data <- cbind(y_raw = test_samples[[current_var]], test_samples)
    
    res_batch_list <- lapply(batchID, function(current_batch) {
      
      train_X_batch <- train_samples[train_samples[[col_batchID]] == current_batch,]
      
      if (targetVal_batchWise) {
        train_y_batch <- Internal.compute_errorRatio(rawVal     = train_X_batch[[current_var]],
                                                     sampleType = train_X_batch[[col_sampleType]],
                                                     targetVal  = targetVal_list[[current_batch]][current_var])
        
        
      } else {
        train_y_batch <- train_y_all[train_samples[[col_batchID]] == current_batch,]
      }
      
      if (selectVar_batchWise) {
        train_X_selected <- train_X_batch[ c(col_order, col_position,
                                             var_selected_list[[current_batch]][[current_var]]) ]
      } else {
        train_X_selected <- train_X_batch[ c(col_order, col_position,
                                             var_selected_list$wholeDataset[[current_var]]) ]
      }
      
      trainSet <- cbind(train_y_batch, train_X_selected)
      testSet  <- test_data[test_data[[col_batchID]] == current_batch,]
      
      var_pred <- Internal.run_ensemble(trainSet = trainSet, testSet = testSet,
                                        mtry_percent = mtry_percent,
                                        nodesize_percent = nodesize_percent,
                                        ... = ..., return_base_res = FALSE)
      
      if (targetVal_batchWise) {
        test_targetVal_all   <- do.call(targetVal_method, list(test_data$y_raw, na.rm = TRUE))
        test_targetVal_batch <- do.call(targetVal_method, list(testSet$y_raw,   na.rm = TRUE))
        var_pred <- var_pred * test_targetVal_all / test_targetVal_batch
      }
      
      names(var_pred) <- testSet$original_idx
      var_pred
    })
    
    res_batch       <- do.call("c", res_batch_list)
    res_batch_order <- res_batch[order(as.numeric(names(res_batch)))]
    res_batch_df    <- data.frame(res_batch_order)
    
    names(res_batch_df) <- current_var
    res_batch_df
    
  }, var_selected_list = var_selected_list, targetVal_list = targetVal_list,
  targetVal_batchWise = targetVal_batchWise, selectVar_batchWise = selectVar_batchWise,
  train_samples = train_samples, test_samples = test_samples, col_sampleID = col_sampleID,
  col_sampleType = col_sampleType, col_batchID = col_batchID, col_order = col_order,
  col_position = col_position, batchID = batchID, mtry_percent = mtry_percent,
  nodesize_percent = nodesize_percent, targetVal_method = targetVal_method, ... = ..., cl = cl)
  
  parallel::stopCluster(cl)
  
  message("  - Merging results...")
  check_order <- sapply(res_var[-1], function(x) {
    any(row.names(x) != row.names(res_var[[1]]))
  })
  if (any(check_order)) stop("    Error occurs when merging data!")
  
  res_var_df <- do.call("cbind", res_var)
  
  test_samples[names(res_var_df)] <- res_var_df
  test_samples$original_idx <- NULL
  
  idx_norm_na <- is.na(test_samples)
  if (any(idx_norm_na)) test_samples[idx_norm_na] <- as.numeric(test_samples_bak[idx_norm_na])
  
  idx_norm_zero <- test_samples < 0
  if (any(idx_norm_zero)) test_samples[idx_norm_zero] <- as.numeric(test_samples_bak[idx_norm_zero])
  
  test_samples[test_samples_bak == 0] <- 0
  if (any(idx_test_na)) test_samples[idx_test_na] <- NA
  
  message("+ Completed.   ", Sys.time())
  test_samples
}

#' compute_RSD
#' 
#' This function computes the relative standard deviation (RSD)
#'  
#' @return return RSD
#' 
compute_RSD <- function(input_data) {
  val_RSD <- sd(input_data, na.rm = TRUE) / mean(input_data, na.rm = TRUE)
  val_RSD
}

#' select_variable
#' 
#' This function provides an advanced option to select metabolite variables from external dataset(s). The selected variables (as a list) can be further passed to argument \code{selectVar_external} in function \code{\link{run_TIGER}} for a customised data correction.
#' 
#' @param train_num a numeric data.frame \strong{only} including the metabolite values of training samples (can be quality control samples). Information such as injection order or well position need to be excluded. Row: sample. Column: metabolite variable. See Examples.
#' @param test_num an optional numeric data.frame including the metabolite values of test samples (can be subject samples). If provided, the column names of \code{test_num} should correspond to the column names of \code{train_num}. Row: sample. Column: metabolite variable. If \code{NULL}, the variables will be selected based on \code{train_num} only. See Examples.
#' @param train_batchID \code{NULL} or a vector corresponding to \code{train_num} to specify the batch of each sample. Ignored if \code{selectVar_batchWise = FALSE}. See Examples.
#' @param test_batchID \code{NULL} or a vector corresponding to \code{test_num} to specify the batch of each sample. Ignored if \code{selectVar_batchWise = FALSE}. See Examples.
#' @param selectVar_corType a character string indicating correlation (\code{"cor"}, default) or partial correlation (\code{"pcor"}) is to be used. Can be abbreviated. See Details. \strong{Note}: computing partial correlations of a large dataset can be very time-consuming.
#' @param selectVar_corMethod a character string indicating which correlation coefficient is to be computed. One of \code{"spearman"} (default) or \code{"pearson"}. Can be abbreviated. See Details.
#' @param selectVar_minNum an integer specifying the minimum number of the selected variables. If \code{NULL}, no limited, but 1 at least. See Details. Default: 5.
#' @param selectVar_maxNum an integer specifying the maximum number of the selected variables. If \code{NULL}, no limited, but \code{ncol(train_num) - 1} at most. See Details. Default: 10.
#' @param selectVar_batchWise (advanced) logical. Specify whether the variable selection should be performed based on each batch. Default: \code{FALSE}. \strong{Note}: if \code{TRUE}, batch ID of each sample are required. The support of batch-wise variable selection is provided for data requiring special processing (for example, data with strong batch effects). But in most case, batch-wise variable selection is not necessary. Setting \code{TRUE} might make the algorithm less robust. See Details.
#' @param coerce_numeric logical. If \code{TRUE}, values in \code{train_num} and  \code{test_num} will be coerced to numeric before the computation. The columns cannot be coerced will be removed (with warnings). See Examples. Default: \code{FALSE}.
#'  
#' @return the function returns a list of length one containing the selected variables computed on the whole dataset.
#' 
select_variable <- function(train_num, test_num = NULL,
                            train_batchID = NULL, test_batchID = NULL,
                            selectVar_corType   = c("cor", "pcor"),
                            selectVar_corMethod = c("spearman", "pearson"),
                            selectVar_minNum = 5, selectVar_maxNum = 10,
                            selectVar_batchWise = FALSE,
                            coerce_numeric = FALSE) {
  
  message("+ Selecting highly-correlated variables...   ", Sys.time())
  
  selectVar_corType   <- match.arg(selectVar_corType)
  selectVar_corMethod <- match.arg(selectVar_corMethod)
  
  if(coerce_numeric) {
    train_num <- as.data.frame(sapply(train_num, as.numeric))
    idx_NA <- sapply(train_num, function(x) {
      all(is.na(x))
    })
    
    if (sum(idx_NA) > 0) {
      train_num <- train_num[,!idx_NA]
      warning("  ", sum(idx_NA), " column(s) in train_num removed due to non-numeric values." )
    }
    
    if (!is.null(test_num)) {
      test_num <- as.data.frame(sapply(test_num, as.numeric))
      idx_NA <- sapply(test_num, function(x) {
        all(is.na(x))
      })
      if (sum(idx_NA) > 0) {
        test_num <- test_num[,!idx_NA]
        warning("  ", sum(idx_NA), " column(s) in test_num removed due to non-numeric values." )
      }
    }
  } else {
    if (!all(sapply(train_num, is.numeric))) stop("  The values in train_num should be numeric!")
    if (!is.null(test_num) & !all(sapply(test_num, is.numeric))) stop("  The values in test_num should be numeric!")
  }
  
  if (!is.null(test_num)) {
    if (!all(names(test_num) %in% names(train_num))) stop("  Variables in training and test data cannot match!")
    train_num <- train_num[names(test_num)]
  }
  
  selectVar_minNum <- ifelse(is.null(selectVar_minNum), 1, selectVar_minNum)
  selectVar_maxNum <- ifelse(is.null(selectVar_maxNum), (ncol(train_num) - 1), selectVar_maxNum)
  
  selectVar_minNum <- max(as.integer(selectVar_minNum), 1)
  selectVar_maxNum <- max(as.integer(selectVar_maxNum), 1)
  
  selectVar_maxNum <- min(selectVar_maxNum, ncol(train_num) - 1)
  selectVar_minNum <- min(selectVar_minNum, selectVar_maxNum)
  
  if (selectVar_batchWise) {
    train_num_list    <- split(train_num, f = train_batchID)
    batch_names_train <- sort(names(train_num_list))
    batch_names_len   <- length(batch_names_train)
    
    if (!is.null(test_num)) {
      test_num_list     <- split(test_num,  f = test_batchID)
      batch_names_test  <- sort(names(test_num_list))
      
      if (any(batch_names_train != batch_names_test) | batch_names_len != length(batch_names_test)) stop("  Batch names of train and test samples cannot match: their column names should be identical!")
    } else test_num_list <- NULL
    
    selected_var_list <- lapply(1:batch_names_len, function(batch_name_idx) {
      
      one_batch_name <- batch_names_train[[batch_name_idx]]
      message("  - Current batch: ", one_batch_name, "  (", batch_name_idx, "/", batch_names_len, ")")
      
      message("      Computing correlation coefficients...")
      cor_info <- Internal.compute_cor(train_num = train_num_list[[one_batch_name]],
                                       test_num = test_num_list[[one_batch_name]],
                                       selectVar_corType = selectVar_corType,
                                       selectVar_corMethod = selectVar_corMethod)
      
      message("      Selecting variables...")
      selected_var <- Internal.select_variable(cor_info = cor_info,
                                               selectVar_minNum = selectVar_minNum,
                                               selectVar_maxNum = selectVar_maxNum)
    })
    names(selected_var_list) <- batch_names_train
  } else {
    message("  - Computing correlation coefficients...")
    cor_info <- Internal.compute_cor(train_num = train_num, test_num = test_num,
                                     selectVar_corType = selectVar_corType,
                                     selectVar_corMethod = selectVar_corMethod)
    message("  - Selecting variables...")
    selected_var <- Internal.select_variable(cor_info = cor_info,
                                             selectVar_minNum = selectVar_minNum,
                                             selectVar_maxNum = selectVar_maxNum)
    selected_var_list <- list(wholeDataset = selected_var)
  }
  
  selected_var_list
}

# Internal.boxplot.stats <- function(x, # borrowed from grDevices::boxplot.stats()
#                                    coef = 1.5, do.conf = TRUE, do.out = TRUE) {
#   if (coef < 0)
#     stop("'coef' must not be negative")
#   nna <- !is.na(x)
#   n <- sum(nna)
#   stats <- stats::fivenum(x, na.rm = TRUE)
#   iqr <- diff(stats[c(2, 4)])
#   names(iqr) <- NULL
#   if (coef == 0)
#     do.out <- FALSE
#   else {
#     lower_limit <- (stats[2L] - coef * iqr)
#     upper_limit <- (stats[4L] + coef * iqr)
#     names(lower_limit) <- NULL
#     names(upper_limit) <- NULL
#     out <- if (!is.na(iqr)) {
#       x < lower_limit | x > upper_limit
#     }
#     else !is.finite(x)
#     if (any(out[nna], na.rm = TRUE))
#       stats[c(1, 5)] <- range(x[!out], na.rm = TRUE)
#   }
#   conf <- if (do.conf)
#     stats[3L] + c(-1.58, 1.58) * iqr/sqrt(n)
#   list(stats = stats, n = n, conf = conf,
#        out = if (do.out) x[out & nna] else numeric(),
#        iqr = iqr, lower_limit = lower_limit, upper_limit = upper_limit)
# }

#' Internal.remove_NA
#' 
#' This function is an internal function to remove NA values
#'
#' @return data with NA values removed
#' 
Internal.remove_NA <- function(input_data_num, data_label = NULL,
                               replace_with_zero = FALSE) {
  data_na_sample_idx <- apply(input_data_num, 1, function(x) any(is.na(x)))
  data_na_sample_sum <- sum(data_na_sample_idx)
  
  if (data_na_sample_sum > 0) {
    if (!is.null(data_label)) warning(paste0("    ", data_na_sample_sum, " sample(s) in ", data_label, " contain(s) NA."))
    
    if (replace_with_zero) {
      input_data_num[is.na(input_data_num)] <- 0
    } else {
      input_data_num <- input_data_num[!data_na_sample_idx,]
      if (nrow(input_data_num) == 0) stop(paste0("    Variable selection failed: ", data_label, " contain too many NA!"))
    }
  }
  
  input_data_num
}
