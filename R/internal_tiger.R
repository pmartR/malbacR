#' Internal.createFolds
#' 
#' This function is an internal function to create folds
#' 
#' @param y description
#' @param k description
#' @param list description
#' @param returnTrain description
#'  
#' @return return created folds
#' 
#' @export
#' 
Internal.createFolds <- function (y, k = 10, list = TRUE, returnTrain = FALSE) {
  # borrowed from caret::createFolds()
  # package caret has been cited in our original paper
  if (is.numeric(y)) {
    cuts <- floor(length(y)/k)
    if (cuts < 2)
      cuts <- 2
    if (cuts > 5)
      cuts <- 5
    breaks <- unique(stats::quantile(y, probs = seq(0, 1, length = cuts)))
    y <- cut(y, breaks, include.lowest = TRUE)
  }
  if (k < length(y)) {
    y <- factor(as.character(y))
    numInClass <- table(y)
    foldVector <- vector(mode = "integer", length(y))
    for (i in 1:length(numInClass)) {
      min_reps <- numInClass[i]%/%k
      if (min_reps > 0) {
        spares <- numInClass[i]%%k
        seqVector <- rep(1:k, min_reps)
        if (spares > 0)
          seqVector <- c(seqVector, sample(1:k, spares))
        foldVector[which(y == names(numInClass)[i])] <- sample(seqVector)
      }
      else {
        foldVector[which(y == names(numInClass)[i])] <- sample(1:k,
                                                               size = numInClass[i])
      }
    }
  }
  else foldVector <- seq(along = y)
  if (list) {
    out <- split(seq(along = y), foldVector)
    names(out) <- paste("Fold", gsub(" ", "0",
                                     format(seq(along = out))), sep = "")
    if (returnTrain)
      out <- lapply(out, function(data, y) y[-data], y = seq(along = y))
  }
  else out <- foldVector
  out
}

#' Internal.compute_cor
#' 
#' This function is an internal function to compute correlation
#' 
#' @param train_num description
#' @param test_num description
#' @param selectVar_corType description
#' @param selectVar_corMethod description
#'  
#' @return return created folds
#' 
#' @export
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

Internal.compute_targetVal <- function(input_col, targetVal_method = c("median", "mean"),
                                       targetVal_removeOutlier = TRUE) {
  
  input_col <- input_col[is.finite(input_col)]
  if (targetVal_removeOutlier) outlier_res <- Internal.boxplot.stats(input_col, coef = 1.5)$out
  
  if (targetVal_method == "mean") {
    if (targetVal_removeOutlier) input_col[input_col %in% outlier_res] <- mean(input_col[!input_col %in% outlier_res], na.rm = TRUE)
    res <- mean(input_col, na.rm = TRUE)
  } else {
    if (targetVal_removeOutlier) input_col[input_col %in% outlier_res] <- median(input_col[!input_col %in% outlier_res], na.rm = TRUE)
    res <- median(input_col, na.rm = TRUE)
  }
  res
}

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

#' Internal.run_ensemble
#' 
#' This function is an internal function to run ensembl learning
#' 
#' @param trainSet description
#' @param testSet description
#' @param mtry_percent description
#' @param nodesize_percent description
#' @param return_base_res description
#'  
#' @return return created folds
#' 
#' @export
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
    
    folds_train <- Internal.createFolds2(1:length(trainSet$y), k = 5, returnTrain = TRUE)
    
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
#' This function is an internal function to run ensembl learning
#' 
#' @param test_samples description
#' @param train_samples description
#' @param col_sampleID description
#' @param col_sampleType description
#' @param col_batchID description
#' @param col_order description
#' @param col_position description
#'  
#' @return return tiger batch corrected results
#' 
#' @export
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
  parallel::clusterExport(cl = cl, varlist = c("Internal.compute_errorRatio", "Internal.run_ensemble","Internal.createFolds"), envir = environment())
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