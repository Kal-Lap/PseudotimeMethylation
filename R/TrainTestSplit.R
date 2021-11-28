#' @title Stratified train-validation-test split by ratio
#'
#' @description   Split data into train-val-test set, stratified by age bin.
#' @param data_list a list of methylation and meta data data frame respectively
#' @param train_size a numeric value of ratio to split into train set
#' @param val_size (optional) a numeric value of ratio to split into val set
#' @param test_size a numeric value of ratio to split into the test set
#' @param age_bin_num number of binnings of chronological age for stratification
#' @return a list of 4 or 6
#' \describe{
#'   \itemize{
#'     \item train_meth : training data
#'     \item test_meth: test data
#'     \item val_meth (optional) : validation data
#'     \item train_meta : training meta data
#'     \item test_meta : test meta data
#'     \item val_meta (optional): validation meta data
#'   }
#' }
#' @export
#' @examples
#' \dontrun{ratio_train_val_test_split(list("meth_data" = meth_data, "meta_data" = meta_data),
#'                            train_size = 0.2,
#'                            test_size = 0.6, val_size = 0.1,
#'                            age_bin_num = 25)
#'}
ratio_train_val_test_split <- function(data_list,
                                       train_size,
                                       test_size,
                                       val_size = 0.0,
                                       age_bin_num = 20){


  if (test_size < 0 | train_size <0 | val_size <0){
    stop("Negative size is invalid.")
  }
  if (test_size > 1 | train_size > 1 | val_size > 1){
    stop("Size larger than 1 is invalid")
  }
  if (train_size + val_size + test_size !=1.0) {
    stop("Sum of train, test, and validation size does not add up to 1.")
  }
  if (test_size == 0.0){
    stop("Please select the size for test set")
  }
  if (train_size == 0.0){
    stop("Please select the size for train set")
  }

  meta <- data_list$meta_data
  meth <- data_list$meth_data


  if (val_size == 0.0){
    inds <- splitTools::partition(meta$chrono_age,
                                  shuffle = TRUE,
                                  n_bins = age_bin_num,
                                  p = c(train = train_size, test = test_size))
    split_data <- list(train_meth = meth[,inds$train],
                       test_meth = meth[,inds$test],
                       train_meta = meta[inds$train,],
                       test_meta = meta[inds$test,])
  } else{
    inds <- splitTools::partition(meta$chrono_age,
                                  shuffle = TRUE,
                                  p = c(train = train_size,
                                        valid = val_size,
                                        test = test_size))

    split_data <- list(train_meth = meth[,inds$train],
                       test_meth = meth[,inds$test],
                       val_meth = meth[,inds$valid],
                       train_meta = meta[inds$train,],
                       test_meta = meta[inds$test,],
                       val_meta = meta[inds$valid,])

  }

  return(split_data)
}


#' @title  Train-validation-test split by experiment label
#'
#' @description   Split data into train-val-test set by given experiment label
#' @param data_list a list of methylation and meta data data frame respectively
#' @param train_exp a character vector of experiments in train set
#' @param val_exp (optional) a character vector of experiments in val set
#' @param test_exp a character vector of experiments in test set
#' @return a list of 4 or 6
#' \describe{
#'   \itemize{
#'     \item train_meth : training data
#'     \item test_meth: test data
#'     \item val_meth (optional) : validation data
#'     \item train_meta : training meta data
#'     \item test_meta : test meta data
#'     \item val_meta (optional): validation meta data
#'   }
#' }
#' @export
#' @examples
#' \dontrun{exp_train_val_test_split(list("meth_data" = meth_data, "meta_data" = meta_data),
#'                          train_exp = c('GSE482504', 'GSE235456'),
#'                          test_exp = c('GSE621336','GSE125673'),
#'                          val_exp = c('GSE34354'))}
#'
exp_train_val_test_split <- function(data_list,
                                     train_exp,
                                     test_exp,
                                     val_exp){
  meta <- data_list$meta_data
  meth <- data_list$meth_data

  if (any(is.na(meta$exp_labels))){
    stop("Missing value in experiment labels in meta data")
  }
  if (missing(val_exp)){
    train_ind = which(meta$exp_labels %in% train_exp)
    test_ind = which(meta$exp_labels %in% test_exp)

    split_data <- list(train_meth = meth[,train_ind],
                       test_meth = meth[,test_ind],
                       train_meta = meta[train_ind,],
                       test_meta = meta[test_ind,])
  }else{
    train_ind = which(meta$exp_labels %in% train_exp)
    test_ind = which(meta$exp_labels %in% test_exp)
    val_ind = which(meta$exp_labels %in% val_exp)

    split_data <- list(train_meth = meth[,train_ind],
                       test_meth = meth[,test_ind],
                       val_meth = meth[,val_ind],
                       train_meta = meta[train_ind,],
                       test_meta = meta[test_ind,],
                       val_meta = meta[val_ind,])
  }

  return(split_data)
}

#' @title  Train-validation-test split by train-test label
#'
#' @description   Split data into train-val-test set by train-test label
#' @param data_list a list of methylation and meta data data frame respectively
#' @param train_label a string keyword of train set label
#' @param val_label (optional) a string keyword of val set label
#' @param test_label a string keyword of test set label
#' @return a list of 4 or 6
#' \describe{
#'   \itemize{
#'     \item train_meth : training data
#'     \item test_meth: test data
#'     \item val_meth (optional) : validation data
#'     \item train_meta : training meta data
#'     \item test_meta : test meta data
#'     \item val_meta (optional): validation meta data
#'   }
#' }
#' @export
#' @examples
#' \dontrun{labels_train_test_val_split(list("meth_data" = meth_data, "meta_data" = meta_data),
#'                          train_label = "train",
#'                          test_label = "test",
#'                          val_label = "val")
#'
#' labels_train_test_val_split(list("meth_data" = meth_data, "meta_data" = meta_data),
#'                          train_label = "train",
#'                          test_label = "test")}
#'
labels_train_test_val_split<- function(data_list,
                                       train_label,
                                       test_label,
                                       val_label){
  meta <- data_list$meta_data
  meth <- data_list$meth_data

  if (any(is.na(meta$train_test_labels))){
    stop("Missing value in train-test labels in meta data")
  }
  if (missing(val_label)){

    train_ind = which(meta$train_test_labels %in% c(train_label))
    test_ind = which(meta$train_test_labels %in% c(test_label))

    split_data <- list(train_meth = meth[,train_ind],
                       test_meth = meth[,test_ind],
                       train_meta = meta[train_ind,],
                       test_meta = meta[test_ind,])
  }else{
    train_ind = which(meta$train_test_labels %in% c(train_label))
    test_ind = which(meta$train_test_labels %in% c(test_label))
    val_ind = which(meta$train_test_labels %in% c(val_label))

    split_data <- list(train_meth = meth[,train_ind],
                       test_meth = meth[,test_ind],
                       val_meth = meth[,val_ind],
                       train_meta = meta[train_ind,],
                       test_meta = meta[test_ind,],
                       val_meta = meta[val_ind,])
  }

  return(split_data)
}
