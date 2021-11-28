#' @title Check valid input
#'
#' @description  Check dimension and data type of input data frames if valid
#' @param meth_data a data matrix
#' \describe{ Each column contains methylation data of each sample.
#' Each row contains methylation data at each site.
#' }
#' @param meta_data a data frame
#' \describe{Each row contains meta data of each sample (NA if not applicable).
#' Each column contains meta data in the following order.
#' \itemize{
#' \item chrono_age : numeric chronological age
#' \item epigen_age : numeric epigenetic age
#' \item exp_labels : experiment label
#' \item train_test_labels : label of whether the sample is in train/test/validation set
#' }}
#' @return a list of data frame such that
#'  list("meth_data" = meth_data, "meta_data" = meta_data)
#' @export
#' @examples
#' \dontrun{check_valid_data(meth_data, meta_data)}
#'
check_valid_data <- function(meth_data, meta_data){
  #check for valid column and row numbers
  if(ncol(meth_data) != nrow(meta_data)){
    stop('Number of methylation data columns does not match meta data rows. \n',
         'Methylation data: ',ncol(meth_data),' columns\n',
         'Meta data: ',nrow(meta_data),' columns')
  }
  if(nrow(meth_data)==0){
    stop('Methylation data has 0 rows.')
  }
  if(nrow(meta_data)==0){
    stop('Meta data has 0 rows.')
  }
  if(ncol(meth_data)==0){
    stop('Methylation data has 0 columns.')
  }
  if(ncol(meta_data)!=5){
    stop('Meta data does not have 5 columns.
         Please enter NA for columns with no input.
         See example dataset for more details')
  }

  #check for missing values
  if(any(is.na(meth_data))){
    stop('Methylation data cannot contain missing values.')
  }
  if(!all(apply(meth_data,2,is.numeric))){
    stop('Methylation data must only contain numeric values.')
  }
  if(!all(apply(meta_data, 2, function(x) {
    if (is.numeric(x) && any(is.na(x))){
      return(FALSE)
    } else {
      return(TRUE)
    }
  }))){
    stop('All numeric columns in meta data cannot have missing data.')
  }


  #check if values are valid
  if(!all(apply(meth_data,2,function(x) {all(ifelse(x<1.0 & x>0.0,TRUE, FALSE))} ))){
    stop('All methylation data must be smaller than 1.0 and larger than 0.0.')
  }
  if(!all(apply(meta_data, 2, function(x) {
    if (is.numeric(x) && any(ifelse(x<0.0,TRUE, FALSE))){
      return(FALSE)
    } else {
      return(TRUE)
    }
  }))){
    stop('All age (numeric entries) in meta data must be larger than 0.0.')
  }

  #check and rename columns and rows
  if(!all(c("sample_labels","chrono_age","epigen_age","exp_labels", "train_test_labels") %in% colnames(meta_data))){
    stop("Meta data columns are not named correctly. Please name them as
          sample_labels, chrono_age, epigen_age, exp_labels,  train_test_labels, respectively.
          See example dataset for more details.")
  }
  if( any(is.na(meta_data$sample_labels))){
    message("Missing label entry. Renaming samples by index.")
    rownames(meta_data) <- paste('sample', seq_len(nrow(meta_data)), sep = '-')
    colnames(meth_data) <- paste('sample', seq_len(nrow(meta_data)), sep = '-')
  }

  return(list("meth_data" = meth_data, "meta_data" = meta_data))

}

#' @title Remove unwanted experiment
#'
#' @description Remove unwanted experiment given a vector of experiment labels
#' @param data_list  a list of methylation and meta data data frame respectively
#' @param bad_exps a character vector of experiment labels
#' @return
#' \describe{
#'   a list of methylation and meta data data frames without unwanted experiment
#'   such that list("meth_data" = meth_data, "meta_data" = meta_data)
#' }
#' @export
#' @examples
#' \dontrun{remove_bad_exp(list("meth_data" = meth_data, "meta_data" = meta_data),
#'                c('GSE45363','GSE34564'))}

remove_bad_exp <- function(data_list, bad_exps){
  if(length(bad_exps) != 0){
    meth <- data_list$meth_data
    meta <- data_list$meta_data
    #bad_ind <-meta$exp_labels %in% bad_exps
    return(list("meth_data" = meth[,!is.element(meta$exp_labels, bad_exps)],
                "meta_data" = meta[!is.element(meta$exp_labels, bad_exps),]))
  }
  return(data_list)
}



