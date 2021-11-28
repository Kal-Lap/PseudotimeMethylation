#' @title Select site by variance
#'
#' @description Select site with larger Pearson's Correlation Coefficient
#' with chronological age than the given threshold
#' @param meth_data a data frame of methylation data
#' @param meta_data a data frame of meta data
#' @param var_thresh a numeric value of variance threshold
#' @return  a vector of index that satisfy the variance threshold
#' @export
#' @examples
#' \dontrun{select_site_var(meth_data, meta_data)
#' select_site_var(meth_data, meta_data, var_thresh = 0.002)}
select_site_var <- function(meth_data, meta_data, var_thresh = 0.001){
  row_var <- apply(meth_data, 1, stats::var)
  return(which(row_var) > var_thresh)
}

#' @title Select site by Pearson's Correlation Coefficient
#'
#' @description Select site with larger Pearson's Correlation Coefficient (PCC)
#' with chronological age than the given threshold
#' @param meth_data a data frame of methylation data
#' @param meta_data a data frame of meta data
#' @param pcc_thresh a numeric value of PCC threshold
#' @return  a vector of index that satisfy the PCC threshold
#' @export
#' @examples
#' \dontrun{select_site_pcc(meth_data, meta_data)
#' select_site_pcc(meth_data, meta_data, pcc_thresh = 0.002)}
select_site_pcc <- function(meth_data, meta_data, pcc_thresh = 0.5){
  row_pcc <- apply(meth_data, 1, stats::cor, meta_data$chrono_age)
  return(which(abs(row_pcc) > pcc_thresh))
}
