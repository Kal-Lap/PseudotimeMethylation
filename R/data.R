#' Methylation matrix of 6251 blood samples after preprocess
#'
#' A dataset containing the methylation levels of 6251 samples at 9447 CpG sites
#'
#' @format A data frame with 9447 rows and 6251 variables:
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/}
"blood_meth_data"

#' Meta data of 6251 methylation blood array samples
#'
#' A dataset containing the meta data of 6251 methylation array samples
#'
#' @format A data frame with 6251 rows and 5 variables:
#' \describe{
#'   \item{sample_labels}{label of sample}
#'   \item{chrono_age}{chronological age of sample}
#'   \item{epigen_age}{epigenetic state of sample}
#'   \item{exp_labels}{experiment or group label of sample}
#'   \item{train_test_labels}{label of whether sample is in train/val/test set}
#' }
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/}
"blood_meta_data"

#' Methylation matrix of 675 brain samples after preprocess
#'
#' A dataset containing the methylation levels of 675 samples at 8510 CpG sites
#'
#' @format A data frame with 8510 rows and 675 variables:
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/}
"brain_meth_data"

#' Meta data of 675 brain methylation array samples
#'
#' A dataset containing the meta data of 675 methylation array samples
#'
#' @format A data frame with 675 rows and 5 variables:
#' \describe{
#'   \item{sample_labels}{label of sample}
#'   \item{chrono_age}{chronological age of sample}
#'   \item{epigen_age}{epigenetic state of sample}
#'   \item{exp_labels}{experiment or group label of sample}
#'   \item{train_test_labels}{label of whether sample is in train/val/test set}
#' }
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/}
"brain_meta_data"
