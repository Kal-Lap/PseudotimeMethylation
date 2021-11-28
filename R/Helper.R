#' @title print sample number in an experiment
#'
#' @description print sample number of given experiment label
#' @param meta_data a data frame of meta data
#' @param exp_label string of experiment label
#' @export
#' @examples
#' \dontrun{print_sample_num_by_train_test(meta_data,'GSE42861')}
print_sample_num_by_exp <- function(meta_data, exp_label)
{
  print(paste(exp_label,":",length(which(meta_data$exp_labels == exp_label))))
}

#' @title print sample number by experiment in train/test set
#'
#' @description print sample number in each experiment of the
#' given train-test label
#' @param meta_data a data frame of meta data
#' @param train_test_label a string for train, val, or test label
#' @export
#' @examples
#' \dontrun{print_sample_num_by_train_test(meta_data,"test")
#' print_sample_num_by_train_test(meta_data,"train")}
#

print_sample_num_by_train_test <- function(meta_data, train_test_label){
  for (i in unique(meta_data[meta_data$train_test_labels == train_test_label, ]$exp_labels)){
    print(paste(i,":",length(which(meta_data[meta_data$train_test_labels == train_test_label, ]$exp_labels == i))))
  }
}
