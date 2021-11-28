## prepare `meth_data` dataset

# retrieve paths to datafiles
blood_meth_file <- system.file(
  "inst",
  "extdata",
  "all_blood_dataset_pseudotime_v.4.csv",
  package = "PseudotimeMethylation"
)


blood_meth_data_df <- read.csv(
  blood_meth_file,
  stringsAsFactors = FALSE,
  encoding = "UTF-8"
)

blood_meth_data <-data.matrix(blood_meth_data_df[,-1])
usethis::use_data(blood_meth_data, overwrite = TRUE)

