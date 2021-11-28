brain_meth_file <- system.file(
  "inst",
  "extdata",
  "brain_pseudotime.csv",
  package = "PseudotimeMethylation"
)

brain_meth_data_df <- read.csv(
  brain_meth_file,
  stringsAsFactors = FALSE,
  encoding = "UTF-8"
)

brain_meth_data <- data.matrix(brain_meth_data_df[,-1])
usethis::use_data(brain_meth_data, overwrite = TRUE)
