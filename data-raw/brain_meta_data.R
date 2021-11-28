brain_meta_file <- system.file(
  "inst",
  "extdata",
  "brain_pseudotime_meta.csv",
  package = "PseudotimeMethylation"
)

brain_meta_data <- read.csv(
  brain_meta_file,
  stringsAsFactors = FALSE,
  encoding = "UTF-8"
)[,-1]
colnames(brain_meta_data) <- c("sample_labels", colnames(brain_meta_data)[-1])
usethis::use_data(brain_meta_data, overwrite = TRUE)
