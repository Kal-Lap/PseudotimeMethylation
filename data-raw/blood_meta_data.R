## Prepare `meta_data` dataset

blood_meta_file <- system.file(
  "inst",
  "extdata",
  "all_blood_dataset_pseudotime_meta_v.4.csv",
  package = "PseudotimeMethylation"
)


blood_meta_data <- read.csv(
  blood_meta_file,
  stringsAsFactors = FALSE,
  encoding = "UTF-8"
)

colnames(blood_meta_data) <- c("sample_labels", colnames(blood_meta_data)[-1])

#Relabel (val -> test, test -> val)
blood_meta_data[blood_meta_data =="val"]<-"t"
blood_meta_data[blood_meta_data == "test"] <-"val"
blood_meta_data[blood_meta_data =="t"] <- "test"

usethis::use_data(blood_meta_data, overwrite = TRUE)

