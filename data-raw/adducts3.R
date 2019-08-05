library(usethis)

adducts3 <- read.table("data-raw/adducts3.tsv",sep="\t",header=T)

usethis::use_data(adducts3,overwrite=TRUE)
