library(usethis)

adducts3 <- read.table("adducts3.tsv",sep="\t",header=T)

usethis::use_data(adducts3,overwrite=TRUE)