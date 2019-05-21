library(usethis)

adducts2 <- read.table("adducts2.tsv",sep="\t",header=T)

usethis::use_data(adducts2,overwrite=TRUE)
