library(usethis)

adducts <- read.table("adducts.tsv",sep="\t",header=T)

usethis::use_data(adducts,overwrite=TRUE)
