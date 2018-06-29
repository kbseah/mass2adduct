context("Basic tests")
library(mass2adduct)

test_that("Included data files are present", {
  expect_true(file.exists(system.file("extdata","msi.csv",package="mass2adduct")))
  expect_true(file.exists(system.file("extdata","mass_example",package="mass2adduct")))
})

d <- msimat(system.file("extdata","msi.csv",package="mass2adduct"),sep=";")
peaklist <- scan(system.file("extdata","mass_example",package="mass2adduct"),what=numeric())

test_that("Import file classes", {
  expect_is(d, "msimat")
  expect_is(peaklist, "numeric")
})

test_that("Conversion of CSV to Triplet input matrix", {
  e <- msimat(system.file("extdata","msi.csv",package="mass2adduct"),
              sep=";",
              crlf=TRUE,
              reformat.large.csv="test")
  expect_is(e,"msimat")
  expect_is(e$mat,"TsparseMatrix")
  expect_equal(d$peaks,e$peaks)
  file.remove(file.path(".","test_rows.list"))
  file.remove(file.path(".","test_cols.list"))
  file.remove(file.path(".","test_vals.list"))
  file.remove(file.path(".","test_peaknames.list"))
  file.remove(file.path(".","test_spotnames.list"))
})

d.diff <- massdiff(d)

test_that("Diff calculation", {
  expect_message(peaklist.diff <- massdiff(peaklist))
  expect_is(d.diff,"massdiff")
  expect_is(peaklist.diff,"massdiff")
  expect_equal(dim(d.diff),c(choose(length(d$peaks),2),3))
})

test_that("Diff histogram and annotation", {
  d.diff.hist <- hist(d.diff)
  expect_is(d.diff.hist,"histogram")
  expect_is(d.diff.hist,"massdiffhist")
  d.diff.hist.annot <- adductMatch(d.diff.hist,add=adducts2)
  expect_is(d.diff.hist.annot,"data.frame")
  expect_equal(dim(d.diff.hist.annot),c(9,5))
})

d.diff.annot <- adductMatch(d.diff,add=adducts2)

test_that("Adduct matching", {
  expect_is(d.diff.annot,"massdiff")
  expect_is(d.diff.annot,"data.frame")
  expect_equal(dim(d.diff.annot),c(11035,5))
})

test_that("Correlation testing", {
  expect_message(d.diff.annot.cor <- corrPairsMSI(d,d.diff.annot))
  expect_is(d.diff.annot.cor,"data.frame")
  expect_equal(sum(d.diff.annot.cor$Significance),3867)
})

test_that("Parallelized code", {
  if (.Platform$OS.type != "unix") {
    skip("Not using Unix")
  }
  if (parallel::detectCores() <= 1) {
    skip("Multiple cores not detected")
  }
  serialversion <- corrPairsMSI(d,d.diff.annot,how="apply")
  parallelversion <- corrPairsMSI(d,d.diff.annot,how="parallel",ncores=2)
  chunkedversion <- corrPairsMSIchunks(d,d.diff.annot,ncores=2,mem.limit=0.01) # Should use two chunks
  expect_equal(serialversion, parallelversion)
  expect_equal(serialversion, chunkedversion)
})

test_that("Platform detection", {
  if (.Platform$OS.type != "windows") {
    skip("Not using Windows")
  }
  expect_error(corrPairsMSI(d,d.diff.annot,how="parallel",ncores=2))
  expect_error(corrPairsMSIchunks(d,d.diff.annot,ncores=2,mem.limit=0.01)) # Should use two chunks
})

# Cleanup
rm(d,peaklist,d.diff,d.diff.annot)