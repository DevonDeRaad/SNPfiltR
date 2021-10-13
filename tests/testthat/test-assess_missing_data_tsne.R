context("assess_missing_data_tsne")
library(RADstackshelpR)

test_that("assess_missing_data_tsne generates output of the appropriate class (list)", {
  #find data in package using CRAN friendly syntax
  f<-assess_missing_data_tsne(vcfR = SNPfiltR::vcfR.example,
                             popmap = SNPfiltR::popmap,
                             thresholds = .6)
  #test that assess_missing_data_tsne returns an object of class "vcfR"
  expect_is(f, "list")
})


test_that("assess_missing_data_tsne generates output of the appropriate class (data.frame), with no specified cutoff", {
  #find data in package using CRAN friendly syntax
  f<-assess_missing_data_tsne(vcfR = SNPfiltR::vcfR.example,
                             popmap = SNPfiltR::popmap)
  #test that assess_missing_data_tsne returns an object of class "vcfR"
  expect_is(f, "data.frame")
})


test_that("assess_missing_data_tsne generates an error if run with a non-vcf file", {
  #expect error trying to read this vector as a vcf file
  expect_error(
    assess_missing_data_tsne(vcfR = SNPfiltR::popmap,
                            popmap = SNPfiltR::popmap)
  )
})


test_that("assess_missing_data_tsne generates an error if run without a popmap", {
  #expect error trying to read this vector as a vcf file
  expect_error(
    assess_missing_data_tsne(vcfR = SNPfiltR::vcfR.example)
  )
})


test_that("assess_missing_data_tsne generates an error if run with an out of bounds cutoff (not 0-1)", {
  #expect error trying to read this vector as a vcf file
  expect_error(
    assess_missing_data_tsne(vcfR = SNPfiltR::vcfR.example,
                            popmap = SNPfiltR::popmap,
                            thresholds=5)
  )
})


test_that("assess_missing_data_tsne generates an error if run with a malformed popmap", {
  #expect error trying to read this vector as a vcf file
  f<-SNPfiltR::popmap
  f<-f[,c(2,1)]
  expect_error(
    assess_missing_data_tsne(vcfR = SNPfiltR::vcfR.example,
                            popmap = f)
  )
})


test_that("assess_missing_data_tsne generates an error if run with a non-numeric iterations input", {
  #expect error trying to read this vector as a vcf file
  expect_error(
    assess_missing_data_tsne(vcfR = SNPfiltR::vcfR.example,
                             popmap = SNPfiltR::popmap,
                             iterations = "k")
  )
})


test_that("assess_missing_data_tsne generates an error if run with a non-numeric initial_dims input", {
  #expect error trying to read this vector as a vcf file
  expect_error(
    assess_missing_data_tsne(vcfR = SNPfiltR::vcfR.example,
                             popmap = SNPfiltR::popmap,
                             initial_dims = "k")
  )
})


test_that("assess_missing_data_tsne generates an error if run with a non-numeric perplecity input", {
  #expect error trying to read this vector as a vcf file
  expect_error(
    assess_missing_data_tsne(vcfR = SNPfiltR::vcfR.example,
                             popmap = SNPfiltR::popmap,
                             perplexity = "k")
  )
})

