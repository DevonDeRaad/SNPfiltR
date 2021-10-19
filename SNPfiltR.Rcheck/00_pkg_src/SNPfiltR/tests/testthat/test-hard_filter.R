context("hard_filter")
library(RADstackshelpR)

test_that("hard_filter generates output of the appropriate class (vcfR)", {
  #find data in package using CRAN friendly syntax
  f<-hard_filter(vcfR = SNPfiltR::vcfR.example,
                 depth = 5,
                 gq = 30)
  #test that hard_filter returns an object of class "vcfR"
  expect_is(f, "vcfR")
})


test_that("hard_filter generates output of the appropriate class (vcfR) even with no specified cutoffs", {
  #find data in package using CRAN friendly syntax
  f<-hard_filter(vcfR = SNPfiltR::vcfR.example)
  #test that hard_filter returns an object of class "vcfR"
  expect_is(f, "vcfR")
})


test_that("hard_filter catches all samples at a given missing data cutoff", {
  #find data in package using CRAN friendly syntax
  f<-hard_filter(vcfR = SNPfiltR::vcfR.example,
                       depth = 5)
  g<-hard_filter(vcfR = f,
                       depth = 5)
  #test that all SNPs are caught at a given cutoff by hard_filter
  expect_equal(nrow(f@gt), nrow(g@gt))
})


test_that("hard_filter generates an error if run with a non-vcf file", {
  #expect error trying to read this vector as a vcf file
  expect_error(
    hard_filter(vcfR = SNPfiltR::popmap)
  )
})


test_that("hard_filter generates an error if run with a non-numeric depth cutoff", {
  #expect error trying to read this vector as a vcf file
  expect_error(
    hard_filter(vcfR = SNPfiltR::vcfR.example,
                      depth = "k")
  )
})


test_that("hard_filter generates an error if run with a non-numeric genotype quality cutoff", {
  #expect error trying to read this vector as a vcf file
  expect_error(
    hard_filter(vcfR = SNPfiltR::vcfR.example,
                gq = "k")
  )
})
