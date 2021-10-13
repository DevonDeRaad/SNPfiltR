context("min_mac")
library(RADstackshelpR)

test_that("min_mac generates output of the appropriate class (vcfR)", {
  #find data in package using CRAN friendly syntax
  f<-min_mac(vcfR = SNPfiltR::vcfR.example,
             min.mac = 2)
  #test that min_mac returns an object of class "vcfR"
  expect_is(f, "vcfR")
})


test_that("min_mac generates output of the appropriate class (vcfR) even with no specified cutoffs", {
  #find data in package using CRAN friendly syntax
  f<-min_mac(vcfR = SNPfiltR::vcfR.example)
  #test that min_mac returns an object of class "vcfR"
  expect_is(f, "vcfR")
})


test_that("min_mac catches all samples at a given missing data cutoff", {
  #find data in package using CRAN friendly syntax
  f<-min_mac(vcfR = SNPfiltR::vcfR.example,
                 min.mac = 5)
  g<-min_mac(vcfR = f,
                 min.mac = 5)
  #test that all SNPs are caught at a given cutoff by min_mac
  expect_equal(nrow(f@gt), nrow(g@gt))
})


test_that("min_mac generates an error if run with a non-vcf file", {
  #expect error trying to read this vector as a vcf file
  expect_error(
    min_mac(vcfR = SNPfiltR::popmap)
  )
})


test_that("min_mac generates an error if run with a non-numeric depth cutoff", {
  #expect error trying to read this vector as a vcf file
  expect_error(
    min_mac(vcfR = SNPfiltR::vcfR.example,
                min.mac = "k")
  )
})
