context("max_depth")
library(RADstackshelpR)

test_that("max_depth generates output of the appropriate class (vcfR)", {
  #find data in package using CRAN friendly syntax
  f<-max_depth(vcfR = SNPfiltR::vcfR.example,
                maxdepth = 100)
  #test that max_depth returns an object of class "vcfR"
  expect_is(f, "vcfR")
})


test_that("max_depth generates output of the appropriate class (vcfR) even with no specified cutoffs", {
  #find data in package using CRAN friendly syntax
  f<-max_depth(vcfR = SNPfiltR::vcfR.example)
  #test that max_depth returns an object of class "vcfR"
  expect_is(f, "character")
})


test_that("max_depth catches all samples at a given missing data cutoff", {
  #find data in package using CRAN friendly syntax
  f<-max_depth(vcfR = SNPfiltR::vcfR.example,
                 maxdepth = 100)
  g<-max_depth(vcfR = f,
                 maxdepth = 100)
  #test that all SNPs are caught at a given cutoff by max_depth
  expect_equal(nrow(f@gt), nrow(g@gt))
})


test_that("max_depth generates an error if run with a non-vcf file", {
  #expect error trying to read this vector as a vcf file
  expect_error(
    max_depth(vcfR = SNPfiltR::popmap)
  )
})


test_that("max_depth generates an error if run with a non-numeric depth cutoff", {
  #expect error trying to read this vector as a vcf file
  expect_error(
    max_depth(vcfR = SNPfiltR::vcfR.example,
                maxdepth = "k")
  )
})

