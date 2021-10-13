context("missing_by_sample")
library(RADstackshelpR)

test_that("missing_by_sample generates output of the appropriate class (vcfR)", {
  #find data in package using CRAN friendly syntax
  f<-missing_by_sample(vcfR = SNPfiltR::vcfR.example,
                    cutoff = .6)
  #test that missing_by_sample returns an object of class "vcfR"
  expect_is(f, "vcfR")
})


test_that("missing_by_sample generates output of the appropriate class (data.frame), with no specified cutoff", {
  #find data in package using CRAN friendly syntax
  f<-missing_by_sample(vcfR = SNPfiltR::vcfR.example)
  #test that missing_by_sample returns an object of class "vcfR"
  expect_is(f, "data.frame")
})


test_that("missing_by_sample catches all samples at a given missing data cutoff", {
  #find data in package using CRAN friendly syntax
  f<-missing_by_sample(vcfR = SNPfiltR::vcfR.example,
                    cutoff = .6)
  g<-missing_by_sample(vcfR = f,
                    cutoff = .6)
  #test that all SNPs are caught at a given cutoff by missing_by_sample
  expect_equal(ncol(f@gt), ncol(g@gt))
})


test_that("missing_by_sample generates an error if run with a non-vcf file", {
  #expect error trying to read this vector as a vcf file
  expect_error(
    missing_by_sample(vcfR = SNPfiltR::popmap,
                   cutoff = .6)
  )
})


test_that("missing_by_sample generates an error if run with an out of bounds cutoff (not 0-1)", {
  #expect error trying to read this vector as a vcf file
  expect_error(
    missing_by_sample(vcfR = SNPfiltR::vcfR.example,
                   cutoff = 5)
  )
})


test_that("missing_by_sample generates an error if run with a malformed popmap", {
  #expect error trying to read this vector as a vcf file
  f<-SNPfiltR::popmap
  f<-f[,c(2,1)]
  expect_error(
    missing_by_sample(vcfR = SNPfiltR::vcfR.example,
                      popmap = f)
  )
})


