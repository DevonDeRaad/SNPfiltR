context("filter_biallelic")
library(RADstackshelpR)

test_that("filter_biallelic generates output of the appropriate class (vcfR)", {
  #find data in package using CRAN friendly syntax
  f<-filter_biallelic(vcfR = SNPfiltR::vcfR.example)
  #test that filter_biallelic returns an object of class "vcfR"
  expect_is(f, "vcfR")
})


test_that("filter_biallelic catches all samples at a given missing data cutoff", {
  #find data in package using CRAN friendly syntax
  f<-filter_biallelic(vcfR = SNPfiltR::vcfR.example)
  g<-filter_biallelic(vcfR = f)
  #test that all SNPs are caught at a given cutoff by filter_biallelic
  expect_equal(nrow(f@gt), nrow(g@gt))
})


test_that("filter_biallelic generates an error if run with a non-vcf file", {
  #expect error trying to read this vector as a vcf file
  expect_error(
    filter_biallelic(vcfR = SNPfiltR::popmap)
  )
})

