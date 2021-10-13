context("filter_allele_balance")
library(RADstackshelpR)

test_that("filter_allele_balance generates output of the appropriate class (vcfR)", {
  #find data in package using CRAN friendly syntax
  f<-filter_allele_balance(vcfR = SNPfiltR::vcfR.example)
  #test that filter_allele_balance returns an object of class "vcfR"
  expect_is(f, "vcfR")
})


test_that("filter_allele_balance catches all samples at a given missing data cutoff", {
  #find data in package using CRAN friendly syntax
  f<-filter_allele_balance(vcfR = SNPfiltR::vcfR.example)
  g<-filter_allele_balance(vcfR = f)
  #test that all SNPs are caught at a given cutoff by filter_allele_balance
  expect_equal(nrow(f@gt), nrow(g@gt))
})


test_that("filter_allele_balance generates an error if run with a non-vcf file", {
  #expect error trying to read this vector as a vcf file
  expect_error(
    filter_allele_balance(vcfR = SNPfiltR::popmap)
  )
})
