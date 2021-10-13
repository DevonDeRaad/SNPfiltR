context("missing_by_snp")
library(RADstackshelpR)

test_that("missing_by_snp generates output of the appropriate class (vcfR)", {
  #find data in package using CRAN friendly syntax
  f<-missing_by_snp(vcfR = SNPfiltR::vcfR.example,
                    cutoff = .6)
  #test that missing_by_snp returns an object of class "vcfR"
  expect_is(f, "vcfR")
})


test_that("missing_by_snp generates output of the appropriate class (data.frame), with no cutoff specified", {
  #find data in package using CRAN friendly syntax
  f<-missing_by_snp(vcfR = SNPfiltR::vcfR.example)
  #test that missing_by_snp returns an object of class "vcfR"
  expect_is(f, "data.frame")
})


test_that("missing_by_snp catches all SNPs at a given missing data cutoff", {
  #find data in package using CRAN friendly syntax
  f<-missing_by_snp(vcfR = SNPfiltR::vcfR.example,
                    cutoff = .6)
  g<-missing_by_snp(vcfR = f,
                    cutoff = .6)
  #test that all SNPs are caught at a given cutoff by missing_by_snp
  expect_equal(nrow(f@gt), nrow(g@gt))
})


test_that("missing_by_snp generates an error if run with a non-vcf file", {
  #expect error trying to read this vector as a vcf file
  expect_error(
    missing_by_snp(vcfR = SNPfiltR::popmap,
                      cutoff = .6)
  )
})


test_that("missing_by_snp generates an error if run with an out of bounds cutoff (not 0-1)", {
  #expect error trying to read this vector as a vcf file
  expect_error(
    missing_by_snp(vcfR = SNPfiltR::vcfR.example,
                   cutoff = 5)
  )
})

