context("distance_thin")
library(RADstackshelpR)

test_that("distance_thin generates output of the appropriate class (vcfR)", {
  #find data in package using CRAN friendly syntax
  f<-distance_thin(vcfR = SNPfiltR::vcfR.example,
                   min.distance = 100)
  #test that distance_thin returns an object of class "vcfR"
  expect_is(f, "vcfR")
})


test_that("distance_thin throws error with no specified cutoffs", {
  #find data in package using CRAN friendly syntax
  expect_error(
    distance_thin(vcfR = SNPfiltR::popmap)
  )
})


test_that("distance_thin catches all samples at a given missing data cutoff", {
  #find data in package using CRAN friendly syntax
  f<-distance_thin(vcfR = SNPfiltR::vcfR.example,
                 min.distance = 100)
  g<-distance_thin(vcfR = f,
                 min.distance = 100)
  #test that all SNPs are caught at a given cutoff by distance_thin
  expect_equal(nrow(f@gt), nrow(g@gt))
})


test_that("distance_thin generates an error if run with a non-vcf file", {
  #expect error trying to read this vector as a vcf file
  expect_error(
    distance_thin(vcfR = SNPfiltR::popmap)
  )
})


test_that("distance_thin generates an error if run with a non-numeric depth cutoff", {
  #expect error trying to read this vector as a vcf file
  expect_error(
    distance_thin(vcfR = SNPfiltR::vcfR.example,
                min.distance = "k")
  )
})

