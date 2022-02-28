### Summary of new changes
I fixed bugs in 'assess_missing_data_pca' and 'assess_missing_data_tsne' functions, which caused the vignette to throw an error and resulted in CRAN removal. I then incremented the package version 0.1.0 -> 0.1.1 for CRAN re-submission.

### Output from rhub::check_for_cran() run on the tarball being submitted here.
## Test environments
- R-hub Windows Server 2022, R-devel, 64 bit
- R-hub Ubuntu Linux 20.04.1 LTS, R-release, GCC
- R-hub Fedora Linux, R-devel, clang, gfortran

## R CMD check results
* checking CRAN incoming feasibility ... NOTE
Maintainer: ‘Devon DeRaad <devonderaad@gmail.com>’

New submission

Package was archived on CRAN

Possibly mis-spelled words in DESCRIPTION:
  Grünwald (18:34)
  Knaus (18:24)
  polymorphism (10:21)
  reproducibly (9:47, 13:2)
  Rmarkdown (13:19)
  vcf (17:32, 20:43)
  vcfR (19:55)

CRAN repository db overrides:
  X-CRAN-Comment: Archived on 2022-01-06 as check issues were not
    corrected in time.