### Summary of new changes
I fixed bugs in the 'distance_thin()' function and added internal checkpoints to other functions like 'hard_filter()' to make sure that user-friendly error messages are returned if the input vcf file does not contain information on the quality metric that the user is attempting to filter based on. I then incremented the package version 0.1.1 -> 1.0.0 for CRAN re-submission. This is now the official release version of the package, corresponding to the publication of the following manuscript: DeRaad, D.A. (2022), SNPfiltR: an R package for interactive and reproducible SNP filtering. Molecular Ecology Resources, 00, 1–15 https://doi.org/10.1111/1755-0998.13618

### Output from rhub::check_for_cran() run on the tarball being submitted here.
## Test environments
- R-hub Windows Server 2022, R-devel, 64 bit
- R-hub Ubuntu Linux 20.04.1 LTS, R-release, GCC
- R-hub Fedora Linux, R-devel, clang, gfortran

## R CMD check results
* checking CRAN incoming feasibility ... NOTE
Maintainer: ‘Devon DeRaad <devonderaad@gmail.com>’

Found the following (possibly) invalid DOIs:
  DOI: 10.1111/1755-0998.12549
    From: DESCRIPTION
    Status: Service Unavailable
    Message: 503
    
# This DOI is valid (checked multiple times), and this NOTE seems to be unavoidable. According to Hadley, these notes can be safely ignored for publishing on CRAN (https://twitter.com/hadleywickham/status/1358170607314235392)