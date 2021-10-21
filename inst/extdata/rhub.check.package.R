#verify that package is ready for CRAN
#install rhub
#install.packages("rhub")

#library rhub
library(rhub)

#make sure to build this tarball using the following code while in the R project for the package:
devtools::build()

#now test package build on the same operating systems that CRAN uses, by pointing this function to the built tarball for the package
cran_prep <- check_for_cran(path = "/Users/devder/Desktop/SNPfiltR_0.0.0.9000.tar.gz")

#now run cran_summary to ensure that the package is ready for CRAN submission
cran_prep$cran_summary()
