
#print package startup message
.onAttach <- function(libname, pkgname) {
  pkg.version <- utils::packageVersion("SNPfiltR")
packageStartupMessage(
  c(paste0("This is SNPfiltR v.", pkg.version, "\n\n"),
  "Detailed usage information is available at: devonderaad.github.io/SNPfiltR/ \n\n",
  "If you use SNPfiltR in your published work, please cite the following papers: \n\n",
  "DeRaad, D.A. 2021. SNPfiltR: an R package for interactive and reproducible SNP filtering. Preprint on Authorea. <http://dx.doi.org/10.22541/au.163976415.53888836/v1> \n\n",
  "Knaus, Brian J., and Niklaus J. Grunwald. 2017. VCFR: a package to manipulate and visualize variant call format data in R. Molecular Ecology Resources 17.1:44-53. http://dx.doi.org/10.1111/1755-0998.12549.")
  )
}
