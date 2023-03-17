#print package startup message
.onAttach <- function(libname, pkgname) {
  pkg.version <- utils::packageVersion("SNPfiltR")
packageStartupMessage(
  c(paste0("This is SNPfiltR v.", pkg.version, "\n\n"),
  "Detailed usage information is available at: devonderaad.github.io/SNPfiltR/ \n\n",
  "If you use SNPfiltR in your published work, please cite the following papers: \n\n",
  "DeRaad, D.A. (2022), SNPfiltR: an R package for interactive and reproducible SNP filtering. Molecular Ecology Resources, 22, 2443–2453. http://doi.org/10.1111/1755-0998.13618 \n\n",
  "Knaus, Brian J., and Niklaus J. Grunwald. 2017. VCFR: a package to manipulate and visualize variant call format data in R. Molecular Ecology Resources, 17.1:44-53. http://doi.org/10.1111/1755-0998.12549")
  )
}
