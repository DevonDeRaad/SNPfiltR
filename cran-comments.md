### Summary of new changes
I fixed multiple bugs that resulted from lack of coverage of some possible idiosyncratic input vcf formats, and added additional internal checkpoints to try and catch issues with input vcf format not conforming to the expectations of the functions. Functions should now be much more likely to throw an informative error if vcf format is different than expected, rather than silently mischaracterizing genotypes internally. Specific details of bug fixes correspond to [SNPfiltR issues #1-3](https://github.com/DevonDeRaad/SNPfiltR/issues). The last thing I did was incremented the package version 1.0.0 -> 1.0.1 for CRAN re-submission.

### Output from rhub::check_for_cran() run on the tarball being submitted here.
## Test environments
- R-hub Windows Server 2022, R-devel, 64 bit
- R-hub Ubuntu Linux 20.04.1 LTS, R-release, GCC
- R-hub Fedora Linux, R-devel, clang, gfortran

## R CMD check results
There are no errors or warnings.
There is one NOTE that is only found on Windows (Server 2022, R-devel 64-bit): 

```
* checking for detritus in the temp directory ... NOTE
Found the following files/directories:
  'lastMiKTeXException'
```
[R-hub issue #503](https://github.com/r-hub/rhub/issues/503) indicates this is likely the result of a bug in MiKTeX, and can safely be ignored.