### Summary of new changes
I fixed a bug in the vignette that caused the build to fail and resulted in the package being removed from CRAN in January 2025. This patched version (1.0.3) should now be suitable for inclusion on CRAN.

### Output from devtools::check() run on this version of the package:

── R CMD check results ─────────────────────────────────────────────────────────────────────────────────────── SNPfiltR 1.0.3 ────
Duration: 42.8s

❯ checking for future file timestamps ... NOTE
  unable to verify current time

0 errors ✔ | 0 warnings ✔ | 1 note ✖

Stack overflow posts suggest that this note can safely be ignored and doesn't reflect a problem with the package.