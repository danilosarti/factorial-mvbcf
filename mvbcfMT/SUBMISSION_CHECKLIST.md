# CRAN submission checklist — mvbcfMT

A maintainer checklist for releasing `mvbcfMT` to CRAN. Every step must be clean
before submitting. This file is listed in `.Rbuildignore` and is not part of the
package tarball.

## 0. One-time: install tooling
```r
install.packages(c("devtools", "roxygen2", "rmarkdown", "knitr", "testthat",
                   "urlchecker", "spelling"))
```

## 1. Regenerate documentation and NAMESPACE
The roxygen tags are the source of truth. Regenerate `man/*.Rd` and `NAMESPACE`:
```r
setwd("mvbcfMT")
Rcpp::compileAttributes()      # refresh RcppExports after any C++ change
devtools::document()           # writes man/ and NAMESPACE from the roxygen tags
```
Confirm `NAMESPACE` exports the intended functions (every exported function
carries an `@export` tag, including `build_components`).

## 2. Fast local check
```r
devtools::check()              # runs R CMD check with tests, examples, vignette
```
Target: 0 errors, 0 warnings, 0 notes (a "New submission" note is expected).

## 3. As-CRAN check
```r
R CMD build .
R CMD check --as-cran mvbcfMT_0.1.0.tar.gz
```

## 4. Cross-platform checks (recommended for a first submission)
```r
devtools::check_win_devel()    # win-builder, R-devel
devtools::check_win_release()
# R-hub v2 (needs a GitHub repo):
# rhub::rhub_setup(); rhub::rhub_check(platforms = c("linux","windows","macos"))
devtools::check_mac_release()  # macOS builder
```

## 5. URLs and spelling
```r
urlchecker::url_check()        # no broken URLs
devtools::spell_check()        # fix genuine typos; add names to inst/WORDLIST
```

## 6. Fill in before submitting
- `cran-comments.md`: record the test environments used and the check summary.
- Optional but recommended: create the public source repository and add
  `URL:` and `BugReports:` to `DESCRIPTION`, then re-run `urlchecker::url_check()`.
- When the accompanying paper has a persistent identifier, add the reference to
  the `Description:` field (as `<doi:...>` or `<arXiv:...>`) and to
  `inst/CITATION`.

## 7. Submit
```r
devtools::release()            # runs final checks and submits to CRAN
```
Respond to the confirmation email from CRAN.

---

## Notes on package hygiene
- Console output from the sampler goes only through `Rcpp::Rcout` (no `printf`
  or `std::cout`).
- Randomness is drawn through R's RNG (`R::runif`), so results respect
  `set.seed()`; there is no separate RNG stream and no OpenMP region calls the R
  API.
- Every exported function has a title, documented arguments, a `@return` and a
  runnable example (heavier examples use tiny MCMC settings inside `\donttest{}`).
- `man/` is generated from the roxygen source by `devtools::document()`; do not
  edit the `.Rd` files by hand.
