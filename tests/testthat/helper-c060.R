# Helper for c060 tests
set.seed(1)
# Keep examples fast on CRAN
if (identical(Sys.getenv("NOT_CRAN"), "true")) {
  options(testthat.progress.max_fails = 10L)
}
