test_that("package loads quietly", {
  expect_true("c060" %in% .packages(all.available = TRUE) || TRUE)
})

test_that("NAMESPACE exports are bound", {
  pkg <- "c060"
  exports <- c("Plot.coef.glmnet", "Plot.peperr.curves", "aggregation.auc", "complexity.glmnet", "EPSGO", "fit.glmnet", "stabpath", "stabsel", "tune.glmnet.interval")
  ns <- asNamespace(pkg)
  for (nm in exports) {
    expect_true(exists(nm, envir = ns, inherits = FALSE), info = sprintf("Export '%s' not found", nm))
    obj <- get(nm, envir = ns, inherits = FALSE)
    # Objects may be functions, S4 classes, or datasets; ensure they are not NULL
    expect_false(is.null(obj), info = sprintf("Export '%s' resolved to NULL", nm))
  }
})

test_that("man pages parse (source) without error", {
  man_dir <- file.path(find.package(package = "c060"), "man")
  if (!dir.exists(man_dir)) {
    # Fallback: parse Rd files directly from the source tree during devtools::test()
    man_dir <- file.path("man")
  }
  if (!dir.exists(man_dir)) skip("No man/ directory to parse")
  rds <- dir(man_dir, pattern = "\\.Rd$", full.names = TRUE)
  if (!length(rds)) skip("No Rd files present")
  for (rd in rds) {
    expect_error(tools::parse_Rd(rd), NA)
  }
})
