test_that("datasets (if any) load successfully", {
  pkg <- "c060"
  data_dir <- system.file("data", package = pkg)
  if (!nzchar(data_dir)) {
    skip("No datasets directory present")
  }
  files <- dir(data_dir, pattern = "\\.(rda|RData)$", ignore.case = TRUE)
  if (!length(files)) {
    skip("No datasets found")
  }
  for (ff in files) {
    ds <- sub("\\.(rda|RData)$", "", ff, ignore.case = TRUE)
    suppressWarnings(suppressMessages(data(list = ds, package = pkg, envir = environment())))
    expect_true(exists(ds, inherits = FALSE), info = sprintf("Dataset '%s' did not load", ds))
    obj <- get(ds, inherits = FALSE)
    ok <- tryCatch({
      if (is.data.frame(obj) || is.matrix(obj)) {
        nrow(obj) > 0 && ncol(as.matrix(obj)) >= 1
      } else if (is.list(obj)) {
        length(obj) > 0
      } else if (is.vector(obj)) {
        length(obj) > 0
      } else { TRUE }  # don't be too strict for exotic objects
    }, error = function(e) FALSE)
    expect_true(ok, info = sprintf("Dataset '%s' appears empty", ds))
    # cleanup
    rm(list = ds, inherits = FALSE)
  }
})
