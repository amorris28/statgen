test_that("package load exposes version and verbosity controls", {
  expect_identical(version(), as.character(utils::packageVersion("statgen")))
  expect_identical(get_verbosity(), "info")

  expect_invisible(set_verbosity("quiet"))
  expect_identical(get_verbosity(), "quiet")

  expect_invisible(set_verbosity("info"))
  expect_identical(get_verbosity(), "info")

  expect_error(set_verbosity("verbose"), "quiet.*info")
})
