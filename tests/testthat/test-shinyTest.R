library(shinytest)

test_that("shinytest", {
  # Don't run these tests on the CRAN build servers
  skip_on_cran()
  
  appdir <- system.file(package = "viewpoly", "tests/shinytest")
  expect_pass(testApp(appdir, compareImages = FALSE))
})