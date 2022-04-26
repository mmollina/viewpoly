library(shinytest)

test_that("shinytest", {
  # Don't run these tests on the CRAN build servers
  skip_on_cran()
  
  expect_pass(testApp(appDir = system.file(package = "viewpoly"),compareImages = FALSE))
})