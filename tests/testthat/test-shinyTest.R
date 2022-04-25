library(shinytest)

test_that("shinytest", {
  # Don't run these tests on the CRAN build servers
  skip_on_cran()
  
  appdir <- gsub("inst", "", system.file(package = "viewpoly"))
  expect_pass(testApp(appDir = appdir, testnames = "shinyTest_all_tabs", compareImages = FALSE))
})