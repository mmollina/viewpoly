library(shinytest)

test_that("shinytest", {
  # Don't run these tests on the CRAN build servers
  skip_on_cran()
  
  expect_pass(testApp(testnames = "shinyTest_all_tabs", compareImages = FALSE))
})