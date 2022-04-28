library(shinytest)

test_that("shinytest", {
  # Don't run these tests on the CRAN build servers
  skip_on_cran()
  
  if(!grepl("inst", system.file(package = "viewpoly"))) 
    appdir <- paste0(system.file(package = "viewpoly"), "/inst")
  else appdir <-  system.file(package = "viewpoly")
  
  expect_pass(testApp(appDir = appdir, testnames = c("github_actions_tests"),compareImages = FALSE))
})