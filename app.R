# Launch the ShinyApp (Do not remove this comment)
# To deploy, run: rsconnect::deployApp()
# Or use the blue button on top of this file

## devtools::install_deps()
## devtools::install_dev_deps()
pkgload::load_all(export_all = FALSE,helpers = FALSE,attach_testthat = FALSE)
options( "golem.app.prod" = TRUE, shiny.autoload.r=FALSE)
viewpoly::run_app() # add parameters here (if any)



