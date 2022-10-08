test_that("Tests uploaded polymapR files",{
  source(system.file("ext/functions4tests.R", package = "viewpoly"))
  library(curl)
  
  # upload MAPpoly
  input.data <- polymapR.map <- list()
  input.data$datapath <- tempfile()
  polymapR.map$datapath <- tempfile()
  
  if(has_internet()){
    download.file("https://www.polyploids.org/sites/default/files/2022-04/tetra_polymapR_dataset.RData", destfile = input.data$datapath)
    download.file("https://www.polyploids.org/sites/default/files/2022-04/tetra_polymapR_map.RData", destfile = polymapR.map$datapath)
    
    viewmap_polymapr <- prepare_polymapR(polymapR.dataset = input.data, polymapR.map, input.type = "disc", ploidy = 4)
    
    expect_equal(check_viewmap(viewmap_obj = viewmap_polymapr),0)
    
    check_viewmap_values(viewmap_polymapr,
                         c(52, 85, 95, 50, 34),
                         c(126, 190), # bases codified as A and B
                         28122)
    
    # VIEWmap tests
    maps <- lapply(viewmap_polymapr$maps, function(x) {
      y <- x$l.dist
      names(y) <- x$mk.names
      y
    })
    
    # Get max size each chromosome
    expect_equal(map_summary(left.lim = 1,
                             right.lim = 50,
                             ch = 3,
                             maps = maps,
                             d.p1 = viewmap_polymapr$d.p1,
                             d.p2 = viewmap_polymapr$d.p2)[[5]], 96.15, tolerance = 0.0001)
    
    # Map summary table
    summary_table <- summary_maps(viewmap_polymapr)
    expect_equal(sum(as.numeric(summary_table$`Map length (cM)`)), 2317.84)
    expect_equal(sum(as.numeric(summary_table$Simplex)), 2028)
    expect_equal(sum(as.numeric(summary_table$`Double-simplex`)), 802)
    expect_equal(sum(as.numeric(summary_table$`Max gap`)), 111.75)
    
    # VIEWgenome tests
    p <- plot_cm_mb(viewmap_polymapr, 1, 1,50)
    
    expect_equal(sum(p$data$l.dist), 28122, tolerance = 0.001)
  } else {
    print("polymapR tests are only executed if internet conection is available.")
  }
})