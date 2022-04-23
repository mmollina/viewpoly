test_that("Tests uploaded MAPpoly files",{
  # upload MAPpoly
  temp <- tempfile()
  download.file("https://www.polyploids.org/sites/default/files/2022-04/tetra_MAPpoly_maps.RData", destfile = temp)
  temp.name <- load(temp)
  input.data <- get(temp.name)
  viewmap_mappoly <- prepare_MAPpoly(input.data)
  
  expect_equal(check_viewmap(viewmap_mappoly),0)
  
  check_viewmap_values(viewmap_mappoly, 
                       c(14, 132, 139, 157, 34), 
                       c(250, 226), # bases codified as A and B
                       50502.07)
  
  # VIEWmap tests 
  maps <- lapply(viewmap_mappoly$maps, function(x) {
    y <- x$l.dist
    names(y) <- x$mk.names
    y
  })
  
  # Get max size each chromosome
  expect_equal(map_summary(left.lim = 1, 
                           right.lim = 50, 
                           ch = 3, 
                           maps = maps, 
                           d.p1 = viewmap_mappoly$d.p1, 
                           d.p2 = viewmap_mappoly$d.p2)[[5]], 134.073, tolerance = 0.0001)
  
  # Map summary table
  summary_table <- summary_maps(viewmap_mappoly)
  expect_equal(sum(as.numeric(summary_table$`Map length (cM)`)), 3259.98)
  expect_equal(sum(as.numeric(summary_table$Simplex)), 2450)
  expect_equal(sum(as.numeric(summary_table$`Double-simplex`)), 1820)
  expect_equal(sum(as.numeric(summary_table$`Max gap`)), 80.51)
  
  # VIEWgenome tests
  p <- plot_cm_mb(viewmap_mappoly, 1, 1,50)
  
  expect_equal(sum(p$data$l.dist), 50502.07, tolerance = 0.001)
  
})