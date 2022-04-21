test_that("Linkage map graphics and tables",{
  
  viewpoly_obj <- prepare_examples("tetra_map")
  
  qtl_profile_plot <- plot_profile(viewpoly_obj$qtl$profile, 
                                   viewpoly_obj$qtl$qtl_info, 
                                   viewpoly_obj$qtl$selected_mks, 
                                   pheno.col = 2, 
                                   lgs.id = 2, 
                                   by_range = TRUE, 
                                   range.min = 30, 
                                   range.max = 120, 
                                   plot=TRUE, 
                                   software = NULL)
  
  expect_equal(sum(qtl_profile_plot$data$SIG, na.rm = TRUE), 43.81917, tolerance = 0.0001)
  
  maps <- lapply(viewpoly_obj$map$maps, function(x) {
    y <- x$l.dist
    names(y) <- x$mk.names
    y
  })
  
  # Get max size each chromosome
  expect_equal(map_summary(left.lim = 1, 
                           right.lim = 50, 
                           ch = 3, 
                           maps = maps, 
                           d.p1 = viewpoly_obj$map$d.p1, 
                           d.p2 = viewpoly_obj$map$d.p2)[[5]], 134.073, tolerance = 0.0001)
  
  # Map summary table
  summary_table <- summary_maps(viewpoly_obj$map)
  expect_equal(sum(as.numeric(summary_table$`Map length (cM)`)), 3259.98)
  expect_equal(sum(as.numeric(summary_table$Simplex)), 2450)
  expect_equal(sum(as.numeric(summary_table$`Double-simplex`)), 1820)
  expect_equal(sum(as.numeric(summary_table$`Max gap`)), 80.51)
  
  p <- plot_cm_mb(viewpoly_obj$map, 1, 1,50)
  
  expect_equal(sum(p$data$l.dist), 50502.07, tolerance = 0.001)
  
})
