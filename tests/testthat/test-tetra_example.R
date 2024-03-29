test_that("tetra example",{
  source(system.file("ext/functions4tests.R", package = "viewpoly"))
  
  # upload examples
  viewpoly_obj <- prepare_examples("tetra_map")
  
  expect_equal(check_viewpoly(viewpoly_obj),0)
  
  check_viewmap_values(viewpoly_obj$map,
                       c(14, 132, 139, 157, 34),
                       c(36, 167, 164, 109),
                       50502.07)
  
  check_viewqtl_qtlpoly_values(viewpoly_obj$qtl,
                               29320.98,
                               0.641104,
                               3.910344e-13,
                               22.44135,
                               0.003149918,
                               -8.511442e-14,
                               1)
  
  # VIEWmap tests
  qtl_profile_plot <- plot_profile(profile = viewpoly_obj$qtl$profile,
                                   qtl_info = viewpoly_obj$qtl$qtl_info,
                                   selected_mks = viewpoly_obj$qtl$selected_mks,
                                   pheno.col = 2:3,
                                   lgs.id = 2,
                                   by_range = TRUE,
                                   range.min = 30,
                                   range.max = 120,
                                   plot=TRUE,
                                   software = NULL)
  
  expect_equal(sum(qtl_profile_plot$data$SIG, na.rm = TRUE), 52.77, tolerance = 0.0001)
  
  maps <- lapply(viewpoly_obj$map$maps, function(x) {
    y <- x$l.dist
    names(y) <- x$mk.names
    y
  })
  
  # Get max size each chromosome
  expect_equal(map_summary(left.lim = 1,
                           right.lim = 50,
                           ch = 2,
                           maps = maps,
                           d.p1 = viewpoly_obj$map$d.p1,
                           d.p2 = viewpoly_obj$map$d.p2)[[2]], 185, tolerance = 0.0001)
  
  # Map summary table
  summary_table <- summary_maps(viewpoly_obj$map, software = "mappoly")
  expect_equal(sum(as.numeric(summary_table$`Map length (cM)`)), 661.98)
  expect_equal(sum(as.numeric(summary_table$Simplex)), 466)
  expect_equal(sum(as.numeric(summary_table$`Double-simplex`)), 522)
  expect_equal(sum(as.numeric(summary_table$`Max gap`)), 20.87)
  
})