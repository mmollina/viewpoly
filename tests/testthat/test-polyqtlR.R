test_that("Tests uploaded polyqtlR files",{
  # upload polyqtlR
  polyqtlR_QTLscan_list <- polyqtlR_qtl_info <- polyqtlR_effects <- list()
  polyqtlR_QTLscan_list$datapath <- tempfile()
  polyqtlR_qtl_info$datapath <- tempfile()
  polyqtlR_effects$datapath <- tempfile()

  download.file("https://www.polyploids.org/sites/default/files/2022-04/tetra_polyqtlR_QTLscan.RData", destfile = polyqtlR_QTLscan_list$datapath)
  download.file("https://www.polyploids.org/sites/default/files/2022-04/tetra_polyqtlR_qtl_info.RData", destfile = polyqtlR_qtl_info$datapath)
  download.file("https://www.polyploids.org/sites/default/files/2022-04/tetra_polyqtlR_effects.RData", destfile = polyqtlR_effects$datapath)

  viewqtl_polyqtlr <- prepare_polyqtlR(polyqtlR_QTLscan_list, polyqtlR_qtl_info, polyqtlR_effects)

  expect_equal(check_viewqtl(viewqtl_obj = viewqtl_polyqtlr),0)

  check_viewqtl_polyqtlr_values(viewqtl_polyqtlr, 287322.2, 103.8376, 210404.7, 6032.891)

  #VIEWqtl tests
  # plotly
  qtl_profile_plot <- plot_profile(viewqtl_polyqtlr$profile,
                                   viewqtl_polyqtlr$qtl_info,
                                   viewqtl_polyqtlr$selected_mks,
                                   pheno.col = 1:5,
                                   lgs.id = 1,
                                   by_range = TRUE,
                                   range.min = 30,
                                   range.max = 120,
                                   plot=TRUE,
                                   software = "polyqtlR")

  expect_equal(sum(qtl_profile_plot$data$SIG, na.rm = TRUE), 2322.896, tolerance = 0.0001)

  # by range
  qtl_profile_data <- plot_profile(viewqtl_polyqtlr$profile,
                                   viewqtl_polyqtlr$qtl_info,
                                   viewqtl_polyqtlr$selected_mks,
                                   pheno.col = 1:5,
                                   lgs.id = 1,
                                   by_range = TRUE,
                                   range.min = 30,
                                   range.max = 120,
                                   plot=FALSE,
                                   software = "polyqtlR")

  expect_equal(sum(qtl_profile_data$lines$SIG, na.rm = TRUE), 2322.896, tolerance = 0.001)
  expect_equal(sum(qtl_profile_data$lines$`Position (cM)`), 252510.3, tolerance = 0.001)
  expect_equal(as.numeric(qtl_profile_data$points$PVAL), 5.566673, tolerance = 0.001)

  # export data
  qtl_profile_data <- plot_profile(viewqtl_polyqtlr$profile,
                                   viewqtl_polyqtlr$qtl_info,
                                   viewqtl_polyqtlr$selected_mks,
                                   pheno.col = 1:5,
                                   lgs.id = 1:5,
                                   by_range = FALSE,
                                   range.min = NULL,
                                   range.max = NULL,
                                   plot=FALSE,
                                   software = "polyqtlR")

  expect_equal(sum(qtl_profile_data$lines$SIG), 17065.12, tolerance = 0.001)
  expect_equal(sum(qtl_profile_data$lines$`Position (cM)`), 728869.4, tolerance = 0.001)
  expect_equal(as.numeric(qtl_profile_data$points$PVAL), 5.566673, tolerance = 0.001)

  # plot exported data
  p <- only_plot_profile(qtl_profile_data)
  expect_equal(sum(p$data$SIG),  17065.12, tolerance = 0.001)

  # effects graphics
  expect_error(data_effects(qtl_info = viewqtl_polyqtlr$qtl_info,
                            effects = viewqtl_polyqtlr$effects,
                            pheno.col = "SG06",
                            p1 = "P1",
                            p2 = "P2",
                            lgs = 2,
                            groups = 2,
                            position = 77,
                            software = "polyqtlR",
                            design = "circle"))


  expect_error(data_effects(qtl_info = viewqtl_polyqtlr$qtl_info,
                            effects = viewqtl_polyqtlr$effects,
                            pheno.col = "SG06",
                            p1 = "P1",
                            p2 = "P2",
                            lgs = 2,
                            groups = 2,
                            position = 77,
                            software = "polyqtlR",
                            design = "digenic"))


  p <- data_effects(qtl_info = viewqtl_polyqtlr$qtl_info,
                    effects = viewqtl_polyqtlr$effects,
                    pheno.col = "SG06",
                    p1 = "P1",
                    p2 = "P2",
                    lgs = 2,
                    groups = 2,
                    position = 77,
                    software = "polyqtlR",
                    design = "bar")

  expect_equal(sum(p[[1]]$data$effect), -0.1016605, tolerance = 0.001)
  expect_equal(names(p[[1]]$data),
               c("pos", "pheno", "LG", "haplo", "effect", "x.axis"),
               tolerance = 0.001)

})