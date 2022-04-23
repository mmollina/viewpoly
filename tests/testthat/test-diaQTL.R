test_that("Tests uploaded diaQTL files",{
  # upload diaQTL
  scan1_list <- scan1_summaries_list <- fitQTL_list <- BayesCI_list <- list()
  scan1_list$datapath <- tempfile()
  scan1_summaries_list$datapath <- tempfile()
  fitQTL_list$datapath <- tempfile()
  BayesCI_list$datapath <- tempfile()

  download.file("https://www.polyploids.org/sites/default/files/2022-04/tetra_diaQTL_scan1_list.RData", destfile = scan1_list$datapath)
  download.file("https://www.polyploids.org/sites/default/files/2022-04/tetra_diaQTL_scan1_summaries_list.RData", destfile = scan1_summaries_list$datapath)
  download.file("https://www.polyploids.org/sites/default/files/2022-04/tetra_diaQTL_fitQTL_list%20%281%29.RData", destfile = fitQTL_list$datapath)
  download.file("https://www.polyploids.org/sites/default/files/2022-04/tetra_diaQTL_BayesCI_list_0.RData", destfile = BayesCI_list$datapath)

  viewqtl_diaqtl <- prepare_diaQTL(scan1_list,
                                   scan1_summaries_list,
                                   fitQTL_list,
                                   BayesCI_list)

  expect_equal(check_viewqtl(viewqtl_obj = viewqtl_diaqtl),0)

  check_viewqtl_diaqtl_values(viewqtl_diaqtl, 402196.1, -2111.243, -10553.94, -1.689251)

  #VIEWqtl tests
  # plotly
  qtl_profile_plot <- plot_profile(viewqtl_diaqtl$profile,
                                   viewqtl_diaqtl$qtl_info,
                                   viewqtl_diaqtl$selected_mks,
                                   pheno.col = 1,
                                   lgs.id = 5,
                                   by_range = TRUE,
                                   range.min = 1,
                                   range.max = 50,
                                   plot=TRUE,
                                   software = "diaQTL")

  expect_equal(sum(qtl_profile_plot$data$SIG, na.rm = TRUE), 815.4399, tolerance = 0.0001)

  # by range
  qtl_profile_data <- plot_profile(viewqtl_diaqtl$profile,
                                   viewqtl_diaqtl$qtl_info,
                                   viewqtl_diaqtl$selected_mks,
                                   pheno.col = 2,
                                   lgs.id = 5,
                                   by_range = TRUE,
                                   range.min = 30,
                                   range.max = 120,
                                   plot=FALSE,
                                   software = "diaQTL")

  expect_equal(sum(qtl_profile_data$lines$SIG, na.rm = TRUE), 4780.668, tolerance = 0.001)
  expect_equal(sum(qtl_profile_data$lines$`Position (cM)`), 21507.75, tolerance = 0.001)
  expect_equal(as.numeric(qtl_profile_data$points$PVAL), -663.8152, tolerance = 0.001)
  expect_equal(as.numeric(qtl_profile_data$points$INF), 24.18, tolerance = 0.001)
  expect_equal(as.numeric(qtl_profile_data$points$SUP), 32.15, tolerance = 0.001)

  # export data
  qtl_profile_data <- plot_profile(viewqtl_diaqtl$profile,
                                   viewqtl_diaqtl$qtl_info,
                                   viewqtl_diaqtl$selected_mks,
                                   pheno.col = 1:2,
                                   lgs.id = 5,
                                   by_range = FALSE,
                                   range.min = NULL,
                                   range.max = NULL,
                                   plot=FALSE,
                                   software = NULL)

  expect_equal(sum(qtl_profile_data$lines$SIG), 17452.68, tolerance = 0.001)
  expect_equal(sum(qtl_profile_data$lines$`Position (cM)`), 43015.5, tolerance = 0.001)
  expect_equal(as.numeric(qtl_profile_data$points$PVAL), c(-615.2460, -663.8152), tolerance = 0.001)
  expect_equal(as.numeric(qtl_profile_data$points$INF), c(22.43, 24.18), tolerance = 0.001)
  expect_equal(as.numeric(qtl_profile_data$points$SUP), c(28.62, 32.15), tolerance = 0.001)

  # plot exported data
  p <- only_plot_profile(qtl_profile_data)
  expect_equal(sum(p$data$SIG), 17452.68, tolerance = 0.001)

  # effects graphics
  p <- data_effects(qtl_info = viewqtl_diaqtl$qtl_info,
                    effects = viewqtl_diaqtl$effects,
                    pheno.col = "FM07",
                    p1 = "P1",
                    p2 = "P2",
                    lgs = 5,
                    groups = 5,
                    position = 26.19,
                    software = "diaQTL",
                    design = "circle")

  expect_equal(sum(p[[1]]$data$Estimates), -0.744064, tolerance = 0.001)
  expect_equal(names(p[[1]]$data),
               c("Estimates", "CI.lower", "CI.upper", "Alleles", "Parent", "Effects", "pheno", "qtl_id", "LG", "Pos", "unique.id"),
               tolerance = 0.001)

  p <- data_effects(qtl_info = viewqtl_diaqtl$qtl_info,
                    effects = viewqtl_diaqtl$effects,
                    pheno.col = "FM07",
                    p1 = "P1",
                    p2 = "P2",
                    lgs = 5,
                    groups = 5,
                    position = 26.19,
                    software = "diaQTL",
                    design = "digenic")

  expect_equal(sum(p[[1]]$data$z), -0.6900337, tolerance = 0.001)
  expect_equal(names(p[[1]]$data),
               c("x", "y", "z"),
               tolerance = 0.001)

  p <- data_effects(qtl_info = viewqtl_diaqtl$qtl_info,
                    effects = viewqtl_diaqtl$effects,
                    pheno.col = "FM07",
                    p1 = "P1",
                    p2 = "P2",
                    lgs = 5,
                    groups = 5,
                    position = 26.19,
                    software = "diaQTL",
                    design = "bar")

  expect_equal(sum(p[[1]]$data$Estimates), -8.881784e-16, tolerance = 0.001)
  expect_equal(names(p[[1]]$data),
               c("Estimates", "CI.lower", "CI.upper", "Alleles", "Parent", "Effects"),
               tolerance = 0.001)

})