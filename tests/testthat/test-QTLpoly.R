test_that("Tests uploaded QTLpoly files",{
  source(system.file("ext/functions4tests.R", package = "viewpoly"))

  # upload QTLpoly
  input.data <- remim.mod <- est.effects <- fitted.mod <- list()
  input.data$datapath <- tempfile()
  remim.mod$datapath <- tempfile()
  est.effects$datapath <- tempfile()
  fitted.mod$datapath <- tempfile()
  
  if(havingIP()){
    options(timeout=200)
    download.file("https://www.polyploids.org/sites/default/files/2022-04/tetra_QTLpoly_effects.RData", destfile = est.effects$datapath)
    download.file("https://www.polyploids.org/sites/default/files/2022-04/tetra_QTLpoly_data.RData", destfile = input.data$datapath)
    download.file("https://www.polyploids.org/sites/default/files/2022-04/tetra_QTLpoly_remim.RData", destfile = remim.mod$datapath)
    download.file("https://www.polyploids.org/sites/default/files/2022-04/tetra_QTLpoly_fitted.RData", destfile = fitted.mod$datapath)
    
    viewqtl_qtlpoly <- prepare_QTLpoly(data = input.data,
                                       remim.mod = remim.mod,
                                       est.effects = est.effects,
                                       fitted.mod = fitted.mod)
    
    expect_equal(check_viewqtl(viewqtl_qtlpoly),0)
    
    check_viewqtl_qtlpoly_values(viewqtl_qtlpoly,
                                 116504,
                                 5.418909,
                                 -1.067457e-11,
                                 299.6155,
                                 0.000160791,
                                 2.340129e-12,
                                 1)
    
    #VIEWqtl tests
    # plotly
    qtl_profile_plot <- plot_profile(viewqtl_qtlpoly$profile,
                                     viewqtl_qtlpoly$qtl_info,
                                     viewqtl_qtlpoly$selected_mks,
                                     pheno.col = 2,
                                     lgs.id = 2,
                                     by_range = TRUE,
                                     range.min = 30,
                                     range.max = 120,
                                     plot=TRUE,
                                     software = NULL)
    
    expect_equal(sum(qtl_profile_plot$data$SIG, na.rm = TRUE), 43.81917, tolerance = 0.0001)
    
    # by range
    qtl_profile_data <- plot_profile(viewqtl_qtlpoly$profile,
                                     viewqtl_qtlpoly$qtl_info,
                                     viewqtl_qtlpoly$selected_mks,
                                     pheno.col = 2,
                                     lgs.id = 2,
                                     by_range = TRUE,
                                     range.min = 30,
                                     range.max = 120,
                                     plot=FALSE,
                                     software = NULL)
    
    expect_equal(sum(qtl_profile_data$lines$SIG, na.rm = TRUE), 43.81917, tolerance = 0.001)
    expect_equal(sum(qtl_profile_data$lines$`Position (cM)`), 8000.109, tolerance = 0.001)
    expect_equal(as.numeric(qtl_profile_data$points$PVAL), 0.000141, tolerance = 0.001)
    expect_equal(as.numeric(qtl_profile_data$points$H2), 0.17, tolerance = 0.001)
    expect_equal(as.numeric(qtl_profile_data$points$INF), 41, tolerance = 0.001)
    expect_equal(as.numeric(qtl_profile_data$points$SUP), 119, tolerance = 0.001)
    
    # export data
    qtl_profile_data <- plot_profile(viewqtl_qtlpoly$profile,
                                     viewqtl_qtlpoly$qtl_info,
                                     viewqtl_qtlpoly$selected_mks,
                                     pheno.col = 2,
                                     lgs.id = 2,
                                     by_range = FALSE,
                                     range.min = NULL,
                                     range.max = NULL,
                                     plot=FALSE,
                                     software = NULL)
    
    expect_equal(sum(qtl_profile_data$lines$SIG), 292.883, tolerance = 0.001)
    expect_equal(sum(qtl_profile_data$lines$`Position (cM)`), 8000.109, tolerance = 0.001)
    expect_equal(as.numeric(qtl_profile_data$points$PVAL), 0.000141, tolerance = 0.001)
    expect_equal(as.numeric(qtl_profile_data$points$H2), 0.17, tolerance = 0.001)
    expect_equal(as.numeric(qtl_profile_data$points$INF), 41, tolerance = 0.001)
    expect_equal(as.numeric(qtl_profile_data$points$SUP), 119, tolerance = 0.001)
    
    # plot exported data
    p <- only_plot_profile(qtl_profile_data)
    expect_equal(sum(p$data$SIG), 292.883, tolerance = 0.001)
    
    # effects graphics
    p <- data_effects(qtl_info = viewqtl_qtlpoly$qtl_info,
                      effects = viewqtl_qtlpoly$effects,
                      pheno.col = "SG06",
                      lgs = 2,
                      groups = 2,
                      position = 77,
                      software = "QTLpoly",
                      design = "circle")
    
    expect_equal(sum(p[[1]]$data$Estimates), -0.0436829, tolerance = 0.001)
    expect_equal(names(p[[1]]$data),
                 c("Estimates", "Alleles", "Parent", "Effects", "pheno", "qtl_id", "LG", "Pos", "unique.id"),
                 tolerance = 0.001)
    
    p <- data_effects(qtl_info = viewqtl_qtlpoly$qtl_info,
                      effects = viewqtl_qtlpoly$effects,
                      pheno.col = "SG06",
                      lgs = 2,
                      groups = 2,
                      position = 77,
                      software = "QTLpoly",
                      design = "digenic")
    
    expect_equal(sum(p[[1]]$data$z), 1.528847e-14, tolerance = 0.001)
    expect_equal(names(p[[1]]$data),
                 c("x", "y", "z"),
                 tolerance = 0.001)
    
    p <-  data_effects(qtl_info = viewqtl_qtlpoly$qtl_info,
                       effects = viewqtl_qtlpoly$effects,
                       pheno.col = "SG06",
                       lgs = 2,
                       groups = 2,
                       position = 77,
                       software = "QTLpoly",
                       design = "bar")
    
    expect_equal(sum(p[[1]]$data$Estimates), 2.184058e-15, tolerance = 0.001)
    expect_equal(names(p[[1]]$data),
                 c("Estimates", "Alleles", "Parent", "Effects"),
                 tolerance = 0.001)
    
    # breeding values table
    pos <- split(viewqtl_qtlpoly$qtl_info[1:3,]$Pos, viewqtl_qtlpoly$qtl_info[1:3,]$pheno)
    breed.values <- breeding_values(viewqtl_qtlpoly$qtl_info,
                                    viewqtl_qtlpoly$probs,
                                    viewqtl_qtlpoly$selected_mks,
                                    viewqtl_qtlpoly$blups,
                                    viewqtl_qtlpoly$beta.hat,
                                    pos)
    
    expect_equal(sum(breed.values$PY06), 155.63)
    expect_equal(sum(breed.values$SG06), 167.16)
    
    # get and plot homologs prob
    data.prob <- calc_homologprob(probs = viewqtl_qtlpoly$probs,
                                  viewqtl_qtlpoly$selected_mks,
                                  1:5)
    
    expect_equal(sum(data.prob$homoprob$probability), 464880, tolerance = 0.001)
    
    input.haplo <- list("Trait:SG06_LG:2_Pos:77_homolog:P1.1", "Trait:FM07_LG:5_Pos:26_homolog:P1.3",
                        "Trait:SG06_LG:2_Pos:77_homolog:P1.3", "Trait:FM07_LG:5_Pos:26_homolog:P2.3")
    
    p1.list <- select_haplo(input.haplo,
                            viewqtl_qtlpoly$probs,
                            viewqtl_qtlpoly$selected_mks,
                            effects.data = p)
    p1 <- p1.list[[1]]
    
    # Test exclude
    input.haplo <- list("Trait:PY06_LG:5_Pos:29_homolog:P1.1")
    exclude.haplo <- list("Trait:FM07_LG:5_Pos:26_homolog:P1.4")
    
    p1.list <- select_haplo(input.haplo = input.haplo,
                            exclude.haplo = exclude.haplo, 
                            probs = viewqtl_qtlpoly$probs, 
                            selected_mks = viewqtl_qtlpoly$selected_mks,
                            effects.data = p)
    p1 <- p1.list[[1]]
    
    expect_equal(sum(p1[[1]]$data$probability), 431.9998, tolerance = 0.0001)
    expect_equal(sum(p1[[2]]$data$probability), 432.0001, tolerance = 0.0001)
    expect_equal(sum(p1[[3]]$data$probability), 432.0005, tolerance = 0.0001)
  } else {
    print("QTLpoly tests are only executed if internet conection is available.")
  }
})
