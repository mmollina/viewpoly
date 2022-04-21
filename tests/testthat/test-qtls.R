test_that("QTL graphics and tables",{
  
  viewpoly_obj <- prepare_examples("tetra_map")
  
  # by range
  qtl_profile_data <- plot_profile(viewpoly_obj$qtl$profile, 
                                   viewpoly_obj$qtl$qtl_info, 
                                   viewpoly_obj$qtl$selected_mks, 
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
  qtl_profile_data <- plot_profile(viewpoly_obj$qtl$profile, 
                                   viewpoly_obj$qtl$qtl_info, 
                                   viewpoly_obj$qtl$selected_mks, 
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
  p <- data_effects(qtl_info = viewpoly_obj$qtl$qtl_info, 
                    effects = viewpoly_obj$qtl$effects, 
                    pheno.col = "SG06", 
                    p1 = "P1", 
                    p2 = "P2",  
                    lgs = 2, 
                    groups = 2, 
                    position = 77, 
                    software = "QTLpoly", 
                    design = "circle")
  
  expect_equal(sum(p[[1]]$data$Estimates), -0.0436829, tolerance = 0.001)
  expect_equal(names(p[[1]]$data), 
               c("Estimates", "Alleles", "Parent", "Effects", "pheno", "qtl_id", "LG", "Pos", "unique.id"), 
               tolerance = 0.001)
  
  p <- data_effects(qtl_info = viewpoly_obj$qtl$qtl_info, 
                    effects = viewpoly_obj$qtl$effects, 
                    pheno.col = "SG06", 
                    p1 = "P1", 
                    p2 = "P2",  
                    lgs = 2, 
                    groups = 2, 
                    position = 77, 
                    software = "QTLpoly", 
                    design = "digenic")
  
  expect_equal(sum(p[[1]]$data$z), 1.528847e-14, tolerance = 0.001)
  expect_equal(names(p[[1]]$data), 
               c("x", "y", "z"), 
               tolerance = 0.001)
  
  p <- data_effects(qtl_info = viewpoly_obj$qtl$qtl_info, 
                    effects = viewpoly_obj$qtl$effects, 
                    pheno.col = "SG06", 
                    p1 = "P1", 
                    p2 = "P2",  
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
  pos <- split(viewpoly_obj$qtl$qtl_info[1:3,]$Pos, viewpoly_obj$qtl$qtl_info[1:3,]$pheno)
  breed.values <- breeding_values(viewpoly_obj$qtl$qtl_info, 
                                  viewpoly_obj$qtl$probs, 
                                  viewpoly_obj$qtl$selected_mks, 
                                  viewpoly_obj$qtl$blups, 
                                  viewpoly_obj$qtl$beta.hat, 
                                  pos)
  
  expect_equal(sum(breed.values$PY06), 5.26)
  expect_equal(sum(breed.values$SG06), 5.36)
  
  # get and plot homologs prob
  data.prob <- calc_homologprob(probs = viewpoly_obj$qtl$probs, 
                                viewpoly_obj$qtl$selected_mks, 
                                1:5)
  
  expect_equal(sum(data.prob$homoprob$probability), 14900, tolerance = 0.001)

  input.haplo <- c("Trait:SG06_LG:2_Pos:77_homolog:P1.1")
  p <- select_haplo(input.haplo, 
                    viewpoly_obj$qtl$probs, 
                    viewpoly_obj$qtl$selected_mks, 
                    effects.data = p)
  expect_equal(sum(p[[1]]$data$probability), 507.9996, tolerance = 0.0001)
  expect_equal(sum(p[[2]]$data$probability), 508.001, tolerance = 0.0001)
  expect_equal(sum(p[[3]]$data$probability), 508.0009, tolerance = 0.0001)
})
