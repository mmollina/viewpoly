# test_that("upload files",{
#   # upload examples
#   viewpoly_obj <- prepare_examples("tetra_map")
#   
#   expect_equal(check_viewpoly(viewpoly_obj),0)
#   
#   check_viewmap_values(viewpoly_obj$map, 
#                        c(14, 132, 139, 157, 34), 
#                        c(36, 167, 164, 109), 
#                        50502.07)
#   
#   check_viewqtl_qtlpoly_values(viewpoly_obj$qtl, 
#                                116504, 
#                                5.418909, 
#                                -1.067457e-11, 
#                                299.6155, 
#                                0.000160791, 
#                                2.340129e-12, 
#                                1)
#   
#   # VIEWmap tests 
#   qtl_profile_plot <- plot_profile(profile = viewpoly_obj$qtl$profile, 
#                                    qtl_info = viewpoly_obj$qtl$qtl_info, 
#                                    selected_mks = viewpoly_obj$qtl$selected_mks, 
#                                    pheno.col = 2:3, 
#                                    lgs.id = 2, 
#                                    by_range = TRUE, 
#                                    range.min = 30, 
#                                    range.max = 120, 
#                                    plot=TRUE, 
#                                    software = NULL)
#   
#   expect_equal(sum(qtl_profile_plot$data$SIG, na.rm = TRUE), 84.46874, tolerance = 0.0001)
#   
#   maps <- lapply(viewpoly_obj$map$maps, function(x) {
#     y <- x$l.dist
#     names(y) <- x$mk.names
#     y
#   })
#   
#   expect_doppelganger("linkage map draw",   draw_map_shiny(left.lim = 1, 
#                                                            right.lim = 50, 
#                                                            ch = 1,
#                                                            d.p1 = viewpoly_obj$map$d.p1,
#                                                            d.p2 = viewpoly_obj$map$d.p2, 
#                                                            maps = maps, 
#                                                            ph.p1 = viewpoly_obj$map$ph.p1, 
#                                                            ph.p2 = viewpoly_obj$map$ph.p2,
#                                                            snp.names = FALSE))
#   
#   expect_doppelganger("linkage map draw names",   draw_map_shiny(left.lim = 1, 
#                                                                  right.lim = 50, 
#                                                                  ch = 1,
#                                                                  d.p1 = viewpoly_obj$map$d.p1,
#                                                                  d.p2 = viewpoly_obj$map$d.p2, 
#                                                                  maps = maps, 
#                                                                  ph.p1 = viewpoly_obj$map$ph.p1, 
#                                                                  ph.p2 = viewpoly_obj$map$ph.p2,
#                                                                  snp.names = TRUE))
#   
#   expect_doppelganger("plot map list", plot_map_list(viewpoly_obj$map))
#   
#   # Get max size each chromosome
#   expect_equal(map_summary(left.lim = 1, 
#                            right.lim = 50, 
#                            ch = 3, 
#                            maps = maps, 
#                            d.p1 = viewpoly_obj$map$d.p1, 
#                            d.p2 = viewpoly_obj$map$d.p2)[[5]], 134.073, tolerance = 0.0001)
#   
#   # Map summary table
#   summary_table <- summary_maps(viewpoly_obj$map)
#   expect_equal(sum(as.numeric(summary_table$`Map length (cM)`)), 3259.98)
#   expect_equal(sum(as.numeric(summary_table$Simplex)), 2450)
#   expect_equal(sum(as.numeric(summary_table$`Double-simplex`)), 1820)
#   expect_equal(sum(as.numeric(summary_table$`Max gap`)), 80.51)
#   
#   #VIEWqtl tests
#   expect_doppelganger("qtl plot",plot_profile(viewpoly_obj$qtl$profile, 
#                                               viewpoly_obj$qtl$qtl_info, 
#                                               viewpoly_obj$qtl$selected_mks, 
#                                               pheno.col = 2, 
#                                               lgs.id = 2, 
#                                               by_range = FALSE, 
#                                               plot=TRUE, 
#                                               software = NULL))
#   
#   # by range
#   qtl_profile_data <- plot_profile(viewpoly_obj$qtl$profile, 
#                                    viewpoly_obj$qtl$qtl_info, 
#                                    viewpoly_obj$qtl$selected_mks, 
#                                    pheno.col = 2, 
#                                    lgs.id = 2, 
#                                    by_range = TRUE, 
#                                    range.min = 30, 
#                                    range.max = 120, 
#                                    plot=FALSE, 
#                                    software = NULL)
#   
#   expect_equal(sum(qtl_profile_data$lines$SIG, na.rm = TRUE), 43.81917, tolerance = 0.001)
#   expect_equal(sum(qtl_profile_data$lines$`Position (cM)`), 8000.109, tolerance = 0.001)
#   expect_equal(as.numeric(qtl_profile_data$points$PVAL), 0.000141, tolerance = 0.001)
#   expect_equal(as.numeric(qtl_profile_data$points$H2), 0.17, tolerance = 0.001)
#   expect_equal(as.numeric(qtl_profile_data$points$INF), 41, tolerance = 0.001)
#   expect_equal(as.numeric(qtl_profile_data$points$SUP), 119, tolerance = 0.001)
#   
#   # export data
#   qtl_profile_data <- plot_profile(viewpoly_obj$qtl$profile, 
#                                    viewpoly_obj$qtl$qtl_info, 
#                                    viewpoly_obj$qtl$selected_mks, 
#                                    pheno.col = 2, 
#                                    lgs.id = 2, 
#                                    by_range = FALSE, 
#                                    range.min = NULL, 
#                                    range.max = NULL, 
#                                    plot=FALSE, 
#                                    software = NULL)
#   
#   expect_equal(sum(qtl_profile_data$lines$SIG), 292.883, tolerance = 0.001)
#   expect_equal(sum(qtl_profile_data$lines$`Position (cM)`), 8000.109, tolerance = 0.001)
#   expect_equal(as.numeric(qtl_profile_data$points$PVAL), 0.000141, tolerance = 0.001)
#   expect_equal(as.numeric(qtl_profile_data$points$H2), 0.17, tolerance = 0.001)
#   expect_equal(as.numeric(qtl_profile_data$points$INF), 41, tolerance = 0.001)
#   expect_equal(as.numeric(qtl_profile_data$points$SUP), 119, tolerance = 0.001)
#   
#   # plot exported data
#   p <- only_plot_profile(qtl_profile_data)
#   expect_equal(sum(p$data$SIG), 292.883, tolerance = 0.001)
#   
#   # effects graphics
#   p <- data_effects(qtl_info = viewpoly_obj$qtl$qtl_info, 
#                     effects = viewpoly_obj$qtl$effects, 
#                     pheno.col = "SG06", 
#                     p1 = "P1", 
#                     p2 = "P2",  
#                     lgs = 2, 
#                     groups = 2, 
#                     position = 77, 
#                     software = "QTLpoly", 
#                     design = "circle")
#   
#   expect_doppelganger("effects circle", plot_effects(p, "QTLpoly", "circle"))
#   
#   expect_equal(sum(p[[1]]$data$Estimates), -0.0436829, tolerance = 0.001)
#   expect_equal(names(p[[1]]$data), 
#                c("Estimates", "Alleles", "Parent", "Effects", "pheno", "qtl_id", "LG", "Pos", "unique.id"), 
#                tolerance = 0.001)
#   
#   p <- data_effects(qtl_info = viewpoly_obj$qtl$qtl_info, 
#                     effects = viewpoly_obj$qtl$effects, 
#                     pheno.col = "SG06", 
#                     p1 = "P1", 
#                     p2 = "P2",  
#                     lgs = 2, 
#                     groups = 2, 
#                     position = 77, 
#                     software = "QTLpoly", 
#                     design = "digenic")
#   
#   expect_equal(sum(p[[1]]$data$z), 1.528847e-14, tolerance = 0.001)
#   expect_equal(names(p[[1]]$data), 
#                c("x", "y", "z"), 
#                tolerance = 0.001)
#   
#   expect_doppelganger("effects digenic", plot_effects(p, "QTLpoly", "digenic"))
#   
#   p <- data_effects(qtl_info = viewpoly_obj$qtl$qtl_info, 
#                     effects = viewpoly_obj$qtl$effects, 
#                     pheno.col = "SG06", 
#                     p1 = "P1", 
#                     p2 = "P2",  
#                     lgs = 2, 
#                     groups = 2, 
#                     position = 77, 
#                     software = "QTLpoly", 
#                     design = "bar")
#   
#   expect_equal(sum(p[[1]]$data$Estimates), 2.184058e-15, tolerance = 0.001)
#   expect_equal(names(p[[1]]$data), 
#                c("Estimates", "Alleles", "Parent", "Effects"), 
#                tolerance = 0.001)
#   
#   expect_doppelganger("effects bar", plot_effects(p, "QTLpoly", "bar"))
#   
#   # breeding values table
#   pos <- split(viewpoly_obj$qtl$qtl_info[1:3,]$Pos, viewpoly_obj$qtl$qtl_info[1:3,]$pheno)
#   breed.values <- breeding_values(viewpoly_obj$qtl$qtl_info, 
#                                   viewpoly_obj$qtl$probs, 
#                                   viewpoly_obj$qtl$selected_mks, 
#                                   viewpoly_obj$qtl$blups, 
#                                   viewpoly_obj$qtl$beta.hat, 
#                                   pos)
#   
#   expect_equal(sum(breed.values$PY06), 5.26)
#   expect_equal(sum(breed.values$SG06), 5.36)
#   
#   # get and plot homologs prob
#   data.prob <- calc_homologprob(probs = viewpoly_obj$qtl$probs, 
#                                 viewpoly_obj$qtl$selected_mks, 
#                                 1:5)
#   
#   expect_equal(sum(data.prob$homoprob$probability), 14900, tolerance = 0.001)
#   
#   input.haplo <- c("Trait:SG06_LG:2_Pos:77_homolog:P1.1")
#   p <- select_haplo(input.haplo, 
#                     viewpoly_obj$qtl$probs, 
#                     viewpoly_obj$qtl$selected_mks, 
#                     effects.data = p)
#   expect_equal(sum(p[[1]]$data$probability), 507.9996, tolerance = 0.0001)
#   expect_equal(sum(p[[2]]$data$probability), 508.001, tolerance = 0.0001)
#   expect_equal(sum(p[[3]]$data$probability), 508.0009, tolerance = 0.0001)
#   
#   # VIEWgenome tests
#   p <- plot_cm_mb(viewpoly_obj$map, 1, 1,50)
#   
#   expect_equal(sum(p$data$l.dist), 50502.07, tolerance = 0.001)
#   
# })
# 
