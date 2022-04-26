# test_that("Tests with custom files",{
# 
#   # upload custom files
#   dosages <- phases <- genetic_map <- mks_pos <- list()
#   dosages$datapath <- system.file("ext/dosage.tsv.gz", package = "viewpoly")
#   phases$datapath <- system.file("ext/phases.tsv.gz", package = "viewpoly")
#   genetic_map$datapath <- system.file("ext/map.tsv.gz", package = "viewpoly")
# 
#   viewmap_obj <- prepare_map_custom_files(dosages, phases, genetic_map)
# 
#   expect_equal(check_viewmap(viewmap_obj),0)
# 
#   check_viewmap_values(viewmap_obj,
#                        c(887, 1063,543,188,54,10),
#                        c(181, 640, 671, 608, 645),
#                        401368.7)
# 
#   selected_mks <- qtl_info <- blups <- beta.hat <- profile <- effects <- probs <- list()
#   selected_mks$datapath = system.file("ext/selected_mks.tsv.gz", package = "viewpoly")
#   qtl_info$datapath = system.file("ext/qtl_info.tsv.gz", package = "viewpoly")
#   blups$datapath = system.file("ext/blups.tsv.gz", package = "viewpoly")
#   beta.hat$datapath = system.file("ext/beta.hat.tsv.gz", package = "viewpoly")
#   profile$datapath = system.file("ext/profile.tsv.gz", package = "viewpoly")
#   effects$datapath = system.file("ext/effects.tsv.gz", package = "viewpoly")
#   probs$datapath = system.file("ext/probs.tsv.gz", package = "viewpoly")
# 
#   viewqtl_obj <- prepare_qtl_custom_files(selected_mks, qtl_info, blups,
#                                           beta.hat, profile, effects, probs)
# 
#   expect_equal(check_viewqtl(viewqtl_obj),0)
# 
#   check_viewqtl_qtlpoly_values(viewqtl_obj,
#                                245183.2,
#                                6.702273,
#                                3.091306e-10,
#                                341.4562,
#                                4.29507e-05,
#                                -7.647564e-11,
#                                1)
# 
# })
