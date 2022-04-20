context("upload functions")

test_that("upload examples",{
  viewpoly_obj <- viewpoly:::prepare_examples("tetra_map")
  
  expect_equal(viewpoly:::check_viewpoly(viewpoly_obj),0)
})

test_that("upload custom files",{
  dosages <- phases <- genetic_map <- mks_pos <- list()
  dosages$datapath <- system.file("ext/dosage.tsv.gz", package = "viewpoly")
  phases$datapath <- system.file("ext/phases.tsv.gz", package = "viewpoly")
  genetic_map$datapath <- system.file("ext/map.tsv.gz", package = "viewpoly")

  viewmap_obj <- viewpoly:::prepare_map_custom_files(dosages, phases, genetic_map)
  
  expect_equal(viewpoly:::check_viewmap(viewmap_obj),0)
  
  selected_mks <- qtl_info <- blups <- beta.hat <- profile <- effects <- probs <- list()
  selected_mks$datapath = system.file("ext/selected_mks.tsv.gz", package = "viewpoly")
  qtl_info$datapath = system.file("ext/qtl_info.tsv.gz", package = "viewpoly")
  blups$datapath = system.file("ext/blups.tsv.gz", package = "viewpoly")
  beta.hat$datapath = system.file("ext/beta.hat.tsv.gz", package = "viewpoly")
  profile$datapath = system.file("ext/profile.tsv.gz", package = "viewpoly")
  effects$datapath = system.file("ext/effects.tsv.gz", package = "viewpoly")
  probs$datapath = system.file("ext/probs.tsv.gz", package = "viewpoly")
  
  viewqtl_obj <- viewpoly:::prepare_qtl_custom_files(selected_mks, qtl_info, blups, 
                                                     beta.hat, profile, effects, probs)
  
  expect_equal(viewpoly:::check_viewqtl(viewqtl_obj),0)
  
})

# test_that("upload MAPpoly",{
#   prepare_MAPpoly()
# })
#  
# test_that("upload polymapR",{
#   prepare_polymapR()
# })
# 
# test_that("upload QTLpoly",{
#   prepare_QTLpoly()
# })
# 
# test_that("upload diaQTL",{
#   prepare_diaQTL()
# })
# 
# test_that("upload polyqtlR",{
#   prepare_polyqtlR()
# })

