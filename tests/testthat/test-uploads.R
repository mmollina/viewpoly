test_that("upload files",{
  check_viewmap_values <- function(viewmap_obj, doses, phases, maps){
    expect_equal(as.vector(table(viewmap_obj$d.p1[[1]])), doses)
    expect_equal(as.vector(table(viewmap_obj$ph.p1[[1]][,3])), phases)
    expect_equal(as.vector(sum(viewmap_obj$maps[[1]][,2])), maps, tolerance = 0.001)
  }
  
  check_viewqtl_values <- function(viewqtl_obj, pos, h2, u.hat, beta.hat, lod, effect, probs){
    expect_equal(sum(viewqtl_obj$selected_mks$pos), pos, tolerance = 0.001)
    expect_equal(sum(viewqtl_obj$qtl_info$h2), h2, tolerance = 0.001)
    expect_equal(sum(viewqtl_obj$blups$u.hat), u.hat, tolerance = 0.001)
    expect_equal(sum(viewqtl_obj$beta.hat$beta.hat), beta.hat, tolerance = 0.001)
    expect_equal(min(viewqtl_obj$profile$LOP), lod, tolerance = 0.001)
    expect_equal(sum(viewqtl_obj$effects$effect), effect, tolerance = 0.001)
    expect_equal(sum(viewqtl_obj$probs[,1,2]), probs, tolerance = 0.001)
  }
  
  # upload examples
  viewpoly_obj <- prepare_examples("tetra_map")
  
  expect_equal(check_viewpoly(viewpoly_obj),0)

  check_viewmap_values(viewpoly_obj$map, 
                       c(14, 132, 139, 157, 34), 
                       c(36, 167, 164, 109), 
                       50502.07)
  
  check_viewqtl_values(viewpoly_obj$qtl, 
                       116504, 
                       5.418909, 
                       -1.067457e-11, 
                       299.6155, 
                       0.000160791, 
                       2.340129e-12, 
                       1)
  
  # upload custom files
  dosages <- phases <- genetic_map <- mks_pos <- list()
  dosages$datapath <- system.file("ext/dosage.tsv.gz", package = "viewpoly")
  phases$datapath <- system.file("ext/phases.tsv.gz", package = "viewpoly")
  genetic_map$datapath <- system.file("ext/map.tsv.gz", package = "viewpoly")

  viewmap_obj <- prepare_map_custom_files(dosages, phases, genetic_map)
  
  expect_equal(check_viewmap(viewmap_obj),0)
  
  check_viewmap_values(viewmap_obj, 
                       c(887, 1063,543,188,54,10), 
                       c(181, 640, 671, 608, 645), 
                       401368.7)
  
  selected_mks <- qtl_info <- blups <- beta.hat <- profile <- effects <- probs <- list()
  selected_mks$datapath = system.file("ext/selected_mks.tsv.gz", package = "viewpoly")
  qtl_info$datapath = system.file("ext/qtl_info.tsv.gz", package = "viewpoly")
  blups$datapath = system.file("ext/blups.tsv.gz", package = "viewpoly")
  beta.hat$datapath = system.file("ext/beta.hat.tsv.gz", package = "viewpoly")
  profile$datapath = system.file("ext/profile.tsv.gz", package = "viewpoly")
  effects$datapath = system.file("ext/effects.tsv.gz", package = "viewpoly")
  probs$datapath = system.file("ext/probs.tsv.gz", package = "viewpoly")
  
  viewqtl_obj <- prepare_qtl_custom_files(selected_mks, qtl_info, blups, 
                                                     beta.hat, profile, effects, probs)
  
  expect_equal(check_viewqtl(viewqtl_obj),0)
  
  check_viewqtl_values(viewqtl_obj, 
                       245183.2, 
                       6.702273, 
                       3.091306e-10, 
                       341.4562, 
                       4.29507e-05, 
                       -7.647564e-11, 
                       1)

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

  # upload QTLpoly
  input.data <- remim.mod <- est.effects <- fitted.mod <- list()
  input.data$datapath <- tempfile()
  remim.mod$datapath <- tempfile()
  est.effects$datapath <- tempfile()
  fitted.mod$datapath <- tempfile()
  
  download.file("https://www.polyploids.org/sites/default/files/2022-04/tetra_QTLpoly_data.RData", destfile = input.data$datapath)
  download.file("https://www.polyploids.org/sites/default/files/2022-04/tetra_QTLpoly_remim.RData", destfile = remim.mod$datapath)
  download.file("https://www.polyploids.org/sites/default/files/2022-04/tetra_QTLpoly_effects.RData", destfile = est.effects$datapath)
  download.file("https://www.polyploids.org/sites/default/files/2022-04/tetra_QTLpoly_fitted.RData", destfile = fitted.mod$datapath)
 
  viewqtl_qtlpoly <- prepare_QTLpoly(data = input.data, 
                                                remim.mod = remim.mod, 
                                                est.effects = est.effects, 
                                                fitted.mod = fitted.mod)
  
  expect_equal(check_viewqtl(viewqtl_qtlpoly),0)
  
  check_viewqtl_values(viewqtl_qtlpoly, 
                       116504, 
                       5.418909, 
                       -1.067457e-11, 
                       299.6155, 
                       0.000160791, 
                       2.340129e-12, 
                       1)
  
})


 
# 
# test_that("upload diaQTL",{
#   prepare_diaQTL()
# })
# 
# test_that("upload polyqtlR",{
#   prepare_polyqtlR()
# })

