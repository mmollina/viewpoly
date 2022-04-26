library(testthat)
library(vdiffr)
library(viewpoly)

#' Viewpoly object sanity check 
#' 
#' 
#' @param viewpoly_obj an object of class \code{viewpoly}
#' 
#' @return if consistent, returns 0. If not consistent, returns a 
#'         vector with a number of tests, where \code{TRUE} indicates
#'         a failed test.
#'         
#' @importFrom testthat expect_equal
#'         
#' @author Cristiane Taniguti, \email{chtaniguti@tamu.edu}
#' @keywords internal
check_viewmap_values <- function(viewmap_obj, doses, phases, maps){
  expect_equal(as.vector(table(viewmap_obj$d.p1[[1]])), doses)
  expect_equal(as.vector(table(viewmap_obj$ph.p1[[1]][,3])), phases)
  expect_equal(as.vector(sum(viewmap_obj$maps[[1]][,2])), maps, tolerance = 0.001)
}

#' Check viewqtl object for QTLpoly uploaded files
#' 
#' @param viewqtl_obj an object of class \code{viewqtl}
#' @param pos sum of the position values to be matched
#' @param h2 sum of the log herdability values to be matched
#' @param u.hat sum of the estimated u values to be matched
#' @param beta.hat sum of the estimated beta values to be matched
#' @param lop sum of the LOP values to be matched
#' @param effect sum of the effects values to be matched
#' @param probs sum of the genotype probabilties values to be matched
#' 
#' @return if consistent, returns 0. If not consistent, returns a 
#'         vector with a number of tests, where \code{TRUE} indicates
#'         a failed test.
#'         
#' @importFrom testthat expect_equal
#' 
#' @author Cristiane Taniguti, \email{chtaniguti@tamu.edu}
#' @keywords internal
check_viewqtl_qtlpoly_values <- function(viewqtl_obj, pos, h2, u.hat, beta.hat, lop, effect, probs){
  expect_equal(sum(viewqtl_obj$selected_mks$pos), pos, tolerance = 0.001)
  expect_equal(sum(viewqtl_obj$qtl_info$h2), h2, tolerance = 0.001)
  expect_equal(sum(viewqtl_obj$blups$u.hat), u.hat, tolerance = 0.001)
  expect_equal(sum(viewqtl_obj$beta.hat$beta.hat), beta.hat, tolerance = 0.001)
  expect_equal(min(viewqtl_obj$profile$LOP), lop, tolerance = 0.001)
  expect_equal(sum(viewqtl_obj$effects$effect), effect, tolerance = 0.001)
  expect_equal(sum(viewqtl_obj$probs[,1,2]), probs, tolerance = 0.001)
}

#' Check viewqtl object for diaQTL uploaded files
#' 
#' @param viewqtl_obj an object of class \code{viewqtl}
#' @param pos sum of the position values to be matched
#' @param ll sum of the log likelihood values to be matched
#' @param DIC sum of the deltaDIC values to be matched
#' @param effect sum of the effects values to be matched
#' 
#' @return if consistent, returns 0. If not consistent, returns a 
#'         vector with a number of tests, where \code{TRUE} indicates
#'         a failed test.
#'         
#' @importFrom testthat expect_equal        
#' 
#' @author Cristiane Taniguti, \email{chtaniguti@tamu.edu}
#' @keywords internal
check_viewqtl_diaqtl_values <- function(viewqtl_obj, pos, ll, DIC, effect){
  expect_equal(sum(viewqtl_obj$selected_mks$pos), pos, tolerance = 0.001)
  expect_equal(sum(viewqtl_obj$qtl_info$LL), ll, tolerance = 0.001)
  expect_equal(sum(viewqtl_obj$profile$deltaDIC), DIC, tolerance = 0.001)
  expect_equal(sum(viewqtl_obj$effects$effect, na.rm = T), effect, tolerance = 0.001)
}

#' Check viewqtl object for polyqtlR uploaded files
#' 
#' @param viewqtl_obj an object of class \code{viewqtl}
#' @param pos sum of the position values to be matched
#' @param thre sum of the threshold values to be matched
#' @param lod sum of the LOD values to be matched
#' @param effect sum of the effects values to be matched
#' 
#' @return if consistent, returns 0. If not consistent, returns a 
#'         vector with a number of tests, where \code{TRUE} indicates
#'         a failed test.
#'         
#' @importFrom testthat expect_equal  
#'         
#' @author Cristiane Taniguti, \email{chtaniguti@tamu.edu}
#' @keywords internal
check_viewqtl_polyqtlr_values <- function(viewqtl_obj, pos, thre, lod, effect){
  expect_equal(sum(viewqtl_obj$selected_mks$pos), pos, tolerance = 0.001)
  expect_equal(sum(viewqtl_obj$qtl_info$thresh), thre, tolerance = 0.001)
  expect_equal(sum(viewqtl_obj$profile$LOD), lod, tolerance = 0.001)
  expect_equal(sum(viewqtl_obj$effects[,4], na.rm = T), effect, tolerance = 0.001)
}

test_check("viewpoly")
