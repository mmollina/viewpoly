
#' Receive upload files and convert to viewmap object
#' 
#' @param mappoly_in mappoly sequence object
#' @param dosages RDS compressed file with both parents dosage information.
#' It should contain four columns: 1) character vector with chromosomes ID; 
#' 2) Character vector with markers ID; 3) Character vector with parent ID; 
#' 4) numerical vector with dosage. 
#' @param phases RDS compressed file with phases information. It should contain:
#' 1) Character vector with chromosome ID; 2) Character vector with marker ID;
#' 3 to (ploidy number)*2 columns with each parents haplotypes.
#' @param genetic_map RDS compressed file with the genetic map information
#' @param example_map if user don´t come in with any other input, use a example files.
#'  
prepare_Mapdata <- function(mappoly_in = NULL,
                            dosages = NULL,
                            phases = NULL,
                            genetic_map = NULL,
                            example_map){
  
  if(!is.null(mappoly_in)){
    viewmap <- read_mappoly_lst(mappoly_prep = mappoly_in)
  } else if(example_map == "bt_map"){
    viewmap <- read_custom_files(system.file("ext/dosage.rds", package = "viewpoly"),
                                 system.file("ext/phases.rds", package = "viewpoly"),
                                 system.file("ext/map.rds", package = "viewpoly"))
  } else if(!(is.null(dosages) | is.null(phases) | is.null(genetic_map))){
    viewmap <- read_custom_files(dosages, phases, genetic_map) 
    return(viewmap)
  } else {
    stop("Please choose one of the option in the previous screen.")
  }
}

#' Read input custom format files
#' 
#' @param dosages RDS compressed file with both parents dosage information.
#' It should contain four columns: 1) character vector with chromosomes ID; 
#' 2) Character vector with markers ID; 3) Character vector with parent ID; 
#' 4) numerical vector with dosage. 
#' @param phases RDS compressed file with phases information. It should contain:
#' 1) Character vector with chromosome ID; 2) Character vector with marker ID;
#' 3 to (ploidy number)*2 columns with each parents haplotypes.
#' @param genetic_map RDS compressed file with the genetic map information
#' 
#' @import dplyr
#' 
read_custom_files <- function(dosages, phases, genetic_map){
  ds <- readRDS(dosages)
  ph <- readRDS(phases)
  map <- readRDS(genetic_map)
  
  parent1 <- unique(ds$parent)[1]
  parent2 <- unique(ds$parent)[2]
  d.p1 <- ds %>% filter(parent == parent1) %>% select(chr, marker, dosages) 
  d.p1 <- split(d.p1$dosages, d.p1$chr)
  
  d.p2 <- ds %>% filter(parent == parent2) %>% select(chr, marker, dosages) 
  d.p2 <- split(d.p2$dosages, d.p2$chr)
  
  maps <- split(map$dist, map$chr)
  maps.names <- split(map$marker, map$chr)
  for(i in 1:length(maps)){
    names(maps[[i]]) <- maps.names[[i]]  
  }
  names(maps) <- NULL
  
  ploidy <- (dim(ph)[2] - 2)/2
  
  ph.p1 <- select(ph, 3:(ploidy +2))
  rownames(ph.p1) <- ph$marker
  ph.p1 <- split(ph.p1, ph$chr)
  ph.p1 <- lapply(ph.p1, as.matrix)
  
  ph.p2 <- select(ph, (ploidy +3):dim(ph)[2])
  rownames(ph.p2) <- ph$marker
  ph.p2 <- split(ph.p2, ph$chr)
  ph.p2 <- lapply(ph.p2, as.matrix)
  
  viewmap_obj <- list(dp = d.p1, 
                      dq = d.p2,
                      ph.p = ph.p1,
                      ph.q = ph.p2,
                      maps = maps)
  return(viewmap_obj)
}


#' Change mappoly:::prepare_map output format
#' 
read_mappoly_lst <- function(mappoly_prep){
  prep <- readRDS(mappoly_prep$datapath)
  viewmap_obj <- list(dp = lapply(prep, "[[", 5),
                      dq = lapply(prep, "[[", 6),
                      ph.p = lapply(prep, "[[", 3),
                      ph.q = lapply(prep, "[[", 4),
                      maps = lapply(prep, "[[", 2))
}

#' take all information needed from qtlpoly objects
#' 
#' @param remim.mod object of class "qtlpoly.model" "qtlpoly.remim".
#' 
#' @param data object of class "qtlpoly.data"
#' 
#' @import largeList
#' @import tidyr
#' 
#' @export
prepare_QTLpoly <- function(data, remim.mod, est.effects, fitted.mod){
  
  # Only selected markers
  lgs.t <- lapply(data$lgs, function(x) data.frame(mk = names(x), pos = x))
  lgs <- data.frame()
  for(i in 1:length(lgs.t)) {
    lgs <- rbind(lgs, cbind(LG = i,lgs.t[[i]]))
  }
  
  rownames(lgs) <- NULL
  qtl_info <- u.hat <- beta.hat <- pvalue <- LOD <- data.frame()
  for(i in 1:length(remim.mod$results)){
    pheno = names(fitted.mod$results)[i]
    lower <- remim.mod$results[[i]]$lower[,1:2]
    upper <- remim.mod$results[[i]]$upper[,1:2]
    qtl <- remim.mod$results[[i]]$qtls[,c(1,2,6)]
    int <- cbind(LG = lower$LG, Pos_lower = lower$Pos_lower, 
                 Pos_upper = upper$Pos_upper, qtl[,2:3])
    int <- cbind(pheno = names(remim.mod$results)[i], int)
    
    h2 <- fitted.mod$results[[i]]$qtls[-dim(fitted.mod$results[[i]]$qtls)[1],c(1:3,7)]
    
    int <- merge(int, h2, by = c("LG", "Pos"))
    
    qtl_info <- rbind(qtl_info, int[order(int$LG, int$Pos),])
    
    u.hat.t <- do.call(cbind, fitted.mod$results[[i]]$fitted$u.hat)
    colnames(u.hat.t) <- names(fitted.mod$results[[i]]$fitted$u.hat)
    u.hat.t <- cbind(haplo = rownames(u.hat.t), pheno , as.data.frame(u.hat.t))
    u.hat.t <- pivot_longer(u.hat.t, cols = c(1:length(u.hat.t))[-c(1:2)], values_to = "u.hat", names_to = "qtl")
    u.hat <- rbind(u.hat, u.hat.t)
    
    beta.hat.t <- data.frame(pheno, beta.hat = fitted.mod$results[[i]]$fitted$beta.hat[,1])
    beta.hat <- rbind(beta.hat, beta.hat.t) 
    
    if(is(remim.mod, "qtlpoly.feim")) SIG <- remim.mod$results[[t]][[3]] else SIG <- -log10(as.numeric(remim.mod$results[[t]][[3]]))
    
    LOD.t <- data.frame(pheno, SIG)
    LOD <- rbind(LOD, LOD.t)
  }
  
  # Rearrange the progeny probabilities into a list
  z <- lapply(seq(dim(data$Z)[2]), function(x) data$Z[ , x, ])
  names(z) <- names(data$Z[1,,1])

  path <- tempfile()
  saveList(z, file = file.path(path, "test.llo"), compress = TRUE) # the example data is a subset of 5 individuals

  result <- list(lgs, qtl_info, u.hat, beta.hat, LOD, probs = file.path(path, "test.llo"))
  return(result)
}

prepare_QTLdata <- function(qtlpoly_in = NULL,
                            p_values = NULL,
                            intervals = NULL,
                            qtls = NULL,
                            example_qtl){
  if(!is.null(qtlpoly_in)){
    qtls <- qtlpoly_in
  } else if(example_qtl == "bt_map"){
    temp <- load(system.file("ext", "qtl_in.rda", package = "viewpoly"))
    qtls <- get(temp)
    qtls[[6]] <- system.file("ext/probs.llo", package = "viewpoly")
  } else if(!(is.null(p_values) | is.null(intervals) | is.null(qtls))){
    qtls <- read_custom_files(dosages, phases, genetic_map) # update here
    return(qtls)
  } else {
    stop("Please choose one of the option in the previous screen.")
  }
}

