
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
    print(mappoly_in)
    viewmap <- read_mappoly_lst(mappoly_prep = mappoly_in)
  } else if(example_map != "none"){
    viewmap <- read_custom_files(system.file("ext/dosage.rds", package = "viewpoly"),
                                 system.file("ext/phases.rds", package = "viewpoly"),
                                 system.file("ext/map.rds", package = "viewpoly"))
  } else if(!(is.null(dosages) | is.null(phases) | is.null(genetic_map))){
    viewmap <- read_custom_files(dosages, phases, genetic_map) 
  } else {
    stop("Please choose one of the option in the previous screen.")
  }
  return(viewmap)
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


#' Change mappoly::: prepare_map output format
#' 
read_mappoly_lst <- function(mappoly_prep){
  cat(mappoly_prep$datapath)
  prep <- readRDS(mappoly_prep$datapath)
  viewmap_obj <- list(dp = lapply(prep, "[[", 5),
                      dq = lapply(prep, "[[", 6),
                      ph.p = lapply(prep, "[[", 3),
                      ph.q = lapply(prep, "[[", 4),
                      maps = lapply(prep, "[[", 2))
}

