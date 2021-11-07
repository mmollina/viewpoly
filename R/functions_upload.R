#' Receive upload files and convert to viewmap.map object
#' 
#' @param mappoly_in mappoly prepare_map RData file
#' @param dosages TSV or TSV.GZ file with both parents dosage information.
#' It should contain four columns: 1) character vector with chromosomes ID; 
#' 2) Character vector with markers ID; 3) Character vector with parent ID; 
#' 4) numerical vector with dosage. 
#' @param phases TSV or TSV.GZ file with phases information. It should contain:
#' 1) Character vector with chromosome ID; 2) Character vector with marker ID;
#' 3 to (ploidy number)*2 columns with each parents haplotypes.
#' @param genetic_map TSV or TSV.GZ file with the genetic map information
#' @param example_map if user don´t come in with any other input, use a example files.
#'
#' @author Cristiane Taniguti, \email{chtaniguti@tamu.edu}
#' 
read_Mapdata <- function(mappoly_in = NULL,
                         dosages = NULL,
                         phases = NULL,
                         genetic_map = NULL,
                         example_map){
  
  if(!is.null(mappoly_in)){
    viewmap <- prepare_MAPpoly(mappoly_prep = mappoly_in)
  } else if(example_map == "bt_map"){
    viewmap <- prepare_map_custom_files(system.file("ext/dosage.tsv.gz", package = "viewpoly"),
                                        system.file("ext/phases.tsv.gz", package = "viewpoly"),
                                        system.file("ext/map.tsv.gz", package = "viewpoly"))
  } else if(!(is.null(dosages) | is.null(phases) | is.null(genetic_map))){
    viewmap <- prepare_map_custom_files(dosages, phases, genetic_map) 
    return(viewmap)
  } else {
    stop("Please choose one of the option in the previous screen to upload genetic map information.")
  }
}

read_QTLdata <- function(qtlpoly_data = NULL,
                         qtlpoly_remim.mod = NULL,
                         qtlpoly_est.effects = NULL,
                         qtlpoly_fitted.mod = NULL,
                         selected_mks = NULL,
                         qtl_info = NULL,
                         blups = NULL,
                         beta.hat = NULL,
                         profile = NULL,
                         effects = NULL,
                         probs = NULL,
                         example_qtl = NULL){
  if(!is.null(qtlpoly_data) | !is.null(qtlpoly_remim.mod) | 
     !is.null(qtlpoly_est.effects) | !is.null(qtlpoly_fitted.mod)){
    qtls <- prepare_QTLpoly(qtlpoly_data, qtlpoly_remim.mod, 
                            qtl_est.effects, qtl_fitted.mod)
  } else if(example_qtl == "bt_map"){
    qtls <- prepare_qtl_custom_files(system.file("ext", "selected_mks.tsv.gz", package = "viewpoly"),
                                     system.file("ext", "qtl_info.tsv.gz", package = "viewpoly"),
                                     system.file("ext", "blups.tsv.gz", package = "viewpoly"),
                                     system.file("ext", "beta.hat.tsv.gz", package = "viewpoly"),
                                     system.file("ext", "profile.tsv.gz", package = "viewpoly"),
                                     system.file("ext", "effects.tsv.gz", package = "viewpoly"),
                                     system.file("ext", "probs.tsv.gz", package = "viewpoly"))
  } else if(!(is.null(qtl_info) | is.null(blups) | 
              is.null(beta.hat) | is.null(profile) | 
              is.null(effects) | is.null(probs) | is.null(selected_mks))){
    qtls <- prepare_qtl_custom_files(selected_mks, qtl_info, blups, beta.hat,
                                     profile, effects, probs)
  } else {
    stop("Please choose one of the option in the previous screen to upload QTL information.")
  }
  return(qtls)
}


#' Read input custom format files
#' 
#' @param dosages TSV or TSV.GZ file with both parents dosage information.
#' It should contain four columns: 1) character vector with chromosomes ID; 
#' 2) Character vector with markers ID; 3) Character vector with parent ID; 
#' 4) numerical vector with dosage. 
#' @param phases TSV or TSV.GZ file with phases information. It should contain:
#' 1) Character vector with chromosome ID; 2) Character vector with marker ID;
#' 3 to (ploidy number)*2 columns with each parents haplotypes.
#' @param genetic_map TSV or TSV.GZ file with the genetic map information
#' 
#' @import dplyr
#' @import vroom
#' 
prepare_map_custom_files <- function(dosages, phases, genetic_map){
  ds <- vroom(dosages)
  ph <- vroom(phases)
  map <- vroom(genetic_map)
  
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
prepare_MAPpoly <- function(mappoly_prep){
  prep <- load(mappoly_prep$datapath) # Check
  viewmap_obj <- list(dp = lapply(prep, "[[", 5),
                      dq = lapply(prep, "[[", 6),
                      ph.p = lapply(prep, "[[", 3),
                      ph.q = lapply(prep, "[[", 4),
                      maps = lapply(prep, "[[", 2))
}

#' take all information needed from qtlpoly objects
#' 
#' @param data object of class "qtlpoly.data"
#' @param remim.mod object of class "qtlpoly.model" "qtlpoly.remim".
#' @param est.effects object of class "qtlpoly.effects"
#' @param fitted.mod object of class "qtlpoly.fitted"
#' 
#' @author Cristiane Taniguti, \email{chtaniguti@tamu.edu}
#' 
#' @import tidyr
#' 
prepare_QTLpoly <- function(data, remim.mod, est.effects, fitted.mod){
  
  # Only selected markers
  lgs.t <- lapply(data$lgs, function(x) data.frame(mk = names(x), pos = x))
  lgs <- data.frame()
  for(i in 1:length(lgs.t)) {
    lgs <- rbind(lgs, cbind(LG = i,lgs.t[[i]]))
  }
  
  rownames(lgs) <- NULL
  qtl_info <- u.hat <- beta.hat <- pvalue <- profile <- effects <- data.frame()
  for(i in 1:length(remim.mod$results)){
    pheno = names(fitted.mod$results)[i]
    lower <- remim.mod$results[[i]]$lower[,1:2]
    upper <- remim.mod$results[[i]]$upper[,1:2]
    qtl <- remim.mod$results[[i]]$qtls[,c(1,2,6)]
    int <- cbind(LG = lower$LG, Pos_lower = lower$Pos_lower, 
                 Pos_upper = upper$Pos_upper, qtl[,2:3])
    int <- cbind(pheno = names(remim.mod$results)[i], int)
    
    if(dim(fitted.mod$results[[i]]$qtls)[1] > 1) {
      h2 <- fitted.mod$results[[i]]$qtls[-dim(fitted.mod$results[[i]]$qtls)[1],c(1:2,7)]
      h2 <- data.frame(apply(h2, 2, unlist))
    }else {
      h2 <- fitted.mod$results[[i]]$qtls[,c(1:2,7)]
    }
    
    int <- merge(int, h2, by = c("LG", "Pos"))
    
    qtl_info <- rbind(qtl_info, int[order(int$LG, int$Pos),])
    
    u.hat.t <- do.call(cbind, fitted.mod$results[[i]]$fitted$u.hat)
    colnames(u.hat.t) <- names(fitted.mod$results[[i]]$fitted$u.hat)
    u.hat.t <- cbind(haplo = rownames(u.hat.t), pheno , as.data.frame(u.hat.t))
    u.hat.t <- pivot_longer(u.hat.t, cols = c(1:length(u.hat.t))[-c(1:2)], values_to = "u.hat", names_to = "qtl")
    u.hat <- rbind(u.hat, u.hat.t)
    u.hat$qtl <- gsub("g", "", u.hat$qtl)
    
    beta.hat.t <- data.frame(pheno, beta.hat = fitted.mod$results[[i]]$fitted$beta.hat[,1])
    beta.hat <- rbind(beta.hat, beta.hat.t) 
    
    if(is(remim.mod, "qtlpoly.feim")) SIG <- remim.mod$results[[i]][[3]] else SIG <- -log10(as.numeric(remim.mod$results[[i]][[3]]))
    
    profile.t <- data.frame(pheno, LOP = SIG)
    profile <- rbind(profile, profile.t)
    
    for(j in 1:length(est.effects$results[[i]]$effects)){
      effects.t <- do.call(rbind, lapply(est.effects$results[[i]]$effects[[j]], function(x) data.frame(haplo = names(x), effect = x)))
      effects.t <- cbind(pheno = pheno, qtl.id= j, effects.t)
      effects <- rbind(effects, effects.t)
    }
  }
  
  # Rearrange the progeny probabilities into a list
  probs <- data$Z
  
  result <- list(selected_mks = lgs, qtl_info = qtl_info, 
                 blups = as.data.frame(u.hat), beta.hat = beta.hat, 
                 profile = profile, effects = effects, probs = probs, software = "QTLpoly")
  return(result)
}

#' 
#' @import vroom
#' @import abind
#' 
prepare_qtl_custom_files <- function(selected_mks, qtl_info, blups, beta.hat,
                                     profile, effects, probs){
  qtls <- list()
  qtls$selected_mks <- as.data.frame(vroom(selected_mks))
  qtls$qtl_info <- as.data.frame(vroom(qtl_info))
  qtls$blups <- as.data.frame(vroom(blups))
  qtls$beta.hat <- as.data.frame(vroom(beta.hat))
  qtls$profile <- as.data.frame(vroom(profile))
  qtls$profile[,2] <- as.numeric(qtls$profile[,2])
  qtls$effects <- as.data.frame(vroom(effects))
  
  probs.t <- vroom(probs)
  ind <- probs.t$ind
  probs.t <- as.data.frame(probs.t[,-1])
  probs.df <- split(probs.t, ind)
  qtls$probs <- abind(probs.df, along = 3)
  
  return(qtls)
}


