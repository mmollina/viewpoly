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
                         polymapR.dataset = NULL,
                         polymapR.map = NULL,
                         dosages = NULL,
                         phases = NULL,
                         genetic_map = NULL,
                         example_map){
  withProgress(message = 'Uploding map data', value = 0, {
    if(!is.null(mappoly_in)){
      incProgress(0.5, detail = paste("Uploading MAPpoly data..."))
      viewmap <- prepare_MAPpoly(mappoly_list = mappoly_in)
    } else if(!is.null(polymapR.dataset) & !is.null(polymapR.map)){ ## Require update
      incProgress(0.5, detail = paste("Uploading polymapR data..."))
      viewmap <- prepare_polymapR(polymapR.dataset, polymapR.map)
    } else if(example_map == "hex_map"){
      incProgress(0.5, detail = paste("Uploading BT example map data..."))
      viewmap <- prepare_map_custom_files(system.file("ext/dosage.tsv.gz", package = "viewpoly"), ## change here for the tetraploid. Parei aqui!!
                                          system.file("ext/phases.tsv.gz", package = "viewpoly"),
                                          system.file("ext/map.tsv.gz", package = "viewpoly")) # change this
    } else if(example_map == "tetra_map"){
      incProgress(0.5, detail = paste("Uploading tetraploid potato example map data..."))
      viewmap <- get(data("viewmap_tetra"))
    } else if(!(is.null(dosages) | is.null(phases) | is.null(genetic_map))){
      incProgress(0.5, detail = paste("Uploading custom format data..."))
      viewmap <- prepare_map_custom_files(dosages, phases, genetic_map) 
      return(viewmap)
    } else {
      stop("Please choose one of the option in the previous screen to upload genetic map information.")
    }
  })
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
  withProgress(message = 'Uploading QTL data', value = 0, {
    if(!is.null(qtlpoly_data) | !is.null(qtlpoly_remim.mod) | 
       !is.null(qtlpoly_est.effects) | !is.null(qtlpoly_fitted.mod)){
      incProgress(0.5, detail = paste("Uploading QTLpoly data..."))
      qtls <- prepare_QTLpoly(qtlpoly_data, qtlpoly_remim.mod, 
                              qtl_est.effects, qtl_fitted.mod)
    } else if(example_qtl == "hex_map"){
      incProgress(0.5, detail = paste("Uploading BT example QTL data..."))
      qtls <- get(data("qtl_bt"))
    } else if(example_qtl == "tetra_map"){
      qtls <- get(data("viewqtl_tetra"))
    } else if(!(is.null(qtl_info) | is.null(blups) | 
                is.null(beta.hat) | is.null(profile) | 
                is.null(effects) | is.null(probs) | is.null(selected_mks))){
      incProgress(0.5, detail = paste("Uploading custom format QTL data..."))
      qtls <- prepare_qtl_custom_files(selected_mks, qtl_info, blups, beta.hat,
                                       profile, effects, probs)
    } else {
      warning("Please choose one of the option in the previous screen to upload QTL information.")
      qtls <- NULL
    }
  })
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
  
  viewmap_obj <- list(d.p1 = d.p1, 
                      d.p2 = d.p2,
                      ph.p1 = ph.p1,
                      ph.p2 = ph.p2,
                      maps = maps)
  return(viewmap_obj)
}

#' convert list of mappoly.map object into viewpoly_map object
#' 
prepare_MAPpoly <- function(mappoly_list){
  
  prep <- lapply(mappoly_list, prepare_map)
  
  viewmap_obj <- list(d.p1 = lapply(prep, "[[", 5),
                      d.p2 = lapply(prep, "[[", 6),
                      ph.p1 = lapply(prep, "[[", 3),
                      ph.p2 = lapply(prep, "[[", 4),
                      maps = lapply(prep, "[[", 2))
}

prepare_polymapR <- function(polymapR.dataset, polymapR.map){ ## Require update
  data <- import_data_from_polymapR(polymapR.dataset, 
                                    ploidy, 
                                    parent1 = "P1", 
                                    parent2 = "P2",
                                    input.type = c("discrete", "probabilistic"),
                                    prob.thres = 0.95,
                                    pardose = NULL, 
                                    offspring = NULL,
                                    filter.non.conforming = TRUE,
                                    verbose = FALSE)
  map_seq <- import_phased_maplist_from_polymapR(data, polymapR.map)
  viewmap <- prepare_MAPpoly(map_seq)
  return(viewmap)
}

# data = tetra_QTLpoly_data
# remim.mod = tetra_QTLpoly_remim
# est.effects = tetra_QTLpoly_effects
# fitted.mod = tetra_QTLpoly_fitted

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
    if(!is.null(dim(fitted.mod$results[[i]]$qtls)[1])){
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
      
      u.hat.t <- do.call(cbind, fitted.mod$results[[i]]$fitted$U)
      colnames(u.hat.t) <- names(fitted.mod$results[[i]]$fitted$U)
      u.hat.t <- cbind(haplo = fitted.mod$results[[i]]$fitted$alleles, pheno , as.data.frame(u.hat.t))
      u.hat.t <- pivot_longer(u.hat.t, cols = c(1:length(u.hat.t))[-c(1:2)], values_to = "u.hat", names_to = "qtl")
      u.hat <- rbind(u.hat, u.hat.t)
      u.hat$qtl <- gsub("g", "", u.hat$qtl)
      
      beta.hat.t <- data.frame(pheno, beta.hat = fitted.mod$results[[i]]$fitted$Beta[,1])
      beta.hat <- rbind(beta.hat, beta.hat.t) 
      
      for(j in 1:length(est.effects$results[[i]]$effects)){
        effects.t <- do.call(rbind, lapply(est.effects$results[[i]]$effects[[j]], function(x) data.frame(haplo = names(x), effect = x)))
        effects.t <- cbind(pheno = pheno, qtl.id= j, effects.t)
        effects <- rbind(effects, effects.t)
      }
    }
    
    if(is(remim.mod, "qtlpoly.feim")) SIG <- remim.mod$results[[i]][[3]] else SIG <- -log10(as.numeric(remim.mod$results[[i]][[3]]))
    
    profile.t <- data.frame(pheno, LOP = SIG)
    profile <- rbind(profile, profile.t)
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

#' Import data from polymapR
#'
#' Function to import datasets from polymapR. 
#' 
#' See examples at \url{https://rpubs.com/mmollin/tetra_mappoly_vignette}.
#'
#' @param input.data  a \code{polymapR} dataset
#' @param ploidy the ploidy level     
#' @param parent1 a character string containing the name (or pattern of genotype IDs) of parent 1
#' @param parent2 a character string containing the name (or pattern of genotype IDs) of parent 2
#' @param input.type Indicates whether the input is discrete ("disc") or probabilistic ("prob") 
#' @param prob.thres threshold probability to assign a dosage to offspring. If the probability 
#'        is smaller than \code{thresh.parent.geno}, the data point is converted to 'NA'.
#' @param pardose matrix of dimensions (n.mrk x 3) containing the name of the markers in the first column, and the 
#'        dosage of parents 1 and 2 in columns 2 and 3. (see polymapR vignette)      
#' @param offspring a character string containing the name (or pattern of genotype IDs) of the offspring 
#'                  individuals. If \code{NULL} (default) it considers all individuals as offsprings, except 
#'                  \code{parent1} and \code{parent2}.  
#' @param filter.non.conforming if \code{TRUE} exclude samples with non 
#'     expected genotypes under no double reduction. Since markers were already filtered in polymapR, the default is 
#'     \code{FALSE}.
#' @param verbose if \code{TRUE} (default), the current progress is shown; if
#'     \code{FALSE}, no output is produced
#'
#' @author Marcelo Mollinari \email{mmollin@ncsu.edu}
#'
#' @references
#'     Bourke PM et al: (2019) PolymapR — linkage analysis and genetic map 
#'     construction from F1 populations of outcrossing polyploids. 
#'     _Bioinformatics_ 34:3496–3502.
#'     \doi{10.1093/bioinformatics/bty1002}
#' 
#'     Mollinari, M., and Garcia, A.  A. F. (2019) Linkage
#'     analysis and haplotype phasing in experimental autopolyploid
#'     populations with high ploidy level using hidden Markov
#'     models, _G3: Genes, Genomes, Genetics_. 
#'     \doi{10.1534/g3.119.400378}
#'     
#' @export import_data_from_polymapR
#' @importFrom reshape2 acast
#' @importFrom dplyr filter arrange
import_data_from_polymapR <- function(input.data, 
                                      ploidy, 
                                      parent1 = "P1", 
                                      parent2 = "P2",
                                      input.type = c("discrete", "probabilistic"),
                                      prob.thres = 0.95,
                                      pardose = NULL, 
                                      offspring = NULL,
                                      filter.non.conforming = TRUE,
                                      verbose = TRUE){
  input.type <- match.arg(input.type)
  if(input.type  ==  "discrete"){
    geno.dose <- input.data[,-match(c(parent1, parent2), colnames(input.data)), drop = FALSE]
    mappoly.data <- structure(list(ploidy = ploidy,
                                   n.ind = ncol(geno.dose),
                                   n.mrk = nrow(geno.dose),
                                   ind.names = colnames(geno.dose),
                                   mrk.names = rownames(geno.dose),
                                   dosage.p1 = input.data[,parent1],
                                   dosage.p2 = input.data[,parent2],
                                   chrom = NA,
                                   genome.pos = NA,
                                   seq.ref = NULL,
                                   seq.alt = NULL,
                                   all.mrk.depth = NULL,
                                   prob.thres = NULL,
                                   geno.dose = geno.dose,
                                   nphen = 0,
                                   phen = NULL,
                                   kept = NULL,
                                   elim.correspondence = NULL),
                              class = "mappoly.data")
  } 
  else {
    if(is.null(pardose)) 
      stop("provide parental dosage.")
    rownames(pardose) <- pardose$MarkerName
    dat <- input.data[,c("MarkerName", "SampleName",paste0("P", 0:ploidy))]
    p1 <- unique(sapply(parent1, function(x) unique(grep(pattern = x, dat[,"SampleName"], value = TRUE))))
    p2 <- unique(sapply(parent2, function(x) unique(grep(pattern = x, dat[,"SampleName"], value = TRUE))))
    if(is.null(offspring)){
      offspring <- setdiff(as.character(unique(dat[,"SampleName"])), c(p1, p2))    
    } else {
      offspring <- unique(grep(pattern = offspring, dat[,"SampleName"], value = TRUE))
    }
    d1 <- input.data[,c("MarkerName", "SampleName", "geno")]
    geno.dose <- reshape2::acast(d1, MarkerName ~ SampleName, value.var = "geno")
    ## get marker names ----------------------
    mrk.names <- rownames(geno.dose)
    ## get number of individuals -------------
    n.ind <- length(offspring)
    ## get number of markers -----------------
    n.mrk <- length(mrk.names)
    ## get individual names ------------------
    ind.names <- offspring
    ## get dosage in parent P ----------------
    dosage.p1 <- as.integer(pardose[mrk.names,"parent1"])
    names(dosage.p1) <- mrk.names
    ## get dosage in parent Q ----------------
    dosage.p2 <- as.integer(pardose[mrk.names,"parent2"])
    names(dosage.p2) <- mrk.names
    ## monomorphic markers
    d.p1 <- abs(abs(dosage.p1-(ploidy/2))-(ploidy/2))
    d.p2 <- abs(abs(dosage.p2-(ploidy/2))-(ploidy/2))
    mrk.names <- names(which(d.p1+d.p2 != 0))
    dosage.p1 <- dosage.p1[mrk.names]
    dosage.p2 <- dosage.p2[mrk.names]
    nphen <- 0
    phen <- NULL
    if (verbose){
      cat("Importing the following data:")
      cat("\n    Ploidy level:", ploidy)
      cat("\n    No. individuals: ", n.ind)
      cat("\n    No. markers: ", n.mrk) 
      cat("\n    No. informative markers:  ", length(mrk.names), " (", round(100*length(mrk.names)/n.mrk,1), "%)", sep = "")
      cat("\n    ...")
    }
    ## get genotypic info --------------------
    MarkerName <- SampleName <- NULL
    geno <- dat %>%
      dplyr::filter(SampleName %in% offspring)  %>%
      dplyr::filter(MarkerName %in% mrk.names) %>%
      dplyr::arrange(SampleName, MarkerName)
    
    colnames(geno) <- c("mrk", "ind", as.character(0:ploidy))
    ind.names <- unique(geno$ind)
    mrk.names <- unique(geno$mrk)
    dosage.p1 <- dosage.p1[mrk.names]
    dosage.p2 <- dosage.p2[mrk.names]
    
    ## transforming na's in expected genotypes using Mendelian segregation
    i.na <- which(apply(geno, 1, function(x) any(is.na(x))))
    if (length(i.na) > 0) {
      m.na <- match(geno[i.na, 1], mrk.names)
      d.p1.na <- dosage.p1[m.na]
      d.p2.na <- dosage.p2[m.na]
      for (i in 1:length(m.na)) geno[i.na[i], -c(1, 2)] <- segreg_poly(ploidy, d.p1.na[i], d.p2.na[i])
    }
    ## dosage info
    if(filter.non.conforming){
      geno.dose <- geno.dose[mrk.names,offspring]  
    } else {
      geno.dose <- dist_prob_to_class(geno = geno, prob.thres = prob.thres)
      if(geno.dose$flag)
      {
        geno <- geno.dose$geno
        geno.dose <- geno.dose$geno.dose
        n.ind <- ncol(geno.dose)
        ind.names <- colnames(geno.dose)
      } else {
        geno.dose <- geno.dose$geno.dose
      }
      geno.dose[is.na(geno.dose)] <- ploidy + 1
    }
    ## returning the 'mappoly.data' object
    if (verbose) cat("\n    Done with reading.\n")
    mappoly.data <- structure(list(ploidy = ploidy,
                                   n.ind = n.ind,
                                   n.mrk = length(mrk.names),
                                   ind.names = ind.names,
                                   mrk.names = mrk.names,
                                   dosage.p1 = dosage.p1,
                                   dosage.p2 = dosage.p2,
                                   chrom = rep(NA, length(mrk.names)),
                                   genome.pos = rep(NA, length(mrk.names)),
                                   seq.ref = NULL,
                                   seq.alt = NULL,
                                   all.mrk.depth = NULL,
                                   prob.thres = prob.thres,
                                   geno = geno,
                                   geno.dose = geno.dose,
                                   nphen = nphen,
                                   phen = phen,
                                   chisq.pval = NULL,
                                   kept = NULL,
                                   elim.correspondence = NULL),
                              class = "mappoly.data")
  }
  if(filter.non.conforming){
    mappoly.data <- filter_non_conforming_classes(mappoly.data)
    Ds <- array(NA, dim = c(ploidy+1, ploidy+1, ploidy+1))
    for(i in 0:ploidy)
      for(j in 0:ploidy)
        Ds[i+1,j+1,] <- segreg_poly(ploidy = ploidy, d.p1 = i, d.p2 = j)
    d.p1op <- cbind(mappoly.data$dosage.p1, mappoly.data$dosage.p2)
    M <- t(apply(d.p1op, 1, function(x) Ds[x[1]+1, x[2]+1,]))
    dimnames(M) <- list(mappoly.data$mrk.names, c(0:ploidy))
    M <- cbind(M, mappoly.data$geno.dose)
    mappoly.data$chisq.pval <- apply(M, 1, mrk_chisq_test, ploidy = ploidy)
  }
  mappoly.data
}

#' Import phased map list from polymapR
#'
#' Function to import phased map lists from polymapR
#' 
#' See examples at \url{https://rpubs.com/mmollin/tetra_mappoly_vignette}.
#' 
#' @param maplist a list of phased maps obtained using function 
#' \code{create_phased_maplist} from package \code{polymapR} 
#' @param mappoly.data a dataset used to obtain \code{maplist}, 
#' converted into class \code{mappoly.data}
#' @param ploidy the ploidy level     
#'
#' @author Marcelo Mollinari \email{mmollin@ncsu.edu}
#'
#' @references
#'     Bourke PM et al: (2019) PolymapR — linkage analysis and genetic map 
#'     construction from F1 populations of outcrossing polyploids. 
#'     _Bioinformatics_ 34:3496–3502.
#'     \doi{10.1093/bioinformatics/bty1002}
#' 
#'     Mollinari, M., and Garcia, A.  A. F. (2019) Linkage
#'     analysis and haplotype phasing in experimental autopolyploid
#'     populations with high ploidy level using hidden Markov
#'     models, _G3: Genes, Genomes, Genetics_. 
#'     \doi{10.1534/g3.119.400378}
#'     
#' @export import_phased_maplist_from_polymapR
import_phased_maplist_from_polymapR <- function(maplist, 
                                                mappoly.data, 
                                                ploidy = NULL){
  input_classes <- c("list")
  if (!inherits(maplist, input_classes)) {
    stop(deparse(substitute(maplist)), " is not a list of phased maps.")
  }
  X <- maplist[[1]]
  if(is.null(ploidy))
    ploidy <- (ncol(X)-2)/2
  MAPs <- vector("list", length(maplist))
  for(i in 1:length(MAPs)){
    X <- maplist[[i]]
    seq.num <- match(X$marker, mappoly.data$mrk.names)
    seq.rf <- mf_h(diff(X$position)) ## Using haldane
    seq.rf[seq.rf <= 1e-05] <- 1e-4
    P = ph_matrix_to_list(X[,3:(ploidy+2)])
    Q = ph_matrix_to_list(X[,3:(ploidy+2) + ploidy])
    names(P) <- names(Q) <- seq.num
    seq.ph <- list(P = P, Q = Q)
    maps <- vector("list", 1)
    maps[[1]] <- list(seq.num = seq.num, seq.rf = seq.rf, seq.ph = seq.ph, loglike = 0)
    MAPs[[i]] <- structure(list(info = list(ploidy = (ncol(X)-2)/2,
                                            n.mrk = nrow(X),
                                            seq.num = seq.num,
                                            mrk.names = as.character(X$marker),
                                            seq.dose.p1 = mappoly.data$dosage.p1[as.character(X$marker)],
                                            seq.dose.p2 = mappoly.data$dosage.p2[as.character(X$marker)],
                                            chrom = rep(i, nrow(X)),
                                            genome.pos = NULL,
                                            seq.ref = NULL,
                                            seq.alt = NULL,
                                            chisq.pval = mappoly.data$chisq.pval[as.character(X$marker)],
                                            data.name = as.character(sys.call())[3], 
                                            ph.thresh = NULL),
                                maps = maps),
                           class = "mappoly.map")
    MAPs[[i]] <- loglike_hmm(MAPs[[i]], mappoly.data)
  }
  MAPs
}

#' prepare maps for plot - from MAPpoly
#' @param void internal function to be documented
#' @keywords internal
prepare_map <- function(input.map, config = "best"){
  if (!inherits(input.map, "mappoly.map")) {
    stop(deparse(substitute(input.map)), " is not an object of class 'mappoly.map'")
  }
  ## Choosing the linkage phase configuration
  LOD.conf <- get_LOD(input.map, sorted = FALSE)
  if(config  ==  "best") {
    i.lpc <- which.min(LOD.conf)
  } else if(config  ==  "all"){
    i.lpc <- seq_along(LOD.conf) } else if (config > length(LOD.conf)) {
      stop("invalid linkage phase configuration")
    } else i.lpc <- config
  ## Gathering marker positions
  map <- data.frame(mk.names = input.map$info$mrk.names,
                    l.dist = cumsum(imf_h(c(0, input.map$maps[[i.lpc]]$seq.rf))),
                    g.chr = input.map$info$chrom,
                    g.dist = input.map$info$genome.pos,
                    alt = if(!is.null(input.map$info$seq.alt)) input.map$info$seq.alt else NA , # get this info from VCF if it is inputted
                    ref = if(!is.null(input.map$info$seq.ref)) input.map$info$seq.ref else NA)
  ## 
  ph.p1 <- ph_list_to_matrix(input.map$maps[[i.lpc]]$seq.ph$P, input.map$info$ploidy)
  ph.p2 <- ph_list_to_matrix(input.map$maps[[i.lpc]]$seq.ph$Q, input.map$info$ploidy)
  dimnames(ph.p1) <- list(map$mk.names, letters[1:input.map$info$ploidy])
  dimnames(ph.p2) <- list(map$mk.names, letters[(1+input.map$info$ploidy):(2*input.map$info$ploidy)])
  if(is.null(input.map$info$seq.alt))
  {
    ph.p1[ph.p1 == 1] <- ph.p2[ph.p2 == 1] <- "A"
    ph.p1[ph.p1 == 0] <- ph.p2[ph.p2 == 0] <- "B"  
  } else {
    for(i in input.map$info$mrk.names){
      ph.p1[i, ph.p1[i,] == 1] <- input.map$info$seq.alt[i]
      ph.p1[i, ph.p1[i,] == 0] <- input.map$info$seq.ref[i]
      ph.p2[i, ph.p2[i,] == 1] <- input.map$info$seq.alt[i]
      ph.p2[i, ph.p2[i,] == 0] <- input.map$info$seq.ref[i]
    }
  }
  d.p1 <- input.map$info$seq.dose.p1
  d.p2 <- input.map$info$seq.dose.p2
  list(ploidy = input.map$info$ploidy, map = map, ph.p1 = ph.p1, ph.p2 = ph.p2, d.p1 = d.p1, d.p2 = d.p2)
}


#' Map functions - from MAPpoly
#'
#' @param void internal function to be documented
#' @keywords internal
#' @export
imf_h <- function(r) {
  r[r >= 0.5] <- 0.5 - 1e-14
  -50 * log(1 - 2 * r)
}

#' Extract the LOD Scores in a \code{'mappoly.map'} object
#' @param x an object of class \code{mappoly.map}
#' @param sorted logical. if \code{TRUE}, the LOD Scores are displayed
#'     in a decreasing order
#' @return a numeric vector containing the LOD Scores
#' @keywords internal
#' @export
get_LOD <- function(x, sorted = TRUE) {
  w <- sapply(x$maps, function(x) x$loglike)
  LOD <- w - max(w)
  if (sorted)
    LOD <- sort(LOD, decreasing = TRUE)
  abs(LOD)
}

#' Linkage phase format conversion: list to matrix
#'
#' This function converts linkage phase configurations from list
#' to matrix form
#'
#' @param L a list of configuration phases
#'
#' @param ploidy ploidy level
#'
#' @return a matrix whose columns represent homologous chromosomes and
#'     the rows represent markers
#'
#' @keywords internal
#' @export
ph_list_to_matrix <- function(L, ploidy) {
  M <- matrix(0, nrow = length(L), ncol = ploidy)
  for (i in 1:nrow(M)) if (all(L[[i]] != 0))
    M[i, L[[i]]] <- 1
  M
}