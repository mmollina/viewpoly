#' Upload example files
#' 
#' @param example character indicating the example dataset selected
#' 
#' @return object of class \code{viewpoly}
#' 
#' 
#' @importFrom utils download.file
#' 
#' @keywords internal
prepare_examples <- function(example){
  viewmap_tetra <- viewqtl_tetra <- NULL
  if(example == "BExMG"){
    temp <- load(system.file("ext/BExMG_viewpoly_chr_updated.RData", package = "viewpoly"))
    obj <- get(temp)
    
  } else if (example == "SWxBE"){
    temp <- load(system.file("ext/SWxBE_viewpoly_chr_updated.RData", package = "viewpoly"))
    obj <- get(temp)
  }
  
  obj$fasta <- c(system.file("ext/OBDH_1.0_formated.fasta.gz",package="viewpoly"))
  
  obj$gff3 <- c(system.file("ext/OBDH_1.0_gene_models_sorted.gff3.gz", package="viewpoly"))
  
  return(obj)
}

#' Converts map information in custom format files to viewmap object
#' 
#' 
#' @param dosages TSV or TSV.GZ file with both parents dosage information.
#' It should contain four columns: 1) character vector with chromosomes ID; 
#' 2) Character vector with markers ID; 3) Character vector with parent ID; 
#' 4) numerical vector with dosage. 
#' @param phases TSV or TSV.GZ file with phases information. It should contain:
#' 1) Character vector with chromosome ID; 2) Character vector with marker ID;
#' 3 to (ploidy number)*2 columns with each parents haplotypes.
#' @param genetic_map TSV or TSV.GZ file with the genetic map information
#' @param mks_pos TSV or TSV.GZ file with table with three columns: 1) marker ID; 
#' 2) genome position; 3) chromosome
#' 
#' @return object of class \code{viewmap}
#' 
#' @import dplyr
#' @import vroom
#' 
#' @keywords internal
prepare_map_custom_files <- function(dosages, phases, genetic_map, mks_pos=NULL){
  parent <- chr <- marker <- NULL
  ds <- vroom(dosages$datapath, progress = FALSE, col_types = cols())
  ph <- vroom(phases$datapath, progress = FALSE, col_types = cols())
  map <- vroom(genetic_map$datapath, progress = FALSE, col_types = cols())
  if(!is.null(mks_pos)) mks_pos <- vroom(mks_pos$datapath, progress = FALSE, col_types = cols()) 
  
  parent1 <- unique(ds$parent)[1]
  parent2 <- unique(ds$parent)[2]
  d.p1 <- ds %>% filter(parent == parent1) %>% select(chr, marker, dosages) 
  d.p1.names <- split(d.p1$marker, d.p1$chr)
  d.p1 <- split(d.p1$dosages, d.p1$chr)
  d.p1 <- Map(function(x,y) {
    names(x) <- y
    return(x)
  }, d.p1, d.p1.names)
  
  d.p2 <- ds %>% filter(parent == parent2) %>% select(chr, marker, dosages) 
  d.p2.names <- split(d.p2$marker, d.p2$chr)
  d.p2 <- split(d.p2$dosages, d.p2$chr)
  d.p2 <- Map(function(x,y) {
    names(x) <- y
    return(x)
  }, d.p2, d.p2.names)
  
  if(!is.null(mks_pos)) pos <- mks_pos[,2][match(map$marker,mks_pos[,1])] else pos <- NA
  
  maps <- data.frame(mk.names = map$marker, 
                     l.dist = map$dist, 
                     g.chr = map$chr, 
                     g.dist = pos,
                     alt = NA,
                     ref= NA)
  
  maps <- split.data.frame(maps, maps$g.chr)
  
  ploidy <- (dim(ph)[2] - 2)/2
  
  ph.p1 <- as.data.frame(select(ph, 3:(ploidy +2)))
  rownames(ph.p1) <- ph$marker
  ph.p1 <- split(ph.p1, ph$chr)
  ph.p1 <- lapply(ph.p1, as.matrix)
  
  ph.p2 <- as.data.frame(select(ph, (ploidy +3):dim(ph)[2]))
  rownames(ph.p2) <- ph$marker
  ph.p2 <- split(ph.p2, ph$chr)
  ph.p2 <- lapply(ph.p2, as.matrix)
  
  structure(list(d.p1 = d.p1, 
                 d.p2 = d.p2,
                 ph.p1 = ph.p1,
                 ph.p2 = ph.p2,
                 maps = maps,
                 software = "custom"), 
            class = "viewmap")
}

#' Converts list of mappoly.map object into viewmap object
#' 
#' @param mappoly_list list with objects of class \code{mappoly.map}
#' 
#' @return object of class \code{viewmap}
#' 
#' 
#' @keywords internal
prepare_MAPpoly <- function(mappoly_list){
  is <- NULL
  
  if(!is(mappoly_list[[1]], "mappoly.map")){
    temp <- load(mappoly_list$datapath)
    mappoly_list <- get(temp)
  }
  prep <- lapply(mappoly_list, prepare_map)
  
  structure(list(d.p1 = lapply(prep, "[[", 5),
                 d.p2 = lapply(prep, "[[", 6),
                 ph.p1 = lapply(prep, "[[", 3),
                 ph.p2 = lapply(prep, "[[", 4),
                 maps = lapply(prep, "[[", 2),
                 software = "MAPpoly"), 
            class = "viewmap")
}

#' Converts polymapR ouputs to viewmap object
#' 
#' @param polymapR.dataset a \code{polymapR} dataset
#' @param polymapR.map output map sequence from polymapR
#' @param input.type indicates whether the input is discrete ("disc") or probabilistic ("prob")
#' @param ploidy ploidy level
#' 
#' @return object of class \code{viewmap}
#' 
#' 
#' @keywords internal
prepare_polymapR <- function(polymapR.dataset, polymapR.map, input.type, ploidy){ 
  
  temp <- load(polymapR.dataset$datapath)
  polymapR.dataset <- get(temp)
  
  temp <- load(polymapR.map$datapath)
  polymapR.map <- get(temp)
  data <- import_data_from_polymapR(input.data = polymapR.dataset, 
                                    ploidy = ploidy, 
                                    parent1 = "P1", 
                                    parent2 = "P2",
                                    input.type = ,
                                    prob.thres = 0.95,
                                    pardose = NULL, 
                                    offspring = NULL,
                                    filter.non.conforming = TRUE,
                                    verbose = FALSE)
  
  map_seq <- import_phased_maplist_from_polymapR(maplist = polymapR.map, 
                                                 mappoly.data = data)
  
  viewmap <- prepare_MAPpoly(mappoly_list = map_seq)
  viewmap$software <- "polymapR"
  
  structure(viewmap, class = "viewmap")
}

#' Converts QTLpoly outputs to viewqtl object
#' 
#' 
#' @param data object of class "qtlpoly.data"
#' @param remim.mod object of class "qtlpoly.model" "qtlpoly.remim".
#' @param est.effects object of class "qtlpoly.effects"
#' @param fitted.mod object of class "qtlpoly.fitted"
#' 
#' @author Cristiane Taniguti, \email{chtaniguti@tamu.edu}
#' 
#' @return object of class \code{viewqtl}
#' 
#' @importFrom tidyr pivot_longer
#' @import dplyr
#' 
#' @keywords internal
prepare_QTLpoly <- function(data, remim.mod, est.effects, fitted.mod){
  is <- NULL
  
  temp <- load(data$datapath)
  data <- get(temp)
  
  temp <- load(remim.mod$datapath)
  remim.mod <- get(temp)
  
  temp <- load(est.effects$datapath)
  est.effects <- get(temp)
  
  temp <- load(fitted.mod$datapath)
  fitted.mod <- get(temp)
  
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
  
  structure(list(selected_mks = lgs, 
                 qtl_info = qtl_info, 
                 blups = as.data.frame(u.hat), 
                 beta.hat = beta.hat, 
                 profile = profile, 
                 effects = effects, 
                 probs = probs, 
                 software = "QTLpoly"), 
            class = "viewqtl")
}

#' Converts diaQTL output to viewqtl object
#' 
#' @param scan1_list list with results from diaQTL \code{scan1} function
#' @param scan1_summaries_list list with results from diaQTL \code{scan1_summaries} function 
#' @param fitQTL_list list with results from diaQTL \code{fitQTL} function
#' @param BayesCI_list list with results from diaQTL \code{BayesCI} function
#' 
#' @return object of class \code{viewqtl}
#' 
#' @importFrom dplyr filter
#' 
#' 
#' @keywords internal
prepare_diaQTL <- function(scan1_list, scan1_summaries_list, fitQTL_list, BayesCI_list){
  marker <- pheno <- NULL
  
  temp <- load(scan1_list$datapath)
  scan1_list <- get(temp)
  
  temp <- load(scan1_summaries_list$datapath)
  scan1_summaries_list <- get(temp)
  
  temp <- load(fitQTL_list$datapath)
  fitQTL_list <- get(temp)
  
  temp <- load(BayesCI_list$datapath)
  BayesCI_list <- get(temp)
  
  selected_mks <- scan1_list[[1]][,c(2,1,3)]
  colnames(selected_mks) <- c("LG", "mk", "pos")
  qtl_info <- data.frame()
  
  for(i in 1:length(scan1_summaries_list)){
    temp <- cbind(pheno = names(scan1_summaries_list)[i],scan1_summaries_list[[i]]$peaks)
    qtl_info <- rbind(qtl_info, temp)
  }
  
  qtls.id <- list()
  qtl_info2 <- data.frame()
  
  profile <- effects <- data.frame()
  for(i in 1:length(fitQTL_list)){
    qtls.id <- colnames(fitQTL_list[[i]]$effects$additive)
    trait <- gsub("Trait: ","",fitQTL_list[[i]]$plots[[1]]$additive$labels$title)
    qtl_temp <-  filter(qtl_info, pheno == trait & marker %in% qtls.id)
    qtl_info2 <- rbind(qtl_info2, qtl_temp)
    # profile
    profile_temp <- data.frame(pheno = trait, deltaDIC = scan1_list[[which(names(scan1_list) == trait)]]$deltaDIC)
    profile <- rbind(profile, profile_temp)
    
    # Sometimes there is a graphic about epistasis that is not described anywhere yet. We ignored it here by now.
    if(any(grepl("epistasis", names(fitQTL_list[[i]]$plots)))) fitQTL_list[[i]]$plots <- fitQTL_list[[i]]$plots[-grep("epistasis", names(fitQTL_list[[i]]$plots))]
    
    for(j in 1:length(fitQTL_list[[i]]$plots)){
      # aditive effect
      temp <- fitQTL_list[[i]]$plots[[j]]$additive$data
      effects.ad.t <- data.frame(pheno = trait, 
                                 haplo = rownames(temp), 
                                 qtl.id = j, 
                                 effect= temp$mean, 
                                 type = "Additive",
                                 CI.lower = temp$CI.lower,
                                 CI.upper = temp$CI.upper)
      
      # digenic effect
      temp <- data.frame(haplo = rownames(fitQTL_list[[i]]$effects$digenic), z = fitQTL_list[[i]]$effects$digenic[,j])
      if(!is.null(temp)){
        effects.di.t <- data.frame(pheno = trait, 
                                   haplo = gsub("[+]", "x", temp$haplo),
                                   qtl.id = j,
                                   effect = as.numeric(temp$z), 
                                   type = "Digenic",
                                   CI.lower = NA,
                                   CI.upper = NA)
        
        effects.t <- rbind(effects.ad.t, effects.di.t)
      } else effects.t <- effects.ad.t
      effects.t <- effects.t[order(effects.t$pheno, effects.t$qtl.id, effects.t$type,effects.t$haplo),]
      effects <- rbind(effects, effects.t)
    }
  }
  
  # Ordering Bayes info according to qtl info
  BayesCI_list_ord <- list()
  for(i in 1:dim(qtl_info2)[1]){
    for(j in 1:length(BayesCI_list)){
      if(any(paste0(BayesCI_list[[j]]$pheno, BayesCI_list[[j]]$marker, BayesCI_list[[j]]$chrom) %in% 
             paste0(qtl_info2$pheno[i], qtl_info2$marker[i], qtl_info2$chrom[i]))){
        BayesCI_list_ord[[i]] <- BayesCI_list[[j]]
      }
    }
  }
  
  if(length(BayesCI_list_ord) != dim(qtl_info2)[1]) BayesCI_list_ord[[length(BayesCI_list_ord) + 1]] <- NULL
  idx <- which(sapply(BayesCI_list_ord, is.null))
  if(length(idx) != 0 | length(BayesCI_list_ord) != dim(qtl_info2)[1]){
    warning(paste0("Bayes confidence interval information (from diaQTL function BayesCI) was not provided for the QTL in chromosome:", qtl_info2[idx, 3], 
                   "; phenotype: ", qtl_info2[idx, 1]))
  }
  
  CI <- lapply(BayesCI_list_ord, function(x) {
    y = c(Pos_lower = x$cM[1], Pos_upper = x$cM[length(x$cM)])
    return(y)
  })
  
  
  idx <- which(sapply(CI, is.null))
  if(length(idx) != 0 | length(BayesCI_list_ord) != dim(qtl_info2)[1]){
    if(length(idx) != 0)
      CI[[idx]] <- c(NA, NA)
    else  CI[[length(CI) + 1]] <- c(NA, NA)
    
  }
  
  CI <- do.call(rbind, CI)
  
  qtl_info <- qtl_info2[,c(3,4,1,6)]
  qtl_info <- cbind(qtl_info, CI)
  qtl_info <- qtl_info[,c(1:3,5,6,4)]
  colnames(qtl_info)[1:2] <- c("LG", "Pos")
  
  structure(list(selected_mks = selected_mks,
                 qtl_info = qtl_info,
                 profile = profile,
                 effects = effects,
                 software = "diaQTL"), 
            class = "viewqtl")
}

#' Converts polyqtlR outputs to viewqtl object
#' 
#' @param polyqtlR_QTLscan_list list containing results from polyqtlR \code{QTLscan_list} function
#' @param polyqtlR_qtl_info data.frame containing the QTL information:LG - group ID; Pos - QTL position (cM); 
#' pheno - phenotype ID; Pos_lower - lower position of confidence interval; Pos_upper - upper position of the confidence interval;
#' thresh - LOD threshold applied
#' @param polyqtlR_effects data.frame with results from \code{visualiseQTLeffects} polyqtlR function
#' 
#' @return object of class \code{viewqtl}
#' 
#' 
#' @keywords internal
prepare_polyqtlR <- function(polyqtlR_QTLscan_list, polyqtlR_qtl_info, polyqtlR_effects){
  
  temp <- load(polyqtlR_QTLscan_list$datapath)
  polyqtlR_QTLscan_list <- get(temp)
  
  temp <- load(polyqtlR_qtl_info$datapath)
  polyqtlR_qtl_info <- get(temp)
  
  temp <- load(polyqtlR_effects$datapath)
  polyqtlR_effects <- get(temp)
  
  # selected markers
  selected_mks <- polyqtlR_QTLscan_list[[1]]$Map
  colnames(selected_mks) <- c("LG", "mk", "pos")
  
  profile <- qtl_info <- effects <- data.frame()
  for(i in 1:length(polyqtlR_QTLscan_list)){
    pheno <- names(polyqtlR_QTLscan_list)[i]
    # profile
    profile_temp <- data.frame(pheno = pheno,
                               LOD = polyqtlR_QTLscan_list[[i]]$QTL.res$LOD)
    profile <- rbind(profile, profile_temp)
  }
  
  structure(list(selected_mks = selected_mks,
                 qtl_info = polyqtlR_qtl_info,
                 profile = profile,
                 effects = polyqtlR_effects,
                 software = "polyqtlR"), 
            class = "viewqtl")
}

#' Converts QTL information in custom files to viewqtl object
#' 
#' @param selected_mks data.frame with: LG - linkage group ID; mk - marker ID; pos - position in linkage map (cM)
#' @param qtl_info data.frame with: LG - linkage group ID; Pos - position in linkage map (cM); 
#' Pheno - phenotype ID; Pos_lower - lower position of confidence interval; 
#' Pos_upper - upper position of the confidence interval; Pval - QTL p-value; h2 - herdability
#' @param blups data.frame with: haplo - haplotype ID; pheno - phenotype ID; qtl - QTL ID; u.hat - QTL estimated BLUPs
#' @param beta.hat data.frame with: pheno - phenotype ID; beta.hat - estimated beta
#' @param profile data.frame with: pheno - phenotype ID; LOP - significance value for the QTL, in this case LOP (can be LOD or DIC depending of the software used)
#' @param effects data.frame with: pheno - phenotype ID; qtl.id - QTL ID; haplo - haplotype ID; effect - haplotype effect value
#' @param probs data.frame with first column (named `ind`) as individuals ID and next columns named with markers ID and containing the genotype probability at each marker
#' 
#' 
#' @return object of class \code{viewqtl}
#' 
#' @import vroom
#' @import abind
#' 
#' @keywords internal
prepare_qtl_custom_files <- function(selected_mks, qtl_info, blups, beta.hat,
                                     profile, effects, probs){
  
  qtls <- list()
  qtls$selected_mks <- as.data.frame(vroom(selected_mks$datapath, progress = FALSE, col_types = cols()))
  qtls$qtl_info <- as.data.frame(vroom(qtl_info$datapath, progress = FALSE, col_types = cols()))
  qtls$blups <- as.data.frame(vroom(blups$datapath, progress = FALSE, col_types = cols()))
  qtls$beta.hat <- as.data.frame(vroom(beta.hat$datapath, progress = FALSE, col_types = cols()))
  qtls$profile <- as.data.frame(vroom(profile$datapath, progress = FALSE, col_types = cols()))
  qtls$profile[,2] <- as.numeric(qtls$profile[,2])
  qtls$effects <- as.data.frame(vroom(effects$datapath, progress = FALSE, col_types = cols()))
  
  probs.t <- vroom(probs$datapath, progress = FALSE, col_types = cols())
  ind <- probs.t$ind
  probs.t <- as.data.frame(probs.t[,-1])
  probs.df <- split(probs.t, ind)
  qtls$probs <- abind(probs.df, along = 3)
  qtls$software <- "custom"
  
  structure(qtls, class = "viewqtl")
}

