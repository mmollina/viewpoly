#' Upload example files
#' 
prepare_examples <- function(example, env.obj= NULL){
  if(example == "hex_map"){
    withProgress(message = 'Working...', value = 0, {
      incProgress(0.5, detail = paste("Uploading BT example map data..."))
      load(system.file("ext/viewmap_hexa.RData", package = "viewpoly"))
    })
    structure(list(map = viewmap_hexa,
                qtl = NULL,
                fasta= system.file("ext/Trifida.Chr01.fa.gz", package="viewpoly"),
                gff3 = system.file("ext/Trifida.Chr01.gff3.gz", package = "viewpoly")),
              class = "viewpoly")
    
  } else if(example == "tetra_map"){
    withProgress(message = 'Working...', value = 0, {
      incProgress(0.5, detail = paste("Uploading tetraploid potato example map data..."))
      load(system.file("ext/viewmap_tetra.RData", package = "viewpoly"))
      load(system.file("ext/viewqtl_tetra.RData", package = "viewpoly"))
    })
    structure(list(map=viewmap_tetra, 
                qtl=viewqtl_tetra, 
                fasta= system.file("ext/Stuberosum.Chr01.fa.gz", package = "viewpoly"),
                gff3 = system.file("ext/Stuberosum.Chr01.gff3.gz", package="viewpoly")),
           class = "viewpoly")
  } else {
    return(env.obj)
  }
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
  
  structure(list(d.p1 = d.p1, 
                 d.p2 = d.p2,
                 ph.p1 = ph.p1,
                 ph.p2 = ph.p2,
                 maps = maps,
                 software = "custom"), 
            class = "viewmap")
}

#' convert list of mappoly.map object into viewpoly_map object
#' 
prepare_MAPpoly <- function(mappoly_list){
  temp <- load(mappoly_list$datapath)
  mappoly_list <- get(temp)
  prep <- lapply(mappoly_list, prepare_map)
  
  structure(list(d.p1 = lapply(prep, "[[", 5),
                 d.p2 = lapply(prep, "[[", 6),
                 ph.p1 = lapply(prep, "[[", 3),
                 ph.p2 = lapply(prep, "[[", 4),
                 maps = lapply(prep, "[[", 2),
                 software = "MAPpoly"), 
            class = "viewmap")
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
  structure(viewmap, class = "viewmap")
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
  withProgress(message = 'Working...', value = 0, {
    incProgress(0.1, detail = paste("Uploading QTLpoly data..."))
    temp <- load(data$datapath)
    data <- get(temp)
    
    temp <- load(remim.mod$datapath)
    remim.mod <- get(temp)
    
    temp <- load(est.effects$datapath)
    est.effects <- get(temp)
    
    temp <- load(fitted.mod$datapath)
    fitted.mod <- get(temp)
    
    # Only selected markers
    incProgress(0.3, detail = paste("Uploading QTLpoly data..."))
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
    incProgress(0.75, detail = paste("Uploading QTLpoly data..."))
    
    # Rearrange the progeny probabilities into a list
    probs <- data$Z
  })
  
  structure(list(selected_mks = lgs, 
                 qtl_info = qtl_info, 
                 blups = as.data.frame(u.hat), 
                 beta.hat = beta.hat, 
                 profile = profile, 
                 effects = effects, 
                 probs = probs, 
                 software = "QTLpoly", class = "viewqtl"))
}

#' Prepare diaQTL output files
#' 
#' 
#' 
prepare_diaQTL <- function(scan1_list, scan1_summaries_list, fitQTL_list, BayesCI_list){
  withProgress(message = 'Working...', value = 0, {
    incProgress(0.1, detail = paste("Uploading diaQTL data..."))
    
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
    incProgress(0.3, detail = paste("Uploading diaQTL data..."))
    
    for(i in 1:length(scan1_summaries_list)){
      temp <- cbind(pheno = names(scan1_summaries_list)[i],scan1_summaries_list[[i]]$peaks)
      qtl_info <- rbind(qtl_info, temp)
    }
    
    qtls.id <- list()
    qtl_info2 <- data.frame()
    incProgress(0.5, detail = paste("Uploading diaQTL data..."))
    
    profile <- effects <- data.frame()
    for(i in 1:length(fitQTL_list)){
      qtls.id <- colnames(fitQTL_list[[i]]$effects$additive)
      trait <- gsub("Trait: ","",fitQTL_list[[i]]$plots[[1]]$additive$labels$title)
      qtl_temp <- qtl_info %>% filter(pheno == trait & marker %in% qtls.id)
      qtl_info2 <- rbind(qtl_info2, qtl_temp)
      # profile
      profile_temp <- data.frame(pheno = trait, deltaDIC = scan1_list[[which(names(scan1_list) == trait)]]$deltaDIC)
      profile <- rbind(profile, profile_temp)
      
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
        temp <- fitQTL_list[[i]]$plots[[j]]$digenic$data
        if(!is.null(temp)){
          effects.di.t <- data.frame(pheno = trait, 
                                     haplo = paste0(temp$x,"x",temp$y),
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
    incProgress(0.75, detail = paste("Uploading diaQTL data..."))
    
    CI <- lapply(BayesCI_list, function(x) {
      y = c(Pos_lower = x$cM[1], Pos_upper = x$cM[length(x$cM)])
      return(y)
    })
    
    CI <- do.call(rbind, CI)
    
    qtl_info <- qtl_info2[,c(3,4,1,6)]
    qtl_info <- cbind(qtl_info, CI)
    qtl_info <- qtl_info[,c(1:3,5,6,4)]
    colnames(qtl_info)[1:2] <- c("LG", "Pos")
    
    incProgress(0.85, detail = paste("Uploading diaQTL data..."))
    
  })
  structure(list(selected_mks = selected_mks,
                 qtl_info = qtl_info,
                 profile = profile,
                 effects = effects,
                 software = "diaQTL"), 
            class = "viewqtl")
}

#' Upload polyqtlR outputs
#' 
prepare_polyqtlR <- function(polyqtlR_QTLscan_list, polyqtlR_IBD, polyqtlR_phenotypes){
  withProgress(message = 'Working...', value = 0, {
    incProgress(0.1, detail = paste("Uploading polyqtlR data..."))
    temp <- load(polyqtlR_QTLscan_list$datapath)
    polyqtlR_QTLscan_list <- get(temp)
    
    temp <- load(polyqtlR_IBD$datapath)
    polyqtlR_IBD <- get(temp)
    
    temp <- load(polyqtlR_phenotypes$datapath)
    polyqtlR_phenotypes <- get(temp)
    
    # selected markers
    selected_mks <- polyqtlR_QTLscan_list[[1]]$Map
    colnames(selected_mks) <- c("LG", "mk", "pos")
    incProgress(0.5, detail = paste("Uploading polyqtlR data..."))
    
    profile <- qtl_info <- effects <- data.frame()
    for(i in 1:length(polyqtlR_QTLscan_list)){
      pheno <- names(polyqtlR_QTLscan_list)[i]
      # profile
      profile_temp <- data.frame(pheno = pheno,
                                 LOD = polyqtlR_QTLscan_list[[i]]$QTL.res$LOD)
      profile <- rbind(profile, profile_temp)
      
      # qtl info
      for(j in 1:length(unique(selected_mks$LG))){
        qtl_info_temp <- find_polyqtlR_Peak(polyqtlR_QTLscan_list[[i]], linkage_group = unique(selected_mks$LG)[j])
        if(!is.null(qtl_info_temp)){
          qtl_info_temp <- data.frame(LG = qtl_info_temp$chromosome, 
                                      Pos = qtl_info_temp$position,
                                      pheno = pheno, 
                                      Pos_lower = NA, # Change here with 1.5 LOD rule ?
                                      Pos_upper = NA,
                                      thresh = qtl_info_temp$thresh)
          qtl_info <- rbind(qtl_info, qtl_info_temp)
        }
        # effects
        effects.t <- polyqtlReffects(IBD_list = polyqtlR_IBD, # input
                                     Phenotype.df = polyqtlR_phenotypes, # input
                                     trait.ID = pheno,
                                     linkage_group = j,
                                     LOD_data = polyqtlR_QTLscan_list[[i]])
        effects <- rbind(effects, effects.t)
      }
    }
    incProgress(0.85, detail = paste("Uploading polyqtlR data..."))
  })
  
  structure(list(selected_mks = selected_mks,
                 qtl_info = qtl_info,
                 profile = profile,
                 effects = effects,
                 software = "polyqtlR"), 
            class = "viewqtl")
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
  qtls$software <- "custom"
  
  structure(qtls, class = "viewqtl")
}

