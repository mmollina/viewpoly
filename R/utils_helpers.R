

#' Addapted from polyqtlR
#' 
find_polyqtlR_Peak <- function(LOD_data, linkage_group) {
  if (is.null(LOD_data$Perm.res)) {
    warning("No significance threshold available to check")
    thresh <- 0
  }
  else {
    thresh <- LOD_data$Perm.res$threshold
  }
  lgdata <- LOD_data$QTL.res[LOD_data$QTL.res$chromosome == 
                               linkage_group, ]
  maxdata <- lgdata[which.max(lgdata$LOD), ]
  if (maxdata$LOD >= thresh) 
    df <- cbind(maxdata, thresh = thresh)
  else df <- NULL
  return(df)
}

#' Addapted from polyqtlR
#' 
#' 
polyqtlReffects <- function(IBD_list, 
                            Phenotype.df, 
                            genotype.ID, 
                            trait.ID.c, 
                            linkage_group, 
                            LOD_data, 
                            cM_range = NULL) {
  geno <- IBD_list[[1]]$offspring
  idx1 <- sum(Phenotype.df[,1] %in% geno)
  idx2 <- sum(Phenotype.df[,2] %in% geno)
  genotype.ID <- which.max(c(idx1, idx2))
  trait.ID <- which(colnames(Phenotype.df) %in% trait.ID.c)
  
  ploidy <- IBD_list[[1]]$ploidy
  ploidy2 <- IBD_list[[1]]$ploidy2
  if (length(linkage_group) != 1 | !is.numeric(linkage_group)) 
    stop("linkage_group should be a single numeric identifier.")
  listdepth <- list.depth(IBD_list)
  if (listdepth != 3) 
    stop("Unexpected input - IBD_list is expected to be a nested list representing 1 or more chromosomes! Please check input.")
  IBDarray <- IBD_list[[linkage_group]]$IBDarray
  IBDtype <- IBD_list[[linkage_group]]$IBDtype
  bivalent_decoding <- IBD_list[[linkage_group]]$biv_dec
  gap <- IBD_list[[linkage_group]]$gap
  if (IBDtype == "genotypeIBD") {
    if (IBD_list[[1]]$method == "hmm_TO") {
      GenotypeCodes <- t(sapply(IBD_list[[1]]$genocodes, function(x) as.numeric(strsplit(as.character(x), "")[[1]])))
      Nstates <- nrow(GenotypeCodes)
      indicatorMatrix <- matrix(0, nrow = Nstates, ncol = ploidy + ploidy2)
      for (r in 1:nrow(indicatorMatrix)) {
        for (h in 1:ncol(GenotypeCodes)) {
          indicatorMatrix[r, GenotypeCodes[r, h]] <- indicatorMatrix[r, GenotypeCodes[r, h]] + 1
        }
      }
    } else {
      GenotypeCodes <- t(sapply(sapply(IBD_list[[1]]$genocodes, function(x) strsplit(x, "")), function(y) sapply(y, function(z) which(letters == z))))
      Nstates <- nrow(GenotypeCodes)
      indicatorMatrix <- matrix(0, nrow = Nstates, ncol = ploidy + ploidy2)
      for (r in 1:nrow(indicatorMatrix)) {
        for (h in 1:ncol(GenotypeCodes)) {
          indicatorMatrix[r, GenotypeCodes[r, h]] <- indicatorMatrix[r, GenotypeCodes[r, h]] + 1
        }
      }
    }
  } else if (IBDtype == "haplotypeIBD") {
    Nstates <- ploidy + ploidy2
    indicatorMatrix <- NULL
  } else {
    stop("Unable to determine IBD type, please check IBD_list elements have a valid $IBDtype tag.")
  }
  if (dim(IBDarray)[2]%%Nstates != 0) 
    stop("Incompatible input detected..")
  if (!is.null(gap)) {
    cMpositions <- as.numeric(unlist(lapply(strsplit(dimnames(IBDarray)[[1]], "cM"), function(x) x[2])))
  }  else {
    cMpositions <- IBD_list[[linkage_group]]$map$position
  }
  subcM <- 1:length(cMpositions)
  cM_range <- cMpositions
  index_range <- dimnames(IBDarray)[[1]]
  if (any(is.na(Phenotype.df[, genotype.ID]))) {
    warning("Missing genotype.ID values detected. Removing and proceeding without...")
    Phenotype.df <- Phenotype.df[!is.na(Phenotype.df[, genotype.ID]), ]
  }
  phenoGeno0 <- intersect(dimnames(IBDarray)[[3]], unique(Phenotype.df[, genotype.ID]))
  pheno <- Phenotype.df[Phenotype.df[, genotype.ID] %in% phenoGeno0 & !is.na(Phenotype.df[, trait.ID]), 
                        c(genotype.ID, trait.ID)]  
  pheno <- pheno[order(pheno[, 1]), ]
  phenoGeno <- as.character(unique(pheno[, 1]))
  popSize <- length(phenoGeno)
  Phenotypes <- pheno[, 2]
  if (popSize != nrow(pheno)) 
    stop("Check input - unequal dimensions of phenotypes and genotype data.\nSuggest to use polyqtlR:::BLUE() to generate consensus phenotypes first.")
  prob_matched <- match(phenoGeno, dimnames(IBDarray)[[3]])
  
  calc.AlEffects <- function(cM.pos = index_range[1], tag) {
    if (tag == "genotypeIBD") {
      haploProbs <- t(IBDarray[cM.pos, , prob_matched]) %*% indicatorMatrix
    } else {
      haploProbs <- t(IBDarray[cM.pos, , prob_matched])
    }
    mean_sd_weights <- function(wghts, min.wghts, zerosum = FALSE) {
      plus.mean <- sum(wghts * Phenotypes)/sum(wghts)
      minus.mean <- sum(min.wghts * Phenotypes)/sum(min.wghts)
      if (zerosum) {
        return(plus.mean - minus.mean)
      } else {
        return(plus.mean - mean(Phenotypes))
      }
    }
    output <- setNames(sapply(1:(ploidy + ploidy2), function(sub.config) {
      weights <- haploProbs[, sub.config]
      min.weights <- 1 - haploProbs[, sub.config]
      mean_sd_weights(weights, min.weights)
    }), paste0("H", 1:(ploidy + ploidy2)))
  }
  AlEffects <- do.call(rbind, lapply(index_range, calc.AlEffects, IBDtype))
  AlEffects <- data.frame(pos = cM_range, pheno = trait.ID.c, LG = linkage_group, AlEffects)
  
  return(AlEffects)
}

#' From polyqtlR
#' 
list.depth <- function (obj, objdepth = 0) {
  if (!is.list(obj)) {
    return(objdepth)
  } else {
    return(max(unlist(lapply(obj, list.depth, objdepth = objdepth + 1))))
  }
}

