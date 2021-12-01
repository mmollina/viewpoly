#' Import data from polymapR
#'
#' Function to import datasets from polymapR. Function from MAPpoly.
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
#' @rdname utils
#' 
#' @keywords internal
#' 
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

#' Filter non-conforming classes in F1, non double reduced population.
#' Function from MAPpoly.
#'
#' @param void internal function to be documented
#' 
#' @rdname utils
#' 
#' @keywords internal
filter_non_conforming_classes <- function(input.data, prob.thres = NULL)
{
  if (!inherits(input.data, "mappoly.data")) {
    stop(deparse(substitute(input.data)), " is not an object of class 'mappoly.data'")
  }
  ploidy <- input.data$ploidy
  dp <- input.data$dosage.p1
  dq <- input.data$dosage.p2
  Ds <- array(NA, dim = c(ploidy+1, ploidy+1, ploidy+1))
  for(i in 0:ploidy)
    for(j in 0:ploidy)
      Ds[i+1,j+1,] <- segreg_poly(ploidy = ploidy, d.p1 = i, d.p2 = j)
  Dpop <- cbind(dp,dq)
  #Gathering segregation probabilities given parental dosages
  M <- t(apply(Dpop, 1, function(x) Ds[x[1]+1, x[2]+1,]))
  M[M != 0] <- 1
  dimnames(M) <- list(input.data$mrk.names, 0:ploidy)
  ##if no prior probabilities
  if(!is.prob.data(input.data)){
    for(i in 1:nrow(M)){
      id0 <- !as.numeric(input.data$geno.dose[i,])%in%(which(M[i,] == 1)-1)
      if(any(id0))
        input.data$geno.dose[i,id0] <- (ploidy+1)     
    }
    return(input.data)
  }
  ## 1 represents conforming classes/ 0 represents non-conforming classes
  dp <- rep(dp, input.data$n.ind)
  dq <- rep(dq, input.data$n.ind)
  M <- M[rep(seq_len(nrow(M)), input.data$n.ind),]
  R <- input.data$geno[,-c(1:2)] - input.data$geno[,-c(1:2)]*M
  id1 <- apply(R, 1, function(x) abs(sum(x))) > 0.3 # if the sum of the excluded classes is greater than 0.3, use segreg_poly
  N <- matrix(NA, sum(id1), input.data$ploidy+1)
  ct <- 1
  for(i in which(id1)){
    N[ct,] <- Ds[dp[i]+1, dq[i]+1, ]
    ct <- ct+1
  }
  input.data$geno[id1,-c(1:2)] <- N
  # if the sum of the excluded classes is greater than zero
  # and smaller than 0.3, assign zero to those classes and normalize the vector
  input.data$geno[,-c(1:2)][R > 0] <- 0
  input.data$geno[,-c(1:2)] <- sweep(input.data$geno[,-c(1:2)], 1, rowSums(input.data$geno[,-c(1:2)]), FUN = "/")
  if(is.null(prob.thres))
    prob.thres <- input.data$prob.thres
  geno.dose <- dist_prob_to_class(geno = input.data$geno, prob.thres = prob.thres)
  if(geno.dose$flag)
  {
    input.data$geno <- geno.dose$geno
    input.data$geno.dose <- geno.dose$geno.dose
  } else {
    input.data$geno.dose <- geno.dose$geno.dose
  }
  input.data$geno.dose[is.na(input.data$geno.dose)] <- ploidy + 1
  input.data$n.ind <- ncol(input.data$geno.dose)
  input.data$ind.names <- colnames(input.data$geno.dose)
  return(input.data)
}

#' Linkage phase format conversion: matrix to list. Function from MAPpoly.
#' 
#' This function converts linkage phase configurations from matrix
#' form to list
#'
#' @param M matrix whose columns represent homologous chromosomes and
#'     the rows represent markers
#' 
#' @return a list of linkage phase configurations
#' 
#' @rdname utils
#' 
#' @keywords internal
ph_matrix_to_list <- function(M) {
  w <- lapply(split(M, seq(NROW(M))), function(x, M) which(x  ==  1))
  w[sapply(w, function(x) length(x)  ==  0)] <- 0
  w
}

#' Is it a probability dataset? Function from MAPpoly.
#'
#' @param void internal function to be documented
#' 
#' @rdname utils
#' 
#' @keywords internal
is.prob.data <- function(x){
  exists('geno', where = x)
}

#' Map functions. From MAPpoly
#'
#' @param void internal function to be documented
#' 
#' @rdname utils
#' 
#' @keywords internal
mf_h <- function(d) 0.5 * (1 - exp(-d/50))

#' Chi-square test. Function from MAPpoly.
#'
#' @param void internal function to be documented
#' 
#' @rdname utils
#' 
#' @keywords internal
mrk_chisq_test <- function(x, ploidy){
  y <- x[-c(1:(ploidy+1))]
  y[y == ploidy+1] <- NA
  y <- table(y, useNA = "always")
  names(y) <- c(names(y)[-length(y)], "NA") 
  seg.exp <- x[0:(ploidy+1)]
  seg.exp <- seg.exp[seg.exp != 0]
  seg.obs <- seg.exp
  seg.obs[names(y)[-length(y)]] <- y[-length(y)]
  pval <- suppressWarnings(stats::chisq.test(x = seg.obs, p = seg.exp[names(seg.obs)])$p.value)
  pval
}

#' Returns the class with the highest probability in 
#' a genotype probability distribution. Function from MAPpoly.
#'
#' @param geno the probabilistic genotypes contained in the object
#'     \code{'mappoly.data'}
#' @param prob.thres probability threshold to select the genotype. 
#'     Values below this genotype are assumed as missing data
#' @return a matrix containing the doses of each genotype and
#'     marker. Markers are disposed in rows and individuals are 
#'     disposed in columns. Missing data are represented by NAs
#'     
#' @importFrom tidyr "%>%"
#' @importFrom reshape2 melt dcast
#' @importFrom dplyr group_by filter arrange
#' 
#' @rdname utils
#' 
#' @keywords internal
dist_prob_to_class <- function(geno, prob.thres = 0.9) {
  a <- reshape2::melt(geno, id.vars = c("mrk", "ind"))
  mrk <- ind <- value <- variable <- NULL # Setting the variables to NULL first
  a$variable <- as.numeric(levels(a$variable))[a$variable]
  b <- a %>%
    dplyr::group_by(mrk, ind) %>%
    dplyr::filter(value > prob.thres) %>%
    dplyr::arrange(mrk, ind, variable)
  z <- reshape2::dcast(data = b[,1:3], formula = mrk ~ ind, value.var = "variable")
  rownames(z) <- z[,"mrk"]
  z <- data.matrix(frame = z[,-1])
  n <- setdiff(unique(geno$mrk), rownames(z))
  if(length(n) > 0)
  {
    ploidy <- matrix(NA, nrow = length(n), ncol = ncol(z), dimnames = list(n, colnames(z)))
    z <- rbind(z,ploidy)
  }
  rm.ind <- setdiff(unique(geno$ind), colnames(z))
  flag <- FALSE
  if(length(rm.ind) > 0){
    flag <- TRUE
    warning("Inividual(s) ", paste(rm.ind, collapse = " "), 
            "\n  did not meet the 'prob.thres' criteria for any of\n  the markers and was (were) removed.")
    geno <- geno %>% dplyr::filter(ind %in% colnames(z))
  }
  z <- z[as.character(unique(geno$mrk)), as.character(unique(geno$ind))]
  list(geno.dose = z, geno = geno, flag = flag)
}

#' Polysomic segregation frequency - Function from MAPpoly
#'
#' Computes the polysomic segregation frequency given a ploidy level
#' and the dosage of the locus in both parents. It does not consider
#' double reduction.
#'
#' @param ploidy the ploidy level
#'
#' @param d.p1 the dosage in parent P
#'
#' @param d.p2 the dosage in parent Q
#'
#' @return a vector containing the expected segregation frequency for
#'     all possible genotypic classes.
#'
#'
#' @author Marcelo Mollinari, \email{mmollin@ncsu.edu}
#'
#' @references
#'     Mollinari, M., and Garcia, A.  A. F. (2019) Linkage
#'     analysis and haplotype phasing in experimental autopolyploid
#'     populations with high ploidy level using hidden Markov
#'     models, _G3: Genes, Genomes, Genetics_. 
#'     \doi{10.1534/g3.119.400378}
#'     
#'     Serang O, Mollinari M, Garcia AAF (2012) Efficient Exact 
#'     Maximum a Posteriori Computation for Bayesian SNP 
#'     Genotyping in Polyploids. _PLoS ONE_ 7(2): 
#'     e30906.
#'     
#'
#' @importFrom stats dhyper
#' 
#' @rdname utils
#' 
#' @keywords internal
segreg_poly <- function(ploidy, d.p1, d.p2) {
  if (ploidy%%2 != 0)
    stop("m must be an even number")
  p.dose <- numeric((ploidy + 1))
  p.names <- character((ploidy + 1))
  seg.p1 <- dhyper(x = c(0:(ploidy + 1)), m = d.p1, n = (ploidy - d.p1), k = ploidy/2)
  seg.p2 <- dhyper(x = c(0:(ploidy + 1)), m = d.p2, n = (ploidy - d.p2), k = ploidy/2)
  M <- tcrossprod(seg.p1, seg.p2)
  for (i in 1:nrow(M)) {
    for (j in 1:ncol(M)) {
      p.dose[i + j - 1] <- p.dose[i + j - 1] + M[i, j]
    }
  }
  p.dose <- p.dose[!is.na(p.dose)]
  for (i in 0:ploidy) p.names[i + 1] <- paste(paste(rep("A", i), collapse = ""), paste(rep("a", (ploidy - i)), collapse = ""), sep = "")
  names(p.dose) <- p.names
  return(p.dose)
}

#' Import phased map list from polymapR
#'
#' Function to import phased map lists from polymapR. Function from MAPpoly.
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
#' @rdname utils
#' 
#' @keywords internal
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
    #MAPs[[i]] <- loglike_hmm(MAPs[[i]], mappoly.data)
  }
  MAPs
}

#' prepare maps for plot - from MAPpoly
#' @param void internal function to be documented
#' 
#' @rdname utils
#' 
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
                    g.dist = if(!is.null(input.map$info$seq.alt)) input.map$info$seq.alt else NA ,
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
#' 
#' @rdname utils
#' 
#' @keywords internal
imf_h <- function(r) {
  r[r >= 0.5] <- 0.5 - 1e-14
  -50 * log(1 - 2 * r)
}

#' Extract the LOD Scores in a \code{'mappoly.map'} object
#' Function from MAPpoly.
#' @param x an object of class \code{mappoly.map}
#' @param sorted logical. if \code{TRUE}, the LOD Scores are displayed
#'     in a decreasing order
#' @return a numeric vector containing the LOD Scores
#' @keywords internal
#' 
#' @rdname utils
#' 
#' @keywords internal
get_LOD <- function(x, sorted = TRUE) {
  w <- sapply(x$maps, function(x) x$loglike)
  LOD <- w - max(w)
  if (sorted)
    LOD <- sort(LOD, decreasing = TRUE)
  abs(LOD)
}

#' Linkage phase format conversion: list to matrix.
#' Function from MAPpoly.
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
#' @rdname utils
#' 
#' @keywords internal
ph_list_to_matrix <- function(L, ploidy) {
  M <- matrix(0, nrow = length(L), ncol = ploidy)
  for (i in 1:nrow(M)) if (all(L[[i]] != 0))
    M[i, L[[i]]] <- 1
  M
}

#' Addapted from polyqtlR
#' 
#' @rdname utils
#' 
#' @keywords internal
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
#' @rdname utils
#' 
#' @keywords internal
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
#' @rdname utils
#' 
#' @keywords internal
list.depth <- function (obj, objdepth = 0) {
  if (!is.list(obj)) {
    return(objdepth)
  } else {
    return(max(unlist(lapply(obj, list.depth, objdepth = objdepth + 1))))
  }
}

