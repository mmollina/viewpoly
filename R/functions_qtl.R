#' Logarithm of \emph{P}-value (LOP) profile plots. Modified version of QTLpoly function.
#'
#' Plots profiled logarithm of score-based \emph{P}-values (LOP) from individual or combined traits.
#'
#' @import ggplot2
#' @import dplyr
#' @import tidyr
#' @importFrom plotly TeX
#' 
plot_profile <- function(profile, qtl_info, selected_mks, pheno.col = NULL, 
                         lgs.id = NULL, by_range = TRUE, range.min = NULL, range.max = NULL, plot=TRUE) {
  
  lgs.size <- selected_mks %>% group_by(LG) %>% group_map(~ tail(.x, 1)) %>% do.call(rbind, .) 
  lgs.size <- lgs.size$pos
  lines <- points <- thre <- map <- data.frame()
  y.dat <- trait.names <- c()
  count <- 0
  
  nphe <- length(pheno.col)
  LGS <- selected_mks$LG
  POS <- selected_mks$pos
  for(p in 1:nphe) { #lines
    TRT <- rep(unique(profile$pheno)[pheno.col[p]], length(LGS))
    SIG <- profile[which(profile$pheno == TRT),2]
    lines <- rbind(lines, data.frame(TRT=as.factor(TRT), LGS=LGS, POS=POS, SIG=SIG))
  }
  
  count <- 0
  y.dat <- c()
  for(p in 1:nphe) { #points
    trait.names <- unique(qtl_info$pheno)[pheno.col[p]]
    if(!is.null(qtl_info)) {
      qtl_info.sub <- qtl_info %>% filter(pheno == trait.names)
      nqtls <-  qtl_info.sub %>% summarize(n()) 
      TRT <- qtl_info.sub$pheno
      LGS <- qtl_info.sub$LG
      POS <- qtl_info.sub$Pos
      INF <- qtl_info.sub$Pos_lower
      SUP <- qtl_info.sub$Pos_upper
      PVAL <- qtl_info.sub$Pval
      H2 <- qtl_info.sub$h2
      points <- rbind(points, data.frame(TRT=TRT, LGS=LGS, POS=POS, INF=INF, SUP=SUP, PVAL = PVAL, H2 = round(H2,2)))
      count <- count+1
      y.dat <- c(y.dat, rep((-0.4*count), nqtls))
    }
  }
  points <- cbind(points, y.dat)
  
  # The axis name change according with software
  y.lab <- colnames(profile)[2]
  if(y.lab == "LOP")  {
    if(by_range){
      y.lab <- "LOP"
    }  else {
      y.lab <- expression(-log[10](italic(P)))
      
    }
  }
  
  # Filter group
  if(!is.null(lgs.id)){
    lines <- lines[which(lines$LGS %in% lgs.id),]
    points <- points[which(points$LGS %in% lgs.id),]
  }
  
  # Interval
  lines$INT <- NA
  for(i in 1:dim(points)[1]){
    idx <- which(lines$POS >= points$INF[i] & 
                   lines$POS <= points$SUP[i] & 
                   lines$LGS == points$LGS[i] &
                   lines$TRT == points$TRT[i])
    lines$INT[idx] <- lines$POS[idx]
  }
  
  # Filter position
  lines$range <- NA
  if(!is.null(range.min)){
    lines$range[which(lines$POS >= range.min & lines$POS <= range.max)] <- lines$SIG[which(lines$POS >= range.min & lines$POS <= range.max)]
    lines$SIG[which(lines$POS > range.min & lines$POS < range.max)] <- NA
  }
  
  dot.height <- data.frame(trt = unique(points$TRT), heigth = unique(points$y.dat))
  y.dat.lines <- dot.height$heigth[match(lines$TRT, dot.height$trt)]
  lines$y.dat <- y.dat.lines
  
  colnames(lines) <- c("Trait", "LG", "Position (cM)", "SIG", "INT", "range","y.dat")
  colnames(points)[1:3] <- c("Trait", "LG", "Position (cM)")
  
  if(max(lgs.size[lgs.id]) > 200) cutx <- 150 else cutx <- 100
  if(length(lgs.size[lgs.id]) > 10) {linesize <- 1} else {cutx <- 50; linesize <- 1.25}
  
  if(plot){
    if(by_range){
      pl <- ggplot(data = lines, aes(x = `Position (cM)`, color = Trait)) +
        facet_grid(.~LG, space = "free") +
        {if(!all(is.na(lines$INT))) geom_path(data=lines, aes(x = INT, y = y.dat), colour = "black")} +
        geom_line(data=lines, aes(y = range, color = Trait), size=linesize, alpha=0.8, lineend = "round", show.legend = F) +
        geom_line(data=lines, aes(y = SIG, shape = Trait),  colour = "gray", size=linesize, alpha=0.8, lineend = "round") +
        scale_x_continuous(breaks=seq(0,max(lgs.size),cutx)) +
        {if(!all(is.na(lines$INT))) geom_point(data=points, aes(y = y.dat, color = Trait), shape = 2, size = 2, stroke = 1, alpha = 0.8)} +
        scale_y_continuous(breaks=seq(0,max(lgs.size, na.rm = T))) +
        guides(color = guide_legend("Trait"), fill = guide_legend("Trait"), shape = guide_legend("Trait")) + 
        labs(y = y.lab, x = "Position (cM)", subtitle="Linkage group") + 
        theme_classic()
    } else {
      pl <- ggplot(data = lines, aes(x = `Position (cM)`, color = Trait, group=1)) +
        facet_grid(.~LG, space = "free") +
        {if(!all(is.na(lines$INT))) geom_path(data=lines, aes(x = INT, y =y.dat), colour = "black")} +
        geom_line(data=lines, aes(y = SIG, color = Trait), size=linesize, alpha=0.8, lineend = "round", show.legend = F) +
        scale_x_continuous(breaks=seq(0,max(lgs.size),cutx)) +
        {if(!all(is.na(lines$INT))) geom_point(data=points, aes(y = y.dat, color = Trait), shape = 2, size = 2, stroke = 1, alpha = 0.8)} +
        scale_y_continuous(breaks=seq(0,max(lgs.size, na.rm = T))) +
        guides(color = guide_legend("Trait"), fill = guide_legend("Trait"), shape = guide_legend("Trait")) + 
        labs(y = y.lab, x = "Position (cM)", subtitle="Linkage group") +
        theme_classic()
    }
  } else {
    pl <- list(lines = lines, points =points, linesize = linesize, 
               cutx = cutx, y.lab = y.lab)
    
    size <- table(pl$lines$Trait)[1]
    pl$lines$x <- rep(1:size, length(table(pl$lines$Trait)))
    pl$lines$x.int <- NA
    pl$lines$x.int[which(!is.na(pl$lines$INT))] <- pl$lines$x[which(!is.na(pl$lines$INT))]
    
    if(dim(points)[1] > 0){
      all <- paste0(pl$lines$Trait, "_", round(pl$lines$`Position (cM)`,2), "_", pl$lines$LG)
      point <- paste0(pl$points$Trait, "_", round(pl$points$`Position (cM)`,2), "_", pl$points$LG)
      pl$points$x <- pl$lines$x[match(point, all)]
    }
    
    pl$lines$SIG[which(pl$lines$SIG == "Inf")] <- NA ## Bugfix!!!
  }
  return(pl)
}

#' Only the plot part of plot_profile function
#' 
only_plot_profile <- function(pl.in){
  
  vlines <- split(pl.in$lines$x, pl.in$lines$LG)
  vlines <- sapply(vlines, function(x) x[1])
  
  pl <- ggplot(data = pl.in$lines, aes(x = x, color = Trait)) +
    {if(!all(is.na(pl.in$lines$INT))) geom_path(data=pl.in$lines, aes(x = x.int, y =y.dat), colour = "black")} +
    geom_line(data=pl.in$lines, aes(y = SIG, color = Trait), size=pl.in$linesize, alpha=0.8, show.legend = F) +
    {if(!all(is.na(pl.in$lines$INT))) geom_point(data=pl.in$points, aes(y = y.dat, color = Trait), shape = 2, size = 2, stroke = 1, alpha = 0.8)} +
    {if(length(vlines) > 1) geom_vline(xintercept=vlines, linetype="dashed", size=.5, alpha=0.8)} +  #threshold
    guides(color = guide_legend("Trait"), fill = guide_legend("Trait"), shape = guide_legend("Trait")) + 
    labs(y = pl.in$y.lab, x = "Linkage group") +
    annotate(x=vlines,y=+Inf,label= paste0("LG", names(vlines)),vjust=1, hjust= -0.1,geom="label") +
    ylim(c(min(pl.in$lines$y.dat),max(pl.in$lines$SIG, na.rm = T) + 3)) +
    theme_classic() + theme(axis.text.x=element_blank(),
                            axis.ticks.x=element_blank())
  
  return(pl)
}

#' Adapted function from QTLpoly
#' 
plot_qtlpoly.effects <- function(qtl_info, effects, pheno.col = NULL, p1 = "P1", p2 = "P2", df.info=NULL, lgs = NULL, position = NULL) {
  if(is.null(pheno.col)) {
    pheno.col <- 1:length(unique(qtl_info$pheno))
  } else {
    pheno.col <- which(unique(qtl_info$pheno) %in% pheno.col)
  }
  
  ploidy <- max(nchar(effects$haplo))
  
  qtl_info.sub <- qtl_info %>% filter(pheno %in% unique(qtl_info$pheno)[pheno.col]) %>%
    filter(Pos %in% position) %>% filter(LG %in% lgs)
  
  total <- split(qtl_info, qtl_info$pheno)
  total <- lapply(total, function(x) paste0(x[,1], "_", x[,2], "_", x[,5]))
  total <- total[match(unique(qtl_info$pheno), names(total))]
  
  sub <- split(qtl_info.sub, qtl_info.sub$pheno)
  sub <- lapply(sub, function(x) paste0(x[,1], "_", x[,2], "_", x[,5]))
  
  group.idx <- list()
  for(i in 1:length(pheno.col)){
    idx <- match(names(sub)[i], names(total))
    group.idx[[idx]] <- match(sub[[i]], total[[idx]])
  }
  
  plots2 <- list()
  for(p in pheno.col) {
    effects.sub <- effects %>% filter(pheno == unique(qtl_info$pheno)[p]) %>% 
      filter(qtl.id %in% group.idx[[p]]) 
    nqtl <- length(unique(effects.sub$qtl.id))
    if(nqtl > 0) {
      plots1 <- list()
      for(q in group.idx[[p]]) {
        data <- effects.sub %>% filter(qtl.id == q)
        if(ploidy == 4) {
          data <- data[1:36,]
          data <- data.frame(Estimates=as.numeric(data$effect), Alleles=data$haplo, Parent=c(rep(p1,4),rep(p2,4),rep(p1,14),rep(p2,14)), Effects=c(rep("Additive",8),rep("Digenic",28)))
          data <- data[-c(12:15,18:21,23:30),]
        }
        if(ploidy == 6) {
          data <- data[-c(18:23,28:33,37:42,45:50,52:63,83:88,92:97,100:105,107:133,137:142,145:150,152:178,181:186,188:214,216:278,299:1763),]
          data <- data.frame(Estimates=as.numeric(data$effect), Alleles=data$haplo, Parent=c(rep(p1,6),rep(p2,6),rep(p1,15),rep(p2,15),rep(p1,20),rep(p2,20)), Effects=c(rep("Additive",12),rep("Digenic",30),rep("Trigenic",40)))
        }
        data$Parent <- factor(data$Parent, levels=unique(data$Parent))
        plot <- ggplot(data[which(data$Effects == "Additive"),], aes(x = Alleles, y = Estimates, fill = Estimates)) +
          geom_bar(stat="identity") +
          scale_fill_gradient2(low = "red", high = "blue", guide = "none") +
          labs(title=unique(qtl_info$pheno)[p], subtitle=paste("QTL", q, "\n")) +
          facet_wrap(. ~ Parent, scales="free_x", ncol = 2, strip.position="bottom") +
          theme_minimal() +
          theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5), axis.text.x.bottom = element_text(hjust = 1, vjust = 0.5))
        plots1[[q]] <- plot
      }
      plots2[[p]] <- plots1
    }
  }
  p.t <- unlist(plots2, recursive = F)
  nulls <- which(sapply(p.t, is.null))
  if(length(nulls) > 0)  p.t <- p.t[-nulls]
  return(p.t)
}


#' Adapted function from QTLpoly
#' 
#' @import dplyr
#' @import tidyr
#' 
breeding_values <- function(qtl_info, probs, selected_mks, blups, beta.hat, pos) {
  pheno.names <- unique(as.character(qtl_info$pheno))
  results <- vector("list", length(pheno.names))
  names(results) <- pheno.names
  
  # possible to index individuals
  phenos <- which(pheno.names %in% names(pos))
  
  for(p in phenos) { # select pheno
    nqtl <- length(pos[[pheno.names[p]]])
    infos <- qtl_info %>% filter(pheno == pheno.names[p])
    infos <- infos[which(infos$Pos %in% pos[[pheno.names[p]]]),]
    markers <- which((round(selected_mks$pos,2) %in% infos$Pos) & (selected_mks$LG %in% infos$LG))
    Z <- probs[,markers,] # select by pos
    u.hat <- blups %>% filter(pheno == pheno.names[p])
    u.hat <- split(u.hat$u.hat, u.hat$qtl)
    
    beta.hat.sub <- beta.hat %>% filter(pheno == pheno.names[p])
    beta.hat.v <- beta.hat.sub$beta.hat
    
    Zu <- vector("list", nqtl)
    if(nqtl > 1) {
      for(m in 1:nqtl) {
        Zu[[m]] <- t(Z[,m,]) %*% u.hat[[m]]
      }
      nind <- dim(Z)[3]
      y.hat <- matrix(rep(beta.hat.v, nind), byrow = FALSE) + Reduce("+", Zu)
    } else if(nqtl == 1) {
      Zu <- t(Z) %*% u.hat[[1]]
      nind <- dim(Z)[2]
      y.hat <- matrix(rep(beta.hat.v, nind), byrow = FALSE) + Zu
    }
    
    colnames(y.hat) <- pheno.names[p]
    results[[p]] <- round(y.hat,2)
  }
  
  id.names <- rownames(results[[which(sapply(results, function(x) !is.null(x)))[1]]])
  results <- as.data.frame(do.call(cbind, results))
  results <- cbind(gen=id.names, results)
  
  return(results)
}

