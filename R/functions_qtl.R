#' Logarithm of \emph{P}-value (LOP) profile plots. Modified version of QTLpoly function.
#'
#' Plots profiled logarithm of score-based \emph{P}-values (LOP) from individual or combined traits.
#'
#' @import ggplot2
#' @import dplyr
#' @import tidyr
#' @importFrom plotly TeX
#' 
#' @rdname viewqtl
#' 
#' @keywords internal
plot_profile <- function(profile, qtl_info, selected_mks, pheno.col = NULL, 
                         lgs.id = NULL, by_range = TRUE, range.min = NULL, 
                         range.max = NULL, plot=TRUE, software = NULL) {
  
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
    trait.names <- unique(profile$pheno)[pheno.col[p]]
    if(!is.null(qtl_info)) {
      qtl_info.sub <- qtl_info %>% filter(pheno == trait.names) %>% filter(LG %in% lgs.id)
      if(dim(qtl_info.sub)[1] > 0){
        nqtls <-  qtl_info.sub %>% summarize(n()) 
        TRT <- qtl_info.sub$pheno
        LGS <- qtl_info.sub$LG
        POS <- qtl_info.sub$Pos
        INF <- qtl_info.sub$Pos_lower
        SUP <- qtl_info.sub$Pos_upper
        PVAL <- qtl_info.sub[,6]
        H2 <- qtl_info.sub$h2
        if(!is.null(H2)){
          points <- rbind(points, data.frame(TRT=TRT, LGS=LGS, POS=POS, INF=INF, SUP=SUP, PVAL = PVAL, H2 = round(H2,2)))
        } else 
          points <- rbind(points, data.frame(TRT=TRT, LGS=LGS, POS=POS, INF=INF, SUP=SUP, PVAL = PVAL))
        count <- count+1
        y.dat <- c(y.dat, rep((-0.5*count), nqtls))
      }
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
  } else if(y.lab == "deltaDIC") {
    lines$SIG <- -lines$SIG
    y.lab <- "-\U0394 DIC"
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
  
  if(dim(points)[1] > 0){
    dot.height <- data.frame(trt = unique(points$TRT), heigth = unique(points$y.dat))
    y.dat.lines <- dot.height$heigth[match(lines$TRT, dot.height$trt)]
    lines$y.dat <- y.dat.lines
    colnames(points)[1:3] <- c("Trait", "LG", "Position (cM)")
  } else lines$y.dat <- NA
  
  colnames(lines) <- c("Trait", "LG", "Position (cM)", "SIG", "INT", "range","y.dat")
  
  
  if(max(lgs.size[lgs.id]) > 200) cutx <- 150 else cutx <- 100
  if(length(lgs.size[lgs.id]) > 10) {linesize <- 1} else {cutx <- 50; linesize <- 1.25}
  
  lines$y.dat <- lines$y.dat + min(lines$SIG, na.rm = T)
  points$y.dat <- points$y.dat + min(lines$SIG, na.rm = T)
  
  scale.max <- round(max(lines$SIG[which(is.finite(lines$SIG))], na.rm = T),0)
  scale.max <- scale.max*1.2
  scale.min <- round(min(lines$SIG[which(is.finite(lines$SIG))], na.rm = T),0)
  
  if(scale.max > 50) {
    lines$y.dat <- lines$y.dat*3
    points$y.dat <- points$y.dat*3
    scale.each <- 10
  } else scale.each = 2 
  
  if(plot){
    if(by_range){
      pl <- ggplot(data = lines, aes(x = `Position (cM)`, color = Trait)) +
        facet_grid(.~LG, space = "free") +
        {if(!all(is.na(lines$INT))) geom_path(data=lines, aes(x = INT, y = y.dat), colour = "black")} +
        geom_line(data=lines, aes(y = range, color = Trait), size=linesize, alpha=0.8, lineend = "round") +
        geom_line(data=lines, aes(y = SIG, shape = Trait),  colour = "gray", size=linesize, alpha=0.8, lineend = "round") +
        scale_x_continuous(breaks=seq(0,max(lgs.size),cutx)) +
        {if(dim(points)[1] > 0) geom_point(data=points, aes(y = y.dat, color = Trait), shape = 2, size = 2, stroke = 1, alpha = 0.8)} +
        scale_y_continuous(breaks=seq(scale.min, scale.max,scale.each)) +
        guides(color = guide_legend("Trait"), color = guide_legend("Trait")) + 
        labs(y = y.lab, x = "Position (cM)", subtitle="Linkage group") + 
        theme_classic()
    } else {
      pl <- ggplot(data = lines, aes(x = `Position (cM)`, color = Trait, group=1)) +
        facet_grid(.~LG, space = "free") +
        {if(!all(is.na(lines$INT))) geom_path(data=lines, aes(x = INT, y =y.dat), colour = "black")} +
        geom_line(data=lines, aes(y = SIG, color = Trait), size=linesize, alpha=0.8, lineend = "round") +
        scale_x_continuous(breaks=seq(0,max(lgs.size),cutx)) +
        {if(dim(points)[1] > 0) geom_point(data=points, aes(y = y.dat, color = Trait), shape = 2, size = 2, stroke = 1, alpha = 0.8)} +
        scale_y_continuous(breaks=seq(scale.min, scale.max, scale.each)) +
        guides(color = guide_legend("Trait"), color = guide_legend("Trait")) + 
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
#' @rdname viewqtl
#' 
#' @keywords internal
only_plot_profile <- function(pl.in){
  
  vlines <- split(pl.in$lines$x, pl.in$lines$LG)
  vlines <- sapply(vlines, function(x) x[1])
  
  pl <- ggplot(data = pl.in$lines, aes(x = x)) +
    {if(!all(is.na(pl.in$lines$INT))) geom_path(data=pl.in$lines, aes(x = x.int, y =y.dat), colour = "black")} +
    geom_line(data=pl.in$lines, aes(y = SIG, color = Trait), size=pl.in$linesize, alpha=0.8) +
    #guides(color = guide_legend("Trait")) + 
    {if(dim(pl.in$points)[1] > 0) geom_point(data=pl.in$points, aes(y = y.dat, color = Trait), shape = 2, size = 2, stroke = 1, alpha = 0.8)} +
    {if(length(vlines) > 1) geom_vline(xintercept=vlines, linetype="dashed", size=.5, alpha=0.8)} +  #threshold
    labs(y = pl.in$y.lab, x = "Linkage group") +
    annotate(x=vlines,y=+Inf,label= paste0("LG", names(vlines)),vjust=1, hjust= -0.1,geom="label") +
    ylim(c(min(pl.in$lines$y.dat),max(pl.in$lines$SIG, na.rm = T) + 3)) +
    theme_classic() + theme(axis.text.x=element_blank(),
                            axis.ticks.x=element_blank())
  return(pl)
}

#' Get effects information
#' 
#' @rdname viewqtl
#' 
#' @keywords internal
data_effects <- function(qtl_info, effects, pheno.col = NULL, 
                         p1 = "P1", p2 = "P2", df.info=NULL, 
                         lgs = NULL, groups = NULL, position = NULL, 
                         software, design = c("bar", "circle", "digenic")) {
  if(is.null(pheno.col)) {
    pheno.col <- 1:length(unique(qtl_info$pheno))
  } else {
    pheno.col <- which(unique(qtl_info$pheno) %in% pheno.col)
  }
  
  if(software == "QTLpoly" | software == "diaQTL"){
    
    if(software == "QTLpoly"){
      ploidy <- max(nchar(effects$haplo))
    } else if(software == "diaQTL") {
      get.size <- effects %>% filter(pheno == unique(qtl_info$pheno)[1] & qtl.id == 1)
      if(nrow(get.size) == 64) ploidy = 4 else ploidy = 2
    } else if(software == "polyqtlR"){
      ploidy <- (dim(effects)[2] - 3)/2
    }
    
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
    
    plots2 <- all.additive <- list()
    count <-  count.p <- 1
    for(p in pheno.col) {
      effects.sub <- effects %>% filter(pheno == unique(qtl_info$pheno)[p]) %>% 
        filter(qtl.id %in% group.idx[[p]]) 
      nqtl <- length(unique(effects.sub$qtl.id))
      if(nqtl > 0) {
        plots1 <- list()
        count.q <- 1
        for(q in group.idx[[p]]) {
          data <- effects.sub %>% filter(qtl.id == q)
          if(ploidy == 4) {
            data <- data[1:36,]
            if(software == "diaQTL"){
              data <- data.frame(Estimates=as.numeric(data$effect), CI.lower = data$CI.lower, CI.upper = data$CI.upper, Alleles=data$haplo, Parent=c(rep(p1,4),rep(p2,4),rep(p1,14),rep(p2,14)), Effects=c(rep("Additive",8),rep("Digenic",28)))
            } else {
              data <- data.frame(Estimates=as.numeric(data$effect), Alleles=data$haplo, Parent=c(rep(p1,4),rep(p2,4),rep(p1,14),rep(p2,14)), Effects=c(rep("Additive",8),rep("Digenic",28)))
              data$Alleles <- gsub("a", paste0(p1,".1_"), data$Alleles)
              data$Alleles <- gsub("b", paste0(p1,".2_"), data$Alleles)
              data$Alleles <- gsub("c", paste0(p1,".3_"), data$Alleles)
              data$Alleles <- gsub("d", paste0(p1,".4_"), data$Alleles)
              data$Alleles <- gsub("e", paste0(p2,".1_"), data$Alleles)
              data$Alleles <- gsub("f", paste0(p2,".2_"), data$Alleles)
              data$Alleles <- gsub("g", paste0(p2,".3_"), data$Alleles)
              data$Alleles <- gsub("h", paste0(p2,".4_"), data$Alleles)
              data$Alleles = substring(data$Alleles,1, nchar(data$Alleles)-1)
            }
          }
          if(ploidy == 6) {
            data <- data[-c(18:23,28:33,37:42,45:50,52:63,83:88,92:97,100:105,107:133,137:142,145:150,152:178,181:186,188:214,216:278,299:1763),] # fix me
            #data <- data[1:298,]
            data <- data.frame(Estimates=as.numeric(data$effect), Alleles=data$haplo, Parent=c(rep(p1,6),rep(p2,6),rep(p1,15),rep(p2,15),rep(p1,20),rep(p2,20)), Effects=c(rep("Additive",12),rep("Digenic",30),rep("Trigenic",40)))
            data$Alleles <- gsub("a", paste0(p1,".1_"), data$Alleles)
            data$Alleles <- gsub("b", paste0(p1,".2_"), data$Alleles)
            data$Alleles <- gsub("c", paste0(p1,".3_"), data$Alleles)
            data$Alleles <- gsub("d", paste0(p1,".4_"), data$Alleles)
            data$Alleles <- gsub("e", paste0(p1,".5_"), data$Alleles)
            data$Alleles <- gsub("f", paste0(p1,".6_"), data$Alleles)
            data$Alleles <- gsub("g", paste0(p2,".1_"), data$Alleles)
            data$Alleles <- gsub("h", paste0(p2,".2_"), data$Alleles)
            data$Alleles <- gsub("i", paste0(p2,".3_"), data$Alleles)
            data$Alleles <- gsub("j", paste0(p2,".4_"), data$Alleles)
            data$Alleles <- gsub("k", paste0(p2,".5_"), data$Alleles)
            data$Alleles <- gsub("l", paste0(p2,".6_"), data$Alleles)
            data$Alleles = substring(data$Alleles,1, nchar(data$Alleles)-1)
          }
          data$Parent <- factor(data$Parent, levels=unique(data$Parent))
          if(design == "bar"){
            plot <- ggplot(data[which(data$Effects == "Additive"),], aes(x = Alleles, y = Estimates, fill = Estimates)) +
              geom_bar(stat="identity") +
              {if(software == "diaQTL") geom_errorbar(aes(ymin=CI.lower, ymax=CI.upper), width=.2, position=position_dodge(.9))} +
              scale_fill_gradient2(low = "red", high = "blue", guide = "none") +
              labs(title=unique(qtl_info$pheno)[p], subtitle=paste("LG:", sapply(strsplit(sub[[count.p]][count.q], "_"), "[",1), 
                                                                   "Pos:", sapply(strsplit(sub[[count.p]][count.q], "_"), "[",2))) +
              facet_wrap(. ~ Parent, scales="free_x", ncol = 2, strip.position="bottom") +
              theme_minimal() +
              theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5), axis.text.x.bottom = element_text(hjust = 1, vjust = 0.5))
            plots1[[q]] <- plot
          } else if(design == "digenic"){
            if(ploidy == 6) data <- data[1:42,]
            temp <- do.call(rbind, strsplit(data$Alleles, "_"))
            data$x <- temp[,1]
            data$y <- temp[,2]
            digenic.effects <- data[which(data$Effects == "Digenic"),]
            additive.effects <- data[which(data$Effects == "Additive"),]
            plot.data <- data.frame(x= c(digenic.effects$x, digenic.effects$y),
                                    y= c(digenic.effects$y, digenic.effects$x),
                                    z= c(digenic.effects$Estimates, digenic.effects$Estimates+
                                           additive.effects$Estimates[match(digenic.effects$x, additive.effects$Alleles)] + 
                                           additive.effects$Estimates[match(digenic.effects$y, additive.effects$Alleles)]))
            
            plot = ggplot(data= plot.data,aes(x= x, y= y, fill= z)) + 
              geom_tile() + scale_fill_gradient2(name="") + 
              labs(title = paste("Trait:", unique(qtl_info$pheno)[p]),
                   subtitle = paste("LG:", sapply(strsplit(sub[[count.p]][count.q], "_"), "[",1), 
                                    "Pos:", sapply(strsplit(sub[[count.p]][count.q], "_"), "[",2))) +
              theme_bw() + xlab("") + ylab("") +
              theme(text = element_text(size=13),axis.text.x = element_text(angle = 90,vjust=0.5,hjust=1)) +
              coord_fixed(ratio=1)
            plots1[[q]] <- plot
          } else if(design == "circle"){
            additive.effects <- data[which(data$Effects == "Additive"),]
            additive.effects$pheno <- unique(qtl_info$pheno)[p]
            additive.effects$qtl_id <- q
            additive.effects$LG <- qtl_info.sub$LG[count]
            additive.effects$Pos <- qtl_info.sub$Pos[count]
            count <- count + 1
            all.additive <- rbind(all.additive, additive.effects)
            plots1 <- NULL
          }
          count.q <- count.q + 1
        }
        plots2[[p]] <- plots1
      }
      count.p <- count.p + 1
    }
    if(design != "circle"){
      p <- unlist(plots2, recursive = F)
      nulls <- which(sapply(p, is.null))
      if(length(nulls) > 0)  p <- p[-nulls]
      return(p)
    } else {
      all.additive$unique.id <- paste0(all.additive$pheno, "/LG:", all.additive$LG, "/Pos:", all.additive$Pos)
      breaks <- seq(round(min(all.additive$Estimates),2), round(max(all.additive$Estimates),2),max(all.additive$Estimates)/2)
      
      lgs <- unique(all.additive$LG)
      p <- list()
      for(i in 1:length(lgs)){
        p[[i]] <- ggplot(data=all.additive[which(all.additive$LG == lgs[i]),], 
                         aes(x=Alleles, y=Estimates, group=unique.id, colour=unique.id)) +
          geom_point() + 
          geom_line() + 
          coord_polar() +
          ylim(min(all.additive$Estimates),max(all.additive$Estimates)) +
          labs(title = paste0("LG", lgs[i])) +
          annotate(x= 0,y=breaks, label= round(breaks,2),geom="text") +
          theme_bw() + 
          theme(axis.title.y=element_blank(),
                axis.text.y=element_blank(),
                axis.ticks.y=element_blank(),
                axis.title.x=element_blank(), 
                legend.title = element_blank()) +
          geom_hline(yintercept=breaks) 
      }
      return(p) 
    }
  } else if(software == "polyqtlR"){
    if(design == "circle" | design == "digenic"){
      stop("Design option not available for: polyqtlR")
    } else {
      effects.df <- effects %>% filter(pheno %in% unique(qtl_info$pheno)[pheno.col]) %>% 
        filter(LG %in% groups) %>% pivot_longer(cols = 4:ncol(.), names_to = "haplo", values_to = "effect") 
      
      effects.df <- effects.df %>% group_by(haplo, pheno) %>% mutate(x.axis = 1:n()) %>% ungroup() %>% as.data.frame()
      
      vlines <- split(effects.df$x.axis, effects.df$LG)
      vlines <- sapply(vlines, function(x) x[1])
      
      p <- list()
      for(i in 1:length(pheno.col)){
        p[[i]] <- effects.df %>% filter(pheno == unique(qtl_info$pheno)[pheno.col][i]) %>%  
          ggplot() +
          geom_path(aes(x=x.axis, y=haplo, col = effect), size = 5)  +
          scale_color_gradient2(low = "purple4", mid = "white",high = "seagreen") +
          {if(length(vlines) > 1) geom_vline(xintercept=vlines, linetype="dashed", size=.5, alpha=0.8)} +
          labs(y = "Haplotype", x = "Linkage group", title = unique(qtl_info$pheno)[pheno.col][i]) +
          annotate(x=vlines,y=+Inf,label= paste0("LG", names(vlines)),vjust= 1, hjust= -0.1,geom="label") +
          coord_cartesian(ylim = c(1,8.5)) +
          theme_classic() + theme(axis.text.x=element_blank(),
                                  axis.ticks.x=element_blank(), legend.title = element_blank())
      }
      return(p)
    }
  }
}

#' Plot effects data
#' 
#' @rdname viewqtl
#' 
#' @keywords internal
plot_effects <- function(data_effects.obj, software, 
                         design = c("bar", "circle", "digenic")){
  if(software == "polyqtlR"){
    p.t <- ggarrange(plotlist = data_effects.obj, common.legend = T, ncol = 1, legend = "right")
  } else {
    if(design == "circle"){
      rows <- ceiling(length(data_effects.obj)/3)
      if(rows == 0) rows <- 1
      p.t <- ggarrange(plotlist = data_effects.obj, nrow = rows, ncol = 3)
    } else {
      rows <- ceiling(length(data_effects.obj)/4)
      if(rows == 0) rows <- 1
      p.t <- ggarrange(plotlist = data_effects.obj, nrow = rows, ncol = 4)
    }
  }
  return(p.t)
}

#' Adapted function from QTLpoly
#' 
#' @import dplyr
#' @import tidyr
#' 
#' @rdname viewqtl
#' 
#' @keywords internal
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

#' Adapted from MAPpoly
#' 
#' @rdname viewqtl
#' 
#' @keywords internal
calc_homologprob  <- function(probs, selected_mks){
  
  input.genoprobs <- probs
  each.split <- vector()
  sizes <- table(selected_mks$LG)
  probs.b <- probs
  input.genoprobs <- list()
  for(i in 1:length(sizes)){
    input.genoprobs[[i]] <- probs.b[,1:sizes[i],]
    probs.b <- probs.b[,-c(1:sizes[i]),]
  }
  
  pos <- split(selected_mks$pos, selected_mks$LG)
  
  df.res <- NULL
  for(j in 1:length(input.genoprobs)){
    stt.names <- dimnames(input.genoprobs[[j]])[[1]] ## state names
    mrk.names <- dimnames(input.genoprobs[[j]])[[2]] ## mrk names
    ind.names <- dimnames(input.genoprobs[[j]])[[3]] ## individual names
    v <- c(2,4,6,8,10,12)
    names(v) <- choose(v,v/2)^2
    ploidy <- v[as.character(length(stt.names))]
    hom.prob <- array(NA, dim = c(ploidy*2, length(mrk.names), length(ind.names)))
    dimnames(hom.prob) <- list(letters[1:(2*ploidy)], mrk.names, ind.names)
    for(i in letters[1:(2*ploidy)])
      hom.prob[i,,] <- apply(input.genoprobs[[j]][grep(stt.names, pattern = i),,], c(2,3), function(x) round(sum(x, na.rm = TRUE),4))
    df.hom <- reshape2::melt(hom.prob)
    map <- data.frame(map.position = pos[[j]], marker = mrk.names)
    colnames(df.hom) <- c("homolog", "marker", "individual", "probability")
    df.hom <- merge(df.hom, map, sort = FALSE)
    df.hom$LG <- names(pos)[j]
    df.res <- rbind(df.res, df.hom)
  }
  if(ploidy == 4){
    df.res$homolog <- gsub("a", paste0("P1.1_"), df.res$homolog)
    df.res$homolog <- gsub("b", paste0("P1.2_"), df.res$homolog)
    df.res$homolog <- gsub("c", paste0("P1.3_"), df.res$homolog)
    df.res$homolog <- gsub("d", paste0("P1.4_"), df.res$homolog)
    df.res$homolog <- gsub("e", paste0("P2.1_"), df.res$homolog)
    df.res$homolog <- gsub("f", paste0("P2.2_"), df.res$homolog)
    df.res$homolog <- gsub("g", paste0("P2.3_"), df.res$homolog)
    df.res$homolog <- gsub("h", paste0("P2.4_"), df.res$homolog)
    df.res$homolog = substring(df.res$homolog,1, nchar(df.res$homolog)-1)
  } else if(ploidy == 6){
    df.res$homolog <- gsub("a", paste0("P1.1_"), df.res$homolog)
    df.res$homolog <- gsub("b", paste0("P1.2_"), df.res$homolog)
    df.res$homolog <- gsub("c", paste0("P1.3_"), df.res$homolog)
    df.res$homolog <- gsub("d", paste0("P1.4_"), df.res$homolog)
    df.res$homolog <- gsub("e", paste0("P1.5_"), df.res$homolog)
    df.res$homolog <- gsub("f", paste0("P1.6_"), df.res$homolog)
    df.res$homolog <- gsub("g", paste0("P2.1_"), df.res$homolog)
    df.res$homolog <- gsub("h", paste0("P2.2_"), df.res$homolog)
    df.res$homolog <- gsub("i", paste0("P2.3_"), df.res$homolog)
    df.res$homolog <- gsub("j", paste0("P2.4_"), df.res$homolog)
    df.res$homolog <- gsub("k", paste0("P2.5_"), df.res$homolog)
    df.res$homolog <- gsub("l", paste0("P2.6_"), df.res$homolog)
    df.res$homolog = substring(df.res$homolog,1, nchar(df.res$homolog)-1)
  }
  
  structure(list(info = list(ploidy = ploidy, 
                             n.ind = length(ind.names)) , 
                 homoprob = df.res), class = "mappoly.homoprob")
}

#' Plots mappoly.homoprob from MAPpoly
#' 
#' @param x an object of class \code{mappoly.homoprob}
#' 
#' @param stack logical. If \code{TRUE}, probability profiles of all homologues
#'              are stacked in the plot (default = FALSE)
#'              
#' @param lg indicates which linkage group should be plotted. If \code{NULL} 
#'           (default), it plots the first linkage group. If 
#'           \code{"all"}, it plots all linkage groups
#'           
#' @param ind indicates which individuals should be plotted. It can be the 
#'            position of the individuals in the dataset or it's name. 
#'            If \code{NULL} (default), the function plots the first 
#'            individual
#'            
#' @param use.plotly if \code{TRUE} (default), it uses plotly interactive graphic
#'
#' @param verbose if \code{TRUE} (default), the current progress is shown; if
#'     \code{FALSE}, no output is produced
#' 
#' @param ... unused arguments
#' @importFrom plotly ggplotly
#' 
#' @rdname viewqtl
#' 
#' @keywords internal
plot.mappoly.homoprob <- function(x, stack = FALSE, lg = NULL, 
                                  ind = NULL, use.plotly = TRUE, 
                                  verbose = TRUE,  ...){
  all.ind <- as.character(unique(x$homoprob$individual))
  #### Individual handling ####
  if(length(ind) > 1){
    if (verbose) message("More than one individual provided: using the first one")
    ind <- ind[1]
  }
  if(is.null(ind)){
    ind <- as.character(all.ind[1])
    df.pr1 <- subset(x$homoprob, individual  ==  ind)  
  } else if(is.numeric(ind)) {
    if(ind > length(all.ind))
      stop("Please chose an individual number between 1 and ", length(all.ind))
    ind <- as.character(all.ind[ind])
    df.pr1 <- subset(x$homoprob, individual  ==  ind)  
  } else if (is.character(ind)){
    if(!ind%in%all.ind)
      stop("Invalid individual name")
  } else stop("Invalid individual name")
  #### LG handling ####
  if(is.null(lg))
    lg <- 1
  if(all(lg == "all"))
    lg <- unique(x$homoprob$LG)
  LG <- individual <- map.position <- probability <- homolog <- NULL
  if(length(lg) > 1 & !stack)
  {
    if (verbose) message("Using 'stack = TRUE' to plot multiple linkage groups")
    stack <- TRUE
  }
  if(stack){
    ##subset linkage group
    if(!is.null(lg)){
      df.pr1 <- subset(x$homoprob, LG%in%lg)
      df.pr1 <- subset(df.pr1, individual  ==  ind)
    } else 
      df.pr1 <- subset(x$homoprob, individual  ==  ind)
    p <- ggplot2::ggplot(df.pr1, ggplot2::aes(x = map.position, y = probability, fill = homolog, color  = homolog)) +
      ggplot2::geom_density(stat = "identity", alpha = 0.7, position = "stack") + 
      ggplot2::ggtitle(ind) + 
      ggplot2::facet_grid(rows = ggplot2::vars(LG)) + 
      ggplot2::ylab(label = "Homologs probabilty") +
      ggplot2::xlab(label = "Map position")
  } else {
    ##subset linkage group
    if(is.null(lg)){
      lg <- 1
      df.pr1 <- subset(x$homoprob, LG %in% lg)
    } else df.pr1 <- subset(x$homoprob, LG %in% lg)
    df.pr1 <- subset(df.pr1, individual  ==  ind)  
    p <- ggplot2::ggplot(df.pr1, ggplot2::aes(x = map.position, y = probability, fill = homolog, color  = homolog)) +
      ggplot2::geom_density(stat = "identity", alpha = 0.7) + 
      ggplot2::ggtitle(paste(ind, "   LG", lg)) + 
      ggplot2::facet_grid(rows = ggplot2::vars(homolog)) + 
      ggplot2::theme_minimal() + 
      ggplot2::ylab(label = "Homologs probabilty") +
      ggplot2::xlab(label = "Map position")
  }
  if(use.plotly)
    p <- plotly::ggplotly(p)
  return(p)
}

#' Plot selected haplotypes
#' 
#' @rdname viewqtl
#' 
#' @keywords internal
select_haplo <- function(input.haplo, probs, selected_mks, effects.data){
  print(input.haplo)
  lgs <- sapply(strsplit(unlist(input.haplo), "_"),function(x) x[grep("LG", x)])
  lgs <- gsub("LG:", "", lgs)
  selec.lg <- selected_mks %>% filter(LG %in% lgs)
  homo.dat <- calc_homologprob(probs = probs, selected_mks = selec.lg)
  pos <- sapply(strsplit(unlist(input.haplo), "_"),function(x) x[grep("Pos", x)])
  pos <- gsub("Pos:", "", pos)
  homo <- sapply(strsplit(unlist(input.haplo), "_"),function(x) x[grep("homolog", x)])
  homo <- gsub("homolog:", "", homo)
  alleles <- effects.data[[1]]$data$Alleles[!grepl("_",effects.data[[1]]$data$Alleles)]
  alleles <- rep(alleles, length(homo))
  idx <- match(homo, alleles)
  
  like.ind.all <- list()
  for(i in 1:length(pos)){
    homoprob_temp <- homo.dat$homoprob %>% 
      filter(map.position %in% pos[i]) %>% filter(LG %in% lgs[i])
    homoprob_temp <- homoprob_temp[order(homoprob_temp$individual, homoprob_temp$homolog),]
    homoprob_temp <- homoprob_temp %>% 
      group_by(map.position, LG, individual) %>% 
      summarise(best = which(probability > 0.5))
    like.ind <-  homoprob_temp$individual[which(homoprob_temp$best %in% idx[i])]
    if(length(like.ind) ==0) like.ind <- NA
    like.ind.all[[i]] <-  like.ind
  }
  like.intersect <- Reduce(intersect, like.ind.all)
  if(length(like.intersect) == 0 | all(is.na(like.intersect))) stop("Any individual contain all the selected homolog/s")
  p <- list()
  for(i in 1:length(like.intersect)){
    p[[i]] <- plot.mappoly.homoprob(homo.dat, 
                                    lg = unique(as.numeric(lgs)), 
                                    ind = as.character(like.intersect)[i], 
                                    use.plotly = FALSE)
  }
  return(p)
}

