draw_map_shiny<-function(left.lim = 0, right.lim = 5, ch = 1,
                         maps, ph.p, ph.q, dp, dq, snp.names=TRUE)
{
  ch <- as.numeric(ch)
  ploidy <- dim(ph.p[[1]])[2]
  # if(is.character(ch))
  #   ch <- as.numeric(strsplit(ch, split = " ")[[1]][3])
  var.col <- RColorBrewer::brewer.pal(n = 4, name = "Set1")
  names(var.col) <- c("A", "T", "C", "G")
  d.col<-c(NA, RColorBrewer::brewer.pal(n = ploidy, name = "Dark2"))
  names(d.col) <- 0:ploidy
  d.col[1]<-NA
  x <- maps[[ch]]
  lab <- names(maps[[ch]])
  zy <- seq(0, 0.5, length.out = ploidy) + 1.5
  pp1 <- ph.p[[ch]]
  pp2 <- ph.q[[ch]]
  dp1<-dq[[ch]]
  dp2<-dp[[ch]]
  x1<-abs(left.lim - x)
  x2<-abs(right.lim - x)
  id.left<-which(x1==min(x1))[1]
  id.right<-rev(which(x2==min(x2)))[1]
  par(mai = c(0.5,0.15,0,0))
  curx<-x[id.left:id.right]
  plot(x = curx,
       y = rep(.5,length(curx)),
       type = "n" , 
       ylim = c(.25, 4.5), 
       #xlim = c(min(curx), max(curx)),
       axes = FALSE)
  lines(c(x[id.left], x[id.right]), c(.5, .5), lwd=15, col = "gray")
  points(x = curx,
         y = rep(.5,length(curx)),
         xlab = "", ylab = "", 
         pch = "|", cex=1.5, 
         ylim = c(0,2))
  axis(side = 1)
  mtext(text = "Distance (cM)", side = 1, adj = 1)
  #Parent TANZANIA
  for(i in 1:ploidy)
  {
    lines(c(x[id.left], x[id.right]), c(zy[i], zy[i]), lwd=10, col = "gray")
    points(x = seq(x[id.left], x[id.right], length.out = length(curx)),
           y = rep(zy[i], length(curx)),
           col = var.col[pp2[id.left:id.right,i]],
           pch = 15,
           cex = 2)
  }
  mtext(text = "Parent 2", side = 2, at = mean(zy), line = -3, font = 4)
  for(i in 1:ploidy)
    mtext(letters[12:7][i], at = zy[i], side = 2,  line = -4, font = 1)
  connect.lines<-seq(x[id.left], x[id.right], length.out = length(curx))
  for(i in 1:length(connect.lines))
    lines(c(curx[i], connect.lines[i]), c(0.575, zy[1]-.05), lwd=0.3)
  points(x = seq(x[id.left], x[id.right], length.out = length(curx)),
         y = zy[ploidy]+0.05+dp2[id.left:id.right]/20,
         col = d.col[as.character(dp2[id.left:id.right])],
         pch = 19, cex = .7)
  corners = par("usr") 
  par(xpd = TRUE) 
  text(x = corners[1]+.5, y = mean(zy[ploidy]+0.05+(1:ploidy/20)), "Doses")
  #Parent BEAUREGARD
  zy<-zy+1.1
  for(i in 1:ploidy)
  {
    lines(c(x[id.left], x[id.right]), c(zy[i], zy[i]), lwd=10, col = "gray")
    points(x = seq(x[id.left], x[id.right], length.out = length(curx)),
           y = rep(zy[i], length(curx)),
           col = var.col[pp1[id.left:id.right,i]],
           pch = 15,
           cex = 2)
  }
  mtext(text = "Parent 1", side = 2, at = mean(zy), line = -3, font = 4)
  points(x = seq(x[id.left], x[id.right], length.out = length(curx)),
         y = zy[ploidy]+0.05+dp1[id.left:id.right]/20,
         col = d.col[as.character(dp1[id.left:id.right])],
         pch = 19, cex = .7)
  corners = par("usr") 
  par(xpd = TRUE) 
  text(x = corners[1]+.5, y = mean(zy[ploidy]+0.05+(1:ploidy/20)), "Doses")
  if(snp.names)
    text(x = seq(x[id.left], x[id.right], length.out = length(curx)),
         y = rep(zy[ploidy]+0.05+.3, length(curx)),
         labels = names(curx),
         srt=90, adj = 0, cex = .7)
  for(i in 1:ploidy)
    mtext(letters[ploidy:1][i], at = zy[i], side = 2,  line = -4, font = 1)
  # legend("bottomleft", legend=c("A", "T", "C", "G", "-"), 
  #        fill =c(var.col, "white"), 
  #        box.lty=0, bg="transparent")
}

## Function from MAPpoly
imf_h <- function(r) {
  r[r >= 0.5] <- 0.5 - 1e-14
  -50 * log(1 - 2 * r)
}

map_summary<-function(left.lim = 0, right.lim = 5, ch = 1,
                      maps, dp, dq){
  ch <- as.numeric(ch)
  # if(is.character(ch))
  #   ch <- as.numeric(strsplit(ch, split = " ")[[1]][3])
  x <- maps[[ch]]
  lab <- names(maps[[ch]])
  ploidy = max(c(dp[[ch]], dq[[ch]])) 
  dp1<-dq[[ch]]
  dp2<-dp[[ch]]
  x1<-abs(left.lim - x)
  x2<-abs(right.lim - x)
  id.left<-which(x1==min(x1))[1]
  id.right<-rev(which(x2==min(x2)))[1]
  curx<-x[id.left:id.right]
  w<-table(paste(dp1[id.left:id.right], dp2[id.left:id.right], sep = "-"))
  M<-matrix(0, nrow = ploidy+1, ncol = ploidy+1, dimnames = list(0:ploidy, 0:ploidy))
  for(i in as.character(0:ploidy))
    for(j in as.character(0:ploidy))
      M[i,j]<-w[paste(i,j,sep = "-")]
  M[is.na(M)]<-0
  return(list(doses = M, number.snps = length(curx), length = diff(range(curx)), cM.per.snp = round(diff(range(curx))/length(curx), 3), full.size = as.numeric(maps[[ch]][length(maps[[ch]])]), number.of.lgs = length(maps)))
}


#' Logarithm of \emph{P}-value (LOP) profile plots. Modified version of QTLpoly function.
#'
#' Plots profiled logarithm of score-based \emph{P}-values (LOP) from individual or combined traits.
#'
#' @param lgs.info level $lgs and $lgs.all of an object of class \code{qtlpoly.data}. 
#' \item{lgs}{a list with selected marker positions per linkage group.}
#' \item{lgs.all}{a list with all marker positions per linkage group.}
#'
#' @param model an object of class \code{qtlpoly.profile} or \code{qtlpoly.remim}.
#'
#' @param pheno.col a numeric vector with the phenotype column numbers to be plotted; if \code{NULL}, all phenotypes from \code{'data'} will be included.
#'
#' @param sup.int if \code{TRUE}, support interval are shown as shaded areas; if \code{FALSE} (default), no support interval is show.
#' 
#' @param main a character string with the main title; if \code{NULL}, no title is shown.
#'
#' @param legend legend position (either "bottom", "top", "left" or "right"); if \code{NULL}, no legend is shown.
#'
#' @param ylim a numeric value pair supplying the limits of y-axis, e.g. c(0,10); if \code{NULL} (default), limits will be provided automatically.
#'
#' @param grid if \code{TRUE}, profiles will be organized in rows (one per trait); if \code{FALSE} (default), profiles will appear superimposed. Only effective when plotting profiles from more than one trait.
#'
#' @return A \pkg{ggplot2} with the LOP profiles for each trait.
#'
#' @seealso \code{\link[qtlpoly]{profile_qtl}},  \code{\link[qtlpoly]{remim}}
#'
#' @examples
#'   \dontrun{
#'   # load raw data
#'   data(maps)
#'   data(pheno)
#'
#'   # estimate conditional probabilities using mappoly package
#'   library(mappoly)
#'   genoprob <- lapply(maps, calc_genoprob)
#'
#'   # prepare data
#'   data <- read_data(ploidy = 6, geno.prob = genoprob, pheno = pheno, step = 1)
#'
#'   # perform remim
#'   remim.mod <- remim(data = data, w.size = 15, sig.fwd = 0.01, sig.bwd = 0.0001,
#'     d.sint = 1.5, n.clusters = 4, plot = "remim")
#'
#'   # plot profiles
#'   for (p in remim.mod$pheno.col) { 
#'     plot_profile(data = data, model = remim.mod, pheno.col = p, ylim = c(0, 10))
#'   } # separate plots
#'     
#'   plot_profile(data = data, model = remim.mod, grid = FALSE) # combined plots
#'   }
#'
#' @author Guilherme da Silva Pereira, \email{gdasilv@@ncsu.edu}
#'
#' @references
#'     Pereira GS, Gemenet DC, Mollinari M, Olukolu BA, Wood JC, Mosquera V, Gruneberg WJ, Khan A, Buell CR, Yencho GC, Zeng ZB (2020) Multiple QTL mapping in autopolyploids: a random-effect model approach with application in a hexaploid sweetpotato full-sib population, \emph{Genetics} 215 (3): 579-595. \url{http://doi.org/10.1534/genetics.120.303080}.
#'
#' @export plot_profile
#' @import ggplot2
#' 
plot_profile <- function(lgs.info, model = model, pheno.col = NULL, sup.int = TRUE, 
                         main = NULL, legend="bottom", ylim = NULL, grid = TRUE, 
                         lgs.id = NULL, range.min = NULL, range.max = NULL, by_range = TRUE, plot=TRUE) {
  
  lgs.size <- sapply(lgs.info[[2]], function(x) x[length(x)])
  lines <- points <- thre <- map <- data.frame()
  y.dat <- trait.names <- c()
  count <- 0
  if(is.null(pheno.col)) pheno.col <- model$pheno.col
  nphe <- length(pheno.col)
  LGS <- c(); for(c in 1:length(lgs.info[[1]])) LGS <- c(LGS, rep(c, length(lgs.info[[1]][[c]])))
  POS <- unlist(lgs.info[[1]])
  for(p in 1:nphe) { #lines
    t <- which(model$pheno.col == pheno.col[p])
    TRT <- rep(names(model$results)[t], length(LGS))
    if(any(class(model) == "qtlpoly.feim")) SIG <- model$results[[t]][[3]] else SIG <- -log10(as.numeric(model$results[[t]][[3]]))
    lines <- rbind(lines, data.frame(TRT=as.factor(TRT), LGS=LGS, POS=POS, SIG=SIG))
  }
  for(p in 1:nphe) { #points
    t <- which(model$pheno.col == pheno.col[p])
    trait.names <- c(trait.names, names(model$results)[t])
    if(!is.null(model$results[[t]]$qtls)) {
      nqtls <- dim(model$results[[t]]$qtls)[1]
      TRT <- rep(names(model$results)[t], nqtls)
      LGS <- model$results[[t]]$qtls[,"LG"]
      POS <- model$results[[t]]$qtls[,"Pos"]
      INF <- model$results[[t]]$lower[,"Pos_lower"]
      SUP <- model$results[[t]]$upper[,"Pos_upper"]
      points <- rbind(points, data.frame(TRT=TRT, LGS=LGS, POS=POS, INF=INF, SUP=SUP))
      count <- count+1
      y.dat <- c(y.dat, rep((-0.3*count), nqtls))
    }
  }
  points$TRT <- factor(points$TRT, levels=trait.names)
  if(any(class(model) == "qtlpoly.feim")) {
    for(p in 1:nphe) { #threshold
      t <- which(model$pheno.col == pheno.col[p])
      LGS <- c(1:length(lgs.info[[1]]))
      TRT <- rep(names(model$results)[t], length(LGS))
      SIG <- rep(model$sig.lod[t], length(LGS))
      thre <- rbind(thre, data.frame(TRT=as.factor(TRT), LGS=LGS, SIG=SIG))
      y.lab <- "LOD"
    }
  } else {
    # y.lab <- "LOP"
    y.lab <- expression(-log[10](italic(P)))
  }
  if(is.null(y.dat)) y.dat <- ylim[1]
  if(grid) y.dat <- -max(lines$SIG[is.finite(lines$SIG)])/7
  
  for(c in 1:length(lgs.info[[1]])) {
    LGS <- rep(c, length(lgs.info[[2]][[c]]))
    map <- rbind(map, data.frame(LGS=LGS, POS=lgs.info[[2]][[c]]))
  }
  if(max(lgs.size) > 200) cutx <- 150 else cutx <- 100
  if(max(lines$SIG[is.finite(lines$SIG)]) < 10) cuty <- 2 else cuty <- 4
  if(length(lgs.info[[1]]) > 10) { addx <- 50; linesize <- 1} else { addx <- 10 ; cutx <- 50; linesize <- 1.25}
  
  # Filter group
  if(!is.null(lgs.id)){
    lines <- lines[which(lines$LGS %in% lgs.id),]
    points <- points[which(points$LGS %in% lgs.id),]
  }
  
  lines$INT <- NA
  
  for(i in 1:dim(points)[1]){
    lines$INT[which(lines$POS >= points$INF[i] & lines$POS <= points$SUP[i] & lines$LGS == points$LGS[i])] <- lines$POS[which(lines$POS >= points$INF[i] & lines$POS <= points$SUP[i] & lines$LGS == points$LGS[i])]
  }
  
  # Filter position
  lines$range <- NA
  if(!is.null(range.min)){
    lines$range[which(lines$POS >= range.min & lines$POS <= range.max)] <- lines$SIG[which(lines$POS >= range.min & lines$POS <= range.max)]
    lines$SIG[which(lines$POS > range.min & lines$POS < range.max)] <- NA
  }
  
  colnames(lines) <- c("Trait", "LG", "Position (cM)", "LOP", "INT", "range")
  colnames(points)[1:3] <- c("Trait", "LG", "Position (cM)")
  points <- cbind(points, y.dat)
  
  if(plot){
    if(by_range){
      pl <- ggplot(data = lines, aes(x = `Position (cM)`, color = Trait)) +
        {if(!all(is.na(lines$INT)) & sup.int) geom_path(data=lines, aes(x = INT, y =y.dat), colour = "black")} +
        geom_line(data=lines, aes(y = range, color = Trait), size=linesize, alpha=0.8, lineend = "round", show.legend = F) +
        geom_line(data=lines, aes(y = LOP, shape = Trait),  colour = "gray", size=linesize, alpha=0.8, lineend = "round") +
        scale_x_continuous(breaks=seq(0,max(lgs.size),cutx)) +
        {if(!all(is.na(lines$INT))) geom_point(data=points, aes(y = y.dat, color = Trait), shape = 2, size = 2, stroke = 1, alpha = 0.8)} +
        scale_y_continuous(breaks=seq(0,max(lgs.size, na.rm = T))) +
        {if(nrow(thre) > 0) geom_hline(data=thre, aes(yintercept=LOP, color=Trait), linetype="dashed", size=.5, alpha=0.8)} +  #threshold
        guides(color = guide_legend("Trait"), fill = guide_legend("Trait"), shape = guide_legend("Trait")) + 
        labs(title=main, y = "LOP", x = "Position (cM)", subtitle="Linkage group") +
        theme_classic()
    } else {
      pl <- ggplot(data = lines, aes(x = `Position (cM)`, color = Trait)) +
        facet_grid(.~LG, space = "free") +
        {if(!all(is.na(lines$INT)) & sup.int) geom_path(data=lines, aes(x = INT, y =y.dat), colour = "black")} +
        geom_line(data=lines, aes(y = LOP, color = Trait), size=linesize, alpha=0.8, lineend = "round", show.legend = F) +
        scale_x_continuous(breaks=seq(0,max(lgs.size),cutx)) +
        {if(!all(is.na(lines$INT))) geom_point(data=points, aes(y = y.dat, color = Trait), shape = 2, size = 2, stroke = 1, alpha = 0.8)} +
        scale_y_continuous(breaks=seq(0,max(lgs.size, na.rm = T))) +
        {if(nrow(thre) > 0) geom_hline(data=thre, aes(yintercept=LOP, color=Trait), linetype="dashed", size=.5, alpha=0.8)} +  #threshold
        guides(color = guide_legend("Trait"), fill = guide_legend("Trait"), shape = guide_legend("Trait")) + 
        labs(title=main, y = "LOP", x = "Position (cM)", subtitle="Linkage group") +
        theme_classic()
    }
  } else {
    pl <- list(lines = lines, points =points, thre =thre, 
               sup.int = sup.int, linesize = linesize, 
               lgs.size = lgs.size, cutx = cutx, 
               main = main, y.dat =y.dat)
  }
  return(pl)
}

#' Only the plot part of plot_profile function
#' 
only_plot_profile <- function(pl.in, by_range=FALSE){
  if(by_range){
    pl <- ggplot(data = pl.in$lines, aes(x = `Position (cM)`, color = Trait)) +
      {if(!all(is.na(pl.in$lines$INT)) & pl.in$sup.int) geom_path(data=pl.in$lines, aes(x = INT, y =pl.in$y.dat), colour = "black")} +
      geom_line(data=pl.in$lines, aes(y = range, color = Trait), size=pl.in$pl.in$linesize, alpha=0.8, lineend = "round", show.legend = F) +
      geom_line(data=pl.in$lines, aes(y = LOP, shape = Trait),  colour = "gray", size=pl.in$pl.in$linesize, alpha=0.8, lineend = "round") +
      scale_x_continuous(breaks=seq(0,max(pl.in$lgs.size),pl.in$cutx)) +
      {if(!all(is.na(pl.in$lines$INT))) geom_point(data=pl.in$points, aes(y = pl.in$y.dat, color = Trait), shape = 2, size = 2, stroke = 1, alpha = 0.8)} +
      scale_y_continuous(breaks=seq(0,max(pl.in$lgs.size, na.rm = T))) +
      {if(nrow(pl.in$thre) > 0) geom_hline(data=pl.in$thre, aes(yintercept=LOP, color=Trait), linetype="dashed", size=.5, alpha=0.8)} +  #threshold
      guides(color = guide_legend("Trait"), fill = guide_legend("Trait"), shape = guide_legend("Trait")) + 
      labs(title=pl.in$main, y = "LOP", x = "Position (cM)", subtitle="Linkage group") +
      theme_classic()
  } else {
    pl <- ggplot(data = pl.in$lines, aes(x = `Position (cM)`, color = Trait)) +
      facet_grid(.~LG, space = "free") +
      {if(!all(is.na(pl.in$lines$INT)) & pl.in$sup.int) geom_path(data=pl.in$lines, aes(x = INT, y =pl.in$y.dat), colour = "black")} +
      geom_line(data=pl.in$lines, aes(y = LOP, color = Trait), size=pl.in$pl.in$linesize, alpha=0.8, lineend = "round", show.legend = F) +
      scale_x_continuous(breaks=seq(0,max(pl.in$lgs.size),pl.in$cutx)) +
      {if(!all(is.na(pl.in$lines$INT))) geom_point(data=pl.in$points, aes(y = pl.in$y.dat, color = Trait), shape = 2, size = 2, stroke = 1, alpha = 0.8)} +
      scale_y_continuous(breaks=seq(0,max(pl.in$lgs.size, na.rm = T))) +
      {if(nrow(pl.in$thre) > 0) geom_hline(data=pl.in$thre, aes(yintercept=LOP, color=Trait), linetype="dashed", size=.5, alpha=0.8)} +  #threshold
      guides(color = guide_legend("Trait"), fill = guide_legend("Trait"), shape = guide_legend("Trait")) + 
      labs(title=pl.in$main, y = "LOP", x = "Position (cM)", subtitle="Linkage group") +
      theme_classic()
  }
  return(pl)
}

#' Adapted function from QTLpoly
#' 
plot_qtlpoly.effects <- function(x, pheno.col = NULL, p1 = "P1", p2 = "P2", df.info=NULL, lgs = NULL, position = NULL) {
  if(is.null(pheno.col)) {
    pheno.col <- 1:length(x$results)
  } else {
    pheno.col <- which(x$pheno.col %in% pheno.col)
  }

  df.info.sub <- df.info %>% filter(pheno %in% unique(df.info$pheno)[pheno.col])
  group.idx <- which(df.info.sub$Pos %in% position & df.info.sub$LG %in% lgs)
  
  plots2 <- list()
  for(p in pheno.col) {
    nqtl <- length(x$results[[p]]$effects[group.idx])
    if(nqtl > 0) {
      plots1 <- list()
      for(q in group.idx) {
        if(x$ploidy == 4) {
          data <- unlist(x$results[[p]]$effects[[q]])[1:36]
          data <- data.frame(Estimates=as.numeric(data), Alleles=names(data), Parent=c(rep(p1,4),rep(p2,4),rep(p1,14),rep(p2,14)), Effects=c(rep("Additive",8),rep("Digenic",28)))
          data <- data[-c(12:15,18:21,23:30),]
        }
        if(x$ploidy == 6) {
          data <- unlist(x$results[[p]]$effects[[q]])[-c(18:23,28:33,37:42,45:50,52:63,83:88,92:97,100:105,107:133,137:142,145:150,152:178,181:186,188:214,216:278,299:1763)]
          data <- data.frame(Estimates=as.numeric(data), Alleles=names(data), Parent=c(rep(p1,6),rep(p2,6),rep(p1,15),rep(p2,15),rep(p1,20),rep(p2,20)), Effects=c(rep("Additive",12),rep("Digenic",30),rep("Trigenic",40)))
        }
        data$Parent <- factor(data$Parent, levels=unique(data$Parent))
        plot <- ggplot(data[which(data$Effects == "Additive"),], aes(x = Alleles, y = Estimates, fill = Estimates)) +
          geom_bar(stat="identity") +
          scale_fill_gradient2(low = "red", high = "blue", guide = "none") +
          labs(title=names(x$results)[p], subtitle=paste("QTL", q, "\n")) +
          facet_wrap(. ~ Parent, scales="free_x", ncol = 2, strip.position="bottom") +
          # facet_grid(Effects ~ Parent, scales="free_x", space="free_x") +
          theme_minimal() +
          theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5), axis.text.x.bottom = element_text(hjust = 1, vjust = 0.5))
        plots1[[q]] <- plot
      }
      plots2[[p]] <- plots1
    }
  }
  return(plots2)
}
