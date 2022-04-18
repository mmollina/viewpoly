
#' Draws linkage map, parents haplotypes and marker doses
#' Adapted from MAPpoly
#' 
#' @param left.lim covered window in the linkage map start position 
#' @param rigth.lim covered window in the linkage map end position
#' @param ch linkage group ID 
#' @param maps list containing a vector for each linkage group markers with marker positions (named with marker names)
#' @param ph.p1 list containing a data.frame for each group with parent 1 estimated phases. The data.frame contain the columns:
#' 1) Character vector with chromosome ID; 2) Character vector with marker ID;
#' 3 to (ploidy number)*2 columns with each parents haplotypes
#' @param ph.p2 list containing a data.frame for each group with parent 2 estimated phases. See ph.p1 parameter description.
#' @param d.p1 list containing a data.frame for each group with parent 1 dosages. The data.frame contain the columns: 
#' 1) character vector with chromosomes ID; 
#' 2) Character vector with markers ID; 3) Character vector with parent ID; 
#' 4) numerical vector with dosage
#' @param d.p2 list containing a data.frame for each group with parent 2 dosages. See d.p1 parameter description
#' @param snp.names logical TRUE/FALSE. If TRUE it includes the marker names in the plot
#' 
#' @return graphic representing selected section of a linkage group
#' 
#' 
#' @keywords internal
draw_map_shiny<-function(left.lim = 0, right.lim = 5, ch = 1,
                         maps, ph.p1, ph.p2, d.p1, d.p2, snp.names=TRUE)
{
  par <- lines <- points <- axis <- mtext <- text <- NULL
  ch <- as.numeric(ch)
  ploidy <- dim(ph.p1[[1]])[2]
  # if(is.character(ch))
  #   ch <- as.numeric(strsplit(ch, split = " ")[[1]][3])
  var.col <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3")
  names(var.col) <- c("A", "T", "C", "G")
  dark2 <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666")
  d.col<-c(NA, dark2[1:ploidy])
  names(d.col) <- 0:ploidy
  d.col[1]<-NA
  x <- maps[[ch]]
  lab <- names(maps[[ch]])
  zy <- seq(0, 0.5, length.out = ploidy) + 1.5
  pp1 <- ph.p1[[ch]]
  pp2 <- ph.p2[[ch]]
  dp1 <- d.p1[[ch]]
  dp2 <- d.p2[[ch]]
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
  #Parent 1
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
  #Parent 2
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

#' Gets summary information from map.
#' Adapted from MAPpoly
#' 
#' @param left.lim covered window in the linkage map start position 
#' @param rigth.lim covered window in the linkage map end position
#' @param ch linkage group ID 
#' @param maps list containing a vector for each linkage group markers with marker positions (named with marker names)
#' @param d.p1 list containing a data.frame for each group with parent 1 dosages. The data.frame contain the columns: 
#' 1) character vector with chromosomes ID; 
#' 2) Character vector with markers ID; 3) Character vector with parent ID; 
#' 4) numerical vector with dosage
#' @param d.p2 list containing a data.frame for each group with parent 2 dosages. See d.p1 parameter description
#' 
#' @return list with linkage map information: doses;  number snps by group; cM per snp; map size; number of linkage groups
#' 
#' 
#' @keywords internal
map_summary<-function(left.lim = 0, right.lim = 5, ch = 1,
                      maps, d.p1, d.p2){
  ch <- as.numeric(ch)
  # if(is.character(ch))
  #   ch <- as.numeric(strsplit(ch, split = " ")[[1]][3])
  x <- maps[[ch]]
  lab <- names(maps[[ch]])
  ploidy = max(c(d.p1[[ch]], d.p2[[ch]])) 
  d.p1<-d.p2[[ch]]
  d.p2<-d.p1[[ch]]
  x1<-abs(left.lim - x)
  x2<-abs(right.lim - x)
  id.left<-which(x1==min(x1))[1]
  id.right<-rev(which(x2==min(x2)))[1]
  curx<-x[id.left:id.right]
  w<-table(paste(d.p1[id.left:id.right], d.p2[id.left:id.right], sep = "-"))
  M<-matrix(0, nrow = ploidy+1, ncol = ploidy+1, dimnames = list(0:ploidy, 0:ploidy))
  for(i in as.character(0:ploidy))
    for(j in as.character(0:ploidy))
      M[i,j]<-w[paste(i,j,sep = "-")]
  M[is.na(M)]<-0
  return(list(doses = M, 
              number.snps = length(curx), 
              length = diff(range(curx)), 
              cM.per.snp = round(diff(range(curx))/length(curx), 3), 
              full.size = as.numeric(maps[[ch]][length(maps[[ch]])]), 
              number.of.lgs = length(maps)))
}

#' Summary maps - adapted from MAPpoly
#'
#' This function generates a brief summary table 
#' 
#' @param viewmap a list of objects of class \code{viewmap}
#' 
#' @return a data frame containing a brief summary of all maps 
#' 
#' @author Gabriel Gesteira, \email{gabrielgesteira@usp.br}
#' @author Cristiane Taniguti, \email{chtaniguti@tamu.edu}
#' 
#' 
#' @keywords internal
summary_maps = function(viewmap){
  
  simplex <- mapply(function(x,y) {
    sum((x == 1 & y == 0) | (x == 0 & y == 1) |
          (x == max(x) & y == (max(y) -1)) | 
          (x == (max(x) -1) & y == max(y)))
  }, viewmap$d.p1, viewmap$d.p2)
  
  double_simplex <- mapply(function(x,y) {
    sum((x == 1 & y == 1) | (x == 3 & y == 3))
  }, viewmap$d.p1, viewmap$d.p2)
  
  max_gap <- sapply(viewmap$maps, function(x) max(diff(x$l.dist)))
  
  results = data.frame("LG" = as.character(seq(1,length(viewmap$maps),1)),
                       "Genomic sequence" = as.character(unlist(lapply(viewmap$maps, function(x) paste(unique(x$g.chr), collapse = "-")))),
                       "Map length (cM)" = round(sapply(viewmap$maps, function(x) x$l.dist[length(x$l.dist)]),2),
                       "Markers/cM" = round(sapply(viewmap$maps, function(x) length(x$l.dist)/x$l.dist[length(x$l.dist)]),2),
                       "Simplex" = simplex,
                       "Double-simplex" = double_simplex,
                       "Multiplex" = sapply(viewmap$maps, function(x) length(x$mk.names)) - (simplex + double_simplex),
                       "Total" = sapply(viewmap$maps, function(x) length(x$mk.names)),
                       "Max gap" = round(max_gap,2),
                       check.names = FALSE, stringsAsFactors = F)
  results = rbind(results, c('Total', NA, sum(as.numeric(results$`Map length (cM)`)), 
                             round(mean(as.numeric(results$`Markers/cM`)),2), 
                             sum(as.numeric(results$Simplex)), 
                             sum(as.numeric(results$`Double-simplex`)), 
                             sum(as.numeric(results$Multiplex)), 
                             sum(as.numeric(results$Total)), 
                             round(mean(as.numeric(results$`Max gap`)),2)))
  return(results)
}


#' Plot a genetic map - Adapted from MAPpoly
#' 
#' This function plots a genetic linkage map(s) 
#'
#' @param  viewmap object of class \code{viewmap}
#'
#' @param horiz logical. If FALSE, the maps are plotted vertically with the first map to the left. 
#'              If TRUE  (default), the maps are plotted horizontally with the first at the bottom
#'
#' @param col a vector of colors for the bars or bar components (default = 'lightgrey')
#'            \code{ggstyle} produces maps using the default \code{ggplot} color palette 
#'            
#' @param title a title (string) for the maps (default = 'Linkage group')
#'
#' @return A \code{data.frame} object containing the name of the markers and their genetic position
#' 
#' @author Marcelo Mollinari, \email{mmollin@ncsu.edu}
#' @author Cristiane Taniguti, \email{chtaniguti@tamu.edu}
#'
#' @references
#'     Mollinari, M., and Garcia, A.  A. F. (2019) Linkage
#'     analysis and haplotype phasing in experimental autopolyploid
#'     populations with high ploidy level using hidden Markov
#'     models, _G3: Genes, Genomes, Genetics_. 
#'     \doi{10.1534/g3.119.400378}
#'
#' 
#' @keywords internal
plot_map_list <- function(viewmap, horiz = TRUE, col = "ggstyle", title = "Linkage group"){
  axis <- NULL
  map.list <- viewmap$maps
  if(all(col  ==  "ggstyle"))
    col  <- gg_color_hue(length(map.list))
  if(length(col) == 1)
    col <- rep(col, length(map.list))
  z <- NULL
  if(is.null(names(map.list)))
    names(map.list) <- 1:length(map.list)
  max.dist <- max(sapply(map.list, function(x) x$l.dist[length(x$l.dist)]))
  if(horiz){
    plot(0, 
         xlim = c(0, max.dist), 
         ylim = c(0,length(map.list)+1), 
         type = "n", axes = FALSE, 
         xlab = "Map position (cM)", 
         ylab = title)
    axis(1)
    for(i in 1:length(map.list)){
      z <- rbind(z, data.frame(mrk = map.list[[i]]$mk.names, 
                               LG = names(map.list)[i], pos = map.list[[i]]$l.dist))
      plot_one_map(map.list[[i]]$l.dist, i = i, horiz = TRUE, col = col[i])   
    }
    axis(2, at = 1:length(map.list), labels = names(map.list), lwd = 0, las = 2)
  } else{
    plot(0, 
         ylim = c(-max.dist, 0), 
         xlim = c(0,length(map.list)+1), 
         type = "n", axes = FALSE, 
         ylab = "Map position (cM)", 
         xlab = title)
    x <- axis(2, labels = FALSE, lwd = 0)
    axis(2, at = x, labels = abs(x))
    for(i in 1:length(map.list)){
      z <- rbind(z, data.frame(mrk = map.list[[i]]$mk.names, 
                               LG = names(map.list)[i],pos = map.list[[i]]$l.dist))
      plot_one_map(map.list[[i]]$l.dist, i = i, horiz = FALSE, col = col[i])  
    }
    axis(3, at = 1:length(map.list), labels = names(map.list), lwd = 0, las = 2)
  }
  invisible(z)
}


#' Color pallet ggplot-like - Adapted from MAPpoly
#'
#' @param n number of colors
#' 
#' @importFrom grDevices hcl col2rgb hsv rgb2hsv
#' 
#' 
#' @keywords internal
gg_color_hue <- function(n) {
  x <- rgb2hsv(col2rgb("steelblue"))[, 1]
  cols = seq(x[1], x[1] + 1, by = 1/n)
  cols = cols[1:n]
  cols[cols > 1] <- cols[cols > 1] - 1
  return(hsv(cols, x[2], x[3]))
}

#' Plot a single linkage group with no phase - from MAPpoly
#'
#' @param x vector of genetic distances
#' @param i margins size
#' @param horiz logical TRUE/FALSE. If TRUE the map is plotted horizontally.
#' @param col color pallete to be used
#' 
#' 
#' @keywords internal
plot_one_map <- function(x, i = 0, horiz = FALSE, col = "lightgray")
{
  rect <- tail <- lines <- NULL
  if(horiz)
  {
    rect(xleft = x[1], ybottom = i-0.25, 
         xright = tail(x,1), ytop = i+0.25,
         col = col)
    for(j in 1:length(x))
      lines(x = c(x[j], x[j]), y = c(i-0.25, i+0.25), lwd = .5)
  } else {
    x <- -rev(x)
    rect(xleft = i-0.25, ybottom = x[1], 
         xright = i+0.25, ytop = tail(x,1),
         col = col)
    for(j in 1:length(x))
      lines(y = c(x[j], x[j]), x = c(i-0.25, i+0.25), lwd = .5)
  }
}
