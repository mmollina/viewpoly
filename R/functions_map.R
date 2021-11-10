draw_map_shiny<-function(left.lim = 0, right.lim = 5, ch = 1,
                         maps, ph.p1, ph.p2, d.p1, d.p2, snp.names=TRUE)
{
  ch <- as.numeric(ch)
  ploidy <- dim(ph.p1[[1]])[2]
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
  pp1 <- ph.p1[[ch]]
  pp2 <- ph.p2[[ch]]
  d.p1<-d.p2[[ch]]
  d.p2<-d.p1[[ch]]
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
         y = zy[ploidy]+0.05+d.p2[id.left:id.right]/20,
         col = d.col[as.character(d.p2[id.left:id.right])],
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
         y = zy[ploidy]+0.05+d.p1[id.left:id.right]/20,
         col = d.col[as.character(d.p1[id.left:id.right])],
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
  return(list(doses = M, number.snps = length(curx), length = diff(range(curx)), cM.per.snp = round(diff(range(curx))/length(curx), 3), full.size = as.numeric(maps[[ch]][length(maps[[ch]])]), number.of.lgs = length(maps)))
}



