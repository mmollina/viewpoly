require(mappoly)
# Loading final map and conditional probabilities
# load("~/repos/BT_map/src/final_map/final_map_and_probs_7_oct_2019.rda")
#source("~/repos/BT_map/src/genome_assisted_ordering/auxiliary_function.R")

# extract_maps<-function(x, dat)
# {
#   w<-round(cumsum(c(0,imf_h(x$maps[[1]]$seq.rf))),1)
#   names(w)<-dat$mrk.names[x$maps[[1]]$seq.num]
#   w
# }
#Variants
#####
# Variants<-NULL
# for(ch in 1:15){
#   cat("chromosome: ", ch, "\n")
#   mappoly.map<-maps.final[[ch]]
#   map<-extract_maps(x = mappoly.map, dat = BT.trifida.triloba.298)
#   Mbr<-ph_list_to_matrix(mappoly.map$maps[[1]]$seq.ph$P, 6)
#   Mtz<-ph_list_to_matrix(mappoly.map$maps[[1]]$seq.ph$Q, 6)
#   dimnames(Mbr)<-list(names(map), letters[1:6])
#   dimnames(Mtz)<-list(names(map), letters[7:12])
#   x<-get_genomeic_position_2(map = map, ref = "Tf", ch = ch)
#   y1<-rep("Tf", length(x))
#   y1[is.na(x)]<-"Tl"
#   y2<-x
#   y2[y1=="Tl"]<-get_genomeic_position_2(map = map, ref = "Tl", ch = ch)[y1=="Tl"]
#   Geno.pos<-data.frame(ref = y1, pos = y2)
#   load(paste0("~/repos/BT_map/src/genotype_calling/trifida/BT_raw_trifida_", ch,".RData"))
#   Atf<-raw.counts[sapply(raw.counts, function(x) !is.null(x))]
#   Btf<-t(sapply(Atf, function(x) colnames(x$par)))
#   load(paste0("~/repos/BT_map/src/genotype_calling/triloba/BT_raw_triloba_", ch,".RData"))
#   Atl<-raw.counts[sapply(raw.counts, function(x) !is.null(x))]
#   Btl<-t(sapply(Atl, function(x) colnames(x$par)))
#   rownames(Btf)<-sapply(strsplit(rownames(Btf), "_"), function(x) x[2])
#   rownames(Btl)<-sapply(strsplit(rownames(Btl), "_"), function(x) x[2])
#   w<-NULL
#   for(i in 1:nrow(Geno.pos))
#   {
#     if(Geno.pos[i,1] == "Tf")
#       w<-rbind(w,Btf[Geno.pos[i, 2] == rownames(Btf)])
#     if(Geno.pos[i,1] == "Tl"){
#       if(sum(Geno.pos[i, 2] == rownames(Btl)) > 1)
#         stop()
#       w<-rbind(w,Btl[Geno.pos[i, 2] == rownames(Btl)])
#     }
#   }
#   Variants<-rbind(Variants, data.frame(cbind(Geno.pos,w), ch = ch))
# }
# save(Variants, file = "~/repos/BT_map/map_shiny/variants.rda")
#####
#load(file = "~/repos/BT_map/map_shiny/variants.rda")

#Maps
#####
# maps<-vector("list", 15)
# for(ch in 1:15){
#   maps[[ch]]<-extract_maps(x = maps.final[[ch]], dat = BT.trifida.triloba.298)
# }
# save(maps, file = "~/repos/BT_map/map_shiny/maps.rda")
#####
#load(file = "~/repos/BT_map/map_shiny/maps.rda")

#Phases
#####
# var.col <- c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3")
# names(var.col) <- c("A", "T", "C", "G")
# ph.tz<-ph.br<-vector("list", 15)
# for(ch in 1:15){
#   M<-ph_list_to_matrix(maps.final[[ch]]$maps[[1]]$seq.ph$P, 6)
#   dimnames(M)<-list(names(maps[[ch]]), letters[1:6])
#   M <- 2 - M
#   var.temp<-NULL
#   for(i in rownames(M))
#     var.temp <- rbind(var.temp, as.character(as.matrix(Variants[i,3:4][M[i,]])))
#   dimnames(var.temp)<-list(names(maps[[ch]]), letters[1:6])
#   ph.br[[ch]]<-var.temp[,6:1]
#   M<-ph_list_to_matrix(maps.final[[ch]]$maps[[1]]$seq.ph$Q, 6)
#   dimnames(M)<-list(names(maps[[ch]]), letters[7:12])
#   M <- 2 - M
#   var.temp<-NULL
#   for(i in rownames(M))
#     var.temp <- rbind(var.temp, as.character(as.matrix(Variants[i,3:4][M[i,]])))
#   dimnames(var.temp)<-list(names(maps[[ch]]), letters[7:12])
#   ph.tz[[ch]]<-var.temp[,6:1]
# }
# save(ph.br, ph.tz, file = "~/repos/BT_map/map_shiny/phs.rda")
#####
#load(file = "~/repos/BT_map/map_shiny/phs.rda")

# doses
#####
# d.tz<-d.br<-vector("list", 15)
# for(ch in 1:15){
#   d.br[[ch]]<-BT.trifida.triloba.298$dosage.p[match(names(maps[[ch]]), BT.trifida.triloba.298$mrk.names)]
#   d.tz[[ch]]<-BT.trifida.triloba.298$dosage.q[match(names(maps[[ch]]), BT.trifida.triloba.298$mrk.names)]
# }
# save(d.br, d.tz, file = "~/repos/BT_map/map_shiny/ds.rda")
#####
#load(file = "~/repos/BT_map/map_shiny/ds.rda")

draw_map_shiny<-function(left.lim = 0, right.lim = 5, ch = 1,
                         maps, ph.br, ph.tz, d.tz, d.br, snp.names=TRUE)
{
  if(is.character(ch))
    ch <- as.numeric(strsplit(ch, split = " ")[[1]][3])
  var.col <- RColorBrewer::brewer.pal(n = 4, name = "Set1")
  names(var.col) <- c("A", "T", "C", "G")
  d.col<-c(NA, RColorBrewer::brewer.pal(n = 6, name = "Dark2"))
  names(d.col) <- 0:6
  d.col[1]<-NA
  x <- maps[[ch]]
  lab <- names(maps[[ch]])
  zy <- seq(0, 0.5, length.out = 6) + 1.5
  pbr <- ph.br[[ch]]
  ptz <- ph.tz[[ch]]
  dbr<-d.br[[ch]]
  dtz<-d.tz[[ch]]
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
  for(i in 1:6)
  {
    lines(c(x[id.left], x[id.right]), c(zy[i], zy[i]), lwd=10, col = "gray")
    points(x = seq(x[id.left], x[id.right], length.out = length(curx)),
           y = rep(zy[i], length(curx)),
           col = var.col[ptz[id.left:id.right,i]],
           pch = 15,
           cex = 2)
  }
  mtext(text = "Tanzania", side = 2, at = mean(zy), line = -3, font = 4)
  for(i in 1:6)
    mtext(letters[12:7][i], at = zy[i], side = 2,  line = -4, font = 1)
  connect.lines<-seq(x[id.left], x[id.right], length.out = length(curx))
  for(i in 1:length(connect.lines))
    lines(c(curx[i], connect.lines[i]), c(0.575, zy[1]-.05), lwd=0.3)
  points(x = seq(x[id.left], x[id.right], length.out = length(curx)),
         y = zy[6]+0.05+dtz[id.left:id.right]/20,
         col = d.col[as.character(dtz[id.left:id.right])],
         pch = 19, cex = .7)
  corners = par("usr") 
  par(xpd = TRUE) 
  text(x = corners[1]+.5, y = mean(zy[6]+0.05+(1:6/20)), "Doses")
  #Parent BEAUREGARD
  zy<-zy+1.1
  for(i in 1:6)
  {
    lines(c(x[id.left], x[id.right]), c(zy[i], zy[i]), lwd=10, col = "gray")
    points(x = seq(x[id.left], x[id.right], length.out = length(curx)),
           y = rep(zy[i], length(curx)),
           col = var.col[pbr[id.left:id.right,i]],
           pch = 15,
           cex = 2)
  }
  mtext(text = "Beauregard", side = 2, at = mean(zy), line = -3, font = 4)
  points(x = seq(x[id.left], x[id.right], length.out = length(curx)),
         y = zy[6]+0.05+dbr[id.left:id.right]/20,
         col = d.col[as.character(dbr[id.left:id.right])],
         pch = 19, cex = .7)
  corners = par("usr") 
  par(xpd = TRUE) 
  text(x = corners[1]+.5, y = mean(zy[6]+0.05+(1:6/20)), "Doses")
  if(snp.names)
    text(x = seq(x[id.left], x[id.right], length.out = length(curx)),
         y = rep(zy[6]+0.05+.3, length(curx)),
         labels = names(curx),
         srt=90, adj = 0, cex = .7)
  for(i in 1:6)
    mtext(letters[6:1][i], at = zy[i], side = 2,  line = -4, font = 1)
  # legend("bottomleft", legend=c("A", "T", "C", "G", "-"), 
  #        fill =c(var.col, "white"), 
  #        box.lty=0, bg="transparent")
}

map_summary<-function(left.lim = 0, right.lim = 5, ch = 1,
                      maps, d.tz, d.br){
  if(is.character(ch))
    ch <- as.numeric(strsplit(ch, split = " ")[[1]][3])
  x <- maps[[ch]]
  lab <- names(maps[[ch]])
  dbr<-d.br[[ch]]
  dtz<-d.tz[[ch]]
  x1<-abs(left.lim - x)
  x2<-abs(right.lim - x)
  id.left<-which(x1==min(x1))[1]
  id.right<-rev(which(x2==min(x2)))[1]
  curx<-x[id.left:id.right]
  w<-table(paste(dbr[id.left:id.right], dtz[id.left:id.right], sep = "-"))
  M<-matrix(0, nrow = 7, ncol = 7, dimnames = list(0:6, 0:6))
  for(i in as.character(0:6))
    for(j in as.character(0:6))
      M[i,j]<-w[paste(i,j,sep = "-")]
  M[is.na(M)]<-0
  return(list(doses = M, number.snps = length(curx), length = diff(range(curx)), cM.per.snp = round(diff(range(curx))/length(curx), 3)))
}
#row: beauregard; column: tanzania
#map_summary(left.lim = 0, right.lim = 10, ch = 1, maps = maps, d.tz = d.tz, d.br = d.br)
  
# a<-sapply(maps, max)
# b<-sapply(maps, length)
# d<-cbind(b,a,round(a/b,2))
# d<-rbind(d, apply(d, 2, function(x) round(sum(x),1)))
# d[16,3]<-round(d[16,2]/d[16,1],2)
# 
# read.csv(file = "~/repos/BT_map/src/genome_assisted_ordering/csv_subdivisions/ch_1_subdivision.csv", row.names = 1)

#draw_map_shiny(left.lim = 0, right.lim = 7, ch = 1, d.tz = d.tz, d.br = d.br,maps = maps, ph.br = ph.br, ph.tz = ph.tz)
