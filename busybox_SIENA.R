setwd("~/busybox-mail-networks")
library("network") # A SIENA dependency; use for graph manipulation
library("sna") # Needed only for gplot function
library("RSiena")
files.folder <- "directed" # Enter the folder within that directory that contains graphs
filenames <- list.files(files.folder, pattern="*.graph", full.names=T)
ldf <- lapply(filenames, read.paj)
# ldf <- list(ldf[[8]],ldf[[9]])
vertices <- NULL
for (i in 1:length(ldf)) {
  delete.vertices(ldf[[i]],which(!has.edges(ldf[[i]]))) # Delete isolates
  vertices <- union(vertices,ldf[[i]]%v%"vertex.names") # Collate all non-isolate actors
}
for (i in 1:length(ldf)) {
  original <- ldf[[i]]%v%"vertex.names" # Need to store this for last line in loop
  to.add <- setdiff(vertices,ldf[[i]]%v%"vertex.names") # Remove already-present from total
  add.vertices(ldf[[i]],length(to.add)) # Add blank vertices
  ldf[[i]]%v%"vertex.names" <- c(original,to.add) # Name the blank vertices
}

rm(filenames,files.folder,i,original,to.add) # Get rid of unneeded objects

array <- array(dim=c(length(vertices),length(vertices),length(ldf))) # SIENA takes in this
for (i in 1:length(ldf)) { # Need to order all matrices the same way
  A <- as.sociomatrix(ldf[[i]]) # Extract a matrix with labeled rows and columns
  A <- A[,order(colnames(A))] # Order by columns
  A <- A[order(rownames(A)),] # Order by rows
  array[,,i] <- A # Write to array
}
rm(A)

# Check aggregate graph: Turn array into list of matrices w/plyr, sum over them, 
# then plot the matrix as a graph with sna's gplot, and save the coordinates
coord <- gplot(Reduce("+",plyr::alply(array,3)))

# Some graph checking measures
apply(array,3,isSymmetric) # Note: if all symmetric, cannot measure reciprocity
triad.census <- triad.census(ldf) # Is there transitivity? 
triad.census # Everything above 030T is a triangle
dyad.census <- dyad.census(ldf)
dyad.census
################################################
################# TRIAD CENSUS #################
################################################
census.triads <- rep(list(network.initialize(3)),16) # length(triads)
names(census.triads) <- colnames(triad.census)
census.triads[["012"]][1,2] <- 1
census.triads[["102"]][1,2] <- 1; census.triads[["102"]][2,1] <- 1
census.triads[["021D"]][2,1:3] <-1
census.triads[["021U"]][1:3,2] <- 1
census.triads[["021C"]][1,2] <- 1;census.triads[["021C"]][2,3] <- 1
census.triads[["111D"]][1:2,3] <- 1; census.triads[["111D"]][3,1] <- 1
census.triads[["111U"]][3,1:2] <- 1; census.triads[["111U"]][1,3] <- 1
census.triads[["030T"]][1,2:3] <- 1; census.triads[["030T"]][3,2] <- 1
census.triads[["030C"]][1,2] <- 1; census.triads[["030C"]][2,3] <- 1; census.triads[["030C"]][3,1] <- 1
census.triads[["201"]][1,2:3] <- 1; census.triads[["201"]][2:3,1] <- 1
census.triads[["120D"]][2:3,1] <- 1; census.triads[["120D"]][1:2,3] <- 1
census.triads[["120U"]][1,2:3] <- 1; census.triads[["120U"]][3,1:2] <- 1
census.triads[["120C"]][1,2:3] <- 1; census.triads[["120C"]][2,3] <- 1; census.triads[["120C"]][3,1] <- 1
census.triads[["210"]][1,2:3] <- 1; census.triads[["210"]][3,1:2] <- 1; census.triads[["210"]][2,3] <- 1
census.triads[["300"]][1,2:3] <- 1; census.triads[["300"]][2,c(1,3)] <- 1; census.triads[["300"]][3,1:2] <- 1

# census.triad.mat <- lapply(census.triads,as.matrix)
# census.triad.ig <- lapply(census.triad.mat,igraph::graph.adjacency,mode="directed")
# 
# for (i in 1:length(census.triad.ig)) {
#   igraph::plot.igraph(census.triad.ig[[i]],layout=coords,asp=1)
# }

coords <- t(matrix(c(c(-1/sqrt(3),0),c(0,1),c(1/sqrt(3),0)),nrow=2,ncol=3))

pdf("triad.census.pdf",width=10,height=7.5)
w <- 15 # Additional width past the first column
h <- 19
rpad <- 1 # Subtracts from width on right side
for (j in 1:nrow(triad.census)) {
  layout(matrix(c(1:h,rep(h+1,h*(w-rpad)),rep(h+2,h*rpad)),nrow=h),
         widths=c(rep(1,h),w-rpad,rpad),
         heights=c(rep(1,h),h,h))
  # par(mfrow=c(1,16))
  par(mar=c(0,0,0,0))
  plot.new()
  for (i in 1:length(census.triads)) {
    # pdf(paste0("triad.census/",names(census.triads)[i],".pdf"),width=2.5,height=2.5)
    gplot(census.triads[[i]],
          # main=names(census.triads)[i],
          coord=coords,
          vertex.col=1,
          jitter=F,
          edge.lwd=1,
          vertex.cex=2,
          arrowhead.cex=3)
    # dev.off()
  }
  plot.new()
  plot.new()
  par(mar=c(4.5,0,1.5,0)) # bottom, left, top, right
  dotchart(log10(rev(triad.census[j,])),xlim=c(.3,log10(10^8.1)),xaxt='n')
  orders <- 9
  seq <- rep(c(1,2,5),orders)*rep(10^(0:(orders-1)),each=3)
  abline(v=log10(seq[-1]),col="lightgray",lty=2)
  plot.new()
}
dev.off()
par(mar=c(5.1,4.1,4.1,2.1))
par(mfrow=c(1,1))

################################################
################################################
################################################

l.out <- lapply(plyr::alply(array,3),rowSums); l.out # Check outdegrees
l.ind <- lapply(plyr::alply(array,3),colSums); l.ind # Check indegrees
l.tot <- l.ind
for (i in 1:length(l.ind)) {
  l.tot[[i]] <- l.ind[[i]] + l.out[[i]]
}
l.tot # Check total degree
total.deg <- Reduce("+",l.tot); total.deg
plot(table(total.deg))
l.act <- lapply(l.out,">",0); l.act # Active in this timestep?
l.red <- Reduce("+",l.act); l.red # How many times active?
table(l.red) # How many nodes were active in only one time step, 2 steps, etc.?
l.iso <- lapply(l.tot,">",0); l.iso # Isolate in this timestep?
l.ext <- Reduce("+",l.iso); l.ext # How many times not an isolate? Should be no zeros.
table(l.ext) # How many nodes were non-isolates in how many timesteps?


# Rough plot of the graph evolution
cumColSums <- apply(simplify2array(l.act),1,cumsum) # List to matrix to cumulative colSums
# Instead of simplify2array, could do do.call&rbind, sapply&unlist, or matrix&unlist
colfun <- colorRampPalette(c("#f7fcf0","#084081")) # Create a gradient
cols <- c("#FFFFFF",colfun(max(cumColSums))) # Never-active nodes will be white
pdf("network-growth.pdf",width=10,height=7.5)
for (i in 1:length(ldf)) {
  gplot(array[,,i],coord=coord,vertex.col=cols[cumColSums[i,]+1])
}
dev.off()

rm(coord,l.out,l.ind,l.tot,l.act,l.iso,l.red,l.ext,cumColSums,colfun,cols,i,j,h,w,orders,rpad,seq)

################################################################################################
################################################################################################
################################################################################################

# BEGIN SIENA
communication <- sienaDependent(array,sparse=F)
highdegree <- coCovar(1*(total.deg>150))

data <- sienaDataCreate(communication,highdegree)
effects <- getEffects(data)
effects <- includeEffects(effects,density,inPop,reciprocity,transTrip,outRateLog)

# effects <- includeEffects(effects,density,inPop,outPop,reciprocity,transTrip)
# effects <- includeEffects(effects,density,inPop,inPopSqrt,outPop,outPopSqrt,reciprocity,transTrip)
# myeffects <- includeTimeDummy(effects,density,reciprocity,timeDummy="5,8")

algo <- sienaAlgorithmCreate(projname='network-growth',cond=FALSE,firstg=.02)
# out <- siena07(algo,data=data,effects=effects)
# out
out1 <- siena07(algo,data=data,effects=effects,verbose=T)
out1
out2 <- siena07(algo,data=data,effects=effects,prevAns=out,returnDeps=T)
out2

# x <- 1:200
# plot(x,-0.1412*x+0.5130*sqrt(x))
# plot(x,0.0787*x-0.1939*sqrt(x))