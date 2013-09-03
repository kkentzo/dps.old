


## ========================================================================================
## cleans up the supplied argument by settings all NaNs and Infs to NA
clean <- function(x, fill.value=NA, clean.na=F, clean.nan=T, clean.inf=T) {

  if (clean.na)
    x[is.na(x)] <- fill.value

  if (clean.nan)
    x[is.nan(as.matrix(x))] <- fill.value

  if (clean.inf)
    x[abs(x) == Inf] <- fill.value

  x
  
}

## ========================================================================================



## calculate and return the moving average using the specified window
## using the filter function
## if N==0, it will just return X
ma <- function(x, window.size=100, sides=2) {
  
  if (window.size > 0) 
    x <- filter(x, rep(1/window.size, window.size), sides=sides)

  x

}





## calculate and return the moving average using the specified window
## by ignores NA values
## if N==0, it will just return X
ma.old <- function(x, window.size) {

  if (window.size > 0)
    x <- rollapply(x, window.size, function(x) mean(x, na.rm=T))

  x
  
}




## manual -- this is slow for lengthy X
ma.alt <- function( x, n=100, sides=2 ) {

  if (n > 0) {

    n <- n - n %% 2
    left <-n %/% 2
    right <- n - left

    x <- sapply(1:length(x), function(i) {
      lo <- i-left;
      hi <- i+right;
      mean(x[ifelse(lo<0,0,lo):ifelse(hi>length(x),length(x),hi)], na.rm=T)
    })
    
  }

  x
  
}




## ========================================================================================
## plots multiple variables in Y (represented by columns in a matrix or data.frame)
## against X 
## if Y is NA, then X becomes Y and X is the index of Y
mplot <- function( x, y=NA, type="l", ltype=1, col=c("blue", "red", "green"), main="",
                   ylab="", xlab="", ylim=NA, xlim=NA, log.take="",
                   xaxt='s', yaxt='s', ... )
{

  if (length(y) == 1 && is.na(y)) {

    y <- as.matrix(x)
    x <- as.matrix(1:nrow(y))
    
  } else {

    y <- as.matrix(y)
    x <- as.matrix(x)
    
  }

  y <- clean(y)

  ## ## filter out NaNs and infinities
  ## y[which(is.nan(y))] <- NA
  ## y[is.infinite(abs(y))] <- NA

  ## turn log into a table of side values
  log.take <- strsplit(log.take, split="")[[1]]
  

  if ('x' %in% log.take) {
    x <- log10(x)
    xaxt <- 'n'
  } 

  if ('y' %in% log.take) {
    y <- log10(y)
    yaxt <- 'n'
  } 

  x <- clean(x)
  y <- clean(y)

  if (length(ylim) != 2) 
    ylim <- range(y, na.rm=TRUE)

  if (length(xlim) != 2)
    xlim <- range(x, na.rm=TRUE)

  ## format plot
  plot(c(0,0), col="white", xlim=xlim, ylim=ylim, 
       main=main, ylab=ylab, xlab=xlab, xaxt=xaxt, yaxt=yaxt, ...)

  
  ## and now plot all lines
  for (i in 1:ncol(y)) {

    ## recycle colors??
    if (length(col) > 1) {
      color <- col[i]
      if (is.na(color))
        color <- col[((i-1) %% length(col)) + 1]
    } else {
      color <- col
    }

    ## recycle types??
    if (length(ltype) > 1) {
      tp <- ltype[i]
      if (is.na(tp))
        tp <- ltype[((i-1) %% length(ltype)) + 1]
    } else {
      tp <- ltype
    }

    if (ncol(x) > 1) 
      xx <- x[,i]
    else 
      xx <- x

    lines(xx, y[,i], col=color, type=type, lty=tp)

  }

  ## collect ... args
  args <- list(...)
  if (length(args[["cex.axis"]]))
    cex.axis <- args$cex.axis
  else
    cex.axis <- 1
  
  decorate.log(paste(log.take, collapse=""), cex=cex.axis)
  
}




## ========================================================================================
## plot the requested error bars
error.bar <- function(x, y, upper, lower=upper, length=0.1, ...) {

  arrows(x, y+upper, x, y-lower, angle=90, code=3, length=length, ...)
  
}


## ========================================================================================
## plot the variables in VAR.MEAN (represented by columns)
## the standard deviation (VAR.STD) is plotted by grey polygons around the actual curves
plot.with.range <- function( var.mean, var.std, x=numeric(0),
                             col=c("blue", "red", "green"),
                             ylab="", xlab="Time", ylim=NA, type="l",
                             log.take="", border=F, alpha.f=0.5,
                             xaxt='s', yaxt='s', ...)
{

  var.mean <- clean(as.matrix(var.mean))
  var.std <- clean(as.matrix(var.std))


  ## first, convert NaNs in the matrices into NAs
  ##var.mean[is.nan(var.mean)] <- NA
  ##var.std[is.nan(var.std)] <- NA

  ## convert all NAs in std to 0s
  var.std[is.na(var.std)] <- 0

  ## figure out number of points
  ## and create a matrix (for filtering out NaNs)
  if (length(x) == 0)
    x <- 1:nrow(var.mean)

  ## turn log into a table of side values
  log.take <- strsplit(log.take, split="")[[1]]

  if ('x' %in% log.take) {
    x <- log10(x)
    xaxt <- 'n'
  } 

  xx <- c(x, rev(x))

  lo <- var.mean - var.std
  hi <- var.mean + var.std

  yy <- rbind(lo, as.matrix(hi[nrow(hi):1,]))

  if ('y' %in% log.take) {
    var.mean <- log10(var.mean)
    yy <- log10(yy)
    yaxt <- 'n'
  } 

  xx <- clean(xx)
  yy <- clean(yy)

  if (length(ylim) != 2) 
    ylim <- range(yy, na.rm=TRUE)

  if (abs(ylim[1]) == Inf || abs(ylim[2]) == Inf) {
    print("plot.with.range : Infinities in ylim. Skipping plotting.")
    return()
  }

  plot(c(0,0), col="white", xlim=range(x, na.rm=TRUE),
       ylim=ylim, ylab=ylab, xlab=xlab, xaxt=xaxt, yaxt=yaxt, ...)

  ## decorate axes logarithmically??
  if (length(log.take) > 0) {
    ## collect ... args
    args <- list(...)
    if (length(args[["cex.axis"]]))
      cex.axis <- args$cex.axis
    else
      cex.axis <- 1
    decorate.log(paste(log.take, collapse=""), cex=cex.axis)
  }

  ## first plot the sd regions (the polygons)
  for (i in 1:ncol(var.mean)) {

    if (border) {
      color <- col[i]
      if (is.na(color))
        color <- col[((i-1) %% length(col)) + 1]
      
      polygon.border <- adjustcolor(color, alpha.f=0.5)
    } else {
      polygon.border <- NA
    }

    ## plot the data -- **filtered** so as to exclude the NAs in var.mean
    polygon(xx[! is.na(yy[,i])], yy[! is.na(yy[,i]),i], border=polygon.border,
            lty="dashed", col=adjustcolor("gray", alpha.f=alpha.f))
    
  }

  ## now plot the lines
  for (i in 1:ncol(var.mean)) {
    color <- col[i]
    if (is.na(color))
      color <- col[i %% length(col)]

    lines(x[! is.na(yy[,i])], var.mean[,i][! is.na(yy[,i])], col=color, type=type)

  }
  
}




## =============================================================================

## Will draw the axes with proper logarithmic ticks
## REQUIREMENTS :
## 1/ You should **NOT** use the log="x|y" option in plot()
##    instead, you should plot the log10(data) in the respective axis
## 2/ You should use xaxt='n' and/or yaxt='n' in plot()
decorate.log <- function( log.take="xy", cex=1 ) {

  ## turn log into a table of side values
  log.take <- strsplit(log.take, split="")[[1]]

  ## decorate axes 
  for (s in log.take) {

    side <- switch(s, x=1, y=2, NA)

    if (! is.na(side)) {

      ## get range of ticks
      ticks.range <- range(axTicks(side))
      ## compute major ticks from range
      major.ticks <- seq(floor(ticks.range[1]), ceiling(ticks.range[2]), 1)

      ## compute intermediate minor ticks 
      minor.ticks <- Reduce(c, sapply(major.ticks[-length(major.ticks)],
                                      function(x) seq(10^x, 10^(x+1), length.out=10)[-10]),
                            NULL)
      ## add the last item (from the range) and take its logarithm as well
      minor.ticks <- log10(c(minor.ticks, 10^major.ticks[length(major.ticks)]))

      ## form labels at major ticks 
      labels <- sapply(major.ticks, function(i) as.expression(bquote(10 ^ .(i))))
      ## draw axis
      axis(side, at=major.ticks, labels=labels, cex.axis=cex)
      axis(side, lwd=0, lwd.ticks=1, at=minor.ticks,
           labels=rep("", length(minor.ticks)), tcl=-0.25)
      
    }
    
  }
  
}





## ========================================================================================
## detect and return the separate clusters (islands) in a 2-d matrix
## clusters are represented as (nx2) matrix elements of a list
## where rows in each matrix indicate entries belonging to the same cluster
## (based on proximity) and columns denote (row,col) indices into the Z matrix
## separation is based on zero elements (i.e. clusters are the non-zero islands in matrix Z)
detect.clusters <- function( z ) {


  ## figure out the indices of non-zero matrix elements
  ##z.pos <- matrix(which(z>0, arr.ind=T), byrow=F, ncol=2)
  ##z.pos <- cbind(z.pos, rep(0, nrow(z.pos)))

  

  ## this function starts from a given POINT and assembles the surrounding cluster
  ## it updates the Z.CIDS matrix so as to mark the cluster points with the CID (cluster id)
  assemble.cluster <- function(cid, pos, z.cids) {

    ## convert position tuple to matrix
    pos <- matrix(pos, ncol=2)

    ## create structure with position and the cid matrix
    res <- list(cluster=pos,
                z.cids=z.cids)

    ## mark current position as belonging to cluster CID
    res$z.cids[pos] <- cid

    ## the following conditions make provisions for avoiding matrix edges

    ## go up??
    if (pos[1] > 1 && z[pos + c(-1,0)] > 0 && res$z.cids[pos + c(-1,0)] == 0) {
      new.res <- assemble.cluster(cid, pos + c(-1,0), res$z.cids)
      res$cluster <- rbind(res$cluster, new.res$cluster)
      res$z.cids <- res$z.cids + new.res$z.cids
    }
    ## go down??
    if (pos[1] < nrow(z) && z[pos + c(1,0)] > 0 && res$z.cids[pos + c(1,0)] == 0) {
      new.res <- assemble.cluster(cid, pos + c(1,0), res$z.cids)
      res$cluster <- rbind(res$cluster, new.res$cluster)
      res$z.cids <- res$z.cids + new.res$z.cids
    }
    ## go right??
    if (pos[2] < ncol(z) && z[pos + c(0, 1)] > 0 && res$z.cids[pos + c(0,1)] == 0) {
      new.res <- assemble.cluster(cid, pos + c(0,1), res$z.cids)
      res$cluster <- rbind(res$cluster, new.res$cluster)
      res$z.cids <- res$z.cids + new.res$z.cids
    }
    ## go left??
    if (pos[2] > 1 && z[pos + c(0, -1)] > 0 && res$z.cids[pos + c(0,-1)] == 0) {
      new.res <- assemble.cluster(cid, pos + c(0,-1), res$z.cids)
      res$cluster <- rbind(res$cluster, new.res$cluster)
      res$z.cids <- res$z.cids + new.res$z.cids
    }

    ## return result
    res
    
  }
  

  ## start with no clusters
  clusters <- list()

  ## build cluster id matrix
  z.cids <- matrix(rep(0, prod(dim(z))), ncol=ncol(z))

  i.cluster <- 1

  for (i.row in 1:nrow(z))
    for (i.col in 1:ncol(z))
      if (z[i.row,i.col] > 0 && z.cids[i.row,i.col] == 0) {
        res <- assemble.cluster(i.cluster, c(i.row,i.col), z.cids)
        clusters[[i.cluster]] <- res$cluster
        z.cids <- res$z.cids
        i.cluster <- i.cluster + 1
      }
                                                  
  

  clusters

}




## =============================================================================
downsample <- function( x, freq, fun=mean ) {

  x <- as.data.frame(x)
  
  n.groups <- ceiling(nrow(x) / freq)

  group.indices <- vector()

  for (i in 1:n.groups)
    group.indices <- c(group.indices, rep(i, freq))

  r <- aggregate(x, by=list(group.indices[1:nrow(x)]), fun)

  r
  
}
