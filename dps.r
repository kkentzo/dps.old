

## needs the following libraries 
library(tseries)
library(hdf5)
library(gplots)
library(zoo)


## source some helper functions
source("rutils.r")



## ============================================================================
##                     LOAD HDF FILES
## ============================================================================



## load the data from file FNAME
dps.load <- function( fname ) {

  print(sprintf("Loading file %s", fname))

  results <- hdf5load(fname, load=F, tidy=T)

  results
                         

}








## ============================================================================
##                     PLOT RESULTS (DYNAMICS/HISTS)
## ============================================================================



## ======================================================================
## calculates and returns the statistical association of the pair
## of variables specified by NAMES
## DYNAMICS should be either the global, inter or intra data frames
## STEPS.RANGE should be a 2-tuple
## NORM.BY specifies how the covariance should be normalized
## options include "" (covariance), "x" or "y" (regression coeffs),
## or "xy" (correlation coefficient)
dps.calc.association <- function( dynamics, names, norm.by="xy", plot=F, steps.range=NA ) {

  if (length(names) == 1)
    names <- strsplit(names, "\\.")[[1]]

  if (length(steps.range) != 2) 
    steps.range <- 1:nrow(dynamics$C)
  else
    steps.range <- steps.range[1]:steps.range[2]

  ## calculate normalization factor
  norm.denom <- Reduce(`*`,
                       lapply(strsplit(norm.by, "")[[1]],
                              function (token) switch(token,
                                                      x=dynamics$V[steps.range,names[1]],
                                                      y=dynamics$V[steps.range,names[2]])),
                       1)


  assoc <- dynamics$C[steps.range,paste(names, collapse=".")] / sqrt(norm.denom)
  ##sqrt(dynamics$V[,names[1]] * dynamics$V[,names[2]])

  assoc[abs(assoc)==Inf] <- NA

  if (plot) {

    rng <- range(assoc, na.rm=T)

    x <- seq(from=rng[1], to=rng[2], by=diff(rng) / 100)
    y <- density(assoc, na.rm=T)
    ##y <- dnorm(x, mean=mean(assoc, na.rm=T), sd=sd(assoc, na.rm=T))

    h <- hist(assoc, plot=F)

    ## store settings
    par.bak <- par(no.readonly=TRUE)
    par(mfrow=c(2,1), cex.lab=1.2)

    ##title <- expression(paste("Histogram of ", rho, "(", bar(beta), ",", f, ")"))


    ## plot assoc per time frame
    breaks <- 5

    indices <- cut(steps.range, breaks=breaks, include.lowest=T, labels=F)


    colors.hist <- colorpanel(breaks, low="blue", mid="green", high="red")

    hists <- lapply(1:breaks, function (i) hist(assoc[which(indices==i)], plot=F))

    y.min <- min(sapply(hists, function (h) min(h$counts, na.rm=T)))
    y.max <- max(sapply(hists, function (h) max(h$counts, na.rm=T)))
    x.min <- min(sapply(hists, function (h) min(h$mids, na.rm=T)))
    x.max <- max(sapply(hists, function (h) max(h$mids, na.rm=T)))

    ## plot histogram
    hist(assoc, freq=F, ylim=range(h$density, y$y), xlim=c(x.min, x.max),
         xlab="Correlation Coefficient")
    ## plot density estimation
    lines(y$x, y$y, col="blue")

    means <- lapply(1:breaks, function(i) mean(assoc[which(indices==i)], na.rm=T))
    sds <- lapply(1:breaks, function(i) sd(assoc[which(indices==i)], na.rm=T))

    ## plot time-progressive histograms
    for (i in 1:breaks) {

      h <- hists[[i]]

      print(sprintf("mean=%.3f | sd=%.3f (break=%d)", means[[i]], sds[[i]], i))

      if (i==1) {

        plot(h$mids, h$counts, t="l", lwd=2, col=colors.hist[i],
             xlim=c(x.min, x.max), ylim=c(y.min, y.max))
        abline(v=0)
        
      } else {

        lines(h$mids, h$counts, lwd=2, col=colors.hist[i])
        
      }
    }

    legend("topleft", c("start", "middle", "end"),
           lwd=2, col=c("blue", "green", "red"))
    

    ## === QQ PLOT ===
    ## qqnorm(assoc, col="blue")
    ## qqline(assoc)

    ## restore settings
    par(par.bak)
    
  }

  assoc

}





## ==========================================================================
## calculates and returns the 3 components of the price equation
## for beta, kappa and alpha
## STEPS.RANGE is a 2-tuple
dps.calc.price <- function(results, window.size=500,
                           steps.range=NA, decompose=F)
{

  if (length(steps.range) != 2) 
    take.steps <- 1:length(results$dynamics$global$M$fitness)
  else
    take.steps <- steps.range[1]:steps.range[2]
  
  dz <- list(fitness=ma(results$dynamics$global$M$fitness[take.steps],
               window.size=window.size),
             window.size=window.size)

  for (z.name in c("beta", "kappa", "alpha")) {

    ## calculate and store transmission bias
    dz[[z.name]] <- data.frame(
                      tbias=ma(results$dynamics$global$M[[sprintf("t%s", z.name)]][take.steps],
                        window.size=window.size) / dz$fitness)

    ## calculate covariances
    for (level in c("global", "inter", "intra")) {

      dz[[z.name]][[sprintf("%s", level)]] <-
        ma(dps.calc.association(results$dynamics[[level]],
                                sprintf("%s.fitness", z.name),
                                norm.by="",
                                steps.range=steps.range),
           window.size=window.size) / dz$fitness

      if (decompose) {

        for (name in c("nr", "ht", "death")) {

            dz[[z.name]][[sprintf("%s.%s", level, name)]] <-
              ma(dps.calc.association(results$dynamics[[level]],
                                      sprintf("%s.%s", z.name, name),
                                      norm.by="",
                                      steps.range=steps.range),
                 window.size=window.size) / dz$fitness
        }
      }
    }

    ## Calculate total selection
    dz[[z.name]][["total"]] <- dz[[z.name]]$global + dz[[z.name]]$tbias
    
  }

  
  dz
  
}




## =========================================================================
## plot all the Price components (SORTED BY SELECTION LEVEL)
## DZ can be either the output of dps.calc.price() or a results objects
## use FITNESS.POSTFIX="1" to use fitness1 and "2" to use fitness2
dps.plot.price <- function( results, dz=NULL, window.size=0, steps.range=NA,
                            fitness.postfix="", xlim=NA, bw.adjust=1 ) {

  par.bak <- par(no.readonly=TRUE)
  layout(rbind(rep(1,4), matrix(2:9, ncol=4)))

  ## calculate price components??
  if (is.null(dz))
    dz <- dps.calc.price(results, window.size, steps.range=NA)

  ## calculate steps.range
  if (length(steps.range) != 2) 
    steps.range <- 1:nrow(dz$beta)
  else
    steps.range <- steps.range[1]:steps.range[2]

  ## plot evolutionary dynamics of beta, kappa
  plot.with.range(cbind(results$dynamics$global$M$beta[steps.range],
                        results$dynamics$global$M$kappa[steps.range],
                        results$dynamics$global$M$alpha[steps.range]),
                  cbind(sqrt(results$dynamics$global$V$beta[steps.range]),
                        sqrt(results$dynamics$global$V$kappa[steps.range]),
                        sqrt(results$dynamics$global$V$alpha[steps.range])),
                  x=steps.range,
                  main="Plasmid Replication Parameters",
                  col=c("blue", "red", "green"),
                  xlab="Time", ylab="")
  legend("topleft", c(expression(beta), expression(kappa), expression(alpha)),
         lwd=1, col=c("blue", "red", "green"))


  ## PLOT SELECTION (TOTAL, WITHIN AND BETWEEN)
  for (level in c("total", "inter", "intra", "tbias")) {

    ## plot time series
    b <- dz$beta[[level]][steps.range]
    k <- dz$kappa[[level]][steps.range]
    a <- dz$alpha[[level]][steps.range]

    print(sprintf("%s | beta=%.3e | kappa=%.3e | alpha=%.3e",
                  level, mean(b,na.rm=T),mean(k,na.rm=T), mean(a,na.rm=T)))

    mplot(steps.range, cbind(b, k, a), xlab="Time",
          main=sprintf("%s%s", level, fitness.postfix))
    abline(h=0)
          

    ## plot densities
    b <- density(b, adjust=bw.adjust, na.rm=T)
    k <- density(k, adjust=bw.adjust, na.rm=T)
    a <- density(a, adjust=bw.adjust, na.rm=T)

    mplot(cbind(b$x, k$x, a$x), cbind(b$y, k$y, a$y),
          main=sprintf("%s%s", level, fitness.postfix))
    abline(v=0)
    
  }


  print("== Discrepancies ==")
  for (name in c("beta", "kappa", "alpha")) {

    ii <- mean(dz[[name]]$inter[steps.range] + 
               dz[[name]]$intra[steps.range], na.rm=T)
    g <- mean(dz[[name]]$global[steps.range], na.rm=T)
    
    print(sprintf("%s : ABS(DIFF)=%.3e", name, abs(g-ii)))

  }

  par(par.bak)

}





## =========================================================================
## plot all the Price components (SORTED BY EVOLUTIONARY VARIABLES)
## DZ can be either the output of dps.calc.price() or a results objects
## use FITNESS.POSTFIX="1" to use fitness1 and "2" to use fitness2
dps.plot.price2 <- function( results, dz=NULL, window.size=0, steps.range=NA,
                             xlim=NA, bw.adjust=1 ) {

  par.bak <- par(no.readonly=TRUE)
  layout(rbind(rep(1,4), matrix(2:9, ncol=4)))

  ## calculate price components??
  if (is.null(dz))
    dz <- dps.calc.price(results, window.size, steps.range=NA)

  ## calculate steps.range
  if (length(steps.range) != 2) 
    steps.range <- 1:nrow(dz$beta)
  else
    steps.range <- steps.range[1]:steps.range[2]

  ## plot evolutionary dynamics of beta, kappa, alpha
  plot.with.range(cbind(results$dynamics$global$M$beta[steps.range],
                        results$dynamics$global$M$kappa[steps.range],
                        results$dynamics$global$M$alpha[steps.range]),
                  cbind(sqrt(results$dynamics$global$V$beta[steps.range]),
                        sqrt(results$dynamics$global$V$kappa[steps.range]),
                        sqrt(results$dynamics$global$V$alpha[steps.range])),
                  x=steps.range,
                  main="Plasmid Replication Parameters",
                  col=c("blue", "red", "green"),
                  xlab="Time", ylab="")
  legend("topleft", c(expression(beta), expression(kappa), expression(alpha)),
         lwd=1, col=c("blue", "red", "green"))



  ## PLOT SELECTION FOR EACH VARIABLE
  for (name in c("beta", "kappa", "alpha", "tbias")) {

    if (name == "tbias") {

      ## plot time series
      b <- dz$beta$tbias[steps.range]
      k <- dz$kappa$tbias[steps.range]
      a <- dz$alpha$tbias[steps.range]
      print(sprintf("tbias | beta=%.3e | kappa=%.3e | alpha=%.3e",
                    mean(b, na.rm=T), mean(k, na.rm=T), mean(a, na.rm=T)))
      
      mplot(steps.range, cbind(b, k, a), xlab="Time",
            main=sprintf("%s", name))
      abline(h=0)
      
      ## plot densities
      b <- density(b, adjust=bw.adjust, na.rm=T)
      k <- density(k, adjust=bw.adjust, na.rm=T)
      a <- density(a, adjust=bw.adjust, na.rm=T)

      mplot(cbind(b$x, k$x, a$x), cbind(b$y, k$y, a$y),
            main=sprintf("%s", name))
      abline(v=0)
      
      
    } else {
      ## plot time series
      total = dz[[name]]$total[steps.range]
      inter = dz[[name]]$inter[steps.range]
      intra = dz[[name]]$intra[steps.range]

      print(sprintf("%s | total=%.3e | inter=%.3e | intra=%.3e",
                    name, mean(total,na.rm=T),mean(inter,na.rm=T), mean(intra,na.rm=T)))

      mplot(steps.range, cbind(inter, intra, total), xlab="Time",
            main=sprintf("%s", name), col=c("blue", "red", "black"))
      abline(h=0)
      legend("bottomleft", c("total", "inter", "intra"),
             lwd=1, col=c("black", "blue", "red"))

      ## plot densities
      total <- density(total, adjust=bw.adjust, na.rm=T)
      inter <- density(inter, adjust=bw.adjust, na.rm=T)
      intra <- density(intra, adjust=bw.adjust, na.rm=T)

      mplot(cbind(total$x, inter$x, intra$x),
            cbind(total$y, inter$y, intra$y),
            main=sprintf("%s", name), col=c("black", "blue", "red"))
      abline(v=0)
      legend("topleft", c("total", "inter", "intra"),
             lwd=1, col=c("black", "blue", "red"))
    }
  }

  print("== Discrepancies ==")
  for (name in c("beta", "kappa", "alpha")) {

    ii <- mean(dz[[name]]$inter[steps.range] + 
               dz[[name]]$intra[steps.range], na.rm=T)
    g <- mean(dz[[name]]$global[steps.range], na.rm=T)
    
    print(sprintf("%s : ABS(DIFF)=%.3e",
                  name, abs(g-ii)))

  }

  par(par.bak)

}









dps.plot.price.prediction <- function(results, steps.range=NA) {

  ## calculate steps.range
  if (length(steps.range) != 2) 
    steps.range <- 1:nrow(results$dynamics$global$M)
  else
    steps.range <- steps.range[1]:steps.range[2]
  

  par.bak <- par(no.readonly=TRUE)
  layout(matrix(1:3, ncol=1, byrow=T))

  ## calculate price components for window size 0
  dz <- dps.calc.price(results, 0, steps.range=NA)

  ## plot dynamics
  for (name in c("beta", "kappa", "alpha")) {

    ## get actual dynamics
    adyn <- results$dynamics$global$M[[name]][steps.range]
    ## calculate one-step-ahead predicted dynamics
    pdyn <- c(adyn[1], adyn[steps.range] + dz[[name]]$total[steps.range])
    pdyn <- pdyn[1:length(pdyn)-1]

    ## plot actual and predicted
    mplot(steps.range, cbind(adyn, pdyn), main=name, xlab="Time")

    ## calculate actual dz
    dz.actual <- diff(adyn)
    d <- dz.actual[steps.range] - dz[[name]]$total[steps.range]
    print(sprintf("%s : DIFF(DZ)=%.3e", name, mean(d, na.rm=T)))
    
  }

  par(par.bak)
  
}







## ============================================================================
## plot some basic features of the supplied RESULTS
dps.analyze <- function(results, window.size=1000) {

  ## store settings
  par.bak <- par(no.readonly=TRUE)
  layout(matrix(c(rep(1,2), 2:5), ncol=2, byrow=T))

  ## plot plasmid replication parameters
  plot.with.range(cbind(results$dynamics$global$M$beta,
                        results$dynamics$global$M$kappa,
                        results$dynamics$global$M$alpha),
                  cbind(sqrt(results$dynamics$global$V$beta),
                        sqrt(results$dynamics$global$V$kappa),
                        sqrt(results$dynamics$global$V$alpha)),
                  main="Plasmid Replication Parameters",
                  col=c("blue", "red", "green"),
                  xlab="Time", ylab="")
  legend("topleft", c(expression(beta), expression(kappa), expression(alpha)),
         lwd=1, col=c("blue", "red", "green"))


  ## plot division rate
  div.rate <- ma(results$dynamics$counters$div.inf / results$dynamics$counters$n,
                 window.size=window.size)
  plot(div.rate, col="blue", t="l", xlab="Time", ylab="",
       main="Average Host Division Rate")


  ## plot ave copy number
  cn <- ma(results$dynamics$counters$cn / results$dynamics$counters$n, window.size=window.size)
  plot(cn, col="blue", t="l", xlab="Time", ylab="",
       main="Average Copy Number")


  ## plot replication rate per plasmid
  rep.rate <- ma(results$dynamics$counters$rep / results$dynamics$counters$cn, window.size=window.size)
  plot(rep.rate, col="blue", t="l", xlab="Time", ylab="",
       main="Average Plasmid Replication Rate")


  ## plot replication rate per plasmid
  ht.rate <- ma(results$dynamics$counters$ht / results$dynamics$counters$cn, window.size=window.size)
  plot(ht.rate, col="blue", t="l", xlab="Time", ylab="",
       main="Average HT Rate")

  par(par.bak)
  
  
}







## ============================================================================
## plot some basic features of the supplied RESULTS (NAME is either "age" or "cn")
dps.plot.histogram <- function( results, name="age" ) {

  x <- results$histograms[[name]]$x
  y <- colSums(results$histograms[[name]]$y)

  x <- x[y!=0]
  y <- y[y!=0]

  plot(x, log10(y), t="b", xlab=name, ylab="Frequency", yaxt="n")
  decorate.log("y")
  
  
}









## ============================================================================
## plots the joint histogram of the requested pair (either "bk", "ba", "ka")
## over time
dps.plot.joint.dist <- function( results, pair="bk", fparts=NA, lv=15 ) {


  get.expression <- function( axis ) {

    switch(substr(pair, axis, axis), b=expression(beta), k=expression(kappa),
           a=expression(alpha), NA)
    
  }

  ## is results a path to results??
  if (is.character(results)) 
    ## load results from that path
    results <- dsps.load(results)

  x <- results$histograms[[pair]]$x
  y <- results$histograms[[pair]]$y

  if (all(is.na(fparts)))
    fparts <- 1:dim(results$histograms[[pair]]$z)[1]
    
  steps.per.fpart <- results$settings$steps %/% dim(results$histograms[[pair]]$z)[1]
  n.rows <- floor(sqrt(length(fparts)))
  n.cols <- length(fparts) %/% n.rows + ifelse(length(fparts) %% n.rows, 1, 0)

  layout(matrix(1:(n.rows*n.cols), nrow=n.rows, ncol=n.cols, byrow=T))

  levels <- pretty(results$histograms[[pair]]$z[fparts,,], lv)

  for (fpart in fparts) {
  
    z <- results$histograms[[pair]]$z[fpart,,] / sum(results$histograms[[pair]]$z[fpart,,])

    image(x, y, z, col=colorpanel(length(levels), low="white", high="red"),
          ##col=colorpanel(length(levels), low="white", mid="grey", high="black"),
          xlab=get.expression(1), ylab=get.expression(2), 
          main=sprintf("Steps %d-%d", (fpart-1)*steps.per.fpart, fpart*steps.per.fpart))
    
    abline(v=x, lty="dotted", lwd=0.5)
    abline(h=y, lty="dotted", lwd=0.5)

  }

}
