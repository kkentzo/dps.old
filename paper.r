## insert all our functions
source('dps.r')


## === use the following files for results ===

## || NO-CONJ results with pconj=0 ||
## for RESULTS.BKA use results.bka <- dps.load('/data/dps/bka/results.2.16.h5')
## for DZ use load("/data/dps/bka/dz.2.16.xdr")


## || CONJ results with pconj=1e-2 ||
## for RESULTS.BKA use dps.load('/data/dpsm/htg.bka/results.2.1.h5')
## for DZ use load("/data/dpsm/htg.bka/dz.full.2.1.xdr")

## || results for varying pconj ||
## for RESULTS use results <- dps.pp.load('/data/dpsm/htg/results.xdr')


## EPS plotting :
## setEPS(); postscript("fname.eps"); plot(); dev.off()




## ===========================================================================
## plots the evolutionary dynamics of beta, kappa and alpha
## to plot use:  png("fig1.png", width=800, height=500)
plot.fig1 <- function( results.bka, highlight=F ) {


  layout(matrix(1))
  ## plot plasmid replication parameters
  plot.with.range(cbind(results.bka$dynamics$global$M$beta,
                        results.bka$dynamics$global$M$kappa,
                        results.bka$dynamics$global$M$alpha),
                  cbind(sqrt(results.bka$dynamics$global$V$beta),
                        sqrt(results.bka$dynamics$global$V$kappa),
                        sqrt(results.bka$dynamics$global$V$alpha)),
                  main="", col=c("blue", "red", "green"), ylim=c(0,1),
                  xlab="Evolutionary Time", ylab="", xaxt='n')

  ## legend("topleft", c(expression(beta), expression(kappa), expression(alpha)),
  ##        lwd=1, col=c("blue", "red", "green"))

  if (highlight) {

    start <- 1e5
    end <- 2e5

    x <- c(start, start, end, end)
    y <- c(0, 1, 1, 0)

    polygon(x, y, border=NA, col=adjustcolor("gray", alpha.f=0.6))
    
  }

  ## mark the time ticks properly
  xticks.at <- 10000 * seq(0, 100, 25)
  axis(1, at=xticks.at)
  
}






## ===========================================================================
## plots the dynamics of the Price equation components for bka
## to plot use: pdf("fig2.pdf", width=15, height=10)
plot.fig2 <- function(results.bka, dz, steps.range=c(1e5, 2e5),
                      price.ylim=c(-1.5e-5, 1.5e-5)) {


  ## calculate steps.range
  if (length(steps.range) != 2) 
    steps.range <- 1:nrow(results.bka$dynamics$global$M)
  else
    steps.range <- steps.range[1]:steps.range[2]

  par.bak <- par(no.readonly=TRUE)
  layout(matrix(1:6, ncol=3, byrow=F))

  alphabet <- c("A", "D", "B", "E", "C", "F")

  i <- 1

  for (name in c("beta", "kappa", "alpha")) {


    legend.loc <- switch(name,
                         beta="topright",
                         kappa="topleft",
                         alpha="topright")

    dyn <- results.bka$dynamics$global$M[[name]][steps.range]
    dyn.sd <- sqrt(results.bka$dynamics$global$V[[name]][steps.range])
    dyn.ylim <- range(c(dyn-dyn.sd, dyn+dyn.sd))

    ## plot the evolutionary dynamics
    plot.with.range(dyn, dyn.sd, x=steps.range, ylim=dyn.ylim,
                    main=switch(name,
                      beta=expression(paste("Evolutionary Dynamics of ", bar(beta))),
                      kappa=expression(paste("Evolutionary Dynamics of ", bar(kappa))),
                      alpha=expression(paste("Evolutionary Dynamics of ", bar(alpha)))),
                    xlab="Evolutionary Time", ylab="", col="black",
                    cex.main=2, cex.axis=1.5, cex.lab=1.5)
    
    ## write figure label
    text(steps.range[1], dyn.ylim[2], labels=alphabet[i], cex=4, adj=c(0, 1))

    i <- i+1


    total = dz[[name]]$total[steps.range]
    inter = dz[[name]]$inter[steps.range]
    intra = dz[[name]]$intra[steps.range]
    tbias = dz[[name]]$tbias[steps.range]

    

    ## plot the Price Equation components
    mplot(steps.range,
          cbind(dz[[name]]$total[steps.range],
                dz[[name]]$inter[steps.range],
                dz[[name]]$intra[steps.range],
                dz[[name]]$tbias[steps.range]),
          main=switch(name,
            beta=expression(paste("Price Equation Components of ",
                Delta,bar(beta), "  (x", 10^-5,")")),
            kappa=expression(paste("Price Equation Components of ",
                Delta,bar(kappa), "  (x", 10^-5,")")),
            alpha=expression(paste("Price Equation Components of ",
                Delta,bar(alpha), "  (x", 10^-5,")"))),
          ylim=price.ylim, yaxt="n",
          xlab="Evolutionary Time", ylab="", 
          col=c("black", "blue", "red", "green"),
          cex.main=2, cex.axis=1.5, cex.lab=1.5)
    ## write figure label
    text(steps.range[1], price.ylim[2], labels=alphabet[i], cex=4, adj=c(0, 1))
    ## draw the Y axis (major ticks)
    axis(2, at=seq(-1.5e-5, 1.5e-5, 0.5e-5),
         labels=seq(-1.5, 1.5, 0.5),
         cex.axis=1.5)
    ## draw the Y axis (minor ticks)
    ##axis(2, at=c(-1.5e-5, -0.5e-5, 0, 0.5e-5, 1.5e-5), labels=NA)
    
    abline(h=0)

    ## plot a frame??
    ## x <- c(rep(169466, 2), rep(177372, 2))
    ## y <- c(-1.5e-5, 1.5e-5, 1.5e-5, -1.5e-5)
    ## polygon(x, y, border=NA, col=adjustcolor("gray", alpha.f=0.6))

    i <- i + 1

  }

  par(par.bak)
  
}



## ===========================================================================
## plots the metrics for hosts as a function of pconj
## for plotting use: pdf("fig3.pdf", width=8, height=3)
plot.fig3 <- function( results ) {

  ## store settings
  par.bak <- par(no.readonly=TRUE)
  layout(matrix(1:3, nrow=1, byrow=T))

  ## plot host performance
  mplot(results$pconj, results$M$custom$cn,
        main="Host Copy Number", ylab="", xlab=expression(p[c]),
        log.take="x")
  ## plot host division rate
  mplot(results$pconj, results$M$custom$div.inf,
        main="Host Division Rate", ylab="", xlab=expression(p[c]),
        log.take="x")
  ## plot host death rate
  mplot(results$pconj, results$M$custom$death,
        main="Host Death Rate", ylab="", xlab=expression(p[c]),
        log.take="x")
  

  par(par.bak)

}





## ===========================================================================
## plot the mean plasmid values as a function of pconj
## to plot use :  pdf("fig4.pdf", width=6, height=4)
plot.fig4 <- function( results ) {

  ## plot beta
  mplot(results$pconj,
        cbind(results$M$global$M$beta,
              results$M$global$M$kappa,
              results$M$global$M$alpha),
        main="", ylim=c(0.4, 1), xlab=expression(p[c]),
        log.take="x")

  legend("left", c(expression(beta), expression(kappa), expression(alpha)),
         lwd=1, box.lwd=0, col=c("blue", "red", "green"))

}












## ===========================================================================
## plots the metrics for plasmids as a function of pconj
## for plotting use: pdf("fig5.pdf", width=8, height=4)
plot.fig5 <- function( results ) {

  ## store settings
  par.bak <- par(no.readonly=TRUE)
  layout(matrix(1:2, nrow=1, byrow=T))

  ## plot plasmid rep rate
  mplot(results$pconj, results$M$global$M$nr,
        main="Plasmid Replication Rate", ylab="", xlab=expression(p[c]),
        log.take="x")
  ## plot plasmid conj rate
  mplot(results$pconj, results$M$global$M$ht,
        main="Plasmid Conjugation Rate", ylab="", xlab=expression(p[c]),
        log.take="x")

  par(par.bak)

}




## ===========================================================================
## plot plasmid parameter variance as a function of pconj
## for plotting use: pdf("fig6.pdf", width=8, height=3)
plot.fig6 <- function(results) {

  par.bak <- par(no.readonly=TRUE)
  layout(matrix(1:3, nrow=1))

  ## plot beta
  mplot(results$pconj.values,
        cbind(sqrt(results$M$global$V$beta),
              sqrt(results$M$inter$V$beta),
              sqrt(results$M$intra$V$beta)),
        ylim=c(log10(1e-3), log10(4e-2)),
        main=expression(sigma(beta)),
        xlab=expression(p[c]), type="l", log="xy")

  legend("bottomright", c("global", "inter", "intra"),
         lwd=1, col=c("blue", "red", "green"))

  mplot(results$pconj.values,
        cbind(sqrt(results$M$global$V$kappa),
              sqrt(results$M$inter$V$kappa),
              sqrt(results$M$intra$V$kappa)),
        ylim=c(log10(1e-3), log10(4e-2)),
        main=expression(sigma(kappa)),
        xlab=expression(p[c]), type="l", log="xy")
  
  mplot(results$pconj.values,
        cbind(sqrt(results$M$global$V$alpha),
              sqrt(results$M$inter$V$alpha),
              sqrt(results$M$intra$V$alpha)),
        ylim=c(log10(1e-3), log10(4e-2)),
        main=expression(sigma(alpha)),
        xlab=expression(p[c]), type="l", log="xy")
        
  par(par.bak)

  
}







## ===========================================================================
## plots the averages of the Price equation components as a function of pconj
## to plot use: pdf("fig7.pdf", width=16, height=6)
plot.fig7 <- function(results) {


  par.bak <- par(no.readonly=TRUE)
  layout(matrix(1:3, nrow=1, byrow=T))


  ## create a 1-row/3-col plot with blue:inter, red:intra and different
  ## line types for each component (n^R, n^H ...)

  for (name in c("beta", "kappa", "alpha")) {

    ylim <- switch(name,
                   beta=c(-2e-5, 2e-5),
                   kappa=c(-2e-5, 2e-5),
                   alpha=c(-2e-6, 2e-6))
    ##ylim <- NA

    main <- switch(name,
                   beta=expression(paste("Selection on ", beta, "  (x", 10^-5, ")")),
                   kappa=expression(paste("Selection on ", kappa, "  (x", 10^-5, ")")),
                   alpha=expression(paste("Selection on ", alpha, "  (x", 10^-6, ")")))

    ##ylim <- c(-2e-5, 2e-5)

    inter <- cbind(sapply(c("fitness", "nr", "ht", "death"),
                          function (fname)
                          results$M$inter$C[[sprintf("%s.%s", name, fname)]],
                          USE.NAMES=F))

    intra <- cbind(sapply(c("fitness", "nr", "ht", "death"),
                          function (fname)
                          results$M$intra$C[[sprintf("%s.%s", name, fname)]],
                          USE.NAMES=F))

    tbias <- results$M$global$M[[sprintf("t%s", name)]]

    ## change signs in DEATH columns
    inter[,4] <- - inter[,4]
    intra[,4] <- - intra[,4]

    mplot(results$pconj.values,
          cbind(inter, intra, tbias), 
          col=c(rep("blue",4), rep("red",4), "green"), ylim=ylim,
          main=main, log.take="x",
          xlab=expression(p[c]), yaxt="n",
          cex.main=2.5, cex.lab=2, cex.axis=2, ltype=1:4)

    ## draw the Y axis (major ticks)
    axis(2,
         at=switch(name,
           beta=seq(-2e-5, 2e-5, 1e-5),
           kappa=seq(-2e-5, 2e-5, 1e-5),
           alpha=seq(-2e-6, 2e-6, 1e-6)),
         labels=switch(name,
           beta=seq(-2, 2, 1),
           kappa=seq(-2, 2, 1),
           alpha=seq(-2, 2, 1)),
         cex.axis=2)
    
    abline(h=0)

    if (name == "beta") {
      legend("topleft",
             c("Total Selection (intra)",
               expression(paste("Selection due to ", n^R, " (intra)")),
               expression(paste("Selection due to ", n^H, " (intra)")),
               expression(paste("Selection due to ", n^D, " (intra)"))),
             lwd=1, lty=1:4, bty="n", cex=2, col="red")
      legend("bottomleft",
             c("Total Selection (inter)",
               expression(paste("Selection due to ", n^R, " (inter)")),
               expression(paste("Selection due to ", n^H, " (inter)")),
               expression(paste("Selection due to ", n^D, " (inter)"))),
             lwd=1, lty=1:4, bty="n", cex=2, col="blue")
      legend("topleft", "Transmission Bias", inset=c(0, 0.3),
             lwd=1, bty="n", cex=2, col="green")
    }
  }

  par(par.bak)    

}


## ==================================================================================
## plot the intra-cellular covariances of (beta,alpha) and (kappa, alpha)
## as a function of pconj
## to plot use :  pdf("fig8.pdf", width=8, height=5)
plot.fig8 <- function( results ) {

  level = "intra"

  layout(matrix(1))

  mplot(results$pconj,
        cbind(results$M[[level]]$C$beta.alpha,
              results$M[[level]]$C$kappa.alpha),
              ##results$M[[level]]$C$alpha.fitness),
        xlab=expression(p[c]), ylab="", type="l",
        main=expression(paste("Average Within-Host Covariances",
            "  (x", 10 ^ -6, ")")),
        log.take="x", yaxt="n",
        ylim=c(-4e-6, 4e-6))
  axis(2, at=seq(-4e-6, 4e-6, 2e-6),
       labels=seq(-4, 4, 2),
       cex.axis=1)
  abline(h=0)

  ## plot error bars
  ## error.bar(log10(results$pconj),
  ##           cbind(results$M[[level]]$C$beta.alpha,
  ##                 results$M[[level]]$C$kappa.alpha),
  ##                 ##results$M[[level]]$C$alpha.fitness),
  ##           cbind(results$S[[level]]$C$beta.alpha,
  ##                 results$S[[level]]$C$kappa.alpha)  / sqrt(results$runs),
  ##           ##results$S[[level]]$C$alpha.fitness),
  ##           col=matrix(rep(c("blue", "red"),
  ##             length(results$pconj)), ncol=2, byrow=T))
  
  ## legend("topleft",
  ##        c(expression(paste("<", E[i], "[cov"[j], "(", beta, ",", alpha, ")]>")),
  ##          expression(paste("<", E[i], "[cov"[j], "(", kappa, ",", alpha, ")]>"))),
  ##        lwd=1, col=c("blue", "red"))
                      
        
}



## ==================================================================================
## ==================================================================================
##                              AD-HOC FUNCTIONS
## ==================================================================================
## ==================================================================================



