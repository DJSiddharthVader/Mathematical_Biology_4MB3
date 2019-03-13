#!/usr/bin/env Rscript
#suppressMessages(library('ggplot2'))
library('WaveletComp')

quiet <- function(x) {
    sink(tempfile())
    on.exit(sink())
    invisible(force(x))
}
read.ymdc <- function(filepath){
    df <- read.csv2(file=filepath,sep=',',skip=6)
    datestrcol <- paste(df[[1]],df[[2]],df[[3]],sep='-')
    df[[5]] <- as.Date(datestrcol)
    df <- as.data.frame(df)
    colnames(df) <- c("year","month","day","cases","date")
    return(df)
}

smoothData <- function(counts,a){
    startvals <- counts[seq(1,a-1)]
    endvals <- counts[seq(length(counts)-a+1,length(counts))]
    tosmooth <- seq(a,length(counts)-a)
    smoothed <- c()
    for (i in tosmooth){
        xt <- (1/(2*a+1))*sum(counts[seq(i-a,i+a)])
        smoothed <- c(smoothed,xt)
    }
    smoothed <- c(startvals,smoothed,endvals)
    return(smoothed)
}

time.plot <- function(df,add=FALSE,smooth=FALSE,a=30,...){
    #to smooth?
    if (smooth){
        #print(paste("smoothing with moving average, a =",a))
        date <- df$date
        cases <- smoothData(df$cases,a)
        #cases <- filter(df$cases,1,method='convolution',sides=2)
    } else { #normal plotting
        date <- df$date
        cases <- df$cases
    }
    #if add
    if (add){
        lines(date,cases,...)
    } else {
        plot(date,cases,...)
    }
}

periodogram <- function(df,method='pgram',timestart=1,timeend=-1,add=FALSE,...){
    #specify time range
    if (timeend == -1){
        timeend <- length(df$date)-1
    }
    #get power spectrum
    cases <- df$cases
    cases <- cases[seq(timestart,timeend)]
    print(length(cases))
    specdata <- spectrum(cases,plot=FALSE,method=method)
    period <- lapply(specdata$freq, function(x) { return(1/x) })
    spec <- specdata$spec
    #if add
    if (add){
        lines(period,spec,...)
    } else {
        plot(period,spec,...)
    }
}

maxPower <- function(df){
    specdata <- spectrum(df$cases,plot=FALSE,method='pgram')
    return(max(specdata$spec))
}

segmentData <- function(widths,len){
    lengths <- c()
    pl <- 0
    for (w in widths){
        pl <- pl + as.integer(w*len)
        lengths <- c(lengths,pl)
    }
    return(lengths)
}
multipanel <- function(df,widths){
    #set vars and layout
    a=30
    cols <- rainbow(length(widths)+1)
    pxlim = c(0,250)
    maxpower <- maxPower(df)
    pylim = c(0,(maxpower*1.8))
    lengths <- segmentData(widths,length(df$date))
    segs <- length(widths)
    lmat <- rbind(1:segs,rep(segs+1,segs),rep(segs+2,segs))
    layout(lmat,heights=1,widths=widths)
    #periodograms
    par(mar=c(4,4,0,0))
    periodogram(df,type='l',timestart=1,timeend=lengths[1],xlim=pxlim,ylab='Power',ylim=pylim,col=cols[1])
    for (i in seq(1,length(widths)-1)){
        par(mar=c(4,0,0,0))
        periodogram(df,type='l',timestart=lengths[i],timeend=lengths[i+1],xlim=pxlim,col=cols[i+1],ylim=pylim,yaxt='n')
    }
    #Time series
    par(mar=c(4,4,0,0))
    time.plot(df,col='red',type='l',xlab='Time',ylab='Cases')
    time.plot(df,a=a,smooth=TRUE,add=TRUE,col='blue',type='l')
    legend("topright", legend=c("Raw",paste("Smoothed,a=",a)), col=c("red",'blue'), box.lty=0, lty=1:1, cex=0.8)
    #waveletData
    wavelet <- quiet(analyze.wavelet(df,
                                   'cases',
                                   loess.span=0,
                                   make.pval=TRUE,
                                   n.sim=10,
                                   verbose=FALSE,
                                   upperPeriod=250))
    #waveletplot
    par(mar=c(4,4,0,0))
    wt.image(wavelet,
             color.key = 'i',
             n.levels = 250,
             show.date=TRUE,
             timelab='Time',
             periodlab='Period',
             plot.ridge=FALSE,
             #plot.legend=FALSE,
             legend.params = list(lab = "Power",
                                  shrink=0.9,
                                  mar = 2.7,
                                  n.ticks = 6,
                                  label.digits=2,
                                  width=1
                                 ),
             )
}
