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

periodogram <- function(df,method='pgram',timestart=1,timerange=-1,add=FALSE,...){
    #specify time range
    if (timerange == -1){
        timerange <- length(df$date)-1
    }
    #get power spectrum
    specdata <- spectrum(df$cases[seq(timestart,timestart+timerange)],plot=FALSE,method=method)
    period <- lapply(specdata$freq, function(x) { return(1/x) })
    spec <- specdata$spec
    #if add
    if (add){
        lines(period,spec,...)
    } else {
        plot(period,spec,...)
    }
}

multipanel <- function(df){
    #waveletData
    wavelet <- quiet(analyze.wavelet(df,
                                   'cases',
                                   loess.span=0,
                                   make.pval=TRUE,
                                   n.sim=10,
                                   verbose=FALSE,
                                   upperPeriod=200))
    #waveletplot
    wt.image(wavelet,
             color.key = "quantile",
             n.levels = 250,
             show.date=TRUE,
             legend.params = list(lab = "wavelet power levels",
                                  mar = 2.7,
                                  n.ticks = 5,
                                  label.digits=5
                                 ),
             )
    par(mfrow = c(2,1))
    time.plot(df,col='red',type='l')
    time.plot(df,smooth=TRUE,add=TRUE,col='blue',type='l')
    periodogram(df,type='l',timestart=1,col='red',xlim=c(0,250))
    periodogram(df,type='l',timestart=1300,timerange=1000,add=TRUE,col='blue')
    periodogram(df,type='l',timestart=2300,timerange=360,add=TRUE,col='green')
}

#if name == main
#if (sys.nframe() == 0){
#    args <-  commandArgs(trailingOnly=TRUE)
#    df <- read.ymdc('/home/sidreed/Documents/4thyear/Math4MB3/assg4/answers/sections/meas_uk__lon_1944-94_wk.csv')
#    multipanel(df)
#}
