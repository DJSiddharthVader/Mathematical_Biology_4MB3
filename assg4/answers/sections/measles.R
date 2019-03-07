#!/usr/bin/env Rscript

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

time.plot <- function(df,add=FALSE,smooth=TRUE,a,...){
    #to smooth?
    if (smooth | add){
        print(paste("smoothing with moving average, a =",a))
        date <- df$date
        cases <- smoothData(df$cases,a)
        #cases <- filter(df$cases,1,method='convolution',sides=2)
    } else { #normal plotting
        time <- df$date
        cases <- df$cases
    }
    #if add
    if (add){
        lines(date,cases,...)
    } else {
        plot(date,cases,...)
    }
}

periodogram <- function(df,method,timestart=1,timerange=-1,add=FALSE,...){
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

multipanel <- function(){
}
#df <- read.ymdc('./sections/meas_uk__lon_1944-94_wk.csv')
