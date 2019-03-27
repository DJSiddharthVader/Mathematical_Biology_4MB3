#!usr/bin/env Rscript
library('deSolve')

singleWaterModel <- function(t,vars,params){
    with(as.list(c(vars,params)),{
        dx <- mu - mu*x - x*y*beta_i - x*w*beta_w # dx/dt
        dy <- x*y*beta_i + x*w*beta_w - gamma*y - mu*y - alpha*y # dy/dt
        dz <- gamma*y - mu*z # dz/dt
        dw <- beta_w*y - sigma*w #dW/dt
        vec.fld <- list(c(dx,dy,dz,dw))
        return(vec.fld)
    })
}

solveSingleModel <- function(func,ic,params,tmax=1,steps=5000){
    with(as.list(c(ic,params)),{
        times <- seq(0,tmax,by=tmax/steps)
        soln <- ode(y=ic,
                    times=times,
                    func=func,
                    parms=params)
        return(soln)
    })
}

plotSoln <- function(soln,...){
    scol <- 'red'
    icol <- 'blue'
    rcol <- 'green'
    #plotting
    plot(x=0,y=0,type='l',...)
    lines(soln[,'x'],type='l',col=scol)
    lines(soln[,'y'],type='l',col=icol)
    lines(soln[,'z'],type='l',col=rcol)
    legend("topright",
           legend=c('S','I','R'),
           col=c('red','blue','green'),
           inset=0.02,
           box.lty=0,
           lty=1,
           cex=0.8)
}

distmat <- function(rows,cols,scale=0.07,method='euclidian'){
    patchpairs <- expand.grid(1:rows,1:cols)
    patches <- dim(patchpairs)[1]
    dmat = matrix(rep(-1,patches*patches),nrow=patches,ncol=patches)
    maxdist <- dist(rbind(c(1,1),c(rows,cols)),method=method)[1][1]
    for (p1 in 1:patches){
        i1 <- patchpairs[p1,1]
        j1 <- patchpairs[p1,2]
        for (p2 in 1:patches){
            i2 <- patchpairs[p2,1]
            j2 <- patchpairs[p2,2]
            p1p2dist <- dist(rbind(c(i1,j1),c(i2,j2)),method=method)[1][1]
            if (p1 == p2){
                 mij <- 1-(p1p2dist/maxdist)
            } else {
                mij <- (1-(p1p2dist/maxdist))*scale
                if (mij > 1){mij <- 1}
            }
            dmat[p1,p2] <- mij
        }
    }
    return(dmat)
}

neighourmat <- function(rows,cols,influence=0.15){
    patchpairs <- expand.grid(1:rows,1:cols)
    patches <- dim(patchpairs)[1]
    dmat = matrix(rep(-1,patches*patches),nrow=patches,ncol=patches)
    for (p1 in 1:patches){
        i1 <- patchpairs[p1,1]
        j1 <- patchpairs[p1,2]
        for (p2 in 1:patches){
            i2 <- patchpairs[p2,1]
            j2 <- patchpairs[p2,2]
            h <- (((i1-1) <= i2) & ((i1+1) >= i2))
            v <- (((j1-1) <= j2) & ((j1+1) >= j2))
            if (p1 == p2){
                dmat[p1,p2] <- 1
            } else if (h & v){
                dmat[p1,p2] <- influence
            } else {
                dmat[p1,p2] <- 0
            }
        }
    }
    return(dmat)
}

initModel <- function(rows,cols,i0,w0){
    p = rows*cols
    x <- rep(1-i0,p)
    y <- rep(i0,p)
    z <- rep(0,p)
    w <- rep(w0,p)
    return(c(x,y,z,w))
}

multiWaterModel <- function(t,state,params){
    with(as.list(c(params)),{
        x <- state[(0*patches+1):(1*patches+0)]
        y <- state[(1*patches+1):(2*patches+0)]
        z <- state[(2*patches+1):(3*patches+0)]
        w <- state[(3*patches+1):(4*patches+0)]
        dx <- mu - mu*x - dmat %*% x*y*beta_i - dmat %*% x*w*beta_w # dx/dt
        dy <- dmat %*% x*y*beta_i + dmat %*% x*w*beta_w - gamma*y - mu*y - alpha*y # dy/dt
        dz <- gamma*y - mu*z # dz/dt
        dw <- dmat %*% y*beta_w - sigma*w #dW/dt
        vec.fld <- list(c(dx,dy,dz,dw))
        return(vec.fld)
    })
}

solveMultiModel <- function(func,ic,params,tmax=1,steps=5000){
    with(as.list(c(ic,params)),{
        ivec <- initModel(rows,cols,i0,w0)
        times <- seq(0,tmax,by=tmax/steps)
        soln <- ode.1D(y=ivec,
                       times=times,
                       func=func,
                       parms=params,
                       nspec=4,
                       dimens=patches,
                       names=c('S','I','R','W'))
        return(soln)
    })
}

multiPlotSoln <- function(rows,cols,soln,...){
    patches <- rows*cols
    scol <- 'red'
    icol <- 'blue'
    rcol <- 'green'
    #get data
    s <- soln[,(0*patches+2):(1*patches+1)]
    i <- soln[,(1*patches+2):(2*patches+1)]
    r <- soln[,(2*patches+2):(3*patches+1)]
    w <- soln[,(3*patches+2):(4*patches+1)]
    #plotting
    par(mfrow=c(rows,cols))
    for (patch in 1:patches){
        initPlot(patch,rows,cols,...)
        lines(s[,patch],type='l',col=scol,lty=1)
        lines(i[,patch],type='l',col=icol,lty=1)
        lines(r[,patch],type='l',col=rcol,lty=1)
    }
}

initPlot <- function(patch,rows,cols,
                     xlab='Time',
                     ylab='Population Fraction',
                     ...){
    left <- (patch > cols*(rows -1))
    bottom <- ((patch %% cols) == 1)
    if (left & bottom){
        par(mar=c(4,4,0,0))
        plot(x=0,y=0,type='l',xlab=xlab,ylab=ylab,...)
    } else if (left) {
        par(mar=c(4,0,0,0))
        plot(x=0,y=0,type='l',xlab=xlab,yaxt='n',...)
    } else if (bottom) {
        par(mar=c(0,4,0,0))
        plot(x=0,y=0,type='l',ylab=ylab,xaxt='n',...)
    } else {
        par(mar=c(0,0,0,0))
        plot(x=0,y=0,type='l',xaxt='n',yaxt='n',...)
    }
    topright <- (patch == cols)
    if (topright){
        legend("topright",
               legend=c('S','I','R'),
               col=c('red','blue','green'),
               box.lty=0,
               lty=1,
               cex=0.8)
    }
}

