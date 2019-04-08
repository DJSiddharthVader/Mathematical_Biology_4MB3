#!usr/bin/env Rscript
library(deSolve)
library(ReacTran)
library(stringr)
library(oce)
library(progress)

#Single Patch
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
    scol <- 'green'
    icol <- 'red'
    rcol <- 'black'
    #plotting
    plot(x=0,y=0,type='l',...)
    lines(soln[,'x'],type='l',col=scol)
    lines(soln[,'y'],type='l',col=icol)
    lines(soln[,'z'],type='l',col=rcol)
    legend("topright",
           legend=c('S','I','R'),
           col=c('green','red','black'),
           inset=0.02,
           box.lty=0,
           lty=1,
           cex=0.8)
}
#Multi Patch
makepb <- function(total){
    pb <- progress_bar$new(format = "[:bar] :percent eta: :eta",total = total, clear = FALSE, width= 75)
    return(pb)
}
initModel <- function(rows,cols,i0,w0){
    p = rows*cols
    i <- matrix(rep(0,p),rows,cols)
    i[sample(1:length(i), as.integer(0.08*rows*cols)+1)] <- i0
    s <- matrix(rep(1,p),rows,cols) - i
    r <- matrix(rep(0,p),rows,cols)
    w <- matrix(rep(0,p),rows,cols)
    w[sample(1:length(w), as.integer(0.05*rows*cols)+1)] <- w0
    return(c(s,i,r,w))
}
multiWaterModel <- function(t,state,params){
    with(as.list(c(params)),{
        #compartment matrices
        patches <- rows*cols
        s <- matrix(state[(0*patches+1):(1*patches+0)],rows,cols)
        i <- matrix(state[(1*patches+1):(2*patches+0)],rows,cols)
        r <- matrix(state[(2*patches+1):(3*patches+0)],rows,cols)
        w <- matrix(state[(3*patches+1):(4*patches+0)],rows,cols)
        #terms
        person_infect <- s*i*beta_i
        water_infect <- s*w*beta_w
        sdisperse <- tran.2D(s,dx=dx,dy=dy,D.x=Dsr,D.y=Dsr)$dC
        idisperse <- tran.2D(i,dx=dx,dy=dy,D.x=Di,D.y=Di)$dC
        rdisperse <- tran.2D(r,dx=dx,dy=dy,D.x=Dsr,D.y=Dsr)$dC
        wdisperse <- tran.2D(w,dx=dx,dy=dy,D.x=Dw,D.y=Dw)$dC
        #derivatives
        ds <- mu - mu*s - person_infect - water_infect + sdisperse
        di <- person_infect + water_infect - i*gamma -i*mu -i*alpha + idisperse
        dr <- gamma*i - mu*r + rdisperse
        dw <- xi*i - sigma*w + wdisperse
        vec.fld <- list(c(ds,di,dr,dw))
        return(vec.fld)
    })
}
sanitationModel <- function(t,state,params){
    with(as.list(c(params)),{
        #compartment matrices
        patches <- rows*cols
        s <- matrix(state[(0*patches+1):(1*patches+0)],rows,cols)
        i <- matrix(state[(1*patches+1):(2*patches+0)],rows,cols)
        r <- matrix(state[(2*patches+1):(3*patches+0)],rows,cols)
        w <- matrix(state[(3*patches+1):(4*patches+0)],rows,cols)
        #terms
        person_infect <- s*i*beta_i
        water_infect <- s*w*beta_w
        sdisperse <- tran.2D(s,dx=dx,dy=dy,D.x=Dsr,D.y=Dsr)$dC
        idisperse <- tran.2D(i,dx=dx,dy=dy,D.x=Di,D.y=Di)$dC
        rdisperse <- tran.2D(r,dx=dx,dy=dy,D.x=Dsr,D.y=Dsr)$dC
        wdisperse <- tran.2D(w,dx=dx,dy=dy,D.x=Dw,D.y=Dw)$dC
        #derivatives
        ds <- mu - mu*s - person_infect - water_infect + sdisperse
        di <- person_infect + water_infect - i*gamma -i*mu -i*alpha + idisperse
        dr <- gamma*i - mu*r + rdisperse
        dw <- xi*i - sigma*w - rho*w + wdisperse
        vec.fld <- list(c(ds,di,dr,dw))
        return(vec.fld)
    })
}
vaccinationModel <- function(t,state,params){
    with(as.list(c(params)),{
        #compartment matrices
        patches <- rows*cols
        s <- matrix(state[(0*patches+1):(1*patches+0)],rows,cols)
        i <- matrix(state[(1*patches+1):(2*patches+0)],rows,cols)
        r <- matrix(state[(2*patches+1):(3*patches+0)],rows,cols)
        w <- matrix(state[(3*patches+1):(4*patches+0)],rows,cols)
        #terms
        person_infect <- s*i*beta_i
        water_infect <- s*w*beta_w
        sdisperse <- tran.2D(s,dx=dx,dy=dy,D.x=Dsr,D.y=Dsr)$dC
        idisperse <- tran.2D(i,dx=dx,dy=dy,D.x=Di,D.y=Di)$dC
        rdisperse <- tran.2D(r,dx=dx,dy=dy,D.x=Dsr,D.y=Dsr)$dC
        wdisperse <- tran.2D(w,dx=dx,dy=dy,D.x=Dw,D.y=Dw)$dC
        #derivatives
        ds <- mu - mu*s - person_infect - water_infect - nu*s + sdisperse
        di <- person_infect + water_infect - i*gamma -i*mu -i*alpha + idisperse
        dr <- gamma*i - mu*r + nu*s + rdisperse
        dw <- xi*i - sigma*w + wdisperse
        vec.fld <- list(c(ds,di,dr,dw))
        return(vec.fld)
    })
}
antibioticModel <- function(t,state,params){
    with(as.list(c(params)),{
        #compartment matrices
        patches <- rows*cols
        s <- matrix(state[(0*patches+1):(1*patches+0)],rows,cols)
        i <- matrix(state[(1*patches+1):(2*patches+0)],rows,cols)
        r <- matrix(state[(2*patches+1):(3*patches+0)],rows,cols)
        w <- matrix(state[(3*patches+1):(4*patches+0)],rows,cols)
        #terms
        person_infect <- s*i*beta_i
        water_infect <- s*w*beta_w
        sdisperse <- tran.2D(s,dx=dx,dy=dy,D.x=Dsr,D.y=Dsr)$dC
        idisperse <- tran.2D(i,dx=dx,dy=dy,D.x=Di,D.y=Di)$dC
        rdisperse <- tran.2D(r,dx=dx,dy=dy,D.x=Dsr,D.y=Dsr)$dC
        wdisperse <- tran.2D(w,dx=dx,dy=dy,D.x=Dw,D.y=Dw)$dC
        #derivatives
        ds <- mu - mu*s - person_infect - water_infect + sdisperse
        di <- person_infect + water_infect - i*gamma -i*mu -i*alpha -eta*i + idisperse
        dr <- gamma*i - mu*r +eta*i + rdisperse
        dw <- xi*i - sigma*w + wdisperse
        vec.fld <- list(c(ds,di,dr,dw))
        return(vec.fld)
    })
}
solveMultiModel <- function(func,ic,params,tmax=1,steps=5000,lrw=99999999){
    with(as.list(c(ic,params)),{
        ivec <- initModel(rows,cols,i0,w0)
        times <- seq(0,tmax,by=tmax/steps)
        soln <- ode.2D(y=ivec,
                       times=times,
                       func=func,
                       parms=params,
                       nspec=4,
                       lrw=lrw,
                       dimens=c(rows,cols))
        return(soln)
    })
}
gifplot <- function(rows,cols,soln,tmax){
    #get data
    patches <- rows*cols
    s <- soln[,(0*patches+2):(1*patches+1)]
    i <- soln[,(1*patches+2):(2*patches+1)]
    r <- soln[,(2*patches+2):(3*patches+1)]
    w <- soln[,(3*patches+2):(4*patches+1)]
    #water bars
    wticks <- seq(0,max(w),by=max(w)/1000)
    wcfnc <- colorRamp(c("blue","green"))
    wlist <- lapply(lapply(wticks,function(x) x/max(w)),wcfnc)
    wmap <- unlist(lapply(wlist,function(x) rgb(x[1], x[2], x[3],maxColorValue=255)))
    #color bars
    ticks <- seq(0,1,by=0.00001)
    cfnc <- colorRamp(c( "red","yellow","blue" ))
    rgblist <- lapply(ticks,cfnc)
    cmap <- unlist(lapply(rgblist,function(x) rgb(x[1], x[2], x[3], maxColorValue=255)))
    imgdir <- 'heat_map_frames'
    dir.create(imgdir,showWarnings=FALSE)
    stepsize <- tmax/dim(i)[1]
    dec <- 4
    pb <- makepb(dim(i)[1])
    for (t in 1:dim(i)[1]){#dim(soln)[2]){
        timestep <- t  #round(t*stepsize,1)
        name <- paste(imgdir,'/','plotframe',str_pad(t,4,pad='0'),'.png',sep='')
        png(name)
        par(mfrow=c(2,2))
        omar <- par('mar')
        #S plot
        drawPalette(zlim=c(0,1),col=cmap)
        image(z=matrix(s[t,],rows,cols),

              main=paste('S,','Mean %=',round(sum(s[t,]/patches),dec)*100,',Var =',round(var(s[t,]),dec)*100),
              zlim=c(0,1),
              col=cmap,
              axes=FALSE,
              ask=FALSE)
        #I plot
        par(mar=omar)
        drawPalette(zlim=c(0,1),col=cmap)
        image(matrix(i[t,],rows,cols),
              main=paste('I,','Mean %=',round(sum(i[t,]/patches),dec)*100,',Var =',round(var(i[t,]),dec)*100),
              zlim=c(0,1),
              col=cmap,
              axes=FALSE,
              ask=FALSE)
        #R plot
        par(mar=omar)
        drawPalette(zlim=c(0,1),col=cmap)
        image(matrix(r[t,],rows,cols),
              main=paste('R,','Mean %=',round(sum(r[t,]/patches),dec)*100,',Var =',round(var(r[t,]),dec)*100),
              zlim=c(0,1),
              col=cmap,
              axes=FALSE,
              ask=FALSE)
        #W plot
        par(mar=omar)
        drawPalette(zlim=c(0,max(w)),col=wmap)
        image(matrix(w[t,],rows,cols),
              main=paste('W,','Mean =',round(mean(w[t,]),dec),',Var =',round(var(w[t,]),dec)),
              zlim=c(0,max(w)),
              col=wmap,
              axes=FALSE,
              ask=FALSE)
        par(oma=c(0,0,2,0))
        title(paste('Compartments at Time =',timestep),outer=TRUE)
        dev.off()
        pb$tick()
    }
    mycommand <- paste('convert ',imgdir,'/*.png -delay 150 filtering.gif',sep='')
    system(mycommand)
}
compareTreatments <- function(rows=80,...){
    set.seed(9)
    rows <- 80
    cols <- rows
    tmax <- 450
    steps <- 900
    #params
    mu = 1/21170 # avg lifespan 58 yrs == 21170
    beta_i = 0.1 #transmission r. S between I
    gamma = 1/8 #rate of recovery (days)
    sigma = 1/22 #rate of water removal
    beta_w = 0.3 #transmission r. I between W
    xi = 1/22 #shedding rate I into W (force it equal to sigma)
    alpha = 0 #death rate by cholera
    #treatment params
    rho = 1/20#parameters for water treatment [day^-1]
    nu = 1/200#parameters for vaccination [day^-1]
    eta = 1/60#parameters for antibiotics [day^-1]
    #spatial params
    dx <- 0.3
    dy <- dx
    Dsr <- 0.001
    Di <- 0
    Dw <- 0.0001
    params <- list(rows=rows,
                   cols=cols,
                   mu=mu,
                   beta_i=beta_i,
                   gamma=gamma,
                   sigma=sigma,
                   beta_w=beta_w,
                   xi=xi,
                   alpha=alpha,
                   rho=rho,
                   nu=nu,
                   eta=eta,
                   dx=dx,
                   dy=dy,
                   Dsr=Dsr,
                   Di=Di,
                   Dw=Dw)
    #ics
    i0 <- 0.01
    w0 <- 0.03
    ic <- list(rows=rows,
               cols=cols,
               i0=i0,
               w0=w0)

    base <- solveMultiModel(multiWaterModel,ic,params,tmax=tmax,steps=steps)
    bi <- base[,(rows*rows+2):(rows*rows*2+1)]
    bim <- apply(bi,1,mean)

    san <- solveMultiModel(sanitationModel,ic,params,tmax=tmax,steps=steps)
    si <- san[,(rows*rows+2):(rows*rows*2+1)]
    sim <- apply(si,1,mean)

    anti <- solveMultiModel(antibioticModel,ic,params,tmax=tmax,steps=steps)
    ai <- anti[,(rows*rows+2):(rows*rows*2+1)]
    aim <- apply(ai,1,mean)

    vacc <- solveMultiModel(vaccinationModel,ic,params,tmax=tmax,steps=steps)
    vi <- vacc[,(rows*rows+2):(rows*rows*2+1)]
    vim <- apply(vi,1,mean)

    overlaymax <- 0.1 #ymax on overlay graph
    treatcol = c("red","cornflowerblue", "orange", "brown") #for overlay graph
    png('spatialTreatments.png')
    plot(x=0, y=0, type = "n",xlim=c(0, tmax), ylim=c(0, overlaymax), xaxs="i", yaxs="i",xlab = "Time (days)", ylab="Proportion (I)")
    lines(bim,col="red",lty=1, lwd=4,)
    lines(sim,col="cornflowerblue",lty=1, lwd=4,)
    lines(vim,col="brown",lty=1, lwd=4,)
    lines(aim,col="orange",lty=1, lwd=4,)
    legend("topright", legend = c("Base", "Sanitation, $\\rho = 0.05$","Vaccination, $\\nu = 0.005$", "Antibiotics, $\\eta = 0.016$"), lty=c(1,1,1,1), col = treatcol)
    dev.off()
}
main3 <- function(){
    set.seed(9)
    rows <- 80
    cols <- rows
    #params
    mu = 1/21170 # avg lifespan 58 yrs == 21170
    beta_i = 0.1 #transmission r. S between I
    gamma = 1/8 #rate of recovery (days)
    sigma = 1/22 #rate of water removal
    beta_w = 0.3 #transmission r. I between W
    xi = 1/22 #shedding rate I into W (force it equal to sigma)
    alpha = 0 #death rate by cholera
    #treatment params
    rho = 1/20#parameters for water treatment [day^-1]
    nu = 1/200#parameters for vaccination [day^-1]
    eta = 1/60#parameters for antibiotics [day^-1]
    #spatial params
    dx <- 0.3
    dy <- dx
    Dsr <- 0.001
    Di <- 0
    Dw <- 0.0001
    params <- list(rows=rows,
                   cols=cols,
                   mu=mu,
                   beta_i=beta_i,
                   gamma=gamma,
                   sigma=sigma,
                   beta_w=beta_w,
                   xi=xi,
                   alpha=alpha,
                   rho=rho,
                   nu=nu,
                   eta=eta,
                   dx=dx,
                   dy=dy,
                   Dsr=Dsr,
                   Di=Di,
                   Dw=Dw)
    #ics
    i0 <- 0.01
    w0 <- 0.03
    ic <- list(rows=rows,
               cols=cols,
               i0=i0,
               w0=w0)
    tmax <- 150
    steps <- 300
    soln <- solveMultiModel(vaccinationModel,ic,params,tmax=tmax,steps=steps)
#    gifplot(rows,cols,soln,tmax)
#    cmd <- "mpv filtering.gif --loop=inf"
#    system(cmd)
    return(soln[,(80*80+2):(80*80*2+1)])
}

#Depreciated
#distmat <- function(rows,cols,scale=0.07,method='euclidian'){
#    patchpairs <- expand.grid(1:rows,1:cols)
#    patches <- dim(patchpairs)[1]
#    dmat = matrix(rep(-1,patches*patches),nrow=patches,ncol=patches)
#    maxdist <- dist(rbind(c(1,1),c(rows,cols)),method=method)[1][1]
#    for (p1 in 1:patches){
#        i1 <- patchpairs[p1,1]
#        j1 <- patchpairs[p1,2]
#        for (p2 in 1:patches){
#            i2 <- patchpairs[p2,1]
#            j2 <- patchpairs[p2,2]
#            p1p2dist <- dist(rbind(c(i1,j1),c(i2,j2)),method=method)[1][1]
#            if (p1 == p2){
#                 mij <- 1-(p1p2dist/maxdist)
#            } else {
#                mij <- (1-(p1p2dist/maxdist))*scale
#                if (mij > 1){mij <- 1}
#            }
#            dmat[p1,p2] <- mij
#        }
#    }
#    return(dmat)
#}
#neighourmat <- function(rows,cols,influence=0.15){
#    patchpairs <- expand.grid(1:rows,1:cols)
#    patches <- dim(patchpairs)[1]
#    dmat = matrix(rep(-1,patches*patches),nrow=patches,ncol=patches)
#    for (p1 in 1:patches){
#        i1 <- patchpairs[p1,1]
#        j1 <- patchpairs[p1,2]
#        for (p2 in 1:patches){
#            i2 <- patchpairs[p2,1]
#            j2 <- patchpairs[p2,2]
#            h <- (((i1-1) <= i2) & ((i1+1) >= i2))
#            v <- (((j1-1) <= j2) & ((j1+1) >= j2))
#            if (p1 == p2){
#                dmat[p1,p2] <- 1
#            } else if (h & v){
#                dmat[p1,p2] <- influence
#            } else {
#                dmat[p1,p2] <- 0
#            }
#        }
#    }
#    return(dmat)
#}
#initModel <- function(rows,cols,i0,w0){
#    p = rows*cols
#    x <- rep(1-i0,p)
#    y <- rep(i0,p)
#    z <- rep(0,p)
#    w <- rep(w0,p)
#    return(c(x,y,z,w))
#}
#multiWaterModel <- function(t,state,params){
#    with(as.list(c(params)),{
#        x <- state[(0*patches+1):(1*patches+0)]
#        y <- state[(1*patches+1):(2*patches+0)]
#        z <- state[(2*patches+1):(3*patches+0)]
#        w <- state[(3*patches+1):(4*patches+0)]
#        dx <- mu - mu*x - nmat %*% x*y*beta_i - wmat %*% x*w*beta_w # dx/dt
#        dy <- nmat %*% x*y*beta_i + wmat %*% x*w*beta_w - gamma*y - mu*y - alpha*y # dy/dt
#        dz <- gamma*y - mu*z # dz/dt
#        dw <- wmat %*% y*beta_w - sigma*w #dW/dt
#        vec.fld <- list(c(dx,dy,dz,dw))
#        return(vec.fld)
#    })
#}
#solveMultiModel <- function(func,ic,params,tmax=1,steps=5000){
#    with(as.list(c(ic,params)),{
#        ivec <- initModel(rows,cols,i0,w0)
#        times <- seq(0,tmax,by=tmax/steps)
#        soln <- ode.1D(y=ivec,
#                       times=times,
#                       func=func,
#                       parms=params,
#                       nspec=4,
#                       dimens=patches,
#                       names=c('S','I','R','W'))
#        return(soln)
#    })
#}
#betamat <- function(N,p){
#    betamat <- p*matrix(rep(1,N^2),N,N) + (1-p)*diag(N)
#    return(betamat)
#}
#checkBoundary <- function(x,y,maxx,maxy,rx,ry){
#    if (x-rx < 1){
#        xmin <- 1
#    } else{
#        xmin <- x-rx
#    }
#    if ((y-ry) < 1){
#        ymin <- 1
#    } else{
#        ymin <- y-ry
#    }
#    if ((x+rx) > max(patches$xpos)){
#        xmax <- max(patches$xpos)
#    } else{
#        xmax <- x+ry
#    }
#    if ((y+ry) > max(patches$ypos)){
#        ymax <- max(patches$ypos)
#    } else{
#        ymax <- y+ry
#    }
#    return(list(xrange=xmin:xmax,yrange=ymin:ymax))
#}
#getNeighbours  <- function(x,y,maxx,maxy,rx=1,ry=1,self=FALSE){
#    ranges <- checkBoundary(x,y,maxx,maxy,rx,ry)
#    patchpairs <- expand.grid(ranges$xrange,ranges$yrange)
#    if (self){
#        return(filtery)
#    } else {
#        filterself <- filtery[(filtery$xpos != x) | (filtery$ypos != y),]
#        return(filterself)
#    }
#}
#spreading <- function(betamat,imat,smat,neighborfnc){
#    #betamat: matrix of betaij for all patche combinations
#    N <- nrow(imat)
#    newmat <- matrix(rep(0,length(imat)),N,N)
#    for (col in 1:N){
#        for (row in 1:N){
#            neighbours <- neighbourfnc(row,col)
#            nsum <- 0
#            for (n in neighbours){nsum <- nsum + betamat[n]*imat[n]}
#            newmat[row,col] <- smat[row,col]*nsum
#        }
#    }
#    return(newmat)
#}
initPlot <- function(patch,rows,cols, xlab='Time', ylab='Population Fraction', ...){
    left <- (patch > cols*(rows -1))
    bottom <- ((patch %% cols) == 1)
    middlex <- as.integer(1:rows)
    middley <- as.integer(1:cols)
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
        legend("right",
               legend=c('S','I','R'),
               col=c('green','red','black'),
               box.lty=0,
               lty=1,
               cex=0.8)
    }
}
multiPlotSoln <- function(rows,cols,soln,...){
    patches <- rows*cols
    scol <- 'green'
    icol <- 'red'
    rcol <- 'black'
    wcol <- 'blue'
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
        lines(w[,patch],type='l',col=wcol,lty=1)
    }
}
