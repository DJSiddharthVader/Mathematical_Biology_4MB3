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
#Base Models
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
        rdisperse <- tran.2D(r,dx=dx,dy=dy,D.x=Dsr,D.y=Dsr)$dC
        wdisperse <- tran.2D(w,dx=dx,dy=dy,D.x=Dw,D.y=Dw)$dC
        #derivatives
        ds <- mu - mu*s - person_infect - water_infect + sdisperse
        di <- person_infect + water_infect - i*gamma -i*mu -i*alpha
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
        rdisperse <- tran.2D(r,dx=dx,dy=dy,D.x=Dsr,D.y=Dsr)$dC
        wdisperse <- tran.2D(w,dx=dx,dy=dy,D.x=Dw,D.y=Dw)$dC
        #derivatives
        ds <- mu - mu*s - person_infect - water_infect + sdisperse
        di <- person_infect + water_infect - i*gamma -i*mu -i*alpha
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
        rdisperse <- tran.2D(r,dx=dx,dy=dy,D.x=Dsr,D.y=Dsr)$dC
        wdisperse <- tran.2D(w,dx=dx,dy=dy,D.x=Dw,D.y=Dw)$dC
        #derivatives
        ds <- mu - mu*s - person_infect - water_infect - nu*s + sdisperse
        di <- person_infect + water_infect - i*gamma -i*mu -i*alpha
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
        rdisperse <- tran.2D(r,dx=dx,dy=dy,D.x=Dsr,D.y=Dsr)$dC
        wdisperse <- tran.2D(w,dx=dx,dy=dy,D.x=Dw,D.y=Dw)$dC
        #derivatives
        ds <- mu - mu*s - person_infect - water_infect + sdisperse
        di <- person_infect + water_infect - i*gamma -i*mu -i*alpha -eta*i
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
compareTreatments <- function(ic,params,tmax,steps,...){
    #Base
    base <- solveMultiModel(multiWaterModel,ic,params,tmax=tmax,steps=steps)
    bi <- base[,(rows*rows+2):(rows*rows*2+1)]
    bim <- apply(bi,1,mean)
    #Sanitation
    san <- solveMultiModel(sanitationModel,ic,params,tmax=tmax,steps=steps)
    si <- san[,(rows*rows+2):(rows*rows*2+1)]
    sim <- apply(si,1,mean)
    #Antibiotics
    anti <- solveMultiModel(antibioticModel,ic,params,tmax=tmax,steps=steps)
    ai <- anti[,(rows*rows+2):(rows*rows*2+1)]
    aim <- apply(ai,1,mean)
    #Vaccination
    vacc <- solveMultiModel(vaccinationModel,ic,params,tmax=tmax,steps=steps)
    vi <- vacc[,(rows*rows+2):(rows*rows*2+1)]
    vim <- apply(vi,1,mean)
    #Legend entries
    sentry <- paste('Sanitation =',params$rho)
    ventry <- paste('Vaccination =',params$nu)
    aentry <- paste('Antibiotics =',params$eta)
    lentries <- c("Base",sentry,ventry,aentry)
    #Plotting
    overlaymax <- 0.1 #ymax on overlay graph
    treatcol = c("red","cornflowerblue", "orange", "brown") #for overlay graph
    plot(x=0, y=0, type = "n",xlim=c(0, tmax), ylim=c(0, overlaymax), xaxs="i", yaxs="i",xlab = "Time (days)", ylab="Proportion (I)")
    lines(bim,col="red",lty=1, lwd=4,)
    lines(sim,col="cornflowerblue",lty=1, lwd=4,)
    lines(vim,col="brown",lty=1, lwd=4,)
    lines(aim,col="orange",lty=1, lwd=4,)
    legend("topright",lty=c(1,1,1,1), col = treatcol,legend=lentries)
}
#Base Plotting
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
main3 <- function(){
    set.seed(9)
    rows <- 80
    cols <- rows
    tmax <- 150
    steps <- 300
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
    soln <- solveMultiModel(vaccinationModel,ic,params,tmax=tmax,steps=steps)
#    gifplot(rows,cols,soln,tmax)
#    cmd <- "mpv filtering.gif --loop=inf"
#    system(cmd)
    return(soln[,(80*80+2):(80*80*2+1)])
}

#High/Low Models
initHLModel <- function(rows,cols,i0,w0){
    p = rows*cols
    l <- matrix(rep(0,p),rows,cols)
    l[sample(1:length(l), as.integer(0.08*rows*cols)+1)] <- i0
    h <- matrix(rep(0,p),rows,cols)
    h[sample(1:length(h), as.integer(0.04*rows*cols)+1)] <- i0
    s <- matrix(rep(1,p),rows,cols) - l - h
    r <- matrix(rep(0,p),rows,cols)
    w <- matrix(rep(0,p),rows,cols)
    w[sample(1:length(w), as.integer(0.05*rows*cols)+1)] <- w0
    return(c(s,l,h,r,w))
}
highLowWaterModel <- function(t,state,params){
    with(as.list(c(params)),{
        #compartment matrices
        patches <- rows*cols
        s <- matrix(state[(0*patches+1):(1*patches+0)],rows,cols)
        l <- matrix(state[(1*patches+1):(2*patches+0)],rows,cols)
        h <- matrix(state[(2*patches+1):(3*patches+0)],rows,cols)
        r <- matrix(state[(3*patches+1):(4*patches+0)],rows,cols)
        w <- matrix(state[(4*patches+1):(5*patches+0)],rows,cols)
        #disperse terms
        sdisperse <- tran.2D(s,dx=dx,dy=dy,D.x=Dsr,D.y=Dsr)$dC
        ldisperse <- tran.2D(l,dx=dx,dy=dy,D.x=Dl,D.y=Dl)$dC
        rdisperse <- tran.2D(r,dx=dx,dy=dy,D.x=Dsr,D.y=Dsr)$dC
        wdisperse <- tran.2D(w,dx=dx,dy=dy,D.x=Dw,D.y=Dw)$dC
        #terms
        lowinf <- beta_l*l*s
        highinf  <- beta_h*h*s
        waterinf <- beta_w*w*s
        #derivatives
        ds <- mu - mu*s - lowinf - highinf - waterinf + sdisperse #dS/dt
        dl <- lowinf + highinf + waterinf - l*(delta + mu + alpha_l) + ldisperse #dI_l/dt
        dh <- delta*l - h*(gamma + mu + alpha_h) #dI_h/dt
        dr <- gamma*h - mu*r + rdisperse #dR/dt
        dw <- xi_l*l + xi_h*h - sigma*w + wdisperse #dW/dt
        vec.fld <- list(c(ds,dl,dh,dr,dw))
        return(vec.fld)
    })
}
sanitationHL <- function(t,state,params){
    with(as.list(c(params)),{
        #compartment matrices
        patches <- rows*cols
        s <- matrix(state[(0*patches+1):(1*patches+0)],rows,cols)
        l <- matrix(state[(1*patches+1):(2*patches+0)],rows,cols)
        h <- matrix(state[(2*patches+1):(3*patches+0)],rows,cols)
        r <- matrix(state[(3*patches+1):(4*patches+0)],rows,cols)
        w <- matrix(state[(4*patches+1):(5*patches+0)],rows,cols)
        #disperse terms
        sdisperse <- tran.2D(s,dx=dx,dy=dy,D.x=Dsr,D.y=Dsr)$dC
        ldisperse <- tran.2D(l,dx=dx,dy=dy,D.x=Dl,D.y=Dl)$dC
        rdisperse <- tran.2D(r,dx=dx,dy=dy,D.x=Dsr,D.y=Dsr)$dC
        wdisperse <- tran.2D(w,dx=dx,dy=dy,D.x=Dw,D.y=Dw)$dC
        #terms
        lowinf <- beta_l*l*s
        highinf  <- beta_h*h*s
        waterinf <- beta_w*w*s
        #derivatives
        ds <- mu - mu*s - lowinf - highinf - waterinf + sdisperse #dS/dt
        dl <- lowinf + highinf + waterinf - l*(delta + mu + alpha_l) + ldisperse #dI_l/dt
        dh <- delta*l - h*(gamma + mu + alpha_h) #dI_h/dt
        dr <- gamma*h - mu*r + rdisperse #dR/dt
        dw <- xi_l*l + xi_h*h - sigma*w - rho*w + wdisperse #dW/dt
        vec.fld <- list(c(ds,dl,dh,dr,dw))
        return(vec.fld)
    })
}
antibioticHL <- function(t,state,params){
    with(as.list(c(params)),{
        #compartment matrices
        patches <- rows*cols
        s <- matrix(state[(0*patches+1):(1*patches+0)],rows,cols)
        l <- matrix(state[(1*patches+1):(2*patches+0)],rows,cols)
        h <- matrix(state[(2*patches+1):(3*patches+0)],rows,cols)
        r <- matrix(state[(3*patches+1):(4*patches+0)],rows,cols)
        w <- matrix(state[(4*patches+1):(5*patches+0)],rows,cols)
        #disperse terms
        sdisperse <- tran.2D(s,dx=dx,dy=dy,D.x=Dsr,D.y=Dsr)$dC
        ldisperse <- tran.2D(l,dx=dx,dy=dy,D.x=Dl,D.y=Dl)$dC
        rdisperse <- tran.2D(r,dx=dx,dy=dy,D.x=Dsr,D.y=Dsr)$dC
        wdisperse <- tran.2D(w,dx=dx,dy=dy,D.x=Dw,D.y=Dw)$dC
        #terms
        lowinf <- beta_l*l*s
        highinf  <- beta_h*h*s
        waterinf <- beta_w*w*s
        #derivatives
        ds <- mu - mu*s - lowinf - highinf - waterinf + sdisperse #dS/dt
        dl <- lowinf + highinf + waterinf - l*(delta + mu + alpha_l) - eta*l + ldisperse #dI_l/dt
        dh <- delta*l - h*(gamma + mu + alpha_h) #dI_h/dt
        dr <- gamma*h - mu*r + rdisperse #dR/dt
        dw <- xi_l*l + xi_h*h - sigma*w + wdisperse #dW/dt
        vec.fld <- list(c(ds,dl,dh,dr,dw))
        return(vec.fld)
    })
}
vaccinationHL <- function(t,state,params){
    with(as.list(c(params)),{
        #compartment matrices
        patches <- rows*cols
        s <- matrix(state[(0*patches+1):(1*patches+0)],rows,cols)
        l <- matrix(state[(1*patches+1):(2*patches+0)],rows,cols)
        h <- matrix(state[(2*patches+1):(3*patches+0)],rows,cols)
        r <- matrix(state[(3*patches+1):(4*patches+0)],rows,cols)
        w <- matrix(state[(4*patches+1):(5*patches+0)],rows,cols)
        #disperse terms
        sdisperse <- tran.2D(s,dx=dx,dy=dy,D.x=Dsr,D.y=Dsr)$dC
        ldisperse <- tran.2D(l,dx=dx,dy=dy,D.x=Dl,D.y=Dl)$dC
        rdisperse <- tran.2D(r,dx=dx,dy=dy,D.x=Dsr,D.y=Dsr)$dC
        wdisperse <- tran.2D(w,dx=dx,dy=dy,D.x=Dw,D.y=Dw)$dC
        #terms
        lowinf <- beta_l*l*s
        highinf  <- beta_h*h*s
        waterinf <- beta_w*w*s
        #derivatives
        ds <- mu - mu*s - lowinf - highinf - waterinf - nu*s + sdisperse #dS/dt
        dl <- lowinf + highinf + waterinf - l*(delta+ mu + alpha_l) + ldisperse #dI_l/dt
        dh <- delta*l - h*(gamma + mu + alpha_h) #dI_h/dt
        dr <- gamma*h - mu*r + nu*s + rdisperse #dR/dt
        dw <- xi_l*l + xi_h*h - sigma*w + wdisperse #dW/dt
        vec.fld <- list(c(ds,dl,dh,dr,dw))
        return(vec.fld)
    })
}
solveHLModel <- function(func,ic,params,tmax=1,steps=500,lrw=9999999){
    with(as.list(c(ic,params)),{
        ivec <- initHLModel(rows,cols,i0,w0)
        times <- seq(0,tmax,by=tmax/steps)
        soln <- ode.2D(y=ivec,
                       times=times,
                       func=func,
                       parms=params,
                       nspec=5,
                       lrw=lrw,
                       dimens=c(rows,cols))
        return(soln)
    })
}
compareHLTreatments <- function(ic,params,tmax,steps...){
    #Base
    base <- solveHLModel(highLowWaterModel,ic,params,tmax=tmax,steps=steps)
    bl <- base[,(1*rows*rows+2):(2*rows*rows+1)]
    blm <- apply(bl,1,mean)
    bh <- base[,(2*rows*rows+2):(3*rows*rows+1)]
    bhm <- apply(bh,1,mean)
    #Sanitation
    sanitation <- solveHLModel(sanitationHL,ic,params,tmax=tmax,steps=steps)
    sl <- sanitation[,(1*rows*rows+2):(2*rows*rows+1)]
    slm <- apply(sl,1,mean)
    sh <- sanitation[,(2*rows*rows+2):(3*rows*rows+1)]
    shm <- apply(sh,1,mean)
    #Antibiotics
    antibiotics <- solveHLModel(antibioticHL,ic,params,tmax=tmax,steps=steps)
    al <- antibiotics[,(1*rows*rows+2):(2*rows*rows+1)]
    alm <- apply(al,1,mean)
    ah <- antibiotics[,(2*rows*rows+2):(3*rows*rows+1)]
    ahm <- apply(ah,1,mean)
    #Vaccination
    vaccination <- solveHLModel(vaccinationHL,ic,params,tmax=tmax,steps=steps)
    vl <- vaccination[,(1*rows*rows+2):(2*rows*rows+1)]
    vlm <- apply(vl,1,mean)
    vh <- vaccination[,(2*rows*rows+2):(3*rows*rows+1)]
    vhm <- apply(vh,1,mean)
    #legend entries
    sentry <- paste('Sanitation =',params$rho)
    ventry <- paste('Vaccination =',params$nu)
    aentry <- paste('Antibiotics =',params$eta)
    lentries <- c("Base",sentry,ventry,aentry)
    #Plotting
    overlaymax <- 0.1 #ymax on overlay graph
    treatcol = c("red","cornflowerblue", "orange", "brown")
    par(mfrow=c(2,1))
    #Low Plotting
    plot(x=0, y=0,type="n",xlim=c(0, tmax),ylim=c(0, overlaymax),xaxs="i",yaxs="i", xlab="Time (days)",ylab="Proportion (I)",main="Low Infectivity")
    lines(blm,col="red",lty=1, lwd=4,)
    lines(slm,col="cornflowerblue",lty=1, lwd=4,)
    lines(alm,col="orange",lty=1, lwd=4,)
    lines(vlm,col="brown",lty=1, lwd=4,)
    legend("topright", lty=c(1,1,1,1), col = treatcol,lentries)
    #High Plotting
    plot(x=0, y=0,type="n",xlim=c(0, tmax),ylim=c(0, overlaymax),xaxs="i",yaxs="i", xlab="Time (days)",ylab="Proportion (I)",main="High Infectivity")
    lines(bhm,col="red",lty=1, lwd=4,)
    lines(shm,col="cornflowerblue",lty=1, lwd=4,)
    lines(vhm,col="brown",lty=1, lwd=4,)
    lines(ahm,col="orange",lty=1, lwd=4,)
    legend("topright", lty=c(1,1,1,1), col = treatcol,lentries)
}


