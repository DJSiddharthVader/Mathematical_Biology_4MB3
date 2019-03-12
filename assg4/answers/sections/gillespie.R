#!/usr/bin/env Rscript
library("deSolve")
SI.Gillespie <- function(beta,N,I0,tmax){
    #i > 0?
    if (I0 == 0){
        stop('initial infected must be > 0 for epidemic to progress')
    }
    #set vars
    i = I0
    incidences = c(i)
    t = 0
    times = c(0)
    s = N - i
    #main loop
    while ( t < tmax ){
        #time till next event
        rnd <- runif(1,0,1)
        dt = (1/beta)*(log((1/(1-rnd))))
        t = t + dt
        #only event is infection
        news = s - (beta/N)*s*i
        newi = i + (beta/N)*s*i
        #update vars
        s = news
        i = newi
        incidences = c(incidences,i)
        times = c(times,t)
    }
    return(list(times,incidences))
}
SI.vector.field <- function(t,vars,params){
    with(as.list(c(params, vars)), {
        dx <- -Beta*x*y
        dy <- Beta*x*y
        vec.fld <- c(dx=dx,dy=dy)
        return(list(vec.fld))
    })
}
SI.Deterministic <- function(beta,N,I0,tmax,timepts=500){
    #i > 0?
    if (I0 == 0){
        stop('initial infected must be > 0 for epidemic to progress')
    }
    #set vars
    ic <- c(x=N - I0,y=I0)
    times <- seq(0,tmax,by=tmax/timepts)
    params <- c(Beta=beta/N)
    #solve ode
    incidences <- ode(ic,times,SI.vector.field,params)
    return(list(times,incidences[,"y"]))
}
multipanel <-  function(realizations,beta,ns,I0,tmax,colors=c('blue')){
    # set color list
    if (length(colors) < realizations){
        colors <- rep(colors,length.out=realizations)
    }
    par(mfrow = c(2,2))
    for (popsize in ns){
        #plot initial stochatic
        result <- SI.Gillespie(beta,popsize,I0,tmax)
        plot(result[[1]],
             result[[2]],
             col=colors[1],
             type="l",
             xlab="Time (t)",
             ylab="Incidence (I(t))",
             main=paste("N =",popsize))
        #plot 30 stochastic realizations
        for (i in seq(0,realizations-1)){
            result <- SI.Gillespie(beta,popsize,I0,tmax)
            lines(result[[1]],result[[2]],col=colors[i],type="l",lty=3)
        }
        #deterministic solution
        result <- SI.Deterministic(beta,popsize,I0,tmax)
        lines(result[[1]],result[[2]],col="red",type="l")
        #legend
        legend("right",
               legend="Deterministic",
               col="red",
               box.lty=0,
               lty=1:1,
               cex=0.8)
    }
}
