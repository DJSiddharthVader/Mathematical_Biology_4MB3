##Single Patch SIRW analysis
library(deSolve)
##Initial Parameters
pars=c(mu = 1/15370, # avg lifespan 58 yrs 
       beta_i = 0.016, #transmission r. S between I
       gamma = 1/7, #rate of recovery (days)
       sigma = 1/14, #rate of water sanitation (or 1/lifetime bacteria)
       beta_w = 0.5, #transmission r. I between W
       xi = 1/14, #shedding rate I into W (force it eaual to sigma?)
       kappa = 0.15, #bacteria concentration that gives 50% of infection
       alpha = 0 #death rate by cholera
)
##Initial Conditions for Phase Portraits##
tmax <- 100
Z0 <-  0 ##recovered class initial
W0 <-  0.5##initial bacterial concentration in the water 

##VECTOR FIELD SIRW BASE MODEL##
  SIRW.vector.field <- function(t, vars, 
                                parms=c(mu = 1/25550, 
                                        beta_i = 0.0016, 
                                        gamma = 1/7, 
                                        sigma = 1/14,
                                        beta_w = 0.5,
                                        xi = 1/14, 
                                        kappa = 0.5, 
                                        alpha = 0 
                                )) {
    with(as.list(c(parms, vars)), {
      dx <- mu*(x+y+z) - mu*x - beta_i*x*y - beta_w*x*w*(1/(kappa+w)) # dS/dt
      dy <- beta_i*y*x + beta_w*x*w - gamma*y - mu*y - alpha*y # dI/dt
      dz <- gamma*y - mu*z # dR/dt
      dw <- xi*y - sigma*w #dW/dt
        vec.fld <- c(dx=dx, dy=dy, dz=dz, dw = dw)
      return(list(vec.fld))
    })
  }
  
##PLOTTING FN for S(x) vs I(y)
  plot.SIRW <- function(ic=c(x=1,y=0,z=0, w=0), tmax=1,
                      times=seq(0,tmax,by=tmax/5000),
                      func, parms, ... ) {
    St <- ode(ic, times, func, parms) 
    lines(St[,"x"], St[,"y"], col="black", lwd=1.5, ... )
  }

#PHASE PORTRAIT for SIRW##
  phase.portrait <- function(S0=seq(0,0.9,by=0.1)[-1],
                             I0=seq(0.9,0, by=-0.1)[-10],
                             Z0=0, W0=0.1){
    plot(x=0, y=0, type = "n", xlim = c(0,1), ylim = c(0,1), #empty subplot
         xlab = "S", ylab= "I", las=1, bty="l")
    lines(x=c(1,0),y=c(0,1),col="grey",lty="dotted") #S + I = 1
    for (j in 1:length(S0)) {
      plot.SIRW(ic=c(x=S0[j],y=I0[j],z=Z0, w = W0), tmax=tmax,
                func=SIRW.vector.field,
                parms = pars 
      )
    }  
  }

###TIME PLOTS##
  #The default initial conditions use a fully susceptible population with the only
  #the W class being "infected"
  
##plot W class over time##  
  plot.Wt <- function(ic=c(x=1,y=0,z=0, w=W0), tmax=100,
                      times=seq(0,tmax,by=tmax/500),
                      func=SIRW.vector.field, parms=pars, ... ) { 
    plot(x=0, y=0, type = "n",xlim=c(0, tmax), ylim=c(0, 1),
         xlab = "Time", ylab= "Bacteria Concentration (W)")
    St <- ode(ic, times, func, parms) 
    lines(times, St[,"w"], col="blue", lwd=1.5, ... )
  } 
  
##Plot I class over time##
  plot.It <- function(ic=c(x=1,y=0,z=0, w=W0), tmax=100,
                     times=seq(0,tmax,by=tmax/500),
                     func=SIRW.vector.field, parms=pars, ... ) { 
    plot(x=0, y=0, type = "n",xlim=c(0, tmax), ylim=c(0, 1),
         xlab = "Time", ylab= "I")
    St <- ode(ic, times, func, parms) 
    lines(times, St[,"y"], col="red", lwd=1.5, ... )
  }

##Plot S class over time##
plot.St <-  function(ic=c(x=1,y=0,z=0, w=W0), tmax=100,
                     times=seq(0,tmax,by=tmax/500),
                     func=SIRW.vector.field, parms=pars, ... ) { 
  plot(x=0, y=0, type = "n",xlim=c(0, tmax), ylim=c(0, 1),
       xlab = "Time", ylab= "I")
  St <- ode(ic, times, func, parms) 
  lines(times, St[,"x"], col="red", lwd=1.5, ... )
}

#plots SIW over time#
plot.dynamic <- function(ic=c(x=0.9,y=0.1,z=0, w=W0), tmax=100,
                         times=seq(0,tmax,by=tmax/500),
                         func=SIRW.vector.field, parms=pars, ... ) { 
  plot(x=0, y=0, type = "n",xlim=c(0, tmax), ylim=c(0, 1),
       xlab = "Time", ylab="Proportion")
  St <- ode(ic, times, func, parms) 
  lines(times, St[,"y"], col="red", lwd=1.5, ... )
  lines(times, St[,"x"], col="green", lwd=1.5, ... )
  lines(times, St[,"w"], col="blue", lwd=1.5, ... )
  legend("topright", legend = c("S", "I", "W"), lty=1, col=c("green", "red", "blue"))
}

##Reproductive Ratio
R_0 <- function(parms=c(mu = 1/25550,
                               beta_i = 0.0016, 
                               gamma = 1/7, 
                               sigma = 1/14, 
                               beta_w = 0.5, 
                               xi = 1/14, 
                               kappa = 0.5, 
                               alpha = 0 
)) {
  with(as.list(parms), {
    return((beta_i+(beta_w/kappa))/(mu + gamma + alpha)) 
 })
}





  