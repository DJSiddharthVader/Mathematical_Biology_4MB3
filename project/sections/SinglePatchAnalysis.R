##Single Patch SIRW analysis
library(deSolve)
library(cholera)

#init time series dataset (with cholera)
#timeSeries <- timeSeries()
#plot(timeSeries())

  ##Initiate (basic with vital dynamics)##
# Assumptions
## 1.4 days before symptoms show (1/para. E class)

  SIRW.vector.field <- function(t, vars, 
                                parms=c(mu = 1/25550, 
                                        beta_i = 1, #transmission r. S between I
                                        gamma = 1/7, #rate of recovery (days)
                                        sigma = 1/14, #rate of water sanitation (or 1/lifetime bacteria)
                                        beta_w = 1, #transmission r. I between W
                                        alpha = 1 #death rate by cholera
                                )) {
    with(as.list(c(parms, vars)), {
      dx <- mu*(x+y+z) - mu*x - beta_i*x*y - beta_w*x*w # dS/dt
      dy <- beta_i*y*x + beta_w*x*w - gamma*y - mu*y - alpha*y # dI/dt
      dz <- gamma*y - mu*z # dR/dt
        dw <- beta_w*y - sigma*w #dW/dt
        vec.fld <- c(dx=dx, dy=dy, dz=dz, dw = dw)
      return(list(vec.fld))
    })
  }
  plot.SIRW <- function(ic=c(x=1,y=0,z=0), tmax=1,
                      times=seq(0,tmax,by=tmax/5000),
                      func, parms, ... ) {
    St <- ode(ic, times, func, parms) 
    lines(St[,"x"], St[,"y"], col="blue", lwd=1.5, ... )
  }
  pars=c(mu = 1/25550, 
          beta_i = 1, #transmission r. S between I
          gamma = 1/7, #rate of recovery (days)
          sigma = 1/14, #rate of water sanitation (or 1/lifetime bacteria)
          beta_w = 1, #transmission r. I between W
          alpha = 1 #death rate by cholera
  )
  #intial conditions
  tmax <- 1000
  S0 <- seq(0,0.9,by=0.1)[-1]
  I0 <- seq(0.9,0, by=-0.1)[-10]
  par(mfrow())
  plot(x=0, y=0, type = "n", xlim = c(0,1), ylim = c(0,1), #empty subplot
         xlab = "Susceptible (S)", ylab= "Infected (I)", main=title)
  lines(x=c(1,0),y=c(0,1),col="grey",lty="dotted")
    for (j in 1:length(S0)) {
      plot.SIRW(ic=c(x=S0[j],y=I0[j],z=0, w = 0.5), tmax=tmax,
              func=SIRW.vector.field,
              parms = pars 
      )
    }
  
  
  