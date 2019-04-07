##Single Patch SIRW analysis##
##############################
library(deSolve)
##Initial Conditions for Phase Portraits##
timemax <- 1200
Z0 <-  0 
W0 <-  0
int <- "low"
y0 <- 0.001
x0 <- 1-y0
xmax <-  150 #xlim of numerical simulations
ymax <- 1 #ylim of numerical simulations
overlaymax <- 0.1 #ymax on overlay graph
treatcol = c("red","cornflowerblue", "orange", "brown") #for overlay graph
##############
##Parameters##
##############
pars=c(mu = 1/21170, # avg lifespan 58 yrs == 21170
       beta_i = 0.1, #transmission r. S between I
       gamma = 1/8, #rate of recovery (days)
       sigma = 1/22, #rate of water removal
       beta_w = 0.3, #transmission r. I between W
       xi = 1/22, #shedding rate I into W (force it equal to sigma)
       alpha = 0, #death rate by cholera
       
       #parameters for intensity model
       beta_l = 0.016, #transmission r. S between I_l
       beta_h = 0.016, #transmission r. S between I_h
       delta = 1/3, #rate of increasing symptom severity
       xi_l = 1/14, #shedding rate I_l into W; maybe equal to sigma
       xi_h = 1/7,
       alpha_l = 0.01, #death rate by cholera with low intensity symptoms
       alpha_h = 0.05, #death rate by cholera with high int. symptoms
       
       #parameters for water treatment [day^-1]
       rho = 1/20, 
       #parameters for vaccination [day^-1]
       nu = 1/200,
       #parameters for antibiotics [day^-1]
       eta = 1/60
)
#################
##VECTOR FIELDS##
#################
#VECTOR FIELD SIRW BASE MODEL
  SIRW.vector.field <- function(t, vars, 
                                parms=pars) {
    with(as.list(c(parms, vars)), {
      dx <- mu - mu*x - beta_i*x*y - beta_w*x*w # dS/dt
      dy <- beta_i*x*y + beta_w*x*w - gamma*y - mu*y - alpha*y # dI/dt
      dz <- gamma*y - mu*z # dR/dt
      dw <- xi*y - sigma*w #dW/dt
        vec.fld <- c(dx=dx, dy=dy, dz=dz, dw = dw)
      return(list(vec.fld))
    })
  }

#VECTOR FIELD treatment 1 MODEL
  treat1.vector.field <- function(t, vars, 
                                  parms=pars) {
    with(as.list(c(parms, vars)), {
      dx <- mu - mu*x - beta_i*x*y - beta_w*x*w # dS/dt
      dy <- beta_i*y*x + beta_w*x*w - gamma*y - mu*y - alpha*y # dI/dt
      dz <- gamma*y - mu*z # dR/dt
      dw <- xi*y - sigma*w - rho*w #dW/dt
      vec.fld <- c(dx=dx, dy=dy, dz=dz, dw = dw)
      return(list(vec.fld))
    })
  }
#treatment 2 (vaccine) model
  treat2.vector.field <- function(t, vars, 
                                  parms=pars) {
    with(as.list(c(parms, vars)), {
      dx <- mu - mu*x - beta_i*x*y - beta_w*x*w - nu*x# dS/dt
      dy <- beta_i*y*x + beta_w*x*w - gamma*y - mu*y - alpha*y # dI/dt
      dz <- gamma*y - mu*z + nu*x # dR/dt
      dw <- xi*y - sigma*w#dW/dt
      vec.fld <- c(dx=dx, dy=dy, dz=dz, dw = dw)
      return(list(vec.fld))
    })
  }
#treatment 3 (antibiotic) model
  treat3.vector.field <- function(t, vars, 
                                  parms=pars) {
    with(as.list(c(parms, vars)), {
      dx <- mu - mu*x - beta_i*x*y - beta_w*x*w# dS/dt
      dy <- beta_i*y*x + beta_w*x*w - gamma*y - mu*y - alpha*y - eta*y # dI/dt
      dz <- gamma*y - mu*z + eta*y # dR/dt
      dw <- xi*y - sigma*w#dW/dt
      vec.fld <- c(dx=dx, dy=dy, dz=dz, dw = dw)
      return(list(vec.fld))
    })
  }
  
  #SEVERITY MODEL VECTOR FIELD
  severity.vector.field <- function(t, vars, 
                                    parms=pars) {
    with(as.list(c(parms, vars)), {
      dx <- mu*(x+l+h+z) - mu*x - beta_l*x*l - beta_h*x*h - beta_w*x*w # dS/dt
      dl <- x*(beta_l*l + beta_h*h) +beta_w*x*w- gamma*l - mu*l - alpha_l*l # dI_l/dt
      dh <- delta*l - h*(gamma + mu+ alpha_h) #dI_h/dt
      dz <- gamma*h - mu*z # dR/dt
      dw <- xi_l*l + xi_h*h - sigma*w #dW/dt
      vec.fld <- c(dx=dx, dl=dl, dh=dh,  dz=dz, dw = dw)
      return(list(vec.fld))
    })
  }

  #DUMPING
  dump.vector.field <- function(t, vars, 
                                parms=pars) {
    with(as.list(c(parms, vars)), {
      dx <- mu - mu*x - beta_i*x*y - beta_w*x*w # dS/dt
      dy <- beta_i*x*y + beta_w*x*w - gamma*y - mu*y - alpha*y # dI/dt
      dz <- gamma*y - mu*z # dR/dt
      dw <- xi*y - 1/500*sigma*w #dW/dt (increased "dumping", set decay rate really low)
      vec.fld <- c(dx=dx, dy=dy, dz=dz, dw = dw)
      return(list(vec.fld))
    })
  }
  
################
##PLOTTING FNs##
################
  plot.SIRW <- function(ic=c(x=1,y=0,z=0, w=0), tmax=1,
                          times=seq(0,tmax,by=tmax/500),
                          func, parms, ... ) {
    St <- ode(ic, times, func, parms) 
    lines(St[,"x"], St[,"y"], col="black", lwd=2, ... )
  }
  plot.SIIRW <- function(which="low", ic=c(x=1,l=0, h=0.1, z=0), tmax=1,
                         times=seq(0,tmax,by=tmax/500),
                         func, parms, ... ) {
    St <- ode(ic, times, func, parms)
    if (int=="low"){
      lines(St[,"x"], St[,"l"], col="black", lwd=1.5, ... )
    }
    else if (int=="high"){
      lines(St[,"x"], St[,"h"], col="black", lwd=1.5, ... )
    }
  }
  #base model  
  plot.base <- function(ic=c(x=x0,y=y0,z=0, w=W0), tmax=timemax,
                        times=seq(0,tmax,by=1),
                        func=SIRW.vector.field, parms=pars, ... ) { 
    plot(x=0, y=0, type = "n",xlim=c(0, xmax), ylim=c(0, ymax),
         xlab = "Time (days)", ylab="Proportion", main = "Base Model")
    St <- ode(ic, times, func, parms) 
    lines(times, St[,"x"], col="green",lty=1, lwd=2, ... )
    lines(times, St[,"y"], col="red",lty=1, lwd=2, ... )
    lines(times, St[,"z"], col="black",lty=1, lwd=2, ... )
    #lines(times, St[,"w"], col="blue",lty=1, lwd=2, ... )
    legend("topright", legend = c("S", "I","R"), lty=c(1,1,1,1), col=c("green", "red", "black"))
  }
  #severity bodel
  plot.severity <- function(ic=c(x=x0,l=y0,h=0,z=0, w=W0), tmax=timemax,
                            times=seq(0,tmax,by=tmax/500),
                            func=severity.vector.field, parms=pars, ... ) { 
    plot(x=0, y=0, type = "n",xlim=c(0, xmax), ylim=c(0, ymax),
         xlab = "Time (days)", ylab="Proportion", main = "Severity Model with Vital Dynamics")
    St <- ode(ic, times, func, parms) 
    lines(times, St[,"l"], col="red",lty=1, lwd=1.5, ... )
    lines(times, St[,"h"], col="red", lty=2, lwd=1.5, ... )
    #lines(times, St[,"w"], col="blue", lty=2, lwd=1.5, ... )
    legend("topright", legend = c( "I_L","I_H", "W"), lty=1, col=c("red", "orange",  "blue"))
  }
  #treat1 model
  plot.treat1 <- function(ic=c(x=x0,y=y0,z=0, w=W0), tmax=timemax,
                          times=seq(0,tmax,by=1),
                          func=treat1.vector.field, parms=pars, ... ) { 
    plot(x=0, y=0, type = "n",xlim=c(0, xmax), ylim=c(0, ymax),
         xlab = "Time (days)", ylab="Proportion", main = "Treatment 1: Water Sanitation")
    St <- ode(ic, times, func, parms) 
    lines(times, St[,"x"], col="green",lty=1, lwd=2, ... )
    lines(times, St[,"y"], col="red",lty=1, lwd=2, ... )
    lines(times, St[,"z"], col="black",lty=1, lwd=2, ... )
    #lines(times, St[,"w"], col="blue",lty=1, lwd=2, ... )
    legend("topright", legend = c("S", "I","R"), lty=c(1,1,1,1), col=c("green", "red", "black"))
  }
  #treat2 model
  plot.treat2 <- function(ic=c(x=x0,y=y0,z=0, w=W0), tmax=timemax,
                          times=seq(0,tmax,by=1),
                          func=treat2.vector.field, parms=pars, ... ) { 
    plot(x=0, y=0, type = "n",xlim=c(0, xmax), ylim=c(0, ymax),
         xlab = "Time (days)", ylab="Proportion", main = "Treatment 2: Vaccination")
    St <- ode(ic, times, func, parms) 
    lines(times, St[,"x"], col="green",lty=1, lwd=2, ... )
    lines(times, St[,"y"], col="red",lty=1, lwd=2, ... )
    lines(times, St[,"z"], col="black",lty=1, lwd=2, ... )
    #lines(times, St[,"w"], col="blue",lty=1, lwd=2, ... )
    legend("topright", legend = c("S", "I","R"), lty=c(1,1,1,1), col=c("green", "red", "black"))
  }
  #treat3 model
  plot.treat3 <- function(ic=c(x=x0,y=y0,z=0, w=W0), tmax=timemax,
                          times=seq(0,tmax,by=1),
                          func=treat3.vector.field, parms=pars, ... ) { 
    plot(x=0, y=0, type = "n",xlim=c(0, xmax), ylim=c(0, ymax),
         xlab = "Time (days)", ylab="Proportion", main = "Treatment 3: Antibiotics")
    St <- ode(ic, times, func, parms) 
    lines(times, St[,"x"], col="green",lty=1, lwd=2, ... )
    lines(times, St[,"y"], col="red",lty=1, lwd=2, ... )
    lines(times, St[,"z"], col="black",lty=1, lwd=2, ... )
    #lines(times, St[,"w"], col="blue",lty=1, lwd=2, ... )
    legend("topright", legend = c("S", "I","R"), lty=c(1,1,1,1), col=c("green", "red", "black"))
  }
  
  plot.dump <- function(ic=c(x=x0,y=y0,z=0, w=W0), tmax=timemax,
                        times=seq(0,tmax,by=1),
                        func=dump.vector.field, parms=pars, ... ) { 
    plot(x=0, y=0, type = "n",xlim=c(0, xmax), ylim=c(0, ymax),
         xlab = "Time (days)", ylab="Proportion", main = "Base Model: 19th Century")
    St <- ode(ic, times, func, parms) 
    lines(times, St[,"x"], col="green",lty=1, lwd=2, ... )
    lines(times, St[,"y"], col="red",lty=1, lwd=2, ... )
    lines(times, St[,"z"], col="black",lty=1, lwd=2, ... )
    #lines(times, St[,"w"], col="blue",lty=1, lwd=2, ... )
    legend("topright", legend = c("S", "I","R"), lty=c(1,1,1,1), col=c("green", "red", "black"))
  }
  
  #overlay of prevalence btwn treatments
 plot.overlay <- function(ic=c(x=x0,y=y0,z=0, w=W0), tmax=timemax,
                                  times=seq(0,tmax,by=1),
                                  func=func, parms=pars, ... ) { 
   Stb <- ode(ic, times, func=SIRW.vector.field, parms)
   St1 <- ode(ic, times, func=treat1.vector.field, parms)
   St2 <- ode(ic, times, func=treat2.vector.field, parms)
   St3 <- ode(ic, times, func=treat3.vector.field, parms)
   plot(x=0, y=0, type = "n",xlim=c(0, xmax), ylim=c(0, overlaymax), xaxs="i", yaxs="i",
        xlab = "Time (days)", ylab="Proportion (I)")
   lines( Stb[,"y"], col="red",lty=1, lwd=4, ... )
   lines( St1[,"y"], col="cornflowerblue",lty=1, lwd=4, ... )
   lines( St2[,"y"], col="orange",lty=1, lwd=4, ... )
   lines( St3[,"y"], col="brown",lty=1, lwd=4, ... )
   #lines(times, St[,"w"], col="blue",lty=3, lwd=1.5, ... )
   legend("topright", legend = c("Base, $\\mathcal R_0 =2.66$", "Sanitation, $\\rho = 0.05$","Vaccination, $\\nu = 0.005$", "Antibiotics, $\\eta = 0.016$"), lty=c(1,1,1,1),
          col = treatcol)
 }
 plot.basedump <- function(ic=c(x=x0,y=y0,z=0, w=W0), tmax=timemax,
                          times=seq(0,tmax,by=1),
                          func=func, parms=pars, ... ) {
   parms1 <- c(mu = 1/21170, # avg lifespan 58 yrs == 21170
               beta_i = 0.1, #transmission r. S between I
               gamma = 1/8, #rate of recovery (days)
               sigma = 1/22, #rate of water removal
               beta_w = 0.6, #transmission r. I between W
               xi = 1/22, #shedding rate I into W (force it equal to sigma)
               alpha = 0) #death rate by cholera
   Stb <- ode(ic, times, func=SIRW.vector.field, parms)
   St1 <- ode(ic, times, func=dump.vector.field, parms)
   St2 <- ode(ic, times, func=SIRW.vector.field, parms1)
   plot(x=0, y=0, type = "n",xlim=c(0, xmax), ylim=c(0, overlaymax), xaxs="i", yaxs="i", las=1,
        xlab = "Time (days)", ylab="Proportion (I)")
   lines( Stb[,"y"], col="red",lty=1, lwd=4, ... )
   lines( St1[,"y"], col="purple",lty=1, lwd=4, ... )
   lines( St2[,"y"], col="black",lty=1, lwd=4, ... )
   legend("topright", legend = c("Base $\\alpha = 0$", "Base $\\alpha = 0.1$","Dumping $\\alpha = 0.1$"), lty=c(1,1), col = c("black","red","purple"))
 }
###################
##PHASE PORTRAITS##
###################
  phase.portrait <- function(S0=seq(0.5,0.9,by=0.1),
                             I0=rep(0.001, times=5),
                             Z0=0, W0=0){
    plot(x=0, y=0, type = "n", xlim = c(0,1), ylim = c(0,0.15), xaxs="i", yaxs="i",las=1, #empty subplot
         xlab = "S", ylab= " ", bty="l", main="Base Model")
    mtext("I", side=2, line=3.2, las=1)
    lines(x=c(1,0),y=c(0,1),col="grey",lty="dotted", lwd=3) #S + I = 1
    for (j in 1:length(S0)) {
      plot.SIRW(ic=c(x=S0[j],y=I0[j],z=Z0, w = W0), tmax=timemax,
                func=SIRW.vector.field,
                parms = pars 
      )
    }
    points(x=S0, y=I0,col="black", pch=21, bg="white", xpd=TRUE, cex=1.2, lwd=2)
  }
  phase.portrait.treat1 <- function(S0=seq(0.5,0.9,by=0.1),
                                    I0=rep(0.001, times=5),
                                    Z0=0, W0=0){
    plot(x=0, y=0, type = "n", xlim = c(0,1), ylim = c(0,.25), xaxs="i", yaxs="i",las=1, #empty subplot
         xlab = "S", ylab=" ", las=1, bty="l", main="Treatment 1: Water Sanitation")
    mtext("I", side=2, line=3.2, las=1)
    lines(x=c(1,0),y=c(0,1),col="grey",lty="dotted", lwd=3) #S + I = 1
    for (j in 1:length(S0)) {
      plot.SIRW(ic=c(x=S0[j],y=I0[j],z=Z0, w = W0), tmax=timemax,
                func=treat1.vector.field,
                parms = pars 
      )
    }
    points(x=S0, y=I0,col="black", pch=21, bg="white", xpd=TRUE, cex=1.2, lwd=2)
  }
  phase.portrait.treat2 <- function(S0=seq(0.5,0.9,by=0.1),
                                    I0=rep(0.001, times=5),
                                    Z0=0, W0=0){
    plot(x=0, y=0, type = "n", xlim = c(0,1), ylim = c(0,.25), xaxs="i", yaxs="i",las=1, #empty subplot
         xlab = "S", ylab=" ",  bty="l", main="Treatment 2: Vaccination")
    mtext("I",  ylab=" ", side=2, line=3.2, las=1)
    lines(x=c(1,0),y=c(0,1),col="grey",lty="dotted", lwd=3) #S + I = 1
    for (j in 1:length(S0)) {
      plot.SIRW(ic=c(x=S0[j],y=I0[j],z=Z0, w = W0), tmax=timemax,
                func=treat2.vector.field,
                parms = pars 
      )
    }
    points(x=S0, y=I0,col="black", pch=21, bg="white", xpd=TRUE, cex=1.2, lwd=2)
  }
  phase.portrait.treat3 <- function(S0=seq(0.5,0.9,by=0.1),
                                    I0=rep(0.001, times=5),
                                    Z0=0, W0=0){
    plot(x=0, y=0, type = "n", xlim = c(0,1), ylim = c(0,.25), xaxs="i", yaxs="i",las=1, #empty subplot
         xlab = "S", ylab=" ", bty="l", main="Treatment 3: Antibiotics")
    mtext("I", side=2, line=3.2, las=1)
    lines(x=c(1,0),y=c(0,1),col="grey",lty="dotted", lwd=3) #S + I = 1
    for (j in 1:length(S0)) {
      plot.SIRW(ic=c(x=S0[j],y=I0[j],z=Z0, w = W0), tmax=timemax,
                func=treat3.vector.field,
                parms = pars 
      )
    }
    points(x=S0, y=I0,col="black", pch=21, bg="white", xpd=TRUE, cex=1.2, lwd=2)
  }
  phase.portrait.sev <-function(which="low", S0=seq(0,0.9,by=0.1)[-1],
                                I0=seq(0.9,0, by=-0.1)[-10],
                                Z0=0, W0=0){
    if (which == "low"){
      plot(x=0, y=0, type = "n", xlim = c(0,1), ylim = c(0,1), #empty subplot
           xlab = "S", ylab= "I_l", las=1, bty="l", main="Severity (Low)")
      lines(x=c(1,0),y=c(0,1),col="grey",lty="dotted", lwd=3) #S + I = 1
      for (j in 1:length(S0)) {
        plot.SIIRW(int=which, ic=c(x=S0[j],l=I0[j], h=0, z=Z0, w = W0), tmax=timemax,
                   func=severity.vector.field,
                   parms = pars)      
      }
    }
    else if (which == "high") {
      plot(x=0, y=0, type = "n", xlim = c(0,1), ylim = c(0,1), #empty subplot
           xlab = "S", ylab= "I_h", las=1, bty="l", main="Severity (High)")
      lines(x=c(1,0),y=c(0,1),col="grey",lty="dotted", lwd=3) #S + I = 1
      for (j in 1:length(S0)) {
        phase.SIIRW(int=which, ic=c(x=S0[j],l=0, h=I0[j],z=Z0, w = W0), tmax=timemax,
                    func=severity.vector.field,
                    parms = pars)      
      }
    }
  }
########
##test##
########
# plot.base()
# plot.dump()
# plot.basedump()
# plot.treat1()
# plot.treat2()
# plot.treat3()
plot.overlay()
# 
# phase.portrait()
# phase.portrait.treat1()
# phase.portrait.treat2()
# phase.portrait.treat3()

############
#Final Size#
############
# fun1 <- function(R_0, Z1){1-exp(-R_0*Z1)}
# Z1 <- seq(1, 8, by=1/20)
# R_0 <- c(1, 1.5, 2)
# R.0lab <-  c(1, 1.5, 2, "exact")
# res1 <- mapply(fun1, list(Z1), R_0)
# tcols <- c("cornflowerblue", "orange", "brown", "black")
# 
# fun <- function(Z){(-1/Z)*log(1-Z)}
# Z <- seq(0, 1, by=1/5000)
# res <- mapply(fun, list(Z))
# matplot(res, Z, type="l", lty=1, lwd=4, xlab="$\\mathcal R_0$", ylab="Final Size", xlim = c(0,8), las=1)
# matlines(Z1, res1, col=tcols, type="l", lty=3, lwd=2, xlab="$\\mathcal R_0$", ylab="Final Size")
# grid(col = "lightgray", lty = "dotted",lwd = par("lwd"), equilogs = TRUE)
# legend("bottomright", legend=R.0lab, title="value of $\\mathcal R_0$", lwd=2, col=tcols)
#####
#fix#
#####

##Initial Conditions Endemic Equilibrium Phase Portraits##
# params that give rise to EE
# timemax <- 500
# Z0 <-  0 
# W0 <-  0
# int <- "low"
# y0 <- 0.01
# x0 <- 1-y0
# 
#################
##Parameters EE##
#################
# pars=c(mu = 0.017, # avg lifespan 58 yrs
#        beta_i = 0.1, #transmission r. S between I
#        gamma = 1/14, #rate of recovery (days)
#        sigma = 1/25, #rate of water removal
#        beta_w = 0.1, #transmission r. I between W
#        xi = 1/25, #shedding rate I into W (force it eaual to sigma?)
#        #kappa = 0.15, #bacteria concentration that gives 50% of infection
#        alpha = 0, #death rate by cholera
#        
#        #parameters for intensity model
#        beta_l = 0.016, #transmission r. S between I_l
#        beta_h = 0.016, #transmission r. S between I_h
#        beta_w = 0.5, #transmission r. I between W
#        delta = 1/3, #rate of increasing symptom severity
#        xi_l = 1/14, #shedding rate I_l into W; maybe equal to sigma
#        xi_h = 1/7,
#        alpha_l = 0.01, #death rate by cholera with low intensity symptoms
#        alpha_h = 0.05, #death rate by cholera with high int. symptoms
#        
#        #parameters for water treatment
#        rho = 1/70,
#        #parameters for vaccination
#        nu = 1/700,
#        #parameters for antibiotics
#        eta = 1/70
# )
#######
##R-0##
#######
##Reproductive Ratio for base model (fix for each treatment)
R_0 <- function(which,parms=pars
) {
  with(as.list(parms), {
    if (which == "base"){
      return((beta_i+beta_w)/(mu + gamma + alpha))
    }
    else if (which == "treat1"){
      return((beta_i + beta_w)/(mu + gamma + alpha))
    }
    else if (which == "treat2"){
      return((beta_i + beta_w)/(mu + gamma + alpha))
    }
    else if (which == "treat3"){
      return((beta_i + beta_w)/(mu + gamma + alpha))
    }
 })
}
R_0(which="base")