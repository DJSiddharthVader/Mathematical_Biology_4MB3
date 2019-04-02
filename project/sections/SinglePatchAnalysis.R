##Single Patch SIRW analysis
library(deSolve)

##Initial Parameters
pars=c(mu = 1/15370, # avg lifespan 58 yrs 
       beta_i = 0.025, #transmission r. S between I
       gamma = 1/7, #rate of recovery (days)
       sigma = 1/25, #rate of water removal
       beta_w = 0.5, #transmission r. I between W
       xi = 5, #shedding rate I into W (force it eaual to sigma?)
       #kappa = 0.15, #bacteria concentration that gives 50% of infection
       alpha = 0, #death rate by cholera
       
       #parameters for intensity model
       beta_l = 0.016, #transmission r. S between I_l
       beta_h = 0.016, #transmission r. S between I_h
       beta_w = 0.5, #transmission r. I between W
       delta = 1/3, #rate of increasing symptom severity
       xi_l = 1/14, #shedding rate I_l into W; maybe equal to sigma
       xi_h = 1/7,
       alpha_l = 0.01, #death rate by cholera with low intensity symptoms
       alpha_h = 0.05, #death rate by cholera with high int. symptoms
       
       #parameters for water treatment
       rho = 0.9,
       #parameters for vaccination
       nu = 0.9,
       #parameters for antibiotics
       eta = 0.9
)
##Initial Conditions for Phase Portraits##
tmax <- 60
Z0 <-  0 
W0 <-  0.2
int <- "low"
x0 <- 0.95
y0 <- 0.05

############################
#######VECTOR FIELDS########
############################

##VECTOR FIELD SIRW BASE MODEL (proportions)##
  SIRW.vector.field <- function(t, vars, 
                                parms=pars) {
    with(as.list(c(parms, vars)), {
      dx <- mu*(x+y+z) - mu*x - beta_i*x*y - beta_w*x*w # dS/dt
      dy <- beta_i*x*y + beta_w*x*w - gamma*y - mu*y - alpha*y # dI/dt
      dz <- gamma*y - mu*z # dR/dt
      dw <- xi*y - sigma*w #dW/dt
        vec.fld <- c(dx=dx, dy=dy, dz=dz, dw = dw)
      return(list(vec.fld))
    })
  }
#SEVERITY MODEL VECTOR FIELD#
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
#VECTOR FIELD treatment 1 MODEL#
  treat1.vector.field <- function(t, vars, 
                                  parms=pars) {
    with(as.list(c(parms, vars)), {
      dx <- mu*(x+y+z) - mu*x - beta_i*x*y - beta_w*x*w # dS/dt
      dy <- beta_i*y*x + beta_w*x*w - gamma*y - mu*y - alpha*y # dI/dt
      dz <- gamma*y - mu*z # dR/dt
      dw <- xi*y - sigma*w - rho*w #dW/dt
      vec.fld <- c(dx=dx, dy=dy, dz=dz, dw = dw)
      return(list(vec.fld))
    })
  }
  #treatment 2 (vaccine) model#
  treat2.vector.field <- function(t, vars, 
                                  parms=pars) {
    with(as.list(c(parms, vars)), {
      dx <- mu*(x+y+z) - mu*x - beta_i*x*y - beta_w*x*w - nu*x# dS/dt
      dy <- beta_i*y*x + beta_w*x*w - gamma*y - mu*y - alpha*y # dI/dt
      dz <- gamma*y - mu*z + nu*x # dR/dt
      dw <- xi*y - sigma*w#dW/dt
      vec.fld <- c(dx=dx, dy=dy, dz=dz, dw = dw)
      return(list(vec.fld))
    })
  }
  #treatment 3 (antibiotic) model#
  treat3.vector.field <- function(t, vars, 
                                  parms=pars) {
    with(as.list(c(parms, vars)), {
      dx <- mu*(x+y+z) - mu*x - beta_i*x*y - beta_w*x*w# dS/dt
      dy <- beta_i*y*x + beta_w*x*w - gamma*y - mu*y - alpha*y - eta*y # dI/dt
      dz <- gamma*y - mu*z + eta*y # dR/dt
      dw <- xi*y - sigma*w#dW/dt
      vec.fld <- c(dx=dx, dy=dy, dz=dz, dw = dw)
      return(list(vec.fld))
    })
  }
  
  #########################
  ####PHASE PORTRAITS######
  #########################
  
  ##PLOTTING FN for PHASE PORTRAITS
  plot.SIRW <- function(ic=c(x=1,y=0,z=0, w=0), tmax=1,
                          times=seq(0,tmax,by=tmax/5000),
                          func, parms, ... ) {
    St <- ode(ic, times, func, parms) 
    lines(St[,"x"], St[,"y"], col="black", lwd=1.5, ... )
  }
  
  plot.SIIRW <- function(which="low", ic=c(x=1,l=0, h=0.1, z=0), tmax=1,
                         times=seq(0,tmax,by=tmax/5000),
                         func, parms, ... ) {
    St <- ode(ic, times, func, parms)
    if (int=="low"){
      lines(St[,"x"], St[,"l"], col="black", lwd=1.5, ... )
    }
    else if (int=="high"){
      lines(St[,"x"], St[,"h"], col="black", lwd=1.5, ... )
    }
  }

  
#PHASE PORTRAIT for SIRW##
  phase.portrait <- function(S0=seq(0,0.9,by=0.1)[-1],
                             I0=seq(0.9,0, by=-0.1)[-10],
                             Z0=0, W0=0.1){
    plot(x=0, y=0, type = "n", xlim = c(0,1), ylim = c(0,1), #empty subplot
         xlab = "S", ylab= "I", las=1, bty="l", main="Base Single Patch")
    lines(x=c(1,0),y=c(0,1),col="grey",lty="dotted") #S + I = 1
    for (j in 1:length(S0)) {
      plot.SIRW(ic=c(x=S0[j],y=I0[j],z=Z0, w = W0), tmax=tmax,
                func=SIRW.vector.field,
                parms = pars 
      )
    }  
  }
  phase.portrait.sev <-function(which="low", S0=seq(0,0.9,by=0.1)[-1],
                                   I0=seq(0.9,0, by=-0.1)[-10],
                                   Z0=0, W0=0.1){
    if (which == "low"){
      plot(x=0, y=0, type = "n", xlim = c(0,1), ylim = c(0,1), #empty subplot
           xlab = "S", ylab= "I_l", las=1, bty="l", main="Severity (Low)")
      lines(x=c(1,0),y=c(0,1),col="grey",lty="dotted") #S + I = 1
      for (j in 1:length(S0)) {
        plot.SIIRW(int=which, ic=c(x=S0[j],l=I0[j], h=0, z=Z0, w = W0), tmax=tmax,
                   func=severity.vector.field,
                   parms = pars)      
      }
    }
    else if (which == "high") {
      plot(x=0, y=0, type = "n", xlim = c(0,1), ylim = c(0,1), #empty subplot
           xlab = "S", ylab= "I_h", las=1, bty="l", main="Severity (High)")
      lines(x=c(1,0),y=c(0,1),col="grey",lty="dotted") #S + I = 1
      for (j in 1:length(S0)) {
        phase.SIIRW(int=which, ic=c(x=S0[j],l=0, h=I0[j],z=Z0, w = W0), tmax=tmax,
                    func=severity.vector.field,
                    parms = pars)      
    }
    }
  }
  
  phase.portrait.treat1 <- function(S0=seq(0,0.9,by=0.1)[-1],
                                    I0=seq(0.9,0, by=-0.1)[-10],
                                    Z0=0, W0=0.1){
    plot(x=0, y=0, type = "n", xlim = c(0,1), ylim = c(0,1), #empty subplot
         xlab = "S", ylab= "I", las=1, bty="l", main="Treatment 1")
    lines(x=c(1,0),y=c(0,1),col="grey",lty="dotted") #S + I = 1
    for (j in 1:length(S0)) {
      plot.SIRW(ic=c(x=S0[j],y=I0[j],z=Z0, w = W0), tmax=tmax,
                func=treat1.vector.field,
                parms = pars 
      )
    }  
  }
  
  phase.portrait.treat2 <- function(S0=seq(0,0.9,by=0.1)[-1],
                                    I0=seq(0.9,0, by=-0.1)[-10],
                                    Z0=0, W0=0.1){
    plot(x=0, y=0, type = "n", xlim = c(0,1), ylim = c(0,1), #empty subplot
         xlab = "S", ylab= "I", las=1, bty="l", main="Treatment 1")
    lines(x=c(1,0),y=c(0,1),col="grey",lty="dotted") #S + I = 1
    for (j in 1:length(S0)) {
      plot.SIRW(ic=c(x=S0[j],y=I0[j],z=Z0, w = W0), tmax=tmax,
                func=treat2.vector.field,
                parms = pars 
      )
    }  
  }
  
  phase.portrait.treat3 <- function(S0=seq(0,0.9,by=0.1)[-1],
                                    I0=seq(0.9,0, by=-0.1)[-10],
                                    Z0=0, W0=0.1){
    plot(x=0, y=0, type = "n", xlim = c(0,1), ylim = c(0,1), #empty subplot
         xlab = "S", ylab= "I", las=1, bty="l", main="Treatment 3")
    lines(x=c(1,0),y=c(0,1),col="grey",lty="dotted") #S + I = 1
    for (j in 1:length(S0)) {
      plot.SIRW(ic=c(x=S0[j],y=I0[j],z=Z0, w = W0), tmax=tmax,
                func=treat3.vector.field,
                parms = pars 
      )
    }  
  }

##plots SIRW over time##
#base model  
plot.base <- function(ic=c(x=x0,y=y0,z=0, w=W0), tmax=60,
                         times=seq(0,tmax,by=tmax/500),
                         func=SIRW.vector.field, parms=pars, ... ) { 
  plot(x=0, y=0, type = "n",xlim=c(0, tmax), ylim=c(0, 1),
       xlab = "Time (days)", ylab="Proportion", main = "Base Model with Vital Dynamics")
  St <- ode(ic, times, func, parms) 
  lines(times, St[,"x"], col="green",lty=3, lwd=1.5, ... )
  lines(times, St[,"y"], col="red",lty=1, lwd=1.5, ... )
  lines(times, St[,"z"], col="black",lty=3, lwd=1.5, ... )
  #lines(times, St[,"w"], col="orange",lty=3, lwd=1.5, ... )
  legend("topright", legend = c("S", "I","R"), lty=c(3,1,3,3), col=c("green", "red", "black"))
}

#severity bodel
plot.severity <- function(ic=c(x=0.95,l=0.05,h=0,z=0, w=W0), tmax=100,
                          times=seq(0,tmax,by=tmax/500),
                          func=severity.vector.field, parms=pars, ... ) { 
  plot(x=0, y=0, type = "n",xlim=c(0, tmax), ylim=c(0, 1),
       xlab = "Time (days)", ylab="Proportion", main = "Severity Model with Vital Dynamics")
  St <- ode(ic, times, func, parms) 
  lines(times, St[,"l"], col="red",lty=1, lwd=1.5, ... )
  lines(times, St[,"h"], col="red", lty=2, lwd=1.5, ... )
  lines(times, St[,"w"], col="blue", lty=2, lwd=1.5, ... )
  legend("topright", legend = c( "I_L","I_H", "W"), lty=1, col=c("red", "orange",  "blue"))
}

#treat1 model
plot.treat1 <- function(ic=c(x=x0,y=y0,z=0, w=W0), tmax=60,
                        times=seq(0,tmax,by=tmax/500),
                        func=treat1.vector.field, parms=pars, ... ) { 
  plot(x=0, y=0, type = "n",xlim=c(0, tmax), ylim=c(0, 1),
       xlab = "Time (days)", ylab="Proportion", main = "Treatment 1 (water sanitation) Model")
  St <- ode(ic, times, func, parms) 
  lines(times, St[,"x"], col="green",lty=3, lwd=1.5, ... )
  lines(times, St[,"y"], col="red",lty=1, lwd=1.5, ... )
  lines(times, St[,"z"], col="black",lty=3, lwd=1.5, ... )
  #lines(times, St[,"w"], col="orange",lty=3, lwd=1.5, ... )
  legend("topright", legend = c("S", "I","R"), lty=c(3,1,3,3), col=c("green", "red", "black"))
}

#treat2 model
plot.treat2 <- function(ic=c(x=x0,y=y0,z=0, w=W0), tmax=60,
                        times=seq(0,tmax,by=tmax/500),
                        func=treat2.vector.field, parms=pars, ... ) { 
  plot(x=0, y=0, type = "n",xlim=c(0, tmax), ylim=c(0, 1),
       xlab = "Time (days)", ylab="Proportion", main = "Treatment 2 (vaccination) Model")
  St <- ode(ic, times, func, parms) 
  lines(times, St[,"x"], col="green",lty=3, lwd=1.5, ... )
  lines(times, St[,"y"], col="red",lty=1, lwd=1.5, ... )
  lines(times, St[,"z"], col="black",lty=3, lwd=1.5, ... )
  #lines(times, St[,"w"], col="orange",lty=3, lwd=1.5, ... )
  legend("topright", legend = c("S", "I","R"), lty=c(3,1,3,3), col=c("green", "red", "black"))
}

#treat3 model
plot.treat3 <- function(ic=c(x=x0,y=y0,z=0, w=W0), tmax=60,
                        times=seq(0,tmax,by=tmax/500),
                        func=treat3.vector.field, parms=pars, ... ) { 
  plot(x=0, y=0, type = "n",xlim=c(0, tmax), ylim=c(0, 1),
       xlab = "Time (days)", ylab="Proportion", main = "Treatment 3 (antibiotics) Model")
  St <- ode(ic, times, func, parms) 
  lines(times, St[,"x"], col="green",lty=3, lwd=1.5, ... )
  lines(times, St[,"y"], col="red",lty=1, lwd=1.5, ... )
  lines(times, St[,"z"], col="black",lty=3, lwd=1.5, ... )
  #lines(times, St[,"w"], col="orange",lty=3, lwd=1.5, ... )
  legend("topright", legend = c("S", "I","R"), lty=c(3,1,3,3), col=c("green", "red", "black"))
}


# ##Reproductive Ratio for base model
# R_0 <- function(parms=c(mu = 1/25550,
#                                beta_i = 0.0016, 
#                                gamma = 1/7, 
#                                sigma = 1/14, 
#                                beta_w = 0.5, 
#                                xi = 1/14, 
#                                alpha = 0 
# )) {
#   with(as.list(parms), {
#     return((beta_i+beta_w)/(mu + gamma + alpha)) 
#  })
# }
  