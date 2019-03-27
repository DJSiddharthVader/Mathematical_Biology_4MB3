#!usr/bin/env Rscript
suppressMessages(library('dplyr'))

ptoi <- function(x,y,ymax){
    return(ymax*(y-1)+x)
}

singleToMatrix <- function(value,rows,cols){
    return(matrix(rep(value,rows*cols), nrow=rows,ncol=cols))
}

initPatches <- function(vp,hp,i0,w0){
    size <- vp*hp
    patches <- data.frame(idx=rep(NA,size),
                          xpos=rep(NA,size),
                          ypos=rep(NA,size),
                          s=rep(NA,size),
                          i=rep(NA,size),
                          r=rep(NA,size),
                          w=rep(NA,size),
                          n=rep(NA,size))
    for (v in 1:vp){
        for (h in 1:hp){
            idx <- ptoi(h,v,vp)
            i0vh <- i0[v,h]
            w0vh <- w0[v,h]
            patches[idx,] <- list(idx,h,v,1-i0vh,i0vh,0,w0vh,1)
        }
    }
    return(patches)
}

checkBoundary <- function(x,y,rx,ry,patches){
    if (x-rx < 1){
        xmin <- 1
    } else{
        xmin <- x-rx
    }
    if ((y-ry) < 1){
        ymin <- 1
    } else{
        ymin <- y-ry
    }
    if ((x+rx) > max(patches$xpos)){
        xmax <- max(patches$xpos)
    } else{
        xmax <- x+ry
    }
    if ((y+ry) > max(patches$ypos)){
        ymax <- max(patches$ypos)
    } else{
        ymax <- y+ry
    }
    return(list(xrange=xmin:xmax,yrange=ymin:ymax))
}

getNeighbours  <- function(x,y,rx,ry,patches,self=FALSE){
    ranges <- checkBoundary(x,y,rx,ry,patches)
    xrange <- ranges$xrange
    filterx <- patches[is.element(patches$xpos,xrange),]
    yrange <- ranges$yrange
    filtery <- filterx[is.element(filterx$ypos,yrange),]
    if (self){
        return(filtery)
    } else {
        filterself <- filtery[(filtery$xpos != x) | (filtery$ypos != y),]
        return(filterself)
    }
}

maxdist <- function(patches){
    xmax <- max(patches$xpos)
    ymax <- max(patches$ypos)
    return(dist(rbind(c(1,1),c(xmax,ymax)),method="manhattan")[1][1])
}

mdist <- function(x,y,x1,y1,norm=FALSE){
    d <- dist(rbind(c(x,y),c(x1,y1)),method="manhattan")[1][1]
    if (typeof(norm) != 'logical'){
        maxdist <- maxdist(norm)
        return(1-(d/maxdist))
    } else {
        return(d)
    }
}

getBetaI <- function(x,y,rx,ry,patches,beta_i){
    neighbours <- getNeighbours(x,y,rx,ry,patches)
    nbeta_i<- 0
    for(i in 1:nrow(neighbours)) {
        row <- neighbours[i,]
        nbeta_i <- nbeta_i + (1-mdist(x,y,row$xpos,row$ypos,norm=patches))*row$i*beta_i
    }
    return(nbeta_i)
}

getWaterLeeching <- function(x,y,rx,ry,patches,beta_w,scaling=0.5){
    neighbours <- getNeighbours(x,y,rx,ry,patches)
    leeching <- 0
    for (n in 1:dim(neighbours)[1]){
        neighbour <- neighbours[n,]
        leeching <- leeching + neighbour$w*(1-mdist(x,y,neighbour$x,neighbour$y,norm=patches))
    }
    return(leeching*scaling)
}

waterModel <- function(vars,params,nbeta_i,nbeta_w){
    with(as.list(c(vars, params)),{
        dx <- mu - mu*x - nbeta_i*x*y - beta_w*x*w # dS/dt
        dy <- nbeta_i*y*x + nbeta_w*x*w - gamma*y - mu*y# - alpha*y # dI/dt
        dz <- gamma*y - mu*z # dR/dt
        dw <- nbeta_w*y - sigma*w #dW/dt
        vec.fld <- as.list(c(s=dx, i=dy, r=dz, w=dw))
        return(vec.fld)
    })
}

patchEvolve <- function(patches,diffsys,params){
    with(as.list(c(params)),{
        size <- dim(patches)[1]
        newpatches <- data.frame(idx=rep(NA,size),
                                 xpos=rep(NA,size),
                                 ypos=rep(NA,size),
                                 s=rep(NA,size),
                                 i=rep(NA,size),
                                 r=rep(NA,size),
                                 w=rep(NA,size),
                                 n=rep(NA,size))
        for (patch in 1:size){
            old <- patches[patch,]
            nbeta_i <- getBetaI(old$xpos,old$ypos,rx,ry,patches,beta_i)
            nbeta_w <- getWaterLeeching(old$xpos,old$ypos,rx,ry,patches,beta_w)
            svars <- c(x=old$s,y=old$i,z=old$r,w=old$w)
            delta <- waterModel(svars,params,nbeta_i,nbeta_w)
            news <- old$s+delta$s
            newi <- old$i+delta$i
            newr <- old$r+delta$r
            neww <- old$w+delta$w
            newpatches[patch,] <- list(old$idx,
                                       old$xpos,
                                       old$ypos,
                                       news,
                                       newi,
                                       newr,
                                       neww,
                                       sum(c(news,newi,newr)))
        }
        return(newpatches)
    })
}

runSimulation <- function(tsteps,diffsys,ic,params){
#   params are c(v: number of vertical patches (integer)
#                h: number of horizontal patches (integer)
#                mu: natural birth/death rate (in [0,1])
#                gamma: cholera recovery rate (in [0,1])
#                betav: rate an infectious person transmits cholera to water (in [0,1])
#                sigma: rate of cholera removal from water (sanitation,antibiotics,death,etc) (in [0,1])
#                kappa: rate at which water infects people (in [0,1])
#                betat: transmission rate of cholera between people in a patch i (in [0,1])
#                alpha: death rate from cholera of those infected
#                rx,ry: neighbourhood in which infections and water spread (x,y ints with x<h, y<v)
#               )
#    ic is c(i0: v*h matrix of initial infective proportions
#            w0: v*h matrix of initial water infective loads
#            )
    with(as.list(c(params,ic)),{
        initialpatches <- initPatches(v,h,i0mat,w0mat)
        timestates <- list(initialpatches)
        times <- 1:tsteps
        for (t in times){
            timestates[[t+1]] <- patchEvolve(timestates[[t]],diffsys,params)
        }
        return(timestates)
    })
}

main <- function(tsteps){
    #params
    set.seed(9)
    v = 5
    h = 5
    mu = 1/25550
    beta_i = 1 #transmission r. S between I
    gamma = 1/7 #rate of recovery (days)
    sigma = 1/14 #rate of water sanitation (or 1/lifetime bacteria)
    beta_w = 1 #transmission r. I between W
    alpha = 0.005 #death rate by cholera
    rx = 1
    ry = 1
    #initial conditions
    i0 = 0.01
    i0mat <- singleToMatrix(i0,h,v)
    w0 = 0.05
    w0mat <- singleToMatrix(w0,h,v)
    ic <- list(i0mat=i0mat,w0mat=w0mat)
    params <- c(v=v,h=h,mu=mu,gamma=gamma,beta_i=beta_i,sigma=sigma,kappa=kappa,beta_w=beta_w,rx=rx,ry=ry,alpha=alpha)
    states <- runSimulation(tsteps,diffsys=waterModel,ic=ic,params=params)
    print(states[[tsteps-1]])
}
warnings()
if (sys.nframe() == 0){
    source('spatialSIR.R')
    main(150)
}

