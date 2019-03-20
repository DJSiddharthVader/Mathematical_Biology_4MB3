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
                          w=rep(NA,size))
    for (v in 1:vp){
        for (h in 1:hp){
            idx <- ptoi(h,v,vp)
            i0vh <- i0[v,h]
            w0vh <- w0[v,h]
            patches[idx,] <- list(idx,h,v,1-i0vh,i0vh,0,w0vh)
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

getBetaI <- function(x,y,rx,ry,patches,betat){
    neighbours <- getNeighbours(x,y,rx,ry,patches)
    betai <- 0
    for(i in 1:nrow(neighbours)) {
        row <- neighbours[i,]
        betai <- betai + (1-mdist(x,y,row$xpos,row$ypos,norm=patches))*row$i*betat
    }
    return(betai)
}

getWaterLeeching <- function(x,y,rx,ry,patches){
    neighbours <- getNeighbours(x,y,rx,ry,patches)
    leeching <- 0
    for (n in 1:dim(neighbours)[1]){
        neighbour <- neighbours[n,]
        leeching <- leeching + neighbour$w*(1-mdist(x,y,neighbour$x,neighbour$y,norm=patches))
    }
    return(leeching)
}

patchEvolve <- function(patches,params){
    size <- dim(patches)[1]
    newpatches <- data.frame(idx=rep(NA,size),
                             xpos=rep(NA,size),
                             ypos=rep(NA,size),
                             s=rep(NA,size),
                             i=rep(NA,size),
                             r=rep(NA,size),
                             w=rep(NA,size))
    with(as.list(params), {
        for (patch in 1:size){
            old <- patches[patch,]
            betai <- getBetaI(old$xpos,old$ypos,rx,ry,patches,betat)
            leeching <- getWaterLeeching(old$xpos,old$ypos,rx,ry,patches)
            deltas <- mu - old$s*(betai+kappa*old$w+mu)
            deltai <- old$s*(betai+kappa*old$w) - old$i*(gamma+alpha+mu)
            deltar <- gamma*old$i - mu*old$r
            deltaw <- betav*old$i + leeching - sigma*old$w
#            deltas <- -old$s*(betai+kappa*old$w)
#            deltai <- old$s*(betai+kappa*old$w) - old$i*gamma
#            deltar <- gamma*old$i
#            deltaw <- betav*old$i + leeching - sigma*old$w
            newpatches[patch,] <- list(old$idx,
                                       old$xpos,
                                       old$ypos,
                                       old$s+deltas,
                                       old$i+deltai,
                                       old$r+deltar,
                                       old$w+deltaw)
        }
        return(newpatches)
    })
}

runSimulation <- function(tsteps,ic,params){
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
            timestates[[t+1]] <- patchEvolve(timestates[[t]],params)
        }
        return(timestates)
    })
}

main <- function(tsteps){
    #params
    set.seed(9)
    v = 5
    h = 5
    mu = 0.06
    gamma = 0.002
    betav = 0.001
    sigma = 0.003
    kappa = 0.002
    betat= 0.001
    alpha = 0.12
    rx = 1
    ry = 1
    #initial conditions
    i0 = 0.01
    i0mat <- singleToMatrix(i0,h,v)
    w0 = 0.05
    w0mat <- singleToMatrix(w0,h,v)
    ic <- list(i0mat=i0mat,w0mat=w0mat)
    params <- c(v=v,h=h,mu=mu,gamma=gamma,betav=betav,sigma=sigma,kappa=kappa,betat=betat,rx=rx,ry=ry,alpha=alpha)
    states <- runSimulation(tsteps,ic=ic,params=params)
    print(states[[tsteps]])
}

if (sys.nframe() == 0){
    source('spatialSIR.R')
    main(11)
}


#DEPRECIATED
#closestWell <- function(x,y,patches){
#    wells <- getWells(patches)
#    mindist = 100000
#    minwell <- c(wells[1,1],wells[1,2])
#    for (well in dim(wells)[1]){
#        wellx = wells[well,1]
#        welly = wells[well,2]
#        md <- mdist(x,y,wellx,welly)
#        if (md < mindist){
#            minwell <- c(wellx,welly)
#        }
#    }
#    return(minwell)
#}
#getWells <- function(patches){
#    welldf <- patches[patches$iswell == TRUE,]
#    #wdf = data.frame(x=rep(NA,dim(welldf)[1]),y=rep(NA,dim(welldf)[1]))
#    numw <- dim(welldf)[1]
#    wdf <- matrix(1:numw*2,ncol=2,nrow=numw)
#    for (i in 1:numw){
#        wdf[i,1] <- welldf$xpos[i]
#        wdf[i,2] <- welldf$ypos[i]
#    }
#    colnames(wdf) <- c('x','y')
#    return(wdf)
#}
#
#initPatches <- function(vp,hp,wells,i0,w0){
#    size <- vp*hp
#    patches <- data.frame(idx=rep(NA,size),
#                          s=rep(NA,size),
#                          i=rep(NA,size),
#                          r=rep(NA,size),
#                          w=rep(NA,size),
#                          xpos=rep(NA,size),
#                          ypos=rep(NA,size),
#                          iswell=rep(NA,size))
#    for (v in 1:vp){
#        for (h in 1:hp){
#            idx <- ptoi(h,v,vp)
#            patches[idx, ] <- list(idx,1-i0,i0,0,w0,v,h,FALSE)
#        }
#    }
#    ywells <- sample(1:vp,wells,replace=F)
#    xwells <- sample(1:hp,wells,replace=F)
#    for (w in 1:wells){
#        widx <- ptoi(xwells[w],ywells[w],vp)
#        patches$iswell[widx] <- TRUE
#    }
#    return(patches)
#}
