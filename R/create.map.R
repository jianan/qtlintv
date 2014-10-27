## from Karl Broman's R/qtl package.
######################################################################
#
# create.map
#
# create a new map with inserted inter-marker locations
#
# Note: map is a vector or a matrix with 2 rows
# 
# stepwidth = "fixed" is the original R/qtl version;
# stepwidth="variable" is for Brian Yandell and the qtlbim package
# stepwidth="max" creates the minimal number of inserted pseudomarkers
#                 to have the maximum stepwidth = step
######################################################################
create.map <-
    function(map, step, off.end, stepwidth = c("fixed", "variable", "max"))
{
  stepwidth <- match.arg(stepwidth)
  if(step<0 || off.end<0) stop("step and off.end must be > 0.")

  if(!is.matrix(map)) { # sex-ave map
    if(stepwidth == "variable") {
      if(off.end > 0) {
        tmp <- names(map)
        ## Append and prepend by off.end value (exact here, no roundoff).
        map <- c(map[1] - off.end, map, map[length(map)] + off.end)
        names(map) <- c("loc000", tmp, "loc999")
      }
      if(step == 0)
          return(unclass(map))

      ## Determine differences and expansion vector.
      dif <- diff(map)
      expand <- pmax(1, floor(dif / step))

      ## Create pseudomarker map.
      a <- min(map) + cumsum(c(0, rep(dif / expand, expand)))

      ## Names are marker names or locNNN.
      namesa <- paste("loc", seq(length(a)), sep = "")
      namesa[cumsum(c(1, expand))] <- names(map)
      names(a) <- namesa

      return(unclass(a))
    }
    if(stepwidth == "max") {
      if(off.end > 0) {
        toadd <- c(map[1] - off.end, map[length(map)]+off.end)

        if(step==0) {
          names(toadd) <- paste("loc", 1:2, sep="")
          map <- sort(c(map, toadd))
          return(unclass(map))
        }

        nmap <- c(map[1] - off.end, map, map[length(map)]+off.end)
      }
      else {
        nmap <- map
        toadd <- NULL
      }

      if(step==0 || (length(map)==1 && off.end==0)) return(unclass(map))
      
      d <- diff(nmap)
      nadd <- ceiling(d/step)-1
      if(sum(nadd) > 0) {
        for(j in 1:(length(nmap)-1)) {
          if(nadd[j]>0)
              toadd <- c(toadd, seq(nmap[j], nmap[j+1], len=nadd[j]+2)[-c(1,nadd[j]+2)])
        }
      }
      if(length(toadd) > 0)  {
        names(toadd) <- paste("loc", 1:length(toadd), sep="")
        map <- sort(c(map, toadd))
      }
      return(unclass(map))
    }

    if(length(map) == 1) { # just one marker!
      if(off.end==0) {
        if(step == 0) step <- 1
        nam <- names(map)
        map <- c(map,map+step)
        names(map) <- c(nam,paste("loc",step,sep=""))
      }
      else {
        if(step==0) m <- c(-off.end,off.end)
        else m <- seq(-off.end,off.end,by=step)
        m <- m[m!=0]
        names(m) <- paste("loc",m,sep="")
        map <- sort(c(m+map,map))
      }
      return(map)
    }

    minloc <- min(map)
    map <- map-minloc

    if(step==0 && off.end==0) return(map+minloc)
    else if(step==0 && off.end > 0) {
      a <- c(floor(min(map)-off.end),ceiling(max(map)+off.end))
      names(a) <- paste("loc", a, sep="")
      return(sort(c(a,map))+minloc)
    }
    else if(step>0 && off.end == 0) {
      a <- seq(floor(min(map)),max(map),
               by = step)
      if(any(is.na(match(a, map)))) {
        a <- a[is.na(match(a,map))]
        names(a) <- paste("loc",a,sep="")
        return(sort(c(a,map))+minloc)
      }
      else return(map+minloc)
    }
    else {
      a <- seq(floor(min(map)-off.end),ceiling(max(map)+off.end+step),
               by = step)
      a <- a[is.na(match(a,map))]

      # no more than one point above max(map)+off.end
      z <- (seq(along=a))[a >= max(map)+off.end]
      if(length(z) > 1) a <- a[-z[-1]]

      names(a) <- paste("loc",a,sep="")
      return(sort(c(a,map))+minloc)
    }
  } # end sex-ave map
  else { # sex-specific map
    if(stepwidth == "variable") {
      if(off.end > 0) {
        tmp <- colnames(map)
        map <- cbind(map[, 1] - off.end, map, map[, ncol(map)] + off.end)
        dimnames(map) <- list(NULL, c("loc000", tmp, "loc999"))
      }
      if(step == 0)
          return(unclass(map))

      ## Determine differences and expansion vector.
      dif <- diff(map[1, ])
      expand <- pmax(1, floor(dif / step))

      ## Create pseudomarker map.
      a <- min(map[1, ]) + cumsum(c(0, rep(dif / expand, expand)))
      b <- min(map[2, ]) + cumsum(c(0, rep(diff(map[2, ]) / expand, expand)))

      namesa <- paste("loc", seq(length(a)), sep = "")
      namesa[cumsum(c(1, expand))] <- dimnames(map)[[2]]
      map <- rbind(a,b)
      dimnames(map) <- list(NULL, namesa)

      return(unclass(map))
    }
    if(stepwidth == "max") {
      if(step==0 && off.end==0) return(unclass(map))
      if(step==0 && off.end>0) {
        if(ncol(map)==1) { # only one marker; assume equal recomb in sexes
          L1 <- L2 <- 1
        }
        else {
          L1 <- diff(range(map[1,]))
          L2 <- diff(range(map[2,]))
        }

        nam <- colnames(map)
        nmap1 <- c(map[1,1]-off.end, map[1,], map[1,ncol(map)]+off.end)
        nmap2 <- c(map[2,1]-off.end*L2/L1, map[2,], map[2,ncol(map)]+off.end*L2/L1)
        map <- rbind(nmap1, nmap2)
        colnames(map) <- c("loc1", nam, "loc2")
        return(unclass(map))
      }

      if(ncol(map)==1) L1 <- L2 <- 1
      else {
        L1 <- diff(range(map[1,]))
        L2 <- diff(range(map[2,]))
      }

      nam <- colnames(map)

      if(off.end > 0) {
        toadd1 <- c(map[1,1] - off.end, map[1,ncol(map)]+off.end)
        toadd2 <- c(map[2,1] + off.end*L2/L1, map[2,ncol(map)]+off.end*L2/L1)

        neword <- order(c(map[1,], toadd1))
        nmap1 <- c(map[1,], toadd1)[neword]
        nmap2 <- c(map[2,], toadd2)[neword]
      }
      else {
        nmap1 <- map[1,]
        nmap2 <- map[2,]
        toadd1 <- toadd2 <- NULL
      }

      d <- diff(nmap1)
      nadd <- ceiling(d/step)-1
      if(sum(nadd) > 0) {
        for(j in 1:(length(nmap1)-1)) {
          if(nadd[j]>0) {
            toadd1 <- c(toadd1, seq(nmap1[j], nmap1[j+1], len=nadd[j]+2)[-c(1,nadd[j]+2)])
            toadd2 <- c(toadd2, seq(nmap2[j], nmap2[j+1], len=nadd[j]+2)[-c(1,nadd[j]+2)])
          }
        }
      }
      newnam <- paste("loc", 1:length(toadd1), sep="")
      
      toadd1 <- sort(toadd1)
      toadd2 <- sort(toadd2)
      neword <- order(c(map[1,], toadd1))
      nmap1 <- c(map[1,], toadd1)[neword]
      nmap2 <- c(map[2,], toadd2)[neword]
      map <- rbind(nmap1, nmap2)
      colnames(map) <- c(nam, newnam)[neword]

      return(unclass(map))
    }

    minloc <- c(min(map[1,]),min(map[2,]))
    map <- unclass(map-minloc)
    markernames <- colnames(map)

    if(step==0 && off.end==0) return(map+minloc)
    else if(step==0 && off.end > 0) {
      map <- map+minloc
      if(ncol(map)==1) { # only one marker; assume equal recomb in sexes
        L1 <- L2 <- 1
      }
      else {
        L1 <- diff(range(map[1,]))
        L2 <- diff(range(map[2,]))
      }

      nam <- colnames(map)
      nmap1 <- c(map[1,1]-off.end, map[1,], map[1,ncol(map)]+off.end)
      nmap2 <- c(map[2,1]-off.end*L2/L1, map[2,], map[2,ncol(map)]+off.end*L2/L1)
      map <- rbind(nmap1, nmap2)
      colnames(map) <- c("loc1", nam, "loc2")
      return(map)
    }
    else if(step>0 && off.end == 0) {

      if(ncol(map)==1) return(map+minloc)

      a <- seq(floor(min(map[1,])),max(map[1,]),
               by = step)
      a <- a[is.na(match(a,map[1,]))]

      if(length(a)==0) return(map+minloc)

      b <- sapply(a,function(x,y,z) {
        ZZ <- min((seq(along=y))[y > x])
        (x-y[ZZ-1])/(y[ZZ]-y[ZZ-1])*(z[ZZ]-z[ZZ-1])+z[ZZ-1] }, map[1,],map[2,])

      m1 <- c(a,map[1,])
      m2 <- c(b,map[2,])

      names(m1) <- names(m2) <- c(paste("loc",a,sep=""),markernames)
      return(rbind(sort(m1),sort(m2))+minloc)
    }
    else {
      a <- seq(floor(min(map[1,])-off.end),ceiling(max(map[1,])+off.end+step),
               by = step)
      a <- a[is.na(match(a,map[1,]))]
      # no more than one point above max(map)+off.end
      z <- (seq(along=a))[a >= max(map[1,])+off.end]
      if(length(z) > 1) a <- a[-z[-1]]

      b <- sapply(a,function(x,y,z,ml) {
        if(x < min(y)) {
          return(min(z) - (min(y)-x)/diff(range(y))*diff(range(z)) - ml)
        }
        else if(x > max(y)) {
          return(max(z) + (x - max(y))/diff(range(y))*diff(range(z)) - ml)
        }
        else {
          ZZ <- min((seq(along=y))[y > x])
          (x-y[ZZ-1])/(y[ZZ]-y[ZZ-1])*(z[ZZ]-z[ZZ-1])+z[ZZ-1]
        }
      }, map[1,],map[2,], minloc[2])
      m1 <- c(a,map[1,])
      m2 <- c(b,map[2,])
      names(m1) <- names(m2) <- c(paste("loc",a,sep=""),markernames)
      return(rbind(sort(m1),sort(m2))+minloc)
    }
  }
}

