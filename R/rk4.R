### Classical Runge-Kutta-fixed-step-integration
### R-Implementation by Th. Petzoldt,
### using some code from Woodrow Setzer
rk4 <- function(y, times, func, parms) {
    if (!is.numeric(y)) stop("`y' must be numeric")
    if (!is.numeric(times)) stop("`times' must be numeric")
    if (!is.function(func)) stop("`func' must be a function")
    if (!is.numeric(parms)) stop("`parms' must be numeric")

    n <- length(y)

    ## Call func once to figure out whether and how many "global"
    ## results it wants to return and some other safety checks
    rho <- environment(func)
    tmp <- eval(func(times[1],y,parms),rho)
    if (!is.list(tmp)) stop("Model function must return a list\n")
    if (length(tmp[[1]]) != length(y))
      stop(paste("The number of derivatives returned by func() (",
                 length(tmp[[1]]),
                 "must equal the length of the initial conditions vector (",
                 length(y),")",sep=""))
     Nglobal <- if (length(tmp) > 1) length(tmp[[2]]) else 0
    
    
    y0 <- y
    out <- c(times[1], y0)
    for (i in 1:(length(times)-1)) {
        t  <- times[i]
        dt <- times[i+1] - times[i]
        F1 <- dt * unlist(func(t,      y0,            parms))
        F2 <- dt * unlist(func(t+dt/2, y0 + 0.5 * F1, parms))
        F3 <- dt * unlist(func(t+dt/2, y0 + 0.5 * F2, parms))
        F4 <- dt * unlist(func(t+dt  , y0 + F3,       parms))
        dy <- (F1 + 2 * F2 + 2 * F3 + F4)/6
        y1 <- y0 + dy
        out<- rbind(out, c(times[i+1], y1))
        y0 <- y1
    }

    nm <- c("time", if (!is.null(attr(y, "names")))
            names(y) else as.character(1:n))
    if (Nglobal > 0) {
        out2 <- matrix( nrow=Nglobal, ncol=ncol(out))
        for (i in 1:ncol(out2))
            out2[,i] <- func(out[1,i],out[-1,i],parms)[[2]]
        out <- rbind(out,out2)
        nm <- c(nm,
                if (!is.null(attr(tmp[[2]],"names"))) names(tmp[[2]])
                else as.character((n+1) : (n + Nglobal)))
    }
    dimnames(out) <- list(NULL, nm)
    out
}
