### lsoda -- solves ordinary differential equation systems defined in
###          func, and outputs values for the times in `times'
###          on input, y contains the initial values for times[1]
###          parms is a vector of parameters for func.  They should not
###          change during the integration. `rtol', and `atol'
###          are, respectively, the relative tolerance parameter, and the
###          absolute tolerance parameter.  `atol' may be scaler or vector.
###          `rtol' is a scaler or a vector.
###
###          The return value is a matrix whose rows correspond to the values
###          in `times', and columns to the elements of `y'.

lsoda <- function(y, times, func, parms, rtol=1e-6, atol=1e-6,
	tcrit = NULL, jacfunc=NULL)
{
    if (!is.numeric(y)) stop("`y' must be numeric")
    n <- length(y)
    if (!is.numeric(times)) stop("`times' must be numeric")
    if (!is.function(func)) stop("`func' must be a function")
    if (!is.numeric(parms)) stop("`parms' must be numeric")
    if (!is.numeric(rtol)) stop("`rtol' must be numeric")
    if (!is.numeric(atol)) stop("`atol' must be numeric")
    if (!is.null(tcrit) & !is.numeric(tcrit)) stop("`tcrit' must be numeric")
    if (!is.null(jacfunc) & !is.function(jacfunc))
        stop("`jacfunc' must be a function")
    if (length(atol) > 1 && length(atol) != n)
        stop("`atol' must either be a scaler, or as long as `y'")
    if (length(rtol) > 1 && length(rtol) != n)
        stop("`rtol' must either be a scaler, or as long as `y'")
    
    rho <- environment(func)
    ## Call func once to figure out whether and how many "global"
    ## results it wants to return
    tmp <- eval(func(times[1],y,parms),rho)
    Nglobal <- if (length(tmp) > 1) length(tmp[[2]]) else 0
    out <- .Call("call_lsoda",as.double(y),as.double(times),func,parms,
                 rtol, atol, rho, tcrit,
                 jacfunc)
    istate <- attr(out,"istate")
    nm <- c("time",
            if (!is.null(attr(y,"names"))) names(y) else as.character(1:n))
    if (Nglobal > 0) {
        out2 <- matrix( nrow=Nglobal, ncol=ncol(out))
        for (i in 1:ncol(out2))
            out2[,i] <- func(out[1,i],out[-1,i],parms)[[2]]
        out <- rbind(out,out2)
        nm <- c(nm,
                if (!is.null(attr(tmp[[2]],"names"))) names(tmp[[2]])
                else as.character((n+1) : (n + Nglobal)))
    }
    attr(out,"istate") <- istate
    dimnames(out) <- list(nm,NULL)
    t(out)
}
    
  


