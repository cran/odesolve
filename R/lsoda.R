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
###
###          'func' may be a character vector instead of an R function.  If
###          so, then if jacfunc is not NULL, it must be a character vector
###          as well.  In these cases, the first value in 'func' is the name
###          of a function to be found in the dll named (without extension)
###          in func[2].  'jacfunc[1]' points to the name of the jacobian.
###          if the function pointed to by func[2] exists in the dll pointed
###          to by func[2], it is extracted and assumed to be the
###          the initializer for the problem.

lsoda <- function(y, times, func, parms, rtol=1e-6, atol=1e-6,
	tcrit = NULL, jacfunc=NULL, verbose=FALSE, dllname=NULL, hmin=0, hmax=Inf)
{
    if (!is.numeric(y)) stop("`y' must be numeric")
    n <- length(y)
    if (!is.numeric(times)) stop("`times' must be numeric")
    if (!is.function(func) && !is.character(func))
      stop("`func' must be a function or character vector")
    if (is.character(func) && (is.null(dllname) || !is.character(dllname)))
      stop("You need to specify the name of the dll or shared library where func can be found (without extension)")
    if (!is.numeric(parms)) stop("`parms' must be numeric")
    if (!is.numeric(rtol)) stop("`rtol' must be numeric")
    if (!is.numeric(atol)) stop("`atol' must be numeric")
    if (!is.null(tcrit) & !is.numeric(tcrit)) stop("`tcrit' must be numeric")
    if (!is.null(jacfunc) && !(is.function(jacfunc) || is.character(jacfunc)))
        stop("`jacfunc' must be a function or character vector")
    if (length(atol) > 1 && length(atol) != n)
        stop("`atol' must either be a scaler, or as long as `y'")
    if (length(rtol) > 1 && length(rtol) != n)
        stop("`rtol' must either be a scaler, or as long as `y'")
    if (!is.numeric(hmin)) stop("`hmin' must be numeric")
    if (hmin < 0) stop ("`hmin' must be a non-negative value")
    if (!is.numeric(hmax)) stop("`hmax' must be numeric")
    if (hmax < 0) stop ("`hmax' must be a non-negative value")
    if (hmax == Inf) hmax <- 0 # in lsoda.f 0 is internally handled as Inf    
    
    ## If func is a character vector, then
    ## copy its value to funcname 
    ## check to make sure it describes
    ## a function in a loaded dll
    ## change the value of funcname to either symbol.C or symbol.For, depending
    ## on which way we find the symbol
    if (is.character(func)) {
      
      funcname <- func
      if (is.loaded(symbol.C(funcname),PACKAGE=dllname)) {
        funcname <- symbol.C(func)
      } else if (is.loaded(symbol.For(funcname), PACKAGE=dllname)) {
        funcname <- symbol.For(funcname)
      } else stop(paste("Unable to find",funcname,"in",dllname))
      ## get the pointer and put it in func
      func <- getNativeSymbolInfo(funcname,PACKAGE=dllname)$address

      ## Now, is there an init function?

      ModelInit <- if (is.loaded(symbol.C(dllname),PACKAGE=dllname)) {
        getNativeSymbolInfo(symbol.C(dllname), PACKAGE=dllname)$address
      } else if (is.loaded(symbol.For(dllname),PACKAGE=dllname)) {
        getNativeSymbolInfo(symbol.For(dllname), PACKAGE=dllname)$address
      } else NULL

      ## Finally, is there a jacobian?

      if (!is.null(jacfunc)) {
        if (!is.character(jacfunc))
          stop("If 'func' is dynloaded, so must 'jacfunc' be")
        jacfuncname <- jacfunc
        if (is.loaded(symbol.C(jacfuncname),PACKAGE=dllname)) {
          jacfuncname <- symbol.C(jacfuncname)
        } else if (is.loaded(symbol.For(jacfuncname),PACKAGE=dllname)) {
          jacfuncname <- symbol.For(jacfuncname)
        } else stop(paste("Unable to find",jacfuncname,"in",dllname))
        jacfunc <- getNativeSymbolInfo(jacfuncname,PACKAGE=dllname)$address
      }

      ## If we go this route, there are no "global" results.
      Nglobal <- 0
      rho <- NULL
    } else {
      ModelInit <- NULL
      rho <- environment(func)
      ## Call func once to figure out whether and how many "global"
      ## results it wants to return and some other safety checks
      tmp <- eval(func(times[1],y,parms),rho)
      if (!is.list(tmp)) stop("Model function must return a list\n")
      if (length(tmp[[1]]) != length(y))
        stop(paste("The number of derivatives returned by func() (",
                   length(tmp[[1]]),
                   "must equal the length of the initial conditions vector (",
                   length(y),")",sep=""))
      Nglobal <- if (length(tmp) > 1) length(tmp[[2]]) else 0
    }
    storage.mode(y) <- storage.mode(times) <- "double"
    out <- .Call("call_lsoda",y,times,func,parms,
                 rtol, atol, rho, tcrit, jacfunc, ModelInit,
                 as.integer(verbose), hmin, hmax, PACKAGE="odesolve")
    istate <- attr(out,"istate")
    nm <- c("time",
            if (!is.null(attr(y,"names"))) names(y) else as.character(1:n))
    if (Nglobal > 0) {
        out2 <- matrix( nrow=Nglobal, ncol=ncol(out))
        for (i in 1:ncol(out2)) {
          y <- out[-1,i]
          names(y) <- nm[-1]
          out2[,i] <- func(out[1,i],y,parms)[[2]]
          }
        out <- rbind(out,out2)
        nm <- c(nm,
                if (!is.null(attr(tmp[[2]],"names"))) names(tmp[[2]])
                else as.character((n+1) : (n + Nglobal)))
    }
    attr(out,"istate") <- istate
    dimnames(out) <- list(nm,NULL)
    t(out)
}
