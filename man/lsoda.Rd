\name{lsoda}
\alias{lsoda}
\title{Solve System of ODE (ordinary differential equation)s.}
\description{Solving initial value problems for stiff or
  non-stiff systems of first-order ordinary differential equations
  (ODEs), The \R function \code{lsoda} provides an interface to the
  Fortran ODE solver of the same name, written by Linda R. Petzold and Alan
  C. Hindmarsh.  The system of ODE's is written as an \R function (which
  may, of course, use \code{\link{.C}}, \code{.Fortran},
  \code{\link{.Call}}, etc., to call foreign code) or be defined in
  compiled code that has been dynamically loaded.  A vector of
  parameters is passed to the ODEs, so the solver may be used as part of
  a modeling package for ODEs, or for parameter estimation using any
  appropriate modeling tool for non-linear models in \R such as
  \code{\link{optim}}, \code{\link[nls]{nls}}, \code{\link{nlm}} or
  \code{\link[nlme]{nlme}}.
}
\usage{
lsoda(y, times, func, parms, rtol, atol, tcrit=NULL, jacfunc=NULL,
verbose=FALSE, dllname=NULL)
}
\arguments{
  \item{y}{the initial values for the ode system.  If \code{y} has a
    name attribute, the names will be used to label the output matrix.}
  \item{times}{times at which explicit estimates for \code{y} are
    desired.  The first value in \code{times} must be the initial time.}
  \item{func}{either a user-supplied function that computes the values of the
    derivatives in the ode system (the \emph{model defininition}) at time
    t, or a character string
    giving the name of a compiled function in a dynamically loaded
    shared library.

    If \code{func}  is a user-supplied function, it must be called as:
    \code{yprime = func(t, y, parms)}.  \code{t} is the current time point
    in the integration, \code{y} is the current estimate of the variables
    in the ode system, and \code{parms} is a vector of parameters (which
    may have a names attribute, desirable in a large system).

    The return value of func should be a list, whose first element is a
    vector containing the derivatives of \code{y} with respect to
    \code{time}, and whose second element is a vector (possibly with a
    \code{\link{names}} attribute) of global values that are required at
    each point in \code{times}.

    If \code{func} is a string, then \code{dllname} must give the name
    of the shared library (without extension) which must be loaded
    before \code{lsoda()} is called.  see Details for more information.
  }
  \item{parms}{any parameters used in \code{func} that should be
    modifiable without rewriting the function.}
  \item{rtol}{relative error tolerance, either a scalar or an array as
    long as \code{y}.  See details.}
  \item{atol}{absolute error tolerance, either a scalar or an array as
    long as \code{y}.  See details.}
  \item{tcrit}{the Fortran routine \code{lsoda} overshoots its targets
    (times points in the vector \code{times}), and interpolates values
    for the desired time points.  If there is a time beyond which
    integration should not proceed (perhaps because of a singularity),
    that should be provided in \code{tcrit}.  Note that it does not
    make sense (though it is not an error) to include times in
    \code{times} past \code{tcrit}, since the solver will stop and
    return at the last point in \code{times} that is earlier than
    \code{tcrit}.}
  \item{jacfunc}{if not \code{NULL}, an \R function that computes
    the jacobian of the system of differential equations
    dydot(i)/dy(j), or a string giving the name of a function or
    subroutine in \file{dllname} that computes the jacobian (see Details
    below for more about this option).  In some circumstances, supplying
    \code{jac} can speed up
    the computations, if the system is stiff.  The \R calling sequence for
    \code{jac} is identical to that of \code{func}.  \code{jac} should
    return a vector whose ((i-1)*length(y) + j)th value is
    dydot(i)/dy(j). That is, return the matrix dydot/dy, where the ith
    row is the derivative of \eqn{dy_i/dt} with respect to \eqn{y_j},
    by columns (the way \R and Fortran store matrices).}
  \item{verbose}{a logical value that, when TRUE, should trigger more
    verbose output from the ode solver.  Currently does not do anything.}
  \item{dllname}{a string giving the name of the shared library (without
    extension) that contains all the compiled function or subroutine
    definitions refered to in \code{func} and \code{jac}.}
}
\details{
  All the hard work is done by the Fortran subroutine \code{lsoda},
  whose documentation should be consulted for details (it is included as
  comments in the source file \file{src/lsoda.f}).  This is based on the
  Feb 24, 1997 version of lsoda, from Netlib.  The following description
  of error control is adapted from that documentation (input arguments
  \code{rtol} and \code{atol}, above):
  
  The input parameters \code{rtol}, and \code{atol} determine the error
  control performed by the solver.  The solver will control the vector
  \bold{e} of estimated local errors in \bold{y}, according to an
  inequality of the form max-norm of ( \bold{e}/\bold{ewt} )
  \eqn{\leq}{<=} 1, where \bold{ewt} is a vector of positive error
  weights.  The values of \code{rtol} and \code{atol} should all be
  non-negative.
  The form of \bold{ewt} is:
  
  \deqn{\mathbf{rtol} \times \mathrm{abs}(\mathbf{y}) + \mathbf{atol}}{\bold{rtol} * abs(\bold{y}) + \bold{atol}}
  
  where multiplication of two vectors is element-by-element.

  If the request for precision exceeds the capabilities of the machine,
  the Fortran subroutine lsoda will return an error code; under some
  circumstances, the \R function \code{lsoda} will attempt a reasonable
  reduction of precision in order to get an answer.  It will write a
  warning if it does so.

  Models may be defined in compiled C or Fortran code, as well as in R. 
  For C, the calling sequence for func must be as in the following
  example:
  
  \code{void myderivs(long int *neq, double *t, double *y, double *ydot)}\cr
  \code{\{}\cr
    \code{ydot[0] = -k1*y[0] + k2*y[1]*y[2];}\cr
    \code{ydot[2] = k3 * y[1]*y[1];}\cr
    \code{ydot[1] = -ydot[0]-ydot[2];}\cr
  \code{\}}\cr

  where \code{*neq} is the number of equations, \code{*t} is the value
  of the independent variable, \code{y} points to a double precision
  array of length \code{*neq} that contains the current value of the
  state variables, and \code{ydot} points to an array that will contain
  the calculated derivatives.
  
  In this example, parameters are kept in a global variable in the C
  code declared as

  \code{static double parms[3];}

  \code{#define} statements are used to make the code more readable, as
  in \code{#define k1 parms[0]}
  
  This is the only way to pass parameters to a compiled C function
  from the calling R code.  Functions that use this mechanism must be
  accompanied by a function with the same name as the shared library
  which has as its sole argument a pointer to a function (see
  declaration below) that fills a double array with double precision
  values, to copy the parameter
  values into the global variable.  In the example here, the library
  is named \file{mymod.so}, a function such as:
  
  \code{void mymod(void (* odeparms)(long int *, double *))}\cr
  \code{\{}\cr
    \code{long int N=3;}\cr
    \code{odeparms(&N, parms);}\cr
  \code{\}}\cr

will be required to initialize the parameter vector.  Here
  \code{mymod} just calls \code{odeparms} with a pointer to a long int
  that contains the dimension of the parameter vector, and a pointer to
  the array that will contain the parameter values.
  
  Models may also be defined in Fortran.  For example:
  
  \code{subroutine myderivs (neq, t, y, ydot)}\cr
  \code{double precision t, y, ydot, parms(3)}\cr
  \code{integer neq}\cr
  \code{dimension y(3), ydot(3)}\cr
  \code{common /myparms/parms}\cr
    
  \code{ydot(1) = -parms(1)*y(1) + parms(2)*y(2)*y(3)}\cr
  \code{ydot(3) = parms(3)*y(2)*y(2)}\cr
  \code{ydot(2) = -ydot(1) - ydot(3)}\cr
    
  \code{return}\cr
  \code{end}\cr
  
  
  In Fortran, parameters may be stored in a common block, in which
  case, the file that contains the model function definition must also
  contain a subroutine, again with the same name as the file which
  contains the model definition:
  
  \code{subroutine mymod(odeparms)}\cr
  \code{external odeparms}\cr
  \code{integer N}\cr
  \code{double precision parms(3)}\cr
  \code{common /myparms/parms}\cr
    
  \code{N = 3}\cr
  \code{call odeparms(N, parms)}\cr
  \code{return}\cr
  \code{end}\cr


  When models are defined in compiled code, there is no provision for
  returning quantities that are not directly solutions of the odes
  (unlike models defined in R code).

  If it is desired to supply a jacobian to the solver, then the jacobian
  must be defined in compiled code if the ode system is.  The C function
  call for such a function must be as in the following example:
  
  \code{void myjac(long int *neq, double *t, double *y, long int *ml,}\cr
  \code{           long int *mu, double *pd, long int *nrowpd)}\cr
  \code{\{}\cr
    \code{pd[0] = -k1;}\cr
    \code{pd[1] = k1;}\cr
    \code{pd[2] = 0.0;}\cr
    \code{pd[(*nrowpd)] = k2*y[2];}\cr
    \code{pd[(*nrowpd) + 1] = -k2*y[2] - 2*k3*y[1];}\cr
    \code{pd[(*nrowpd) + 2] = 2*k3*y[1];}\cr
    \code{pd[(*nrowpd)*2] = k2*y[1];}\cr
    \code{pd[2*(*nrowpd) + 1] = -k2 * y[1];}\cr
    \code{pd[2*(*nrowpd) + 2] = 0.0;}\cr
  \code{\}}\cr
      
The corresponding subroutine in Fortran is:

\code{subroutine myjac (neq, t, y, ml, mu, pd, nrowpd)}\cr
\code{integer neq, ml, mu, nrowpd}\cr
\code{double precision y(*), pd(nrowpd,*), t, parms(3)}\cr
\code{common /myparms/parms}\cr

\code{pd(1,1) = -parms(1)}\cr
\code{pd(2,1) =  parms(1)}\cr
\code{pd(3,1) = 0.0}\cr
\code{pd(1,2) = parms(2)*y(3)}\cr
\code{pd(2,2) = -parms(2)*y(3) - 2*parms(3)*y(2)}\cr
\code{pd(3,2) = 2*parms(3)*y(2)}\cr
\code{pd(1,3) = parms(2)*y(2)}\cr
\code{pd(2,3) = -parms(2)*y(2)}\cr
\code{pd(3,3) = 0.0}\cr

\code{return}\cr
\code{end}\cr

Examples in both C and Fortran are in the \file{dynload} subdirectory of
the \code{odesolve} package directory.
  

}
\value{
  A matrix with up to as many rows as elements in \code{times} and as
  many columns as elements in \code{y} plus the number of "global"
  values returned in the second element of the return from \code{func},
  plus and additional column for the time value.  There will be a row
  for each element in \code{times} unless the Fortran routine `lsoda'
  returns with an unrecoverable error.  If \code{y} has a names
  attribute, it will be used to label the columns of the output value.
  The output will have the attribute \code{istate} which returns the
  conditions under which the last call to lsoda returned.  See the
  source code for an explanation of those values: normal is
  \code{istate = 2}.
}
\references{
  Hindmarsh, Alan C. (1983) ODEPACK, A Systematized Collection of ODE
  Solvers; in p.55--64 of  Stepleman, R.W. et al.[ed.] (1983)
  \emph{Scientific Computing}, North-Holland, Amsterdam.
  
  Petzold, Linda R. (1983)  Automatic Selection of Methods for Solving
  Stiff and Nonstiff Systems of Ordinary Differential Equations.
  \emph{Siam J. Sci. Stat. Comput.} \bold{4}, 136--148.
  
  Netlib: \url{http://www.netlib.org}
}
\note{
  The \file{demo} directory contains some examples of using
  \code{\link[nlme]{gnls}} to estimate parameters in a
  dynamic model.
  }

\author{R. Woodrow Setzer \email{setzer.woodrow@epa.gov}}

\examples{
### lsexamp -- example from lsoda source code

## names makes this easier to read, but may slow down execution.
parms <- c(k1=0.04, k2=1e4, k3=3e7)
my.atol <- c(1e-6,  1e-10,  1e-6)
times <- c(0,4 * 10^(-1:10))
lsexamp <- function(t, y, p)
  {
    yd1 <- -p["k1"] * y[1] + p["k2"] * y[2]*y[3]
    yd3 <- p["k3"] * y[2]^2
    list(c(yd1,-yd1-yd3,yd3),c(massbalance=sum(y)))
  }
exampjac <- function(t, y, p)
  {
    c(-p["k1"],	 p["k1"],  0,
  
        p["k2"]*y[3],
      - p["k2"]*y[3] - 2*p["k3"]*y[2],
                       2*p["k3"]*y[2],
  
      p["k2"]*y[2],  -p["k2"]*y[2],  0
      )
  }
  
require(odesolve)
## measure speed (here and below)
system.time( 
out <- lsoda(c(1,0,0),times,lsexamp, parms, rtol=1e-4, atol= my.atol)
)  
out

## This is what the authors of lsoda got for the example:

## the output of this program (on a cdc-7600 in single precision)
## is as follows..
##
##   at t =  4.0000e-01   y =  9.851712e-01  3.386380e-05  1.479493e-02
##   at t =  4.0000e+00   y =  9.055333e-01  2.240655e-05  9.444430e-02
##   at t =  4.0000e+01   y =  7.158403e-01  9.186334e-06  2.841505e-01
##   at t =  4.0000e+02   y =  4.505250e-01  3.222964e-06  5.494717e-01
##   at t =  4.0000e+03   y =  1.831975e-01  8.941774e-07  8.168016e-01
##   at t =  4.0000e+04   y =  3.898730e-02  1.621940e-07  9.610125e-01
##   at t =  4.0000e+05   y =  4.936363e-03  1.984221e-08  9.950636e-01
##   at t =  4.0000e+06   y =  5.161831e-04  2.065786e-09  9.994838e-01
##   at t =  4.0000e+07   y =  5.179817e-05  2.072032e-10  9.999482e-01
##   at t =  4.0000e+08   y =  5.283401e-06  2.113371e-11  9.999947e-01
##   at t =  4.0000e+09   y =  4.659031e-07  1.863613e-12  9.999995e-01
##   at t =  4.0000e+10   y =  1.404280e-08  5.617126e-14  1.000000e+00

## Using the analytic jacobian speeds up execution a little :

system.time( 
outJ <- lsoda(c(1,0,0),times,lsexamp, parms, rtol=1e-4, atol= my.atol,
              jac = exampjac)
)
  
all.equal(out, outJ) # TRUE

}
\keyword{math}
