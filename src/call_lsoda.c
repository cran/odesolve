#include <time.h>
#include <string.h>

#include "odesolve.h"

void F77_NAME(lsoda)(void (*)(int *, double *, double *, double *),
		     int *, double *, double *, double *,
		     int *, double *, double *, int *, int *,
		     int *, double *,int *,int *, int *,
		     void (*)(int *, double *, double *, int *,
			      int *, double *, int *),
		     int *);


static void lsoda_derivs (int *neq, double *t, double *y, double *ydot)
{
  int i;
  SEXP R_fcall, ans, Time, Y;

  PROTECT(Time = NEW_NUMERIC(1));incr_N_Protect();
  REAL(Time)[0] = *t;
  PROTECT(Y = allocVector(REALSXP,(*neq)));incr_N_Protect();
  setAttrib(Y, R_NamesSymbol, odesolve_Y_names);
  for (i = 0; i < *neq; i++)
    {
      REAL(Y)[i] = y[i];
    }
  PROTECT(R_fcall = lang4(odesolve_deriv_func,Time,Y,odesolve_gparms));
  incr_N_Protect();
  PROTECT(ans = eval(R_fcall, odesolve_envir));
  incr_N_Protect();
  for (i = 0; i < *neq; i++)
    {
	ydot[i] = REAL(VECTOR_ELT(ans,0))[i];
    }
  my_unprotect(4);
}

static void lsoda_jac (int *neq, double *t, double *y, int *ml,
		    int *mu, double *pd, int *nrowpd)
{
  int i, j;
  SEXP R_fcall, ans, Time, Y;

  PROTECT(Time = NEW_NUMERIC(1));incr_N_Protect();
  REAL(Time)[0] = *t;
  PROTECT(Y = allocVector(REALSXP,(*neq)));incr_N_Protect();
  setAttrib(Y, R_NamesSymbol, odesolve_Y_names);
  for (i = 0; i < *neq; i++)
    {
      REAL(Y)[i] = y[i];
    }
  PROTECT(R_fcall = lang4(odesolve_jac_func,Time,Y,odesolve_gparms));
  incr_N_Protect();
  PROTECT(ans = eval(R_fcall, odesolve_envir));
  incr_N_Protect();
  for (i = 0; i < *neq; i++)
    for (j = 0; j < *neq; j++)
    {
      pd[i * (*nrowpd) + j] = REAL(ans)[i * (*neq) + j];
    }
  my_unprotect(4);
}

typedef void deriv_func(int *, double *, double *,double *);
typedef void jac_func(int *, double *, double *, int *,
		      int *, double *, int *);
typedef void init_func(void (*)(int *, double *));

SEXP call_lsoda(SEXP y, SEXP times, SEXP func, SEXP parms, SEXP rtol,
		SEXP atol, SEXP rho, SEXP tcrit, SEXP jacfunc, SEXP initfunc,
		SEXP verbose, SEXP hmin, SEXP hmax)
{
  SEXP yout, yout2, ISTATE;
  int i, j, k, ny, nt, repcount, latol, lrtol, nprot;
  double *xt, *xytmp, *rwork, tin, tout, *Atol, *Rtol;
  int neq, itol, itask, istate, iopt, lrw, liw, *iwork, jt, lrn, lrs,
    mflag, lstamp, lfnm, lunit;
  /* void (*derivs)(int *, double *, double *,double *);
  void (*jac)(int *, double *, double *, int *,
	      int *, double *, int *);
  */
  deriv_func *derivs;
  jac_func *jac;
  init_func *initializer;
  time_t result;

  init_N_Protect();

  ny = LENGTH(y);
  neq = ny;
  mflag = INTEGER(verbose)[0];
  /*  F77_CALL(xsetf)(&mflag); We don't use this anymore */
  xytmp = (double *) R_alloc(neq, sizeof(double));
  for (j = 0; j < ny; j++) xytmp[j] = REAL(y)[j];
  nt = LENGTH(times);
  nprot = 0;
  PROTECT(yout = allocMatrix(REALSXP,ny+1,nt));incr_N_Protect();
  PROTECT(odesolve_gparms = parms); incr_N_Protect();
  PROTECT(odesolve_Y_names = getAttrib(y, R_NamesSymbol)); incr_N_Protect();

  if (inherits(func, "NativeSymbol")) 
    {
      derivs = (deriv_func *) R_ExternalPtrAddr(func);
      /* If there is an initializer, use it here */
      if (!isNull(initfunc))
	{
	  initializer = (init_func *) R_ExternalPtrAddr(initfunc);
	  initializer(Initodeparms);
	}
	  
    } else {
      derivs = (deriv_func *) lsoda_derivs;
      PROTECT(odesolve_deriv_func = func); incr_N_Protect();
      PROTECT(odesolve_envir = rho);incr_N_Protect();
    }
  xt = NUMERIC_POINTER(times);
  latol = LENGTH(atol);
  Atol = (double *) R_alloc((int) latol, sizeof(double));
  lrtol = LENGTH(rtol);
  Rtol = (double *) R_alloc((int) lrtol, sizeof(double));
  if (latol == 1 && lrtol == 1 ) itol = 1;
  if (latol  > 1 && lrtol == 1 ) itol = 2;
  if (latol == 1 && lrtol  > 1 ) itol = 3;
  if (latol  > 1 && lrtol  > 1 ) itol = 4;
  if (tcrit != R_NilValue)
    {
      itask = 4;
    }
  else
    {
      itask = 1;
    }
  istate = 1;
  iopt = 0;
  lrn = 20 + 16 * neq;
  lrs = 22 + 9 * neq + neq * neq;
  if (lrn > lrs) lrw = lrn;
  else lrw = lrs;
  rwork = (double *) R_alloc(lrw, sizeof(double));
  if (itask == 4) rwork[0] = REAL(tcrit)[0];
  liw = 20 + neq;
  iwork = (int *) R_alloc(liw, sizeof(int));

  for (i=4; i<10; i++) {
    rwork[i]=0.0;
    iwork[i]=0;
  }
  rwork[5] = REAL(hmax)[0];
  rwork[6] = REAL(hmin)[0];
  if (rwork[5] >0 || rwork[6] >0) iopt = 1;

  if (!isNull(jacfunc))
    {
      if (inherits(jacfunc,"NativeSymbol"))
	{
	  jac = (jac_func *) R_ExternalPtrAddr(jacfunc);
	}
      else
	{
	  odesolve_jac_func = jacfunc;
	  jac = lsoda_jac;
	}
      jt = 1;
    }
  else
    jt = 2;
  REAL(yout)[0] = REAL(times)[0];
  for (j = 0; j < ny; j++)
    {
      REAL(yout)[j+1] = REAL(y)[j];
    }
  for (i = 0; i < nt-1; i++)
    {
      /* I still need to trap possible error returns of lsoda */
      /* based on return values in istate */
      tin = REAL(times)[i];
      tout = REAL(times)[i+1];
      repcount = 0;
      for (j = 0; j < lrtol; j++) Rtol[j] = REAL(rtol)[j];
      for (j = 0; j < latol; j++) Atol[j] = REAL(atol)[j];
      do
	{
	  if (istate == -1) istate = 2;
	  if (istate == -2)
	    {
	      for (j = 0; j < lrtol; j++) Rtol[j] *= 10.0;
	      for (j = 0; j < latol; j++) Atol[j] *= 10.0;
	      warning("Excessive precision requested.  `rtol' and `atol' have been scaled upwards by the factor %g\n",10.0);
	      istate = 3;
	      
	    }
	  F77_CALL(lsoda) (derivs, &neq, xytmp, &tin, &tout,
			   &itol, NUMERIC_POINTER(rtol), NUMERIC_POINTER(atol), &itask, &istate, &iopt, rwork,
			   &lrw, iwork, &liw, jac, &jt);
	  repcount ++;
	} while (tin < tout && istate >= -2 && repcount < 20);
      if (istate == -3)
	{
	  unprotect_all();
	  error("Illegal input to lsoda\n");
	}
      else
	{
	  REAL(yout)[(i+1)*(ny+1)] = tin;
	  for (j = 0; j < ny; j++)
	    REAL(yout)[(i+1)*(ny + 1) + j + 1] = xytmp[j];
	}
      if (istate < 0) {
	warning("Returning early from lsoda.  Results are accurate, as far as they go\n");
	/* need to redimension yout here, and add the attribute "istate" for */
	/* the most recent value of `istate' from lsoda */
	PROTECT(yout2 = allocMatrix(REALSXP,ny+1,(i+2)));incr_N_Protect();
	for (k = 0; k < i+2; k++)
	  for (j = 0; j < ny+1; j++)
	    REAL(yout2)[k*(ny+1) + j] = REAL(yout)[k*(ny+1) + j];
	break;
      }
    }
  PROTECT(ISTATE = allocVector(INTSXP, 1));incr_N_Protect();
  INTEGER(ISTATE)[0] = istate;
  if (istate > 0)
    {
      setAttrib(yout, install("istate"), ISTATE);
    }
  else
    {
      setAttrib(yout2, install("istate"), ISTATE);
    }
      
  unprotect_all();
  if (istate > 0)
    return(yout);
  else
    return(yout2);
}

