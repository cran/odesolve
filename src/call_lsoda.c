#include <time.h>
#include <string.h>

#include <R.h>
#include <Rdefines.h>

#include "odesolve.h"

/* Globals :*/ 
SEXP deriv_func;
SEXP jac_func;
SEXP envir;
SEXP gparms;

void lsoda_derivs (long int *neq, double *t, double *y, double *ydot)
{
  int i;
  SEXP R_fcall, ans, Time, Y;

  PROTECT(Time = NEW_NUMERIC(1));incr_N_Protect();
  REAL(Time)[0] = *t;
  PROTECT(Y = NEW_NUMERIC(*neq));incr_N_Protect();
  for (i = 0; i < *neq; i++)
    {
      REAL(Y)[i] = y[i];
    }
  PROTECT(R_fcall = lang4(deriv_func,Time,Y,gparms));incr_N_Protect();
  PROTECT(ans = eval(R_fcall, envir));incr_N_Protect();
  for (i = 0; i < *neq; i++)
    {
#if R_VERSION >= R_Version(1, 2, 0)
	ydot[i] = REAL(VECTOR_ELT(ans,0))[i];
#else
	ydot[i] = REAL(VECTOR(ans)[0])[i];
#endif
    }
  my_unprotect(4);
}

void lsoda_jac (long int *neq, double *t, double *y, long int *ml,
		    long int *mu, double *pd, long int *nrowpd)
{
  int i, j;
  SEXP R_fcall, ans, Time, Y;

  PROTECT(Time = NEW_NUMERIC(1));incr_N_Protect();
  REAL(Time)[0] = *t;
  PROTECT(Y = NEW_NUMERIC(*neq));incr_N_Protect();
  for (i = 0; i < *neq; i++)
    {
      REAL(Y)[i] = y[i];
    }
  PROTECT(R_fcall = lang4(jac_func,Time,Y,gparms));incr_N_Protect();
  PROTECT(ans = eval(R_fcall, envir));incr_N_Protect();
  for (i = 0; i < *neq; i++)
    for (j = 0; j < *neq; j++)
    {
      pd[i * (*nrowpd) + j] = REAL(ans)[i * (*neq) + j];
    }
  my_unprotect(4);
}


SEXP call_lsoda(SEXP y, SEXP times, SEXP func, SEXP parms, SEXP rtol,
		SEXP atol, SEXP rho, SEXP tcrit, SEXP jacfunc)
{
  SEXP yout, yout2, ISTATE;
  int i, j, k, ny, nt, repcount, latol, lrtol, nprot;
  double *xt, *xytmp, *rwork, tin, tout;
  long int neq, itol, itask, istate, iopt, lrw, liw, *iwork, jt, lrn, lrs,
    mflag, lstamp, lfnm, lunit;
  void (*jac)(long int *, double *, double *, long int *,
	      long int *, double *, long int *);
  char *timestamp, *fnm="lsoda.log";
  time_t result;

  init_N_Protect();
  lunit = 10;
  lfnm = strlen(fnm);
  time(&result);
  timestamp = asctime(localtime(&result));
  lstamp = strlen(timestamp);
  F77_CALL(openlog)(&lunit, fnm, &lfnm, timestamp, &lstamp);

  ny = LENGTH(y);
  neq = ny;
  mflag = 0;
  F77_CALL(xsetf)(&mflag);
  xytmp = (double *) R_alloc(neq, sizeof(double));
  for (j = 0; j < ny; j++) xytmp[j] = REAL(y)[j];
  nt = LENGTH(times);
  nprot = 0;
  PROTECT(deriv_func = func); incr_N_Protect();
  PROTECT(envir = rho);incr_N_Protect();
  PROTECT(yout = allocMatrix(REALSXP,ny+1,nt));incr_N_Protect();
  PROTECT(gparms = parms); incr_N_Protect();
  xt = NUMERIC_POINTER(times);
  latol = LENGTH(atol);
  lrtol = LENGTH(rtol);
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
  iwork = (long int *) R_alloc(liw, sizeof(long int));
  if (!isNull(jacfunc))
    {
      jac_func = jacfunc;
      jac = lsoda_jac;
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
      do
	{
	  if (istate == -1) istate = 2;
	  if (istate == -2)
	    {
	      for (j = 0; j < lrtol; j++) REAL(rtol)[j] *= rwork[13];
	      for (j = 0; j < latol; j++) REAL(atol)[j] *= rwork[13];
	      warning("Excessive precision requested.  `rtol' and `atol' have been scaled upwards by the factor %g\n",rwork[13]);
	      istate = 3;
	      
	    }
	  F77_CALL(lsoda) (lsoda_derivs, &neq, xytmp, &tin, &tout,
			   &itol, NUMERIC_POINTER(rtol), NUMERIC_POINTER(atol), &itask, &istate, &iopt, rwork,
			   &lrw, iwork, &liw, jac, &jt);
	  repcount ++;
	} while (tin < tout && istate >= -2 && repcount < 20);
      if (istate == -3)
	{
	  F77_CALL(closelog)(&lunit);
	  unprotect_all();
	  error("Illegal input to lsoda.  `lsoda.log' might have more information.\n");
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
  F77_CALL(closelog)(&lunit);
      
  unprotect_all();
  if (istate > 0)
    return(yout);
  else
    return(yout2);
}

