/* odesolve_utils.c */
/* Define some global variables and functions that operate on some of them */

#include <R.h>
#include <Rdefines.h>

/* some functions for keeping track of how many SEXPs 
 * 	are PROTECTed, and UNPROTECTing them in the case of a fortran stop.
 */
long int N_Protected;

void init_N_Protect(void) { N_Protected = 0; }

void incr_N_Protect(void) { N_Protected++; }

void unprotect_all(void) { UNPROTECT((int) N_Protected); }

void my_unprotect(int n)
{
    UNPROTECT(n);
    N_Protected -= n;
}


/* Globals :*/ 
SEXP odesolve_deriv_func;
SEXP odesolve_jac_func;
SEXP odesolve_envir;
SEXP odesolve_gparms;
SEXP odesolve_Y_names;

void Initodeparms(long int *N, double *parms)
{
  long int i, Nparms;

  Nparms = LENGTH(odesolve_gparms);
  if ((*N) != Nparms)
    {
      PROBLEM "Confusion over the length of parms"
      ERROR;
    } 
  else
    {
      for (i = 0; i < *N; i++) parms[i] = REAL(odesolve_gparms)[i];
    }
}
  

