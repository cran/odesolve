/* R-stop -- implement some functions for keeping track of how many SEXPs 
 * 	are PROTECTed, and UNPROTECTing them in the case of a fortran stop.
 * Also, implement Rstop, which calls error, as a replacement
 * for fortran STOP (call Rstop(message) */
#include <R.h>
#include <Rdefines.h>

#include "odesolve.h"

long int N_Protected;

void init_N_Protect(void) { N_Protected = 0; }

void incr_N_Protect(void) { N_Protected++; }

void unprotect_all(void) { UNPROTECT((int) N_Protected); }

void my_unprotect(int n)
{
    UNPROTECT(n);
    N_Protected -= n;
}

void F77_SUB(rstop)(char *message) { error(message); }
