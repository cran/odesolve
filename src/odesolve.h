#include <R.h>
#include <Rdefines.h>

/* odesolve_utils.c globals */
extern SEXP odesolve_deriv_func;
extern SEXP odesolve_jac_func;
extern SEXP odesolve_envir;
extern SEXP odesolve_gparms;
extern SEXP odesolve_Y_names;

/* odesolve_utils.c utilities */
void init_N_Protect(void);
void incr_N_Protect(void);
void unprotect_all(void);
void my_unprotect(int);
void Initodeparms(long int *, double *);

/* Fortran "top-level" functions: */
void F77_NAME(xsetf)(long int *);



