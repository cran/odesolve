
/* ./R-stop.c utilities */
void init_N_Protect(void);
void incr_N_Protect(void);
void unprotect_all(void);
void my_unprotect(int);

void F77_NAME(rstop)(char *message);

/* ./call_lsoda.c  Exports : */
void lsoda_derivs(long int *neq, double *t, double *y, double *ydot);
void lsoda_jac (long int *neq, double *t, double *y, long int *ml,
		long int *mu, double *pd, long int *nrowpd);
SEXP call_lsoda(SEXP y, SEXP times, SEXP func, SEXP parms, SEXP rtol,
		SEXP atol, SEXP rho, SEXP tcrit, SEXP jacfunc);

/* Fortran "top-level" functions: */
void F77_NAME(openlog)(long int *,char *,long int *,char *,long int *);
void F77_NAME(lsoda)(void (*)(long int *, double *, double *, double *),
		     long int *, double *, double *, double *,
		     long int *, double *, double *, long int *, long int *,
		     long int *, double *,long int *,long int *, long int *,
		     void (*)(long int *, double *, double *, long int *,
			      long int *, double *, long int *),
		     long int *);
void F77_NAME(closelog)(long int *);
void F77_NAME(xsetf)(long int *);



