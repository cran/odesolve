/* compile with R CMD SHLIB mymod.c or (on Windows) Rcmd SHLIB mymod.c */
/* Example from lsoda documentation */
#include <R.h> /* gives F77_CALL through R_ext/RS.h */

static double parms[3];

/* A trick to keep up with the parameters */
#define k1 parms[0]
#define k2 parms[1]
#define k3 parms[2]

/* initializer */
void mymod(void (* odeparms)(int *, double *))
{
    int N=3;
    odeparms(&N, parms);
}

/* Derivatives */
void myderivs(int *neq, double *t, double *y, double *ydot)
{
    ydot[0] = -k1*y[0] + k2*y[1]*y[2];
    ydot[2] = k3 * y[1]*y[1];
    ydot[1] = -ydot[0]-ydot[2];
}

/* Jacobian */
void myjac(int *neq, double *t, double *y, int *ml,
                  int *mu, double *pd, int *nrowpd)
{
    pd[0] = -k1;
    pd[1] = k1;
    pd[2] = 0.0;
    pd[(*nrowpd)] = k2*y[2];
    pd[(*nrowpd) + 1] = -k2*y[2] - 2*k3*y[1];
    pd[(*nrowpd) + 2] = 2*k3*y[1];
    pd[(*nrowpd)*2] = k2*y[1];
    pd[2*(*nrowpd) + 1] = -k2 * y[1];
    pd[2*(*nrowpd) + 2] = 0.0;
}
/* End of example */
 
