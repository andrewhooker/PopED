/* file tmdd_qss_one_target.c */
#include <R.h>
static double parms[3];
#define CL parms[0]
#define V parms[1]
#define KA parms[2]

/* initializer  */
void initmod(void (* odeparms)(int *, double *))
{
  int N=3;
  odeparms(&N, parms);
}

/* Derivatives and 1 output variable */
void derivs (int *neq, double *t, double *y, double *ydot,
	     double *yout, int *ip)
{
    
  if (ip[0] <1) error("nout should be at least 1");
    
  ydot[0] = -KA*y[0];
  ydot[1] = KA*y[0] - CL/V*y[1];
  yout[0] = y[0]+y[1];
}

/* END file tmdd_qss_one_target.c */
