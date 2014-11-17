/* file tmdd_qss_one_target.c */
#include <R.h>
static double parms[8];
#define CL parms[0]
#define V1 parms[1]
#define KA parms[2]
#define Q parms[3]
#define V2 parms[4]
#define Favail parms[5]
#define DOSE parms[6]
#define TAU parms[7]

/* initializer  */
void initmod(void (* odeparms)(int *, double *))
{
  int N=8;
  odeparms(&N, parms);
}

/* Derivatives and 1 output variable */
void derivs (int *neq, double *t, double *y, double *ydot,
	     double *yout, int *ip)
{
    
  if (ip[0] <1) error("nout should be at least 1");
    
  ydot[0] = -KA*y[0];
  ydot[1] = KA*y[0] + Q/V2*y[2]- (CL/V1+Q/V1)*y[1];
  ydot[2] = Q/V1*y[1]-Q/V2*y[2];
  yout[0] = y[0]+y[1]+y[2];
  
}

/* END file tmdd_qss_one_target.c */
