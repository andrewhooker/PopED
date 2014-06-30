/* file tmdd_qss_one_target.c */
#include <R.h>
static double parms[14];
#define CL parms[0]
#define V1 parms[1]
#define Q parms[2]
#define V2 parms[3]
#define FAVAIL parms[4]
#define KA parms[5]
#define VMAX parms[6]
#define KMSS parms[7]
#define R0 parms[8]
#define KSSS parms[9]
#define KDEG parms[10]
#define KINT parms[11]
#define DOSE parms[12]
#define SC_FLAG parms[13]

/* initializer  */
void initmod(void (* odeparms)(int *, double *))
{
  int N=14;
  odeparms(&N, parms);
}

/* Derivatives and 1 output variable */
void derivs (int *neq, double *t, double *y, double *ydot,
	     double *yout, int *ip)
{
    
  if (ip[0] <1) error("nout should be at least 1");
    
  double RTOT, CTOT, CFREE;
  RTOT = y[3];
  CTOT= y[1]/V1;
  CFREE = 0.5*((CTOT-RTOT-KSSS)+sqrt((CTOT-RTOT-KSSS)*(CTOT-RTOT-KSSS)+4*KSSS*CTOT));
    
  ydot[0] = -KA*y[0];
  ydot[1] = FAVAIL*KA*y[0]+(Q/V2)*y[2]-(CL/V1+Q/V1)*CFREE*V1-RTOT*KINT*CFREE*V1/(KSSS+CFREE);
  ydot[2] = (Q/V1)*CFREE*V1 - (Q/V2)*y[2];
  ydot[3] = R0*KDEG - KDEG*RTOT - (KINT-KDEG)*(RTOT*CFREE/(KSSS+CFREE));

  yout[0] = y[0]+y[1]+y[2]+y[3];
}

/* The Jacobian matrix */
void jac(int *neq, double *t, double *y, int *ml, int *mu,
	 double *pd, int *nrowpd, double *yout, int *ip)
{
  
  double RTOT, CTOT, CFREE;
  RTOT = y[3];
  CTOT= y[1]/V1;
  CFREE = 0.5*((CTOT-RTOT-KSSS)+sqrt((CTOT-RTOT-KSSS)*(CTOT-RTOT-KSSS)+4*KSSS*CTOT));

  double CFREE_y1,CFREE_y3;
  CFREE_y1 = 0.5*(1/V1+0.5*1/sqrt((CTOT-RTOT-KSSS)*(CTOT-RTOT-KSSS)+4*KSSS*CTOT)*(2*(CTOT-RTOT-KSSS)*(1/V1)+4*KSSS/V1));
  CFREE_y3 = 0.5*(-1+0.5*1/sqrt((CTOT-RTOT-KSSS)*(CTOT-RTOT-KSSS)+4*KSSS*CTOT)*(2*(CTOT-RTOT-KSSS)*(-1)));
  
  pd[0]               = -KA;
  pd[1]               = FAVAIL*KA;
  pd[2]               = 0.0;
  pd[3]               = 0.0;
  pd[(*nrowpd)]       = 0.0;
  pd[(*nrowpd)+1]=-(CL/V1+Q/V1)*CFREE_y1*V1-RTOT*KINT*CFREE_y1/(KSSS+CFREE)+RTOT*KINT*CFREE*V1/((KSSS+CFREE)*(KSSS+CFREE))*CFREE_y1;
  pd[(*nrowpd) + 2]   = (Q/V1)*CFREE_y1*V1;
  pd[(*nrowpd) + 3]   = -(KINT-KDEG)*(RTOT*CFREE_y1/(KSSS+CFREE)-RTOT*CFREE/((KSSS+CFREE)*(KSSS+CFREE))*CFREE_y1);
  pd[(*nrowpd)*2]     = 0.0;
  pd[2*(*nrowpd) + 1] = Q/V2;
  pd[2*(*nrowpd) + 2] = - Q/V2;
  pd[2*(*nrowpd) + 3] = 0.0;
  pd[(*nrowpd)*3]     = 0.0;
  pd[3*(*nrowpd)+1]=(CL/V1+Q/V1)*CFREE_y3*V1-KINT*CFREE*V1/(KSSS+CFREE)-RTOT*KINT*CFREE_y3*V1/(KSSS+CFREE)+RTOT*KINT*CFREE*V1/((KSSS+CFREE)*(KSSS+CFREE))*CFREE_y3;
  pd[3*(*nrowpd) + 2] =  (Q/V1)*CFREE_y3*V1;
  pd[3*(*nrowpd) + 3] = -KDEG-(KINT-KDEG)*(CFREE/(KSSS+CFREE)+RTOT*CFREE_y3/(KSSS+CFREE)-RTOT*CFREE/((KSSS+CFREE)*(KSSS+CFREE))*CFREE_y3);
}
/* END file tmdd_qss_one_target.c */
