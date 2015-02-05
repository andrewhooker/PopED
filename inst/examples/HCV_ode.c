/* file model.c */
#include <R.h>
static double parms[14];
#define p parms[0]
#define d parms[1]
#define e parms[2]
#define s parms[3]
#define KA parms[4]
#define KE parms[5]
#define VD parms[6]
#define EC50 parms[7]
#define n parms[8]
#define delta parms[9]
#define c parms[10]
#define DOSE parms[11]
#define TINF parms[12]
#define TAU parms[13]

/* initializer  */
void initmod(void (* odeparms)(int *, double *))
{
  int N=14;
  odeparms(&N, parms);
}

/* Derivatives */
void derivs (int *neq, double *t, double *y, double *ydot,
	     double *yout, int *ip)
{    

  double In;
  
  if ( fmod(*t,TAU) < TINF ){
    In =  DOSE/TINF;
  }
  else {
    In  =  0;
  }

  ydot[0] = -KA*y[0] + In;
  ydot[1] = KA*y[0] - KE*y[1];
  ydot[2] = s - y[2]*(e*y[4]+d);
  ydot[3] = e*y[2]*y[4]-delta*y[3];
  ydot[4] = p*(1-(pow(y[1]/VD,n)/(pow(y[1]/VD,n)+pow(EC50,n))))*y[3]-c*y[4];

}

/* END file model.c */
