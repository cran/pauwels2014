/* Encodes a simplified network dynamics and structure(gene 6,7,8)*/
/* Although they are set to 1 initially, some experiments require to divide their value by 10. */
   
   
#include <R.h>
#include <math.h>
static double parms[16];

#define mrna6_degradation_rate parms[0]
#define mrna7_degradation_rate parms[1]
#define mrna8_degradation_rate parms[2]
#define p_degradation_rate parms[3]

#define r6_Kd parms[4]
#define r6_h parms[5]
#define r11_Kd parms[6]
#define r11_h parms[7]
#define r12_Kd parms[8]
#define r12_h parms[9]

#define pro6_strength parms[10]
#define pro7_strength parms[11]
#define pro9_strength parms[12]

#define rbs6_strength parms[13]
#define rbs7_strength parms[14]
#define rbs8_strength parms[15]


/* initializer */
void initmod(void (* odeparms)(int *, double *))
{
    int N=16;
    odeparms(&N, parms);
}

/* Derivatives and 1 output variable */
void derivs (int *neq, double *t, double *y, double *ydot,
             double *yout, int *ip)
{   
    if (ip[0] <1) error("nout should be at least 1");
    
    double as7, as9;
    double rs7;
    double g6, g7, g9;
    
    /* Rules (skipping the pro* and rbs* variables which are constant and do not enter */
    /* the equations anyway */
    /* y[0] = g6=cte  y[1]=p6  y[2]=p7  y[3]=p8,  */
    /* y[4] = v6_mrna  y[5] = v7_mrna ... y[6]=v8_mrna  */

    as7 = pow( y[1]/r12_Kd , r12_h ) /(1+ pow( y[1]/r12_Kd , r12_h ) );
    as9 =  pow( y[1]/r11_Kd , r11_h ) /(1 +  pow( y[1]/r11_Kd , r11_h ) );

    rs7 = 1/(1+ pow( y[2]/r6_Kd , r6_h ) );
    
    g6 = y[0];         /* cte */
    g7 = as7 * rs7;
    g9 = as9;

    /* mrna update (bravo pour la numerotation des promoteurs !!) */

    ydot[4] = pro6_strength * g6 - mrna6_degradation_rate * y[4]; /* d(v6) */
    ydot[5] = pro7_strength * g9 - mrna7_degradation_rate * y[5]; /* d(v7) */
    ydot[6] = pro9_strength * g7 - mrna8_degradation_rate * y[6]; /* d(v8) */

    /* protein update (g9=y[0]=cte). ribosome and promoteur numbering  */
    /* according to graph                                              */

    ydot[0] = 0.;
    ydot[1] = rbs6_strength * y[4] - p_degradation_rate * y[1];  /* d(p6) */
    ydot[2] = rbs8_strength * y[5] - p_degradation_rate * y[2];  /* d(p7) */
    ydot[3] = rbs7_strength * y[6] - p_degradation_rate * y[3];  /* d(p8) */

    yout[0] = 1;
}

