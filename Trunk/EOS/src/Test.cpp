// Example C program written by:
// Chris Muzny
// N.I.S.T.
// Chemical Science and Technology Laboratory
// Physical and Chemical Properties of Fluids Division
// (303) 497-5549 
// chris.muzny@nist.gov

//This program demonstrates explicitly linking the subroutines available in
// refprop.dll.  In order to link this code refprop1.h 
// must be available in the current directory.  When executing refprop.dll must be in the dll
// search path (current directory and $PATH).

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "NISTThermo_Extension.h"
#include "EOS.h"

// Some constants...

const int refpropcharlength = 255;
const int filepathlength = 255;
const int lengthofreference = 3;
const int errormessagelength = 255;
const int ncmax = 20; // Note: ncmax is the max number of components
const int numparams = 72;
const int maxcoefs = 50;

int main(int argc, char* argv[]) {
    // Now use the functions.

    // Refprop variables that need to be defined
    //
    // nc = Number of components in the mixture
    // x[NumberOfComponentsInMixtures] = Mole fraction of each component
    // ierr =  An integer flag defining an error
    // hf[] = a character array defining the fluids in a mixture
    // hrf[] = a character array denoting the reference state
    // herr[] = a character array for storing a string - Error message
    // hfmix[] a character array defining the path to the mixture file

    double x[ncmax], xliq[ncmax], xvap[ncmax];
    int ierr;
    char hf[refpropcharlength * ncmax], hrf[lengthofreference],
            herr[errormessagelength], hfmix[refpropcharlength];

    EOS_Init();

    //...initialize the program and set the pure fluid component name
    strcpy(hf,"nitrogen.fld");
    strcpy(hfmix,"hmx.bnc");
    strcpy(hrf,"DEF");
    strcpy(herr,"Ok");
    
    EOS_Set();
    EOS_Set_Reference_Properties(2395800.0, 119.5302, 1.0);
    EOS_Print_Reference_Properties();
    
    double Q[5], Pro[16];
    // RUP
    Q[0] = 0.5;
    Q[1] = 1.0;
    Q[2] = 0.0;
    Q[3] = 0.0;
    Q[4] = 1.0;
    EOS_Get_Properties(2, Q, Pro);
    for(int k = 0; k < 16; k++)
        printf("Property[%2d] = %10.6f\n", k, Pro[k]);
    printf("========\n");
    Q[4] = Pro[2];
    EOS_Get_Properties(4, Q, Pro);
    for(int k = 0; k < 16; k++)
        printf("Property[%2d] = %10.6f\n", k, Pro[k]);
    exit(0);
    
    double t, p, dl, dv;
    double d, q, e, h, s, cv, cp, w;
    double thermvp[29], thermvt[29];
    printf("================================================================\n");
    printf("Critical Point Test\n");
    printf("================================================================\n");
    p = 3395.80;
    d = 11.2839;
    PDEPDFLSH2(&p, &d, x, xliq, xvap, thermvp, &ierr, herr, errormessagelength);
    p = 0.0;
    t = thermvp[0];
    TDETDFLSH2(&t, &d, x, xliq, xvap, thermvt, &ierr, herr, errormessagelength);
    for (int k = 0; k < 29; k++)
        printf("THERMVP[%2d] = %10.4f \t THERMVT[%2d] = %10.4f \n", k, thermvp[k], k, thermvt[k]);
    p = 3395.80;
    d = 11.2839;
    PDEPDFLSH(&p, &d, x, &t, &dl, &dv, xliq, xvap, &q, &e, &h, &s, &cv, &cp, &w, &ierr, herr, errormessagelength);
    printf("Pressure        = %10.4f\n", p);
    printf("Density         = %10.4f\n", d);
    printf("Temperature     = %10.4f\n", t);
    printf("Density Liq     = %10.4f\n", dl);
    printf("Density Vap     = %10.4f\n", dv);
    printf("Composi Liq     = %10.4f\n", xliq[0]);
    printf("Composi Vap     = %10.4f\n", xvap[0]);
    printf("Quality         = %10.4f\n", q);
    printf("Internal Energy = %10.4f\n", e);
    printf("Enthalpy        = %10.4f\n", h);
    printf("Entropy         = %10.4f\n", s);
    printf("Cv              = %10.4f\n", cv);
    printf("Cp              = %10.4f\n", cp);
    printf("Speed Of Sound  = %10.4f\n", w);
    printf("Error Value     = %d\n", ierr);
    printf("Error String    = %s\n", herr);
    p = 0.0;
    TDETDFLSH(&t, &d, x, &p, &dl, &dv, xliq, xvap, &q, &e, &h, &s, &cv, &cp, &w, &ierr, herr, errormessagelength);
    printf("Pressure        = %10.4f\n", p);
    printf("Density         = %10.4f\n", d);
    printf("Temperature     = %10.4f\n", t);
    printf("Density Liq     = %10.4f\n", dl);
    printf("Density Vap     = %10.4f\n", dv);
    printf("Composi Liq     = %10.4f\n", xliq[0]);
    printf("Composi Vap     = %10.4f\n", xvap[0]);
    printf("Quality         = %10.4f\n", q);
    printf("Internal Energy = %10.4f\n", e);
    printf("Enthalpy        = %10.4f\n", h);
    printf("Entropy         = %10.4f\n", s);
    printf("Cv              = %10.4f\n", cv);
    printf("Cp              = %10.4f\n", cp);
    printf("Speed Of Sound  = %10.4f\n", w);
    printf("Error Value     = %d\n", ierr);
    printf("Error String    = %s\n", herr);
    printf("================================================================\n");
    printf("================================================================\n");
    printf("Multiphase Test\n");
    printf("================================================================\n");
    p = 2395.8110;
    d = 11.1839;
    PDEPDFLSH2(&p, &d, x, xliq, xvap, thermvp, &ierr, herr, errormessagelength);
    p = 0.0;
    t = thermvp[0];
    TDETDFLSH2(&t, &d, x, xliq, xvap, thermvt, &ierr, herr, errormessagelength);
    for (int k = 0; k < 29; k++)
        printf("THERMVP[%2d] = %10.4f \t THERMVT[%2d] = %10.4f \n", k, thermvp[k], k, thermvt[k]);
    p = 2395.8110;
    d = 11.1839;
    PDEPDFLSH(&p, &d, x, &t, &dl, &dv, xliq, xvap, &q, &e, &h, &s, &cv, &cp, &w, &ierr, herr, errormessagelength);
    printf("Pressure        = %10.4f\n", p);
    printf("Density         = %10.4f\n", d);
    printf("Temperature     = %10.4f\n", t);
    printf("Density Liq     = %10.4f\n", dl);
    printf("Density Vap     = %10.4f\n", dv);
    printf("Composi Liq     = %10.4f\n", xliq[0]);
    printf("Composi Vap     = %10.4f\n", xvap[0]);
    printf("Quality         = %10.4f\n", q);
    printf("Internal Energy = %10.4f\n", e);
    printf("Enthalpy        = %10.4f\n", h);
    printf("Entropy         = %10.4f\n", s);
    printf("Cv              = %10.4f\n", cv);
    printf("Cp              = %10.4f\n", cp);
    printf("Speed Of Sound  = %10.4f\n", w);
    printf("Error Value     = %d\n", ierr);
    printf("Error String    = %s\n", herr);
    p = 0.0;
    TDETDFLSH(&t, &d, x, &p, &dl, &dv, xliq, xvap, &q, &e, &h, &s, &cv, &cp, &w, &ierr, herr, errormessagelength);
    printf("Pressure        = %10.4f\n", p);
    printf("Density         = %10.4f\n", d);
    printf("Temperature     = %10.4f\n", t);
    printf("Density Liq     = %10.4f\n", dl);
    printf("Density Vap     = %10.4f\n", dv);
    printf("Composi Liq     = %10.4f\n", xliq[0]);
    printf("Composi Vap     = %10.4f\n", xvap[0]);
    printf("Quality         = %10.4f\n", q);
    printf("Internal Energy = %10.4f\n", e);
    printf("Enthalpy        = %10.4f\n", h);
    printf("Entropy         = %10.4f\n", s);
    printf("Cv              = %10.4f\n", cv);
    printf("Cp              = %10.4f\n", cp);
    printf("Speed Of Sound  = %10.4f\n", w);
    printf("Error Value     = %d\n", ierr);
    printf("Error String    = %s\n", herr);
    printf("================================================================\n");
    printf("================================================================\n");
    printf("Multiphase Vapor Saturated Line Test\n");
    printf("================================================================\n");
    p = 2395.8110;
    d = 4.1598;
    PDEPDFLSH2(&p, &d, x, xliq, xvap, thermvp, &ierr, herr, errormessagelength);
    p = 0.0;
    t = thermvp[0];
    TDETDFLSH2(&t, &d, x, xliq, xvap, thermvt, &ierr, herr, errormessagelength);
    for (int k = 0; k < 29; k++)
        printf("THERMVP[%2d] = %10.4f \t THERMVT[%2d] = %10.4f \n", k, thermvp[k], k, thermvt[k]);
    p = 2395.8110;
    d = 4.1598;
    PDEPDFLSH(&p, &d, x, &t, &dl, &dv, xliq, xvap, &q, &e, &h, &s, &cv, &cp, &w, &ierr, herr, errormessagelength);
    printf("Pressure        = %10.4f\n", p);
    printf("Density         = %10.4f\n", d);
    printf("Temperature     = %10.4f\n", t);
    printf("Density Liq     = %10.4f\n", dl);
    printf("Density Vap     = %10.4f\n", dv);
    printf("Composi Liq     = %10.4f\n", xliq[0]);
    printf("Composi Vap     = %10.4f\n", xvap[0]);
    printf("Quality         = %10.4f\n", q);
    printf("Internal Energy = %10.4f\n", e);
    printf("Enthalpy        = %10.4f\n", h);
    printf("Entropy         = %10.4f\n", s);
    printf("Cv              = %10.4f\n", cv);
    printf("Cp              = %10.4f\n", cp);
    printf("Speed Of Sound  = %10.4f\n", w);
    printf("Error Value     = %d\n", ierr);
    printf("Error String    = %s\n", herr);
    p = 0.0;
    TDETDFLSH(&t, &d, x, &p, &dl, &dv, xliq, xvap, &q, &e, &h, &s, &cv, &cp, &w, &ierr, herr, errormessagelength);
    printf("Pressure        = %10.4f\n", p);
    printf("Density         = %10.4f\n", d);
    printf("Temperature     = %10.4f\n", t);
    printf("Density Liq     = %10.4f\n", dl);
    printf("Density Vap     = %10.4f\n", dv);
    printf("Composi Liq     = %10.4f\n", xliq[0]);
    printf("Composi Vap     = %10.4f\n", xvap[0]);
    printf("Quality         = %10.4f\n", q);
    printf("Internal Energy = %10.4f\n", e);
    printf("Enthalpy        = %10.4f\n", h);
    printf("Entropy         = %10.4f\n", s);
    printf("Cv              = %10.4f\n", cv);
    printf("Cp              = %10.4f\n", cp);
    printf("Speed Of Sound  = %10.4f\n", w);
    printf("Error Value     = %d\n", ierr);
    printf("Error String    = %s\n", herr);
    printf("================================================================\n");
    printf("================================================================\n");
    printf("Multiphase Liquid Saturated Line Test\n");
    printf("================================================================\n");
    p = 2395.8110;
    d = 19.1024;
    PDEPDFLSH2(&p, &d, x, xliq, xvap, thermvp, &ierr, herr, errormessagelength);
    p = 0.0;
    t = thermvp[0];
    TDETDFLSH2(&t, &d, x, xliq, xvap, thermvt, &ierr, herr, errormessagelength);
    for (int k = 0; k < 29; k++)
        printf("THERMVP[%2d] = %10.4f \t THERMVT[%2d] = %10.4f \n", k, thermvp[k], k, thermvt[k]);
    p = 2395.8110;
    d = 19.1024;
    PDEPDFLSH(&p, &d, x, &t, &dl, &dv, xliq, xvap, &q, &e, &h, &s, &cv, &cp, &w, &ierr, herr, errormessagelength);
    printf("Pressure        = %10.4f\n", p);
    printf("Density         = %10.4f\n", d);
    printf("Temperature     = %10.4f\n", t);
    printf("Density Liq     = %10.4f\n", dl);
    printf("Density Vap     = %10.4f\n", dv);
    printf("Composi Liq     = %10.4f\n", xliq[0]);
    printf("Composi Vap     = %10.4f\n", xvap[0]);
    printf("Quality         = %10.4f\n", q);
    printf("Internal Energy = %10.4f\n", e);
    printf("Enthalpy        = %10.4f\n", h);
    printf("Entropy         = %10.4f\n", s);
    printf("Cv              = %10.4f\n", cv);
    printf("Cp              = %10.4f\n", cp);
    printf("Speed Of Sound  = %10.4f\n", w);
    printf("Error Value     = %d\n", ierr);
    printf("Error String    = %s\n", herr);
    p = 0.0;
    TDETDFLSH(&t, &d, x, &p, &dl, &dv, xliq, xvap, &q, &e, &h, &s, &cv, &cp, &w, &ierr, herr, errormessagelength);
    printf("Pressure        = %10.4f\n", p);
    printf("Density         = %10.4f\n", d);
    printf("Temperature     = %10.4f\n", t);
    printf("Density Liq     = %10.4f\n", dl);
    printf("Density Vap     = %10.4f\n", dv);
    printf("Composi Liq     = %10.4f\n", xliq[0]);
    printf("Composi Vap     = %10.4f\n", xvap[0]);
    printf("Quality         = %10.4f\n", q);
    printf("Internal Energy = %10.4f\n", e);
    printf("Enthalpy        = %10.4f\n", h);
    printf("Entropy         = %10.4f\n", s);
    printf("Cv              = %10.4f\n", cv);
    printf("Cp              = %10.4f\n", cp);
    printf("Speed Of Sound  = %10.4f\n", w);
    printf("Error Value     = %d\n", ierr);
    printf("Error String    = %s\n", herr);
    printf("================================================================\n");
    printf("================================================================\n");
    printf("Super Heated Vapor Test\n");
    printf("================================================================\n");
    p = 2395.8110;
    d = 4.0598;
    PDEPDFLSH2(&p, &d, x, xliq, xvap, thermvp, &ierr, herr, errormessagelength);
    p = 0.0;
    t = thermvp[0];
    TDETDFLSH2(&t, &d, x, xliq, xvap, thermvt, &ierr, herr, errormessagelength);
    for (int k = 0; k < 29; k++)
        printf("THERMVP[%2d] = %10.4f \t THERMVT[%2d] = %10.4f \n", k, thermvp[k], k, thermvt[k]);
    p = 2395.8110;
    d = 4.0598;
    PDEPDFLSH(&p, &d, x, &t, &dl, &dv, xliq, xvap, &q, &e, &h, &s, &cv, &cp, &w, &ierr, herr, errormessagelength);
    printf("Pressure        = %10.4f\n", p);
    printf("Density         = %10.4f\n", d);
    printf("Temperature     = %10.4f\n", t);
    printf("Density Liq     = %10.4f\n", dl);
    printf("Density Vap     = %10.4f\n", dv);
    printf("Composi Liq     = %10.4f\n", xliq[0]);
    printf("Composi Vap     = %10.4f\n", xvap[0]);
    printf("Quality         = %10.4f\n", q);
    printf("Internal Energy = %10.4f\n", e);
    printf("Enthalpy        = %10.4f\n", h);
    printf("Entropy         = %10.4f\n", s);
    printf("Cv              = %10.4f\n", cv);
    printf("Cp              = %10.4f\n", cp);
    printf("Speed Of Sound  = %10.4f\n", w);
    printf("Error Value     = %d\n", ierr);
    printf("Error String    = %s\n", herr);
    p = 0.0;
    TDETDFLSH(&t, &d, x, &p, &dl, &dv, xliq, xvap, &q, &e, &h, &s, &cv, &cp, &w, &ierr, herr, errormessagelength);
    printf("Pressure        = %10.4f\n", p);
    printf("Density         = %10.4f\n", d);
    printf("Temperature     = %10.4f\n", t);
    printf("Density Liq     = %10.4f\n", dl);
    printf("Density Vap     = %10.4f\n", dv);
    printf("Composi Liq     = %10.4f\n", xliq[0]);
    printf("Composi Vap     = %10.4f\n", xvap[0]);
    printf("Quality         = %10.4f\n", q);
    printf("Internal Energy = %10.4f\n", e);
    printf("Enthalpy        = %10.4f\n", h);
    printf("Entropy         = %10.4f\n", s);
    printf("Cv              = %10.4f\n", cv);
    printf("Cp              = %10.4f\n", cp);
    printf("Speed Of Sound  = %10.4f\n", w);
    printf("Error Value     = %d\n", ierr);
    printf("Error String    = %s\n", herr);
    printf("================================================================\n");
    printf("================================================================\n");
    printf("Subcooled Compressed Liquid Test\n");
    printf("================================================================\n");
    p = 2395.8110;
    d = 20.1024;
    PDEPDFLSH2(&p, &d, x, xliq, xvap, thermvp, &ierr, herr, errormessagelength);
    p = 0.0;
    t = thermvp[0];
    TDETDFLSH2(&t, &d, x, xliq, xvap, thermvt, &ierr, herr, errormessagelength);
    for (int k = 0; k < 29; k++)
        printf("THERMVP[%2d] = %10.4f \t THERMVT[%2d] = %10.4f \n", k, thermvp[k], k, thermvt[k]);
    p = 2395.8110;
    d = 20.1024;
    PDEPDFLSH(&p, &d, x, &t, &dl, &dv, xliq, xvap, &q, &e, &h, &s, &cv, &cp, &w, &ierr, herr, errormessagelength);
    printf("Pressure        = %10.4f\n", p);
    printf("Density         = %10.4f\n", d);
    printf("Temperature     = %10.4f\n", t);
    printf("Density Liq     = %10.4f\n", dl);
    printf("Density Vap     = %10.4f\n", dv);
    printf("Composi Liq     = %10.4f\n", xliq[0]);
    printf("Composi Vap     = %10.4f\n", xvap[0]);
    printf("Quality         = %10.4f\n", q);
    printf("Internal Energy = %10.4f\n", e);
    printf("Enthalpy        = %10.4f\n", h);
    printf("Entropy         = %10.4f\n", s);
    printf("Cv              = %10.4f\n", cv);
    printf("Cp              = %10.4f\n", cp);
    printf("Speed Of Sound  = %10.4f\n", w);
    printf("Error Value     = %d\n", ierr);
    printf("Error String    = %s\n", herr);
    p = 0.0;
    TDETDFLSH(&t, &d, x, &p, &dl, &dv, xliq, xvap, &q, &e, &h, &s, &cv, &cp, &w, &ierr, herr, errormessagelength);
    printf("Pressure        = %10.4f\n", p);
    printf("Density         = %10.4f\n", d);
    printf("Temperature     = %10.4f\n", t);
    printf("Density Liq     = %10.4f\n", dl);
    printf("Density Vap     = %10.4f\n", dv);
    printf("Composi Liq     = %10.4f\n", xliq[0]);
    printf("Composi Vap     = %10.4f\n", xvap[0]);
    printf("Quality         = %10.4f\n", q);
    printf("Internal Energy = %10.4f\n", e);
    printf("Enthalpy        = %10.4f\n", h);
    printf("Entropy         = %10.4f\n", s);
    printf("Cv              = %10.4f\n", cv);
    printf("Cp              = %10.4f\n", cp);
    printf("Speed Of Sound  = %10.4f\n", w);
    printf("Error Value     = %d\n", ierr);
    printf("Error String    = %s\n", herr);
    printf("================================================================\n");
    printf("================================================================\n");
    printf("Super Critical State Test\n");
    printf("================================================================\n");
    p = 4395.8110;
    d = 11.1839;
    PDEPDFLSH2(&p, &d, x, xliq, xvap, thermvp, &ierr, herr, errormessagelength);
    p = 0.0;
    t = thermvp[0];
    TDETDFLSH2(&t, &d, x, xliq, xvap, thermvt, &ierr, herr, errormessagelength);
    for (int k = 0; k < 29; k++)
        printf("THERMVP[%2d] = %10.4f \t THERMVT[%2d] = %10.4f \n", k, thermvp[k], k, thermvt[k]);
    p = 4395.8110;
    d = 11.1839;
    PDEPDFLSH(&p, &d, x, &t, &dl, &dv, xliq, xvap, &q, &e, &h, &s, &cv, &cp, &w, &ierr, herr, errormessagelength);
    printf("Pressure        = %10.4f\n", p);
    printf("Density         = %10.4f\n", d);
    printf("Temperature     = %10.4f\n", t);
    printf("Density Liq     = %10.4f\n", dl);
    printf("Density Vap     = %10.4f\n", dv);
    printf("Composi Liq     = %10.4f\n", xliq[0]);
    printf("Composi Vap     = %10.4f\n", xvap[0]);
    printf("Quality         = %10.4f\n", q);
    printf("Internal Energy = %10.4f\n", e);
    printf("Enthalpy        = %10.4f\n", h);
    printf("Entropy         = %10.4f\n", s);
    printf("Cv              = %10.4f\n", cv);
    printf("Cp              = %10.4f\n", cp);
    printf("Speed Of Sound  = %10.4f\n", w);
    printf("Error Value     = %d\n", ierr);
    printf("Error String    = %s\n", herr);
    p = 0.0;
    TDETDFLSH(&t, &d, x, &p, &dl, &dv, xliq, xvap, &q, &e, &h, &s, &cv, &cp, &w, &ierr, herr, errormessagelength);
    printf("Pressure        = %10.4f\n", p);
    printf("Density         = %10.4f\n", d);
    printf("Temperature     = %10.4f\n", t);
    printf("Density Liq     = %10.4f\n", dl);
    printf("Density Vap     = %10.4f\n", dv);
    printf("Composi Liq     = %10.4f\n", xliq[0]);
    printf("Composi Vap     = %10.4f\n", xvap[0]);
    printf("Quality         = %10.4f\n", q);
    printf("Internal Energy = %10.4f\n", e);
    printf("Enthalpy        = %10.4f\n", h);
    printf("Entropy         = %10.4f\n", s);
    printf("Cv              = %10.4f\n", cv);
    printf("Cp              = %10.4f\n", cp);
    printf("Speed Of Sound  = %10.4f\n", w);
    printf("Error Value     = %d\n", ierr);
    printf("Error String    = %s\n", herr);
    printf("================================================================\n");
    return 0;
}

int test2(int argc, char* argv[]) {
    // Now use the functions.

    // Refprop variables that need to be defined
    //
    // nc = Number of components in the mixture
    // x[NumberOfComponentsInMixtures] = Mole fraction of each component
    // ierr =  An integer flag defining an error
    // hf[] = a character array defining the fluids in a mixture
    // hrf[] = a character array denoting the reference state
    // herr[] = a character array for storing a string - Error message
    // hfmix[] a character array defining the path to the mixture file

    double x[ncmax], xliq[ncmax], xvap[ncmax], f[ncmax];

    int i, ierr;
    char hf[refpropcharlength * ncmax], hrf[lengthofreference],
            herr[errormessagelength], hfmix[refpropcharlength];

    EOS_Init();
    
    //Exlicitely set the fluid file PATH
    //char *FLD_PATH;  
    //FLD_PATH = "C:\\Program Files\\REFPROP\\fluids\\";
    //	  strcpy(hf,FLD_PATH);
    //	  strcpy(hfmix,FLD_PATH);

    //...initialize the program and set the pure fluid component name
    i=1;
    strcpy(hf,"nitrogen.fld");
    strcpy(hfmix,"hmx.bnc");
    strcpy(hrf,"DEF");
    strcpy(herr,"Ok");
    
    //...For a mixture, use the following setup instead of the lines above.
    // Use "|" as the file name delimiter for mixtures
//    i = 3;
//    strcpy(hf, "nitrogen.fld");
//    for (int j=12; j < 255; j++)
//        hf[j] = ' ';
//    strcat(hf, "argon.fld");
//    for (int j=255+9; j < 2*255; j++)
//        hf[j] = ' ';
//    strcat(hf, "oxygen.fld");
//    for (int j=2*255+10; j < 3*255; j++)
//        hf[j] = ' ';
//    for (int j=3*255; j < refpropcharlength*ncmax; j++)
//        hf[j] = ' ';
//    strcpy(hfmix, "hmx.bnc");
//    strcpy(hrf, "DEF");
//    strcpy(herr, "Ok");
//    x[0] = .7812; //Air composition
//    x[1] = .0092;
//    x[2] = .2096;
    
//    //...Call SETUP to initialize the program
//    SETUP(&i, hf, hfmix, hrf, &ierr, herr,
//            refpropcharlength*ncmax, refpropcharlength,
//            lengthofreference, errormessagelength);
//    if (ierr != 0) printf("%s\n", herr);
    
    EOS_Set();
    EOS_Set_Reference_Properties(2495.8110, 119.0734, 1.0);
    EOS_Print_Reference_Properties();
    
    
    double wm, ttp, tnbp, tc, pc, dc, zc, acf, dip, rgas;
    int info_index = 1;
    printf("INFO:\n");
    INFO(&info_index, &wm, &ttp, &tnbp, &tc, &pc, &dc, &zc, &acf, &dip, &rgas);
    printf("WM,ACF,DIP,TTP,TNBP   %10.4f,%10.4f,%10.4f,%10.4f,%10.4f\n", wm, acf, dip, ttp, tnbp);
    printf("TC,PC,DC,RGAS         %10.4f,%10.4f,%10.4f,%10.4f\n", tc, pc, dc, rgas);
    //...Calculate molecular weight of a mixture
    //     wm=WMOL(x)

    
    double t, p, dl, dv;
    //... Get the Critical Pressure, Temperature and Density
    CRITP(x, &t, &p, &dl, &ierr, herr, errormessagelength);
    printf("Critical: T, P, D     %10.4f,%10.4f,%10.4f\n", t, p, dl);
    
    //...Get saturation properties given t,x; for i=1: x is liquid phase
    //.....                                   for i=2: x is vapor phase
    t = 110.0;
    printf("Liquid, SATT @ T=120K:\n");
    SATT(&t, x, &i, &p, &dl, &dv, xliq, xvap, &ierr, herr, errormessagelength);
    printf("P,Dl,Dv,xl[0],xv[0]   %10.4f,%10.4f,%10.4f,%10.4f,%10.4f\n", p, dl, dv, xliq[0], xvap[0]);
    printf("============: %d\n", ierr);
    printf("Vapor, SATT @ T=120K:\n");
    i = 2;
    SATT(&t, x, &i, &p, &dl, &dv, xliq, xvap, &ierr, herr, errormessagelength);
    printf("P,Dl,Dv,xl[0],xv[0]   %10.4f,%10.4f,%10.4f,%10.4f,%10.4f\n", p, dl, dv, xliq[0], xvap[0]);
    printf("============: %d\n", ierr);
    
    //...Calculate saturation properties at a given p. i is same as SATT
    i = 2;
    SATP(&p, x, &i, &t, &dl, &dv, xliq, xvap, &ierr, herr, errormessagelength);
    printf("T,Dl,Dv,xl(1),xv(1)   %10.4f,%10.4f,%10.4f,%10.4f,%10.4f\n", t, dl, dv, xliq[0], xvap[0]);
    printf("============: %d\n", ierr);

    //...Other saturation routines are given in SAT_SUB.FOR

    int j = 1;
    double d, q, e, h, s, cv, cp, w, b, c,
            dpdrho, d2pdd2, dpdt, dhdt_d, dhdt_p, dhdp_t, dhdp_d,
            sigma, dhdd_t, dhdd_p, eta, tcx, pp, tt, hjt, h1, dd;
    int tmp_int = 0;
    
    printf("A===============================================================\n");
    p = 3395.80;
    d = 11.2839;
    double thermvp[29], thermvt[29];
    PDEPDFLSH2(&p, &d, x, xliq, xvap, thermvp, &ierr, herr, errormessagelength);
    p = 0.0;
    t = thermvp[0];
    TDETDFLSH2(&t, &d, x, xliq, xvap, thermvt, &ierr, herr, errormessagelength);
    for (int k = 0; k < 29; k++)
        printf("THERMVP[%2d] = %10.4f \t THERMVT[%2d] = %10.4f \n", k, thermvp[k], k, thermvt[k]);
    p = 3395.80;
    d = 11.2839;
    PDEPDFLSH(&p, &d, x, &t, &dl, &dv, xliq, xvap, &q, &e, &h, &s, &cv, &cp, &w, &ierr, herr, errormessagelength);
    printf("Pressure        = %10.4f\n", p);
    printf("Density         = %10.4f\n", d);
    printf("Temperature     = %10.4f\n", t);
    printf("Density Liq     = %10.4f\n", dl);
    printf("Density Vap     = %10.4f\n", dv);
    printf("Composi Liq     = %10.4f\n", xliq[0]);
    printf("Composi Vap     = %10.4f\n", xvap[0]);
    printf("Quality         = %10.4f\n", q);
    printf("Internal Energy = %10.4f\n", e);
    printf("Enthalpy        = %10.4f\n", h);
    printf("Entropy         = %10.4f\n", s);
    printf("Cv              = %10.4f\n", cv);
    printf("Cp              = %10.4f\n", cp);
    printf("Speed Of Sound  = %10.4f\n", w);
    printf("Error Value     = %d\n", ierr);
    printf("Error String    = %s\n", herr);
    p = 0.0;
    TDETDFLSH(&t, &d, x, &p, &dl, &dv, xliq, xvap, &q, &e, &h, &s, &cv, &cp, &w, &ierr, herr, errormessagelength);
    printf("Pressure        = %10.4f\n", p);
    printf("Density         = %10.4f\n", d);
    printf("Temperature     = %10.4f\n", t);
    printf("Density Liq     = %10.4f\n", dl);
    printf("Density Vap     = %10.4f\n", dv);
    printf("Composi Liq     = %10.4f\n", xliq[0]);
    printf("Composi Vap     = %10.4f\n", xvap[0]);
    printf("Quality         = %10.4f\n", q);
    printf("Internal Energy = %10.4f\n", e);
    printf("Enthalpy        = %10.4f\n", h);
    printf("Entropy         = %10.4f\n", s);
    printf("Cv              = %10.4f\n", cv);
    printf("Cp              = %10.4f\n", cp);
    printf("Speed Of Sound  = %10.4f\n", w);
    printf("Error Value     = %d\n", ierr);
    printf("Error String    = %s\n", herr);
    printf("A===============================================================\n");
    printf("B===============================================================\n");
    p = 2395.8110;
    d = 18.1024;
    PDEPDFLSH2(&p, &d, x, xliq, xvap, thermvp, &ierr, herr, errormessagelength);
    p = 0.0;
    t = thermvp[0];
    TDETDFLSH2(&t, &d, x, xliq, xvap, thermvt, &ierr, herr, errormessagelength);
    for (int k = 0; k < 29; k++)
        printf("THERMVP[%2d] = %10.4f \t THERMVT[%2d] = %10.4f \n", k, thermvp[k], k, thermvt[k]);
    p = 2395.8110;
    d = 18.1024;
    PDEPDFLSH(&p, &d, x, &t, &dl, &dv, xliq, xvap, &q, &e, &h, &s, &cv, &cp, &w, &ierr, herr, errormessagelength);
    printf("Pressure        = %10.4f\n", p);
    printf("Density         = %10.4f\n", d);
    printf("Temperature     = %10.4f\n", t);
    printf("Density Liq     = %10.4f\n", dl);
    printf("Density Vap     = %10.4f\n", dv);
    printf("Composi Liq     = %10.4f\n", xliq[0]);
    printf("Composi Vap     = %10.4f\n", xvap[0]);
    printf("Quality         = %10.4f\n", q);
    printf("Internal Energy = %10.4f\n", e);
    printf("Enthalpy        = %10.4f\n", h);
    printf("Entropy         = %10.4f\n", s);
    printf("Cv              = %10.4f\n", cv);
    printf("Cp              = %10.4f\n", cp);
    printf("Speed Of Sound  = %10.4f\n", w);
    printf("Error Value     = %d\n", ierr);
    printf("Error String    = %s\n", herr);
    p = 0.0;
    TDETDFLSH(&t, &d, x, &p, &dl, &dv, xliq, xvap, &q, &e, &h, &s, &cv, &cp, &w, &ierr, herr, errormessagelength);
    printf("Pressure        = %10.4f\n", p);
    printf("Density         = %10.4f\n", d);
    printf("Temperature     = %10.4f\n", t);
    printf("Density Liq     = %10.4f\n", dl);
    printf("Density Vap     = %10.4f\n", dv);
    printf("Composi Liq     = %10.4f\n", xliq[0]);
    printf("Composi Vap     = %10.4f\n", xvap[0]);
    printf("Quality         = %10.4f\n", q);
    printf("Internal Energy = %10.4f\n", e);
    printf("Enthalpy        = %10.4f\n", h);
    printf("Entropy         = %10.4f\n", s);
    printf("Cv              = %10.4f\n", cv);
    printf("Cp              = %10.4f\n", cp);
    printf("Speed Of Sound  = %10.4f\n", w);
    printf("Error Value     = %d\n", ierr);
    printf("Error String    = %s\n", herr);
    printf("B===============================================================\n");
    p = 0.0;
    d = dl;
    THERM(&t, &d, x, &p, &e, &h, &s, &cv, &cp, &w, &hjt);
    printf("Pressure        = %10.4f\n", p);
    printf("Density         = %10.4f\n", d);
    printf("Temperature     = %10.4f\n", t);
    printf("Internal Energy = %10.4f\n", e);
    printf("Enthalpy        = %10.4f\n", h);
    printf("Entropy         = %10.4f\n", s);
    printf("Cv              = %10.4f\n", cv);
    printf("Cp              = %10.4f\n", cp);
    printf("Speed Of Sound  = %10.4f\n", w);
    p = 0.0;
    d = dv;
    THERM(&t, &d, x, &p, &e, &h, &s, &cv, &cp, &w, &hjt);
    printf("Pressure        = %10.4f\n", p);
    printf("Density         = %10.4f\n", d);
    printf("Temperature     = %10.4f\n", t);
    printf("Internal Energy = %10.4f\n", e);
    printf("Enthalpy        = %10.4f\n", h);
    printf("Entropy         = %10.4f\n", s);
    printf("Cv              = %10.4f\n", cv);
    printf("Cp              = %10.4f\n", cp);
    printf("Speed Of Sound  = %10.4f\n", w);
    printf("A===============================================================\n");
    
    t = 300.0;
    p = 20000.0;
    //...Calculate d from t,p,x
    //...If phase is known: (j=1: Liquid, j=2: Vapor)
    TPRHO(&t, &p, x, &j, &tmp_int, &d, &ierr, herr, errormessagelength);
    printf("T,P,D                 %10.4f,%10.4f,%10.4f\n", t, p, d);

    t = 110;
    p = 2000.0;
    //...If phase is not known, call TPFLSH
    //...Calls to TPFLSH are much slower than TPRHO since SATT must be called first.
    //.....(If two phase, quality is returned as q)
    q = 1.0;
    TPFLSH(&t, &p, x, &d, &dl, &dv, xliq, xvap, &q, &e, &h, &s, &cv, &cp, &w, &ierr, herr, errormessagelength);
    printf("T,P,D,H,CP            %10.4f,%10.4f,%10.4f,%10.4f,%10.4f\n", t, p, d, h, cp);
    printf("T,P,D,L,Va,q          %10.4f,%10.4f,%10.4f,%10.4f,%10.4f,%10.4f\n", t, p, d, xliq[0], xvap[0], q);
    printf("===============================================================\n");
    
    //...Calculate pressure (p), internal energy (e), enthalpy (h), entropy (s),
    //.....isochoric (cv) and isobaric (cp) heat capacities, speed of sound (w),
    //.....and Joule-Thomson coefficient (hjt) from t,d,x
    //.....(subroutines THERM2 and THERM3 contain more properties, see PROP_SUB.FOR)
    THERM(&t, &d, x, &p, &e, &h, &s, &cv, &cp, &w, &hjt);

    //...Calculate pressure
    PRESS(&t, &d, x, &p);

    //...Calculate fugacity
    FGCTY(&t, &d, x, f);

    //...Calculate second and third virial coefficients
    VIRB(&t, x, &b);
    VIRC(&t, x, &c);
    printf("F,B,C                 %10.4f,%10.4f,%10.4f\n", f[0], b, c);

    //...Calculate the derivatives: dP/dD, d^2P/dD^2, dP/dT  (D indicates density)
    //...(dD/dP, dD/dT, and dB/dT are also available, see PROP_SUB.FOR)
    DPDD(&t, &d, x, &dpdrho);
    DPDD2(&t, &d, x, &d2pdd2);
    DPDT(&t, &d, x, &dpdt);
    printf("dP/dD,d2P/dD2,dP/dT   %10.4f,%10.4f,%10.4f\n", dpdrho, d2pdd2, dpdt);


    //...Calculate derivatives of enthalpy with respect to T, P, and D
    DHD1(&t, &d, x, &dhdt_d, &dhdt_p, &dhdd_t, &dhdd_p, &dhdp_t, &dhdp_d);
    printf("Enthalpy derivatives  %10.4f,%10.4f,%10.4f,%10.4f,%10.4f\n",
            dhdt_d, dhdt_p, dhdd_t, dhdd_p / 1000.0, dhdp_t);
    //...Calculate surface tension
    SURFT(&t, &dl, x, &sigma, &ierr, herr, errormessagelength);
    printf("T,SURF. TN.           %10.4f,%10.4f\n", t, sigma);

    //...Calculate viscosity (eta) and thermal conductivity (tcx)
    TRNPRP(&t, &d, x, &eta, &tcx, &ierr, herr, errormessagelength);
    printf("VIS.,TH.CND.          %10.4f,%10.4f\n", eta, tcx * 1000.0);

    //...General property calculation with inputs of t,d,x
    TDFLSH(&t, &d, x, &pp, &dl, &dv, xliq, xvap, &q, &e, &h1, &s, &cv, &cp, &w, &ierr, herr, errormessagelength);
    printf("T, D, P from TDFLSH   %10.4f,%10.4f,%10.4f\n", t, d, pp / 1000.0);

    //...General property calculation with inputs of p,d,x
    PDFLSH(&p, &d, x, &tt, &dl, &dv, xliq, xvap, &q, &e, &h1, &s, &cv, &cp, &w, &ierr, herr, errormessagelength);
    printf("T, D, P from PDFLSH   %10.4f,%10.4f,%10.4f\n", tt, d, p / 1000.0);

    //...General property calculation with inputs of p,h,x
    PHFLSH(&p, &h, x, &tt, &dd, &dl, &dv, xliq, xvap, &q, &e, &s, &cv, &cp, &w, &ierr, herr, errormessagelength);
    printf("T, D, P from PHFLSH   %10.4f,%10.4f,%10.4f\n", tt, dd, p / 1000.0);

    //...General property calculation with inputs of p,s,x
    PSFLSH(&p, &s, x, &tt, &dd, &dl, &dv, xliq, xvap, &q, &e, &h1, &cv, &cp, &w, &ierr, herr, errormessagelength);
    printf("T, D, P from PSFLSH   %10.4f,%10.4f,%10.4f\n", tt, dd, p / 1000.0);

    //...General property calculation with inputs of d,h,x
    DHFLSH(&d, &h, x, &tt, &pp, &dl, &dv, xliq, xvap, &q, &e, &s, &cv, &cp, &w, &ierr, herr, errormessagelength);
    printf("T, D, P from DHFLSH   %10.4f,%10.4f,%10.4f\n", tt, d, pp / 1000.0);

    //...General property calculation with inputs of t,h,x
    //     kr--flag specifying desired root for multi-valued inputs:
    //         1=return lower density root
    //         2=return higher density root
    int kr = 1;
    THFLSH(&t, &h, x, &kr, &pp, &dd, &dl, &dv, xliq, xvap, &q, &e, &s, &cv, &cp, &w, &ierr, herr, errormessagelength);
    printf("T, D, P from THFLSH   %10.4f,%10.4f,%10.4f\n", t, dd, pp / 1000.0);

    //...Other general property calculation routines are given in FLSH_SUB.FOR
    //...and FLASH2.FOR

    //...Calculate melting pressure
    t = 100.0;
    MELTT(&t, x, &p, &ierr, herr, errormessagelength);
    printf("Melting pressure(MPa) %10.4f,%10.4f\n", p / 1000.0, t);

    //...Calculate melting temperature
    MELTP(&p, x, &tt, &ierr, herr, errormessagelength);
    printf("Melting temperature(K)%10.4f,%10.4f\n", tt, p / 1000.0);

    //...Calculate sublimation pressure
    t = 200.0;
    SUBLT(&t, x, &p, &ierr, herr, errormessagelength);
    printf("Sublimation pr.(kPa)  %10.4f,%10.4f\n", p, t);

    //...Calculate sublimation temperature
    SUBLP(&p, x, &tt, &ierr, herr, errormessagelength);
    printf("Sublimation temp.(K)  %10.4f,%10.4f\n", tt, p);

    //...Get limits of the equations and check if t,d,p is a valid point
    //...Equation of state
    //     call LIMITK ('EOS',1,t,d,p,tmin,tmax,Dmax,pmax,ierr,herr)
    //...Viscosity equation
    //     call LIMITK ('ETA',1,t,d,p,tmin,tmax,Dmax,pmax,ierr,herr)
    //...Thermal conductivity equation
    //     call LIMITK ('TCX',1,t,d,p,tmin,tmax,Dmax,pmax,ierr,herr)

    //...Other routines are given in UTILITY.FOR


    return 0;
}
