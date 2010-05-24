/*******************************************************************************
 * File:        Machine_Parameters.c
 * Author:      Ashish Gupta
 * Revision:    1
 ******************************************************************************/

#include <math.h>
#include <stdio.h>
#include <stdarg.h>

#include "NDM_TypeDefs.h"
#include "Utils.h"
#include "Machine_Parameters.h"

#define CONV(i) ((NDMDouble)(i))

MachineParameters mp;

/*------------------------------------------------------------------------------
|
------------------------------------------------------------------------------*/
static void machar(void);
static void init_machine_parameters(void);

/*******************************************************************************
 *
 ******************************************************************************/
MachineParameters *return_machine_parameters(void) {
    init_machine_parameters();
    machar();
    return (&mp);
} /** return_machine_parameters() **/

/*******************************************************************************
 *
 ******************************************************************************/
void show_machine_parameters(void) {
    init_machine_parameters();
    machar();
    print_machine_parameters();
} /** show_machine_parameters() **/

/*******************************************************************************
 * this routine essentially borrowed from Numerical Recipes. I dont feel
 * too bad about this because these is so much open source stuff of this
 * nature and it is easy to write anyway.
 ******************************************************************************/
static void machar(void) {
    int i, itemp, iz, j, k, mx, nxres;
    NDMDouble a, b, beta, betah, betain, one, t, temp, temp1, tempa, two, y, z, zero;

    /*--------------------------------------------------------------------------
    |
    --------------------------------------------------------------------------*/
    one = CONV(1);
    two = one + one;
    zero = one - one;
    a = one;

    /*--------------------------------------------------------------------------
    | return ibeta (floating point radix) by method of
    | M.Malcolm
    --------------------------------------------------------------------------*/
    do {
        a += a;
        temp = a + one;
        temp1 = temp - a;
    } while (EQ(temp1, one));

    b = one;
    do {
        b += b;
        temp = a + b;
        itemp = (int) (temp - a);
    } while (itemp == 0);

    mp.ibeta = itemp;
    beta = CONV(mp.ibeta);

    mp.it = 0;
    /*--------------------------------------------------------------------------
    | determine number of base-ibeta digits in fp mantissa - it
    --------------------------------------------------------------------------*/
    b = one;
    do {
        ++(mp.it);
        b *= beta;
        temp = b + one;
        temp1 = temp - b;
    } while (EQ(temp1, one));

    /*--------------------------------------------------------------------------
    | irnd - rounding check
    --------------------------------------------------------------------------*/
    mp.irnd = 0;
    betah = beta / two;
    temp = a + betah;
    if (!EQ(temp, a))
        mp.irnd = 1;
    tempa = a + beta;
    temp = tempa + betah;
    if (EQ0(mp.irnd) && !EQ(temp, tempa))
        mp.irnd = 2;
    /*--------------------------------------------------------------------------
    | return negep and epsneg
    --------------------------------------------------------------------------*/
    mp.negep = (mp.it) + 3;
    betain = one / beta;
    a = one;
    for (i = 1; i <= mp.negep; i++)
        a *= betain;
    b = a;
    for (;;) {
        temp = one - a;
        if (!EQ(temp, one))
            break;
        a *= beta;
        --(mp.negep);
    }
    mp.negep = -(mp.negep);

    /*--------------------------------------------------------------------------
    | determine eps and machep
    --------------------------------------------------------------------------*/
    mp.epsneg = a;
    mp.machep = -mp.it - 3;
    a = b;
    for (;;) {
        temp = one + a;
        if (!EQ(temp, one))
            break;
        a *= beta;
        ++(mp.machep);
    }
    mp.eps = a;

    /*--------------------------------------------------------------------------
    | determine ngrd
    --------------------------------------------------------------------------*/
    mp.ngrd = 0;
    temp = one + mp.eps;
    if (EQ0(mp.irnd) && !EQ(temp * one, one))
        mp.ngrd = 1;
    i = 0;
    k = 1;
    z = betain;
    t = one + mp.eps;
    nxres = 0;

    /*--------------------------------------------------------------------------
    | loop till underflow
    --------------------------------------------------------------------------*/
    for (;;) {
        y = z;
        z = y*y;
        a = z*one;
        temp = z*t;
        if (EQ0(a + a) || GE(fabs(z), y))
            break;
        temp1 = temp*betain;
        if (EQ(temp1 * beta, z))
            break;
        ++i;
        k += k;
    }

    /*--------------------------------------------------------------------------
    | for non decimal machines (base 10)
    --------------------------------------------------------------------------*/
    if (mp.ibeta != 10) {
        mp.iexp = i + 1;
        mx = k + k;
    } else
        /*----------------------------------------------------------------------
        | for decimal machines
        ----------------------------------------------------------------------*/
    {
        mp.iexp = 2;
        iz = mp.ibeta;
        while (k >= iz) {
            iz *= mp.ibeta;
            ++(mp.iexp);
        }
        mx = iz + iz - 1;
    }

    /*--------------------------------------------------------------------------
    | determine minexp, xmin by looping till underflow
    --------------------------------------------------------------------------*/
    for (;;) {
        mp.xmin = y;
        y *= betain;
        a = y*one;
        temp = y*t;
        if (!EQ0(a + a) && LT(fabs(y), mp.xmin)) {
            ++k;
            temp1 = temp*betain;
            if (EQ(temp1 * beta, y) && !EQ(temp, y)) {
                nxres = 3;
                mp.xmin = y;
                break;
            }
        } else
            break;
    }
    mp.minexp = -k;
    if (mx <= k + k - 3 && mp.ibeta != 10) {
        mx += mx;
        ++(mp.iexp);
    }
    /*--------------------------------------------------------------------------
    | determine maxexp, xmax
    --------------------------------------------------------------------------*/
    mp.maxexp = mx + (mp.minexp);
    mp.irnd += nxres;
    if (mp.irnd >= 2)
        mp.maxexp -= 2;
    i = (mp.maxexp)+(mp.minexp);
    /*--------------------------------------------------------------------------
    | adjust for machines with implicit leding bit in mantissa, and
    | machines with radix  point at extreme right of mantissa
    --------------------------------------------------------------------------*/
    if (mp.ibeta == 2 && !i)
        --(mp.maxexp);
    if (i > 20)
        --(mp.maxexp);
    if (!EQ(a, y))
        mp.maxexp -= 2;
    mp.xmax = one - (mp.epsneg);
    if (!EQ((mp.xmax) * one, mp.xmax))
        mp.xmax = one - beta * (mp.epsneg);
    mp.xmax /= (mp.xmin * beta * beta * beta);
    i = (mp.maxexp)+(mp.minexp) + 3;
    for (j = 1; j <= i; j++) {
        if (mp.ibeta == 2)
            mp.xmax += mp.xmax;
        else
            mp.xmax *= beta;
    }
} /** machar() **/

/*******************************************************************************
 *
 ******************************************************************************/
void print_machine_parameters(void) {
    char interval[] = "**********************************************************";
    init_machine_parameters();
    machar();

    ndm_msg("%s\n", interval);
    ndm_msg("%s\n", interval);
    ndm_msg("A floating point number is given by s x M x B^(e-E)\n");
    ndm_msg("where s is a sign bit (+/-)\n");
    ndm_msg("  and M is an exact positive integer Mantissa\n");
    ndm_msg("  and B is the base of the representation (usually Base 2)\n");
    ndm_msg("  and e is an exact integer exponent\n");
    ndm_msg("  and E is the bias of the exponent\n");
    ndm_msg("%s\n", interval);
    ndm_msg("%s\n", interval);
    ndm_msg("YOUR MACHINE PROVIDES THE FOLLOWING\n");
    ndm_msg("%s\n", interval);
    ndm_msg("base of representation:                                %d\n", mp.ibeta);
    ndm_msg("%s\n", interval);
    ndm_msg("number of base digits in mantissa:                     %d\n", mp.it);
    ndm_msg("%s\n", interval);
    ndm_msg("smallest exponent p  so that B^p+1 != 1:               %d\n", mp.machep);
    ndm_msg("%s\n", interval);
    ndm_msg("B^p or floating point precision (fpp):                 %g\n", mp.eps);
    ndm_msg("%s\n", interval);
    ndm_msg("smallest p1 so that B^p1-1 != 1                        %d\n", mp.negep);
    ndm_msg("%s\n", interval);
    ndm_msg("B^p1 or another measure of fpp:                        %g\n", mp.epsneg);
    ndm_msg("%s\n", interval);
    ndm_msg("number of bits in exponent:                            %d\n", mp.iexp);
    ndm_msg("%s\n", interval);
    ndm_msg("smallest p2 such that B^p2 has no lead 0 in mantissa:  %d\n", mp.minexp);
    ndm_msg("%s\n", interval);
    ndm_msg("B^p2 (smallest usable floating point number:           %g\n", mp.xmin);
    ndm_msg("%s\n", interval);
    ndm_msg("smallest p3 such B^p3 causes overflow:                 %d\n", mp.maxexp);
    ndm_msg("%s\n", interval);
    ndm_msg("largest usable floating point value:                   %g\n", mp.xmax);
    ndm_msg("%s\n", interval);
    ndm_msg("Number of guard digits when trucating prod. of 2 mant. %d\n", mp.ngrd);
    ndm_msg("%s\n", interval);
    switch (mp.irnd) {
        case (2):
            ndm_msg("Your machine rounds to IEEE standard\n");
            break;
        case (5):
            ndm_msg("Your machine rounds to IEEE standard\n");
            break;
        case (1):
            ndm_msg("Your machine rounds, but not to IEEE standard\n");
            break;
        case (4):
            ndm_msg("Your machine rounds, but not to IEEE standard\n");
            break;
        case (0):
            ndm_msg("BEWARE - YOUR MACHINE TRUNCATES RESULTS\n");
            break;
        case (3):
            ndm_msg("BEWARE - YOUR MACHINE TRUNCATES RESULTS\n");
            break;
        default:
            ndm_msg("Unrecognised result of roud off test - see NDM developers\n");
    }
    ndm_msg("%s\n", interval);

} /** print_machine_parameters() **/

/*******************************************************************************
 *
 ******************************************************************************/
static void init_machine_parameters(void) {
    mp.ibeta = -1;
    mp.it = -1;
    mp.machep = -1;
    mp.eps = -1;
    mp.negep = -1;
    mp.epsneg = -1;
    mp.iexp = -1;
    mp.minexp = -1;
    mp.xmin = -1;
    mp.maxexp = -1;
    mp.xmax = -1;
    mp.irnd = -1;
    mp.ngrd = -1;

} /** init_machine_parameters() **/

#undef CONV

