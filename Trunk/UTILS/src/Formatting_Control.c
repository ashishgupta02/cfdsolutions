/*******************************************************************************
 * File:        Formatting_Control.c
 * Author:      Ashish Gupta
 * Revision:    1
 ******************************************************************************/

#include <math.h>
#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>

#include "NDM_TypeDefs.h"
#include "Formatting_Control.h"

static int length[2];

/*******************************************************************************
 * Compute the word length required for displaying all signiificant figures.
 ******************************************************************************/
void determine_significant_digits(NDMDouble val) {
    /*--------------------------------------------------------------------------
    | Minimum length for ct1==3 in original version
    --------------------------------------------------------------------------*/
    int ct1 = 3;

    NDMDouble c;
    NDMDouble b = modf(val, &c);

    while (c >= 10.0) {
        c *= 0.1;
        ct1++;
    }

    if (b > 0.0) {
        b *= 1.0e+3;

        while (b < 1.0) {
            b *= 10.0;
            ct1++;
        }
    }

    length[1] = ct1;
    length[0] = ct1 + 6;

} /** determine_significant_digits() **/

/*******************************************************************************
 *
 ******************************************************************************/
int *return_significant_digits_for_time(void) {
    return (length);
}

/*******************************************************************************
 *
 ******************************************************************************/
void ndm_sprintf(char *c, const char *format, ...) {
    va_list arglist;
    va_start(arglist, format);
    vsprintf(c, format, arglist);
    va_end(arglist);
}

