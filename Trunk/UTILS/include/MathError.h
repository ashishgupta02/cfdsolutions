/*******************************************************************************
 * File:        MathError.h
 * Author:      Ashish Gupta
 * Revision:    1
 ******************************************************************************/

#ifndef _MATHERROR_H
#define _MATHERROR_H

#include <math.h>
#include <stdio.h>

#include "NDM_TypeDefs.h"

#if !defined(UNDERFLOW) || !defined(OVERFLOW) || !defined(PLOSS)|| \
    !defined(DOMAIN) ||  !defined(SING) || !defined(TLOSS)
#ifndef  HAS_NO_EXCEPTION
#define HAS_NO_EXCEPTION 1
#endif
#endif

#ifdef NEEDS_IEEEFP_H
#include <ieeefp.h>
#endif

#ifdef _HPUX_SOURCE
#define finite(x) isfinite(x)
#endif


/*******************************************************************************
* Keep C++ compilers from getting confused
*******************************************************************************/
#if defined __cplusplus
extern "C" {
#endif

extern int math_err_msg_level;
extern int math_err_count;
extern int math_err_max_msg;
#ifndef HAS_NO_EXCEPTION
#if defined(__cplusplus)
struct exception {
    int type;
    char *name;
    NDMDouble arg1;
    NDMDouble arg2;
    NDMDouble retval;
};
#endif
    extern struct exception math_err_exception;
#endif
extern int math_err_check_finite;
extern int math_err_nan;

/*******************************************************************************
* print result of math exception
*******************************************************************************/
void fprint_math_err(void);

/*******************************************************************************
* check if value is not finite
*******************************************************************************/
#ifndef HAS_NO_FINITE
#define value_not_finite(d) (finite(d) ? 0 : 1)
#else
#define value_not_finite(d) 0
#endif

/*******************************************************************************
* check if values in array are not finite
* when math_err_check_finite is set
*******************************************************************************/
int array_not_finite(int npoints, NDMDouble d[]);


/*******************************************************************************
* check if values in array are not finite
* when math_err_check_finite is set print location
*******************************************************************************/
#define check_array_is_finite(npoints, array) private_check_array_is_finite( \
            (npoints), (array) , #array, __FILE__, __LINE__)


/*******************************************************************************
* check if values in array are not finite
* when math_err_check_finite is set print location
*******************************************************************************/
void private_check_array_is_finite(int npoints, NDMDouble d[],
                                   const char *var, char *file, int line);

/*******************************************************************************
* Keep C++ compilers from getting confused
*******************************************************************************/
#if defined __cplusplus
}
#endif

#endif /* _MATHERROR_H */
