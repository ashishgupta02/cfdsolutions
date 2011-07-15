/*******************************************************************************
 * File:        MathError.c
 * Author:      Ashish Gupta
 * Revision:    4
 ******************************************************************************/

#include "License.h"

#include <stdio.h>
#include <math.h>

#include "MathError.h"
#include "Utils.h"
#include "Log.h"

int math_err_msg_level = 0;
int math_err_check_finite = 0;
int math_err_max_msg = 100;
int math_err_count = 0;
int math_err_nan = 0;

#ifndef HAS_NO_EXCEPTION
/*******************************************************************************
 * mathematical error handling -> found in math.h
 * struct exception
 *  {
 *   int type;       -) which means either DOMAIN, SING, OVERFLOW, UNDERFLOW,...
 *   char *name;
 *   NDMDouble arg1;
 *   NDMDouble arg2;
 *   NDMDouble retval;
 *  };
 *******************************************************************************/
struct exception math_err_exception = {0, NULL, 0.0, 0.0, 0.0};

/*******************************************************************************
 * report math exceptions depending on value of math_err_msg_level:
 * 0: do nothing
 * 1: save first occourance of exception
 * 2: print every exceptions (up to a number of math_err_max_msg)
 *
 * do _not_ change the name of this function, it's defined in libm
 *******************************************************************************/
int matherr(struct exception *e) {
    const int type = e->type;
    const int exit_status = 1;

    if (math_err_msg_level == 0)
        return exit_status;

    switch (type) {
        case 0: return exit_status;
        case UNDERFLOW: return exit_status;
        case OVERFLOW: return exit_status;
        case PLOSS: return exit_status;
    }

    if (math_err_count == 0)
        math_err_exception = *e;

    math_err_count++;

#ifndef VECTORCOMPUTER
    if (math_err_msg_level == 1)
        return exit_status;

    if (math_err_count <= math_err_max_msg)
        switch (type) {
            case DOMAIN:
                ndm_msg("math exception DOMAIN in %s arg1 %e arg2 %e ret %e\n",
                        e->name, e->arg1, e->arg2, e->retval);
                break;
            case SING:
                ndm_msg("math exception SING in %s arg1 %e arg2 %e ret %e\n",
                        e->name, e->arg1, e->arg2, e->retval);
                break;
            case TLOSS:
                ndm_msg("math exception TLOSS in %s arg1 %e arg2 %e ret %e\n",
                        e->name, e->arg1, e->arg2, e->retval);
                break;
            default:
                ndm_msg("math exception in %s arg1 %e arg2 %e ret %e\n",
                        e->name, e->arg1, e->arg2, e->retval);
        }
#endif

    return exit_status;

} /** matherr() **/

/*******************************************************************************
 * print result of math exception
 *******************************************************************************/
void fprint_math_err(void) {
    const int type = math_err_exception.type;

    if (math_err_count)
        switch (type) {
            case DOMAIN:
                ndm_errmsg("math exception DOMAIN in %s arg1 %e arg2 %e ret %e\n",
                        math_err_exception.name, math_err_exception.arg1,
                        math_err_exception.arg2, math_err_exception.retval);
                break;
            case SING:
                ndm_errmsg("math exception SING in %s arg1 %e arg2 %e ret %e\n",
                        math_err_exception.name, math_err_exception.arg1,
                        math_err_exception.arg2, math_err_exception.retval);
                break;
            case TLOSS:
                ndm_errmsg("math exception TLOSS in %s arg1 %e arg2 %e ret %e\n",
                        math_err_exception.name, math_err_exception.arg1,
                        math_err_exception.arg2, math_err_exception.retval);
                break;
            default:
                ndm_errmsg("math exception in %s arg1 %e arg2 %e ret %e\n",
                        math_err_exception.name, math_err_exception.arg1,
                        math_err_exception.arg2, math_err_exception.retval);
        }

    return;

} /** fprint_math_err() **/

#else /* HAS_NO_EXCEPTION */

void fprint_math_err(void) {
    return;
}

#endif

/*******************************************************************************
 * check if values in array are finite
 * when math_err_check_finite is set
 *******************************************************************************/
int array_not_finite(int npoints, NDMDouble d[]) {
    int i, count = 0, res = 0;

    if (math_err_check_finite == 0)
        return 0;

    if (math_err_check_finite == 1) {
#ifndef VECTORCOMPUTER
        {
            for (i = 0; i < npoints; i++)
                count += value_not_finite(d[i]);
            res = count ? 1 : 0;
        }
#else
        {
            NDMDouble sum = 0;
#include "nodep.h"

            for (i = 0; i < npoints; i++)
                sum += d[i];
            res = value_not_finite(sum);
        }
#endif
    }

    if (math_err_check_finite == 2) {
        for (i = 0; i < npoints; i++) {
            if (value_not_finite(d[i]) && count <= math_err_max_msg) {
                ndm_errmsg(
                        "MATH WARNING: point %d of %d points has a not finite value\n",
                        i, npoints);
                count++;
            }
        }

        if (count)
            ndm_errmsg("MATH WARNING: %d points have not finite values\n", count);

        res = count ? 1 : 0;
    }

    math_err_nan += res;
    return res;

} /**  array_not_finite() **/

/*******************************************************************************
 * check if values in array are not finite
 * when math_err_check_finite is set print location
 *******************************************************************************/
void private_check_array_is_finite(int npoints, NDMDouble d[],
        const char *var, char *file, int line) {
    if (array_not_finite(npoints, d))
        ndm_errmsg("\nFILE '%s', LINE %d:\n"
            "MATH WARNING: %s has not finite value\n",
            file, line, var);

    return;
} /**  private_check_array_is_finite() **/

