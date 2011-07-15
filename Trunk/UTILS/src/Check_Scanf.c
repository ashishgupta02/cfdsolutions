/*******************************************************************************
 * File:        Check_Scanf.c
 * Author:      Ashish Gupta
 * Revision:    4
 ******************************************************************************/

#include "License.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>

#include "Check_String.h"
#include "NDM_Dimensions.h"
#include "Utils.h"
#include "Log.h"
#include "Check_Scanf.h"

/*******************************************************************************
 * an alternative to the standard scanf library function:
 *  -  this function use fgets to read the line and vsscanf to convert it
 *  -  the given input line has to match the format list
 *  -  the routine does not stop until the requred input is given
 * ATTENTION
 *  if ISOC89 is set only %d, %g %e %f and %s is supported by this function
 *******************************************************************************/
int check_scanf(char *fmt, ...) {
    char buffer[NDM_BUF_SZ];
    char *f = fmt;
    char *cdum;

    int got_it = 0;
    int no_fmt_spec = 0;

    va_list args;
    va_start(args, fmt);

    /*--------------------------------------------------------------------------
    | count the number of conversion items in the format string
    --------------------------------------------------------------------------*/
    while ((f = strstr(f, "%")) != NULL) {
        f += 2;
        no_fmt_spec++;
    }

    /*--------------------------------------------------------------------------
    | read the line until the input is correct and complete,
    | the ISOC99 function vsscanf is used if the ISOC99 extensions are available
    --------------------------------------------------------------------------*/
#ifdef __USE_ISOC99
    while (!got_it) {
        cdum = fgets(buffer, NDM_BUF_SZ - 1, stdin);

        if (vsscanf(buffer, fmt, args) == no_fmt_spec)
            got_it++;
        else
            ndm_msg("\nWrong or incomplete input, retry: ");
    }
#else
    /*--------------------------------------------------------------------------
    | read the line until the input is correct and complete,
    | the ISOC99 function vsscanf is replaced by a rough implementation
    | with reduced functionality
    --------------------------------------------------------------------------*/
    while (!got_it) {
        char tok_buf[NDM_BUF_SZ];
        char tok_fmt[NDM_BUF_SZ];
        char *s, *t;

        int no_read = 0;
        int i;

        /*----------------------------------------------------------------------
        | read the line and cut leading blanks and the first % in format
        ----------------------------------------------------------------------*/
        fgets(buffer, NDM_BUF_SZ - 1, stdin);

        s = buffer;
        s += strspn(s, " ");

        t = fmt;
        t += strspn(t, " ");
        t = strbftok(t, tok_fmt, NDM_BUF_SZ, "%");

        /*----------------------------------------------------------------------
        | convert all string using sscanf
        ----------------------------------------------------------------------*/
        for (i = 0; i < no_fmt_spec; i++) {
            s = strbftok(s, tok_buf, NDM_BUF_SZ, " ,");
            t = strbftok(t, tok_fmt, NDM_BUF_SZ, "%");

            if (strpbrk(tok_fmt, "d"))
                no_read += sscanf(tok_buf, "%d", va_arg(args, int *));
            else if (strpbrk(tok_fmt, "fge"))
                no_read += sscanf(tok_buf, "%lf", va_arg(args, double *));
            else if (strpbrk(tok_fmt, "s"))
                no_read += sscanf(tok_buf, "%s", va_arg(args, char *));
        }

        /*----------------------------------------------------------------------
        | check for valid input,
        | the variable argument list has to be reopened in case of invalid input
        ----------------------------------------------------------------------*/
        if (no_read == no_fmt_spec)
            got_it++;
        else {
            ndm_msg("\nWrong or incomplete input, retry: ");
            va_end(args);
            va_start(args, fmt);
        }
    }
#endif

    /*--------------------------------------------------------------------------
    | close the va list and return the number of variables read
    --------------------------------------------------------------------------*/
    va_end(args);

    return no_fmt_spec;

} /** check_scanf(char *fmt, ...) **/

