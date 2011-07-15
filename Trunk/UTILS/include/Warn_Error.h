/*******************************************************************************
 * File:        Warn_Error.h
 * Author:      Ashish Gupta
 * Revision:    4
 ******************************************************************************/

#include "License.h"

#ifndef _WARN_ERROR_H
#define _WARN_ERROR_H

#include <stdio.h>

/*******************************************************************************
* Keep C++ compilers from getting confused
*******************************************************************************/
#if defined __cplusplus
extern "C" {
#endif

/*******************************************************************************
* Equivalent to fprintf(stderr, fmt, ...) + exit(1).
* Print to stderr and terminate the program.
*******************************************************************************/
int ndm_exit(const char *fmt, ...);
extern int terminate(const char *fmt, ...);

/*******************************************************************************
* Error message and set flag to terminate later
* needed in parallel if error occurs in conditional constructs not visted
* on all domains in order to avoid hanging connections
*******************************************************************************/
extern int schedule_terminate(const char *fmt, ...);

/*******************************************************************************
* check error flag and exit: to be called outside conditional constructs
* in order to terminate all processes
*******************************************************************************/
extern int fatal_error_check(void);

/*------------------------------------------------------------------------------
| 'fatal_error()' and 'warning()' are called exactly like 'printf()'.
| They both print the filename and linenumber from which they are called in the
| source code as well as the text that comes from the arguments to 'stderr'.
| 'warning()' afterwards returns, 'fatal_error()' terminates the running program
------------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/

#define fatal_error _print_position_and_return(__FILE__, __LINE__, terminate)

#define fatal_error_flag \
        _print_position_and_return(__FILE__, __LINE__, schedule_terminate)

/*----------------------------------------------------------------------------*/

#define warning _print_position_and_return(__FILE__, __LINE__, _warning)

/*----------------------------------------------------------------------------*/

#ifdef DEBUGGING_OUTPUT
# define I_AM_HERE \
    (printf("-> %s: %d\n", __FILE__, __LINE__), fflush(stdout))
#else
# define I_AM_HERE
#endif


/*------------------------------------------------------------------------------
| The following definitions/declarations are private and should only be used by
| the macros 'fatal_error' and 'warning'.
------------------------------------------------------------------------------*/
typedef int (*_message_fun)(const char *fmt, ...);

extern _message_fun _print_position_and_return(const char *file,
                                               int line,
                                               _message_fun call);

extern int _warning(const char *fmt, ...);


/*******************************************************************************
* Keep C++ compilers from getting confused
*******************************************************************************/
#if defined __cplusplus
}
#endif

#endif /* _WARN_ERROR_H */
