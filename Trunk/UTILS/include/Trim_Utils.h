/*******************************************************************************
 * File:        Trim_Utils.h
 * Author:      Ashish Gupta
 * Revision:    4
 ******************************************************************************/

#include "License.h"

#ifndef _TRIM_UTILS_H
#define	_TRIM_UTILS_H

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>
#include <float.h>

#ifdef HAVE_MPI
#include <mpi.h>
#endif

/*******************************************************************************
 * Keep C++ compilers from getting confused
 *******************************************************************************/
#if defined __cplusplus
extern "C" {
#endif

/*------------------------------------------------------------------------------
| Definition to avoid the use of FILENAME_MAX which is possibly too small
------------------------------------------------------------------------------*/
#define MAX_FILENAME_LEN 1024

/*------------------------------------------------------------------------------
| define TRUE and FALSE conditionals
------------------------------------------------------------------------------*/
#ifndef TRUE
#define TRUE 1
#endif

#ifndef FALSE
#define FALSE 0
#endif

/*------------------------------------------------------------------------------
| define OK and ERROR conditionals
------------------------------------------------------------------------------*/
#ifndef OK
#define OK 0
#endif

#ifndef ERROR
#define ERROR 1
#endif

/*------------------------------------------------------------------------------
| define MAXIMUM AND MINIMUM FLOATS/DOUBLES TO AVOID OVERFLOW/UNDERFLOW
| these are system dependent constants
------------------------------------------------------------------------------*/
#ifndef FLT_MIN
#define FLT_MIN 1.e-20
#endif

#ifndef FLT_MAX
#define FLT_MAX 1.e+20
#endif

#ifndef DBL_MIN
#define DBL_MIN 1.e-200
#endif

#ifndef DBL_MAX
#define DBL_MAX 1.e+200
#endif

#ifndef FLT_EPSILON
#define FLT_EPSILON 1.e-7
#endif

#ifndef DBL_EPSILON
#define DBL_EPSILON 1.e-16
#endif

/*------------------------------------------------------------------------------
|
------------------------------------------------------------------------------*/
#ifndef SQR
#define SQR(x) ((x) * (x))
#endif

/*------------------------------------------------------------------------------
|
------------------------------------------------------------------------------*/
#ifndef MIN
#define MIN(a, b) ((a) < (b) ? (a) : (b))
#endif

/*------------------------------------------------------------------------------
|
------------------------------------------------------------------------------*/
#ifndef MAX
#define MAX(a, b) ( (a) > (b) ? (a) : (b))
#endif

/*------------------------------------------------------------------------------
|
------------------------------------------------------------------------------*/
#ifndef MAX0
#define MAX0(a) ((a) > 0 ? (a) : 0)
#endif

/*------------------------------------------------------------------------------
|
------------------------------------------------------------------------------*/
#ifndef MIN0
#define MIN0(a) ((a) < 0 ? (a) : 0)
#endif

/*------------------------------------------------------------------------------
|
------------------------------------------------------------------------------*/
#ifndef MIN3
#define MIN3(a,b,c) ((a < b) ? ((a < c) ? a : c) : ((b < c) ? b : c))
#endif

/*------------------------------------------------------------------------------
|
------------------------------------------------------------------------------*/
#ifndef MAX3
#define MAX3(a,b,c) ((a > b) ? ((a > c) ? a : c) : ((b > c) ? b : c))
#endif

/*------------------------------------------------------------------------------
|
------------------------------------------------------------------------------*/
#ifndef ABS
#define ABS(a) ((a) < 0 ? -1 * (a) : (a) )
#endif

/*------------------------------------------------------------------------------
|
------------------------------------------------------------------------------*/
#ifndef ISIGN
#define ISIGN(a) ((a) < 0 ? -1  : 1 )
#endif

/*------------------------------------------------------------------------------
|
------------------------------------------------------------------------------*/
#ifndef SIGN
#define SIGN(a,b)     ( ((b) >= 0.0) ? fabs(a) : -fabs(a))
#endif

/*------------------------------------------------------------------------------
| Vectoriseable power function x^y
------------------------------------------------------------------------------*/
#ifndef POW
#define POW(x, y) exp( (y) * log(x) )
#endif

/*------------------------------------------------------------------------------
| Vectoriseable power function x^y without overflow or underflow depending on x
| for powers 1e-30 <= y <= 10 and x <= 1e30
------------------------------------------------------------------------------*/
#ifndef POW_LIM
#define POW_LIM(x, y) exp( (y) * log(MAX((x), 1e-30)) )
#endif

/*------------------------------------------------------------------------------
|
------------------------------------------------------------------------------*/
#ifndef atanh
#define atanh(x) (0.5 * log((1. + (x)) / (1. - (x))))
#endif

/*------------------------------------------------------------------------------
| Some vendors define M_PI in 'math.h' or 'float.h', others do not!
------------------------------------------------------------------------------*/
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/*------------------------------------------------------------------------------
| Convert form RAD TO DEG and DEG TO RAD
------------------------------------------------------------------------------*/
#ifndef RAD_TO_DEG
#define RAD_TO_DEG(x) ((x) * 57.29577951308232)
#endif

#ifndef DEG_TO_RAD
#define DEG_TO_RAD(x) ((x) * 0.017453292519943295)
#endif

/*------------------------------------------------------------------------------
|
------------------------------------------------------------------------------*/
#ifndef ACOS
#define ACOS(x) acos(MIN(MAX(x, -1.0), 1.0))
#endif

/*------------------------------------------------------------------------------
| Base comparison macros
------------------------------------------------------------------------------*/
#define EQUAL(a,b)   (ABS((a) - (b)) <= (SYSTEM_EPS * MAX(ABS(a), ABS(b))) \
                                      ? TRUE : FALSE )
#define GREATER(a,b) ( (a) - (b)      > (SYSTEM_EPS * MAX(ABS(a), ABS(b))) \
                                      ? TRUE : FALSE )

/*------------------------------------------------------------------------------
| Since (a==0) occurs often we build a special macro for this case
------------------------------------------------------------------------------*/
#define EQ0(a)    (ABS((a)) < (10.0 * SYSTEM_MINNUM)     ? TRUE : FALSE )

/*-----------------------------------------------------------------------------
| Derived comparison MACROS
------------------------------------------------------------------------------*/
#define EQ(a,b) (EQUAL((a),(b)))
#define GE(a,b) (EQUAL((a),(b))    ||  GREATER((a),(b)))
#define LE(a,b) (EQUAL((a),(b))    || !GREATER((a),(b)))
#define GT(a,b) (GREATER((a),(b))  && !EQUAL((a),(b)))
#define LT(a,b) (!GREATER((a),(b)) && !EQUAL((a),(b)))

/*------------------------------------------------------------------------------
|
------------------------------------------------------------------------------*/
#define NDM_ZERO          0.0e+0
#define DBL_ZERO          1.0E-15
#define DBL_TOLERANCE     1.0E-7
#define DBL_RES_TOLERANCE 1.0E-10

#ifndef _ASSERT
#define _ASSERT(x) {assert(x);}
#endif

/* Define Function Libraries */
/* Print error message to stderr and exit */
void error(const char *fmt, ...);
/* Print warning message to stderr and continue */
void warn(const char *fmt, ...);
/* Print information message to stdout and continue */
void info(const char *fmt, ...);
/* Returns the size of file on disk */
double file_size(const char *fname);
/* Check if a file exists */
int file_exists(char *file);
/* Checks if pathname exists and is executable */
int is_executable(char *file);
/* Locate and build pathname to executable */
char *find_executable(char *exename);
/* Get pathname to a file */
char *find_file(char *filename, char *exename);
/* Check if 2 files are the same */
int same_file(char *file1, char *file2);
/* Create a temporary file */
int temporary_file(char *name);
/* Removes the link name */
int unlink_name(char *name);
/* Make a copy of a file */
void copy_file(char *oldfile, char *newfile);
/* Blank the string */
void str_blank(char *str);
/* Sorting Functions */
void sort_low2high(int n, int list[], double f[]);
void sort_high2low(int n, int list[], double f[]);

/*******************************************************************************
 * Keep C++ compilers from getting confused
 *******************************************************************************/
#if defined __cplusplus
}
#endif

#endif	/* _TRIM_UTILS_H */

