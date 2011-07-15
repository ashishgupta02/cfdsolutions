/*******************************************************************************
 * File:        Utils.h
 * Author:      Ashish Gupta
 * Revision:    4
 ******************************************************************************/

#include "License.h"

#ifndef _UTILS_H
#define	_UTILS_H

#include "NDM_TypeDefs.h"
#include "Trim_Utils.h"

/*******************************************************************************
 * Keep C++ compilers from getting confused
 *******************************************************************************/
#if defined __cplusplus
extern "C" {
#endif

/*------------------------------------------------------------------------------
| Maximum markers allowed for each boundary part
------------------------------------------------------------------------------*/
#define MAX_MARKERS 4096

/*------------------------------------------------------------------------------
| set the final constants according to the setting of NDM_USE_DOUBLE,
| the constant SYSTEM_MINNUM itself is used in parts of the code, keep
| this definition temporally
------------------------------------------------------------------------------*/
#ifndef SYSTEM_MAXNUM
#ifdef NDM_USE_DOUBLE
#define SYSTEM_MAXNUM DBL_MAX
#else
#define SYSTEM_MAXNUM FLT_MAX
#endif
#endif

#ifndef SYSTEM_MINNUM
#ifdef NDM_USE_DOUBLE
#define SYSTEM_MINNUM DBL_MIN
#else
#define SYSTEM_MINNUM FLT_MIN
#endif
#endif

#ifndef SYSTEM_EPS
#ifdef NDM_USE_DOUBLE
#define SYSTEM_EPS DBL_EPSILON
#else
#define SYSTEM_EPS FLT_EPSILON
#endif
#endif

/*------------------------------------------------------------------------------
| Finds the other point when the face and one point are known
------------------------------------------------------------------------------*/
#define OPNT(fpoint, face, thisp) (fpoint[face][0] + fpoint[face][1] - (thisp))

/*******************************************************************************
 * Keep C++ compilers from getting confused
 *******************************************************************************/
#if defined __cplusplus
}
#endif

#endif	/* _UTILS_H */

