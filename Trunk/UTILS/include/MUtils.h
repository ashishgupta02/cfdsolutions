/*******************************************************************************
 * File:        MUtils.h
 * Author:      Ashish Gupta
 * Revision:    4
 ******************************************************************************/

#include "License.h"

#ifndef _MUTILS_H
#define	_MUTILS_H

#ifdef	__cplusplus
extern "C" {
#endif

double DetMat2x2(double a, double b, double c, double d);
void   InvMat2x2(double* M, double* N);
void   MulMat2x2(double* M, double* N, double* O);
double DetMat3X3(double m[3][3]);

#ifdef	__cplusplus
}
#endif

#endif	/* _MUTILS_H */

