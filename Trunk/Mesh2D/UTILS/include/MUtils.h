/* 
 * File:   MUtils.h
 * Author: Ashish Gupta
 *
 * Created on February 5, 2010, 7:22 PM
 */

#ifndef _MUTILS_H
#define	_MUTILS_H

#ifdef	__cplusplus
extern "C" {
#endif

double DetMat2x2(double a, double b, double c, double d);
void InvMat2x2(double* M, double* N);
void MulMat2x2(double* M, double* N, double* O);

#ifdef	__cplusplus
}
#endif

#endif	/* _MUTILS_H */

