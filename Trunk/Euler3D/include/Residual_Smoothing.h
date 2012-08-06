/*******************************************************************************
 * File:        Residual_Smoothing.h
 * Author:      Ashish Gupta
 * Revision:    4
 ******************************************************************************/

#include "License.h"

#ifndef _RESIDUAL_SMOOTHING_H
#define	_RESIDUAL_SMOOTHING_H

extern double *RM_MaxEigenValue;
extern double *RM_SumMaxEigenValue;

// Residual Smoothing Functions
void Residual_Smoothing_Init(void);
void Residual_Smoothing_Finalize(void);
void Residual_Smoothing_Reset(void);
void Residual_Smoothing(void);

#endif	/* _RESIDUAL_SMOOTHING_H */

