/*******************************************************************************
 * File:        AUSM_Fluxes.h
 * Author:      Ashish Gupta
 * Revision:    4
 ******************************************************************************/

#include "License.h"

#ifndef _AUSM_FLUXES_H
#define	_AUSM_FLUXES_H

#include "Vector3D.h"

// AUSM Scheme Functions
void AUSM_Init(void);
void AUSM_Finalize(void);
void AUSM_Reset(void);
void Compute_AUSMFlux(int node_L, int node_R, Vector3D areavec, double *Flux_AUSM, int AddTime);
void Compute_Residual_AUSM(int AddTime);

#endif	/* _AUSM_FLUXES_H */

