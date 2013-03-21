/*******************************************************************************
 * File:        HLLC_Fluxes.h
 * Author:      Ashish Gupta
 * Revision:    4
 ******************************************************************************/

#include "License.h"

#ifndef _HLLC_FLUXES_H
#define	_HLLC_FLUXES_H

#include "Vector3D.h"

// HLLC Scheme Functions
void HLLC_Init(void);
void HLLC_Finalize(void);
void HLLC_Reset(void);
void Compute_HLLCFlux(int node_L, int node_R, Vector3D areavec, double *Flux_HLLC, int AddTime);
void Compute_HLLCFlux_LMFix(int node_L, int node_R, Vector3D areavec, double *Flux_HLLC, int AddTime);
void Compute_Residual_HLLC(int AddTime);

#endif	/* _HLLC_FLUXES_H */

