/*******************************************************************************
 * File:        StegerWarming_Fluxes.h
 * Author:      Ashish Gupta
 * Revision:    1
 * Package:     Multi-Phase Computational Fluid Dynamics Flow Solver
 ******************************************************************************/

#include "License.h"

#ifndef _STEGERWARMING_FLUXES_H
#define	_STEGERWARMING_FLUXES_H

#include "Vector3D.h"

// Steger Warming Scheme Functions
void StegerWarming_Init(void);
void StegerWarming_Finalize(void);
void StegerWarming_Reset(void);
void Compute_StegerWarmingFlux(int node_L, int node_R, Vector3D areavec, double *Flux_StegerWarming, int AddTime);
void Compute_Residual_StegerWarming(int AddTime);

#endif	/* _STEGERWARMING_FLUXES_H */

